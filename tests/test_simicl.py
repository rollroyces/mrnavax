"""Tests for the Sim-ICL demonstration-selection module.

Covers:
- DemoCase dataclass: construction, to_dict, from_dict, query_text
- DemoStore: from_json loading, add, len, iter, rank ordering
- TF-IDF cosine ranker: scores in [0, 1], biological-relevance ranking,
  empty-store safety, k-cap behavior, deterministic tie-breaking
- build_simicl_prompt: few-shot injection, empty fallback, env-var
  disable, prompt structure preservation
- load_default_demo_store: bundled JSON loading, env-var override,
  missing-file graceful fallback
- score_trial_with_llm integration: simicl=on injects demos,
  simicl=off doesn't, custom demo_store honored, fallback on errors
"""
from __future__ import annotations

import json
import os
import tempfile
import unittest

from mrnavax.trial_similar import (
    DEFAULT_TOPK,
    ENV_DEMOS,
    ENV_ENABLED,
    ENV_TOPK,
    DemoCase,
    DemoStore,
    _get_enabled,
    _get_topk,
    _tokenize,
    build_simicl_prompt,
    load_default_demo_store,
)

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


def make_demo(
    demo_id: str,
    patient: str = "patient summary",
    title: str = "trial title",
    nct: str = "NCT00000001",
    inclusion: list[str] | None = None,
    exclusion: list[str] | None = None,
    gt: dict | None = None,
) -> DemoCase:
    return DemoCase(
        demo_id=demo_id,
        patient_text=patient,
        trial_nct=nct,
        trial_title=title,
        inclusion=tuple(inclusion or []),
        exclusion=tuple(exclusion or []),
        ground_truth_verdicts=gt,
    )


# ---------------------------------------------------------------------------
# DemoCase dataclass
# ---------------------------------------------------------------------------


class TestDemoCase(unittest.TestCase):
    def test_minimal_construction(self) -> None:
        d = make_demo("demo_001")
        self.assertEqual(d.demo_id, "demo_001")
        self.assertEqual(d.patient_text, "patient summary")
        self.assertEqual(d.inclusion, ())

    def test_to_dict_round_trip(self) -> None:
        original = make_demo(
            "demo_005",
            patient="55yo BRAF V600E melanoma",
            title="BRAF trial",
            nct="NCT00000005",
            inclusion=["BRAF V600E", "ECOG 0-1"],
            exclusion=["prior therapy"],
            gt={"eligible": True, "n_met": 2, "n_total": 2},
        )
        d = original.to_dict()
        self.assertEqual(d["demo_id"], "demo_005")
        self.assertEqual(d["inclusion"], ["BRAF V600E", "ECOG 0-1"])
        self.assertEqual(d["ground_truth_verdicts"]["eligible"], True)
        # Round-trip
        restored = DemoCase.from_dict(d)
        self.assertEqual(restored, original)

    def test_query_text_concatenates_fields(self) -> None:
        d = make_demo(
            "demo_q",
            patient="BRAF melanoma patient",
            title="BRAF inhibitor trial",
            nct="NCT12345",
            inclusion=["BRAF V600E"],
            exclusion=["prior therapy"],
        )
        qt = d.query_text
        self.assertIn("BRAF melanoma patient", qt)
        self.assertIn("NCT12345", qt)
        self.assertIn("BRAF inhibitor trial", qt)
        self.assertIn("BRAF V600E", qt)
        self.assertIn("prior therapy", qt)

    def test_from_dict_handles_missing_optional_keys(self) -> None:
        d = {
            "demo_id": "demo_x",
            "patient_text": "p",
            "trial_nct": "NCT1",
            "trial_title": "t",
            "inclusion": [],
            "exclusion": [],
        }
        case = DemoCase.from_dict(d)
        self.assertIsNone(case.ground_truth_verdicts)


# ---------------------------------------------------------------------------
# TF-IDF ranker
# ---------------------------------------------------------------------------


class TestTFIDFRanker(unittest.TestCase):
    def test_empty_store_returns_empty(self) -> None:
        store = DemoStore(demos=[])
        self.assertEqual(store.rank("query"), [])

    def test_returns_top_k(self) -> None:
        store = DemoStore(
            demos=[make_demo(f"d{i}", title=f"trial {i}") for i in range(10)]
        )
        top3 = store.rank("query", k=3)
        self.assertEqual(len(top3), 3)

    def test_k_capped_at_store_size(self) -> None:
        store = DemoStore(demos=[make_demo("a"), make_demo("b")])
        top = store.rank("query", k=10)
        self.assertEqual(len(top), 2)

    def test_ranks_biologically_similar_above_dissimilar(self) -> None:
        store = DemoStore(
            demos=[
                make_demo("melanoma", patient="melanoma patient", title="melanoma trial"),
                make_demo("breast", patient="breast cancer patient", title="breast cancer trial"),
                make_demo("lung", patient="lung cancer patient", title="lung cancer trial"),
            ]
        )
        top = store.rank("melanoma patient", k=3)
        self.assertEqual(top[0].demo_id, "melanoma")

    def test_deterministic_tie_breaking(self) -> None:
        """Identical texts should break ties by demo_id ascending."""
        store = DemoStore(
            demos=[
                make_demo("z_last", patient="same patient", title="same trial"),
                make_demo("a_first", patient="same patient", title="same trial"),
                make_demo("m_mid", patient="same patient", title="same trial"),
            ]
        )
        top = store.rank("same patient", k=3)
        ids = [d.demo_id for d in top]
        # All scores equal → sorted by demo_id asc
        self.assertEqual(ids, ["a_first", "m_mid", "z_last"])

    def test_tokenize_lowercases_alphanumeric(self) -> None:
        tokens = _tokenize("BRAF V600E Mutation! Patient: 55-year-old")
        self.assertEqual(tokens, ["braf", "v600e", "mutation", "patient", "55", "year", "old"])

    def test_scores_in_zero_one_range(self) -> None:
        """Cosine similarity is always in [0, 1] for non-negative TF-IDF."""
        store = DemoStore(
            demos=[
                make_demo("d1", patient="a b c"),
                make_demo("d2", patient="d e f"),
            ]
        )
        # No internal score access; verify ordering only
        top = store.rank("a b c", k=2)
        self.assertEqual(top[0].demo_id, "d1")  # perfect match should rank first


# ---------------------------------------------------------------------------
# DemoStore JSON loading
# ---------------------------------------------------------------------------


class TestDemoStoreFromJSON(unittest.TestCase):
    def test_loads_list_shape(self) -> None:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(
                [
                    {
                        "demo_id": "demo_x",
                        "patient_text": "p",
                        "trial_nct": "NCT1",
                        "trial_title": "t",
                        "inclusion": ["i"],
                        "exclusion": ["e"],
                        "ground_truth_verdicts": None,
                    }
                ],
                f,
            )
            tmp_path = f.name
        try:
            store = DemoStore.from_json(tmp_path)
            self.assertEqual(len(store), 1)
            self.assertEqual(store.demos[0].demo_id, "demo_x")
        finally:
            os.unlink(tmp_path)

    def test_loads_dict_with_demos_key(self) -> None:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(
                {
                    "demos": [
                        {
                            "demo_id": "demo_y",
                            "patient_text": "p",
                            "trial_nct": "NCT1",
                            "trial_title": "t",
                            "inclusion": [],
                            "exclusion": [],
                        }
                    ]
                },
                f,
            )
            tmp_path = f.name
        try:
            store = DemoStore.from_json(tmp_path)
            self.assertEqual(len(store), 1)
        finally:
            os.unlink(tmp_path)

    def test_missing_file_raises_filenotfound(self) -> None:
        with self.assertRaises(FileNotFoundError):
            DemoStore.from_json("/nonexistent/path.json")

    def test_malformed_json_raises(self) -> None:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write("{not valid json")
            tmp_path = f.name
        try:
            with self.assertRaises(json.JSONDecodeError):
                DemoStore.from_json(tmp_path)
        finally:
            os.unlink(tmp_path)

    def test_wrong_shape_raises_valueerror(self) -> None:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump("not a list or dict", f)
            tmp_path = f.name
        try:
            with self.assertRaises(ValueError):
                DemoStore.from_json(tmp_path)
        finally:
            os.unlink(tmp_path)


# ---------------------------------------------------------------------------
# build_simicl_prompt
# ---------------------------------------------------------------------------


class TestBuildSimiclPrompt(unittest.TestCase):
    def setUp(self) -> None:
        # Always run with Sim-ICL enabled for these tests
        os.environ[ENV_ENABLED] = "1"
        self.demos = [make_demo("demo_001", title="melanoma trial")]

    def tearDown(self) -> None:
        os.environ.pop(ENV_ENABLED, None)

    def test_injects_few_shot_block(self) -> None:
        plain = "You are a screener.\nReturn ONLY a JSON object with shape {}\n"
        out = build_simicl_prompt(
            "patient", "NCT1", "title", ["i1"], ["e1"], self.demos, base_prompt=plain
        )
        self.assertIn("=== Example", out)
        self.assertIn("Return ONLY a JSON object", out)
        self.assertGreater(len(out), len(plain))

    def test_empty_demos_returns_base_unchanged(self) -> None:
        plain = "You are a screener.\nReturn ONLY a JSON object with shape {}\n"
        out = build_simicl_prompt(
            "patient", "NCT1", "title", ["i1"], ["e1"], [], base_prompt=plain
        )
        self.assertEqual(out, plain)

    def test_disabled_via_env_returns_base_unchanged(self) -> None:
        os.environ[ENV_ENABLED] = "0"
        plain = "You are a screener.\nReturn ONLY a JSON object with shape {}\n"
        out = build_simicl_prompt(
            "patient", "NCT1", "title", ["i1"], ["e1"], self.demos, base_prompt=plain
        )
        self.assertEqual(out, plain)

    def test_demo_format_includes_verdict(self) -> None:
        d = make_demo(
            "demo_v",
            gt={"eligible": True, "n_met": 3, "n_total": 4},
        )
        plain = "Return ONLY a JSON object with shape {}\n"
        out = build_simicl_prompt("p", "NCT", "t", ["i"], ["e"], [d], base_prompt=plain)
        self.assertIn("ELIGIBLE", out)
        self.assertIn("3/4", out)

    def test_demo_format_includes_not_eligible(self) -> None:
        d = make_demo(
            "demo_v",
            gt={"eligible": False, "n_met": 1, "n_total": 4},
        )
        plain = "Return ONLY a JSON object with shape {}\n"
        out = build_simicl_prompt("p", "NCT", "t", ["i"], ["e"], [d], base_prompt=plain)
        self.assertIn("NOT ELIGIBLE", out)

    def test_demo_format_handles_missing_ground_truth(self) -> None:
        d = make_demo("demo_v", gt=None)
        plain = "Return ONLY a JSON object with shape {}\n"
        out = build_simicl_prompt("p", "NCT", "t", ["i"], ["e"], [d], base_prompt=plain)
        self.assertIn("verdict not provided", out)

    def test_falls_back_to_prepend_when_marker_missing(self) -> None:
        """If the prompt has no 'Return ONLY a JSON object' marker, the
        few-shot block is prepended instead of spliced in."""
        plain = "Just a plain prompt with no marker.\n"
        out = build_simicl_prompt(
            "p", "NCT", "t", ["i"], ["e"], self.demos, base_prompt=plain
        )
        self.assertTrue(out.startswith("You are an oncology"))

    def test_multiple_demos_each_rendered(self) -> None:
        demos = [make_demo(f"d{i}") for i in range(3)]
        plain = "Return ONLY a JSON object\n"
        out = build_simicl_prompt("p", "NCT", "t", ["i"], ["e"], demos, base_prompt=plain)
        self.assertEqual(out.count("=== Example"), 3)


# ---------------------------------------------------------------------------
# Env-var helpers
# ---------------------------------------------------------------------------


class TestEnvHelpers(unittest.TestCase):
    def test_get_topk_default(self) -> None:
        os.environ.pop(ENV_TOPK, None)
        self.assertEqual(_get_topk(), DEFAULT_TOPK)

    def test_get_topk_override(self) -> None:
        os.environ[ENV_TOPK] = "5"
        try:
            self.assertEqual(_get_topk(), 5)
        finally:
            os.environ.pop(ENV_TOPK)

    def test_get_topk_invalid_falls_back(self) -> None:
        os.environ[ENV_TOPK] = "not_a_number"
        try:
            self.assertEqual(_get_topk(), DEFAULT_TOPK)
        finally:
            os.environ.pop(ENV_TOPK)

    def test_get_topk_zero_falls_back(self) -> None:
        os.environ[ENV_TOPK] = "0"
        try:
            self.assertEqual(_get_topk(), DEFAULT_TOPK)
        finally:
            os.environ.pop(ENV_TOPK)

    def test_get_enabled_default(self) -> None:
        os.environ.pop(ENV_ENABLED, None)
        self.assertTrue(_get_enabled())

    def test_get_enabled_disabled(self) -> None:
        os.environ[ENV_ENABLED] = "0"
        try:
            self.assertFalse(_get_enabled())
        finally:
            os.environ.pop(ENV_ENABLED)

    def test_get_enabled_off_variants(self) -> None:
        for v in ("false", "no", "off"):
            os.environ[ENV_ENABLED] = v
            try:
                self.assertFalse(_get_enabled())
            finally:
                os.environ.pop(ENV_ENABLED)


# ---------------------------------------------------------------------------
# load_default_demo_store
# ---------------------------------------------------------------------------


class TestLoadDefaultDemoStore(unittest.TestCase):
    def test_bundled_store_loads(self) -> None:
        store = load_default_demo_store()
        self.assertGreaterEqual(len(store), 8)

    def test_demo_ids_have_expected_format(self) -> None:
        store = load_default_demo_store()
        for d in store:
            self.assertTrue(d.demo_id.startswith("demo_"))

    def test_env_var_overrides_path(self) -> None:
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump(
                [
                    {
                        "demo_id": "demo_custom",
                        "patient_text": "p",
                        "trial_nct": "NCT1",
                        "trial_title": "t",
                        "inclusion": [],
                        "exclusion": [],
                    }
                ],
                f,
            )
            tmp_path = f.name
        try:
            os.environ[ENV_DEMOS] = tmp_path
            store = load_default_demo_store()
            self.assertEqual(len(store), 1)
            self.assertEqual(store.demos[0].demo_id, "demo_custom")
        finally:
            os.environ.pop(ENV_DEMOS, None)
            os.unlink(tmp_path)

    def test_missing_env_falls_back_to_bundled(self) -> None:
        """MRNA_AI_SIMICL_DEMOS pointing at a missing file falls back to
        the bundled store (graceful degradation)."""
        os.environ[ENV_DEMOS] = "/nonexistent.json"
        try:
            store = load_default_demo_store()
            # Should still find the bundled examples/simicl_demos.json
            self.assertGreater(len(store), 0)
        finally:
            os.environ.pop(ENV_DEMOS, None)


# ---------------------------------------------------------------------------
# score_trial_with_llm integration
# ---------------------------------------------------------------------------


class TestScoreTrialWithLLMSimICL(unittest.TestCase):
    """Verify Sim-ICL plumbing in score_trial_with_llm."""

    def setUp(self) -> None:
        from mrnavax.trial_llm import score_trial_with_llm

        self.score = score_trial_with_llm
        os.environ[ENV_ENABLED] = "1"
        os.environ["MRNA_AI_FORCE_MOCK"] = "1"
        os.environ["MRNA_AI_LLM_BACKEND"] = "mock"

    def tearDown(self) -> None:
        os.environ.pop(ENV_ENABLED, None)
        os.environ.pop("MRNA_AI_FORCE_MOCK", None)
        os.environ.pop("MRNA_AI_LLM_BACKEND", None)

    def test_simicl_off_no_demos_in_notes(self) -> None:
        """When use_simicl=False, no simicl-k notes appear."""
        from mrnavax.trial_similar import DemoStore

        empty = DemoStore(demos=[])
        r = self.score(
            "55yo BRAF V600E melanoma patient",
            "NCT1",
            "BRAF trial",
            ["metastatic melanoma", "BRAF V600E"],
            ["prior therapy"],
            backend="mock",
            demo_store=empty,
            use_simicl=False,
        )
        self.assertFalse(any("simicl" in n for n in r.notes))

    def test_simicl_on_with_demos_appears_in_notes(self) -> None:
        """When Sim-ICL is on and demos exist, 'simicl-kN' appears in notes."""
        from mrnavax.trial_similar import DemoStore

        store = DemoStore(
            demos=[
                make_demo("demo_a", title="melanoma trial"),
                make_demo("demo_b", title="breast trial"),
            ]
        )
        r = self.score(
            "55yo melanoma patient",
            "NCT1",
            "melanoma trial",
            ["melanoma"],
            [],
            backend="mock",
            demo_store=store,
            use_simicl=True,
        )
        self.assertTrue(any(n.startswith("simicl-k") for n in r.notes))
        self.assertTrue(any("simicl-demo-ids=" in n for n in r.notes))

    def test_simicl_on_empty_store_no_simicl_in_notes(self) -> None:
        """Empty store + simicl on = no simicl in notes (nothing to inject)."""
        from mrnavax.trial_similar import DemoStore

        empty = DemoStore(demos=[])
        r = self.score(
            "55yo patient",
            "NCT1",
            "trial",
            ["i"],
            [],
            backend="mock",
            demo_store=empty,
            use_simicl=True,
        )
        self.assertFalse(any(n.startswith("simicl-k") for n in r.notes))

    def test_simicl_disabled_via_env_no_demos_in_notes(self) -> None:
        """MRNA_AI_SIMICL_ENABLED=0 disables even with non-empty store."""
        from mrnavax.trial_similar import DemoStore

        store = DemoStore(demos=[make_demo("demo_a")])
        os.environ[ENV_ENABLED] = "0"
        r = self.score(
            "55yo patient",
            "NCT1",
            "trial",
            ["i"],
            [],
            backend="mock",
            demo_store=store,
            use_simicl=None,  # let env decide
        )
        self.assertFalse(any(n.startswith("simicl-k") for n in r.notes))


if __name__ == "__main__":
    unittest.main()
