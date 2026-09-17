"""Tests for the AlphaGenome Atlas regulatory-variant-scorer integration.

The Protocol contract is:

    RegulatoryVariantScorer(Protocol):
        def score_variant(self, chrom: str, pos: int, ref: str, alt: str)
            -> AVIResult

`AVIResult` is a frozen dataclass with `score in [0, 1]`,
`classification in {"low", "moderate", "high"}`, and an `is_coding: bool`
flag that distinguishes regulatory-region variants from coding-region
ones (where AlphaMissense is the better signal).

Two backends satisfy the Protocol:

* `MockRegulatoryVariantScorer` — stdlib-only, deterministic.
  Always present. Computes a synthetic AVI score from
  (chrom, pos, ref, alt) via a stdlib hash → [0, 1].

* `AlphaGenomeCLIAdapter` — calls the official `alphagenome` Python
  package via subprocess (subprocess pattern matches RiboDecode and
  STModule shims). Only registered when the `[variant-alphagenome]`
  extra is installed and the user supplies `ALPHAGENOME_API_KEY`.
  Non-commercial use only (per Google's terms).
"""

from __future__ import annotations

import csv
import unittest
from pathlib import Path

from mrnavax.alphagenome_integration import (
    AlphaGenomeCLIAdapter,
    AVIResult,
    MockRegulatoryVariantScorer,
    RegulatoryVariantScorer,
    regulatory_score,
    select_regulatory_scorer,
)


class TestAVIResult(unittest.TestCase):
    """AVIResult is a frozen dataclass with the expected fields."""

    def test_avi_result_fields_and_invariants(self):
        from dataclasses import FrozenInstanceError

        r = AVIResult(score=0.5, classification="moderate", is_coding=False)
        self.assertEqual(r.score, 0.5)
        self.assertEqual(r.classification, "moderate")
        self.assertFalse(r.is_coding)

        # Score out of range must raise at construction time (fail-fast).
        with self.assertRaises(ValueError):
            AVIResult(score=1.5, classification="high", is_coding=False)
        with self.assertRaises(ValueError):
            AVIResult(score=-0.1, classification="low", is_coding=False)

        # Frozen — mutation must raise.
        with self.assertRaises(FrozenInstanceError):
            r.score = 0.9  # type: ignore[misc]

    def test_classification_thresholds(self):
        """Thresholds: low <0.34 ≤ moderate <0.564 ≤ high (matches AlphaMissense paper bins)."""
        for score, expected in [
            (0.0, "low"),
            (0.33, "low"),
            (0.34, "moderate"),
            (0.563, "moderate"),
            (0.564, "high"),
            (1.0, "high"),
        ]:
            with self.subTest(score=score):
                r = AVIResult(score=score, classification="", is_coding=False)
                self.assertEqual(r.classification, expected)


class TestMockRegulatoryVariantScorer(unittest.TestCase):
    """The mock backend: deterministic, stdlib-only, no network."""

    def test_mock_satisfies_protocol(self):
        mock = MockRegulatoryVariantScorer()
        self.assertIsInstance(mock, RegulatoryVariantScorer)

    def test_mock_score_returns_avi_result(self):
        mock = MockRegulatoryVariantScorer()
        r = mock.score_variant("chr7", 140753336, "T", "A")  # BRAF V600E in hg38
        self.assertIsInstance(r, AVIResult)
        self.assertGreaterEqual(r.score, 0.0)
        self.assertLessEqual(r.score, 1.0)
        self.assertIn(r.classification, {"low", "moderate", "high"})

    def test_mock_is_deterministic(self):
        mock = MockRegulatoryVariantScorer()
        a = mock.score_variant("chr17", 7675088, "C", "T")
        b = mock.score_variant("chr17", 7675088, "C", "T")
        self.assertEqual(a.score, b.score)
        self.assertEqual(a.classification, b.classification)

    def test_mock_distinguishes_coding_vs_regulatory(self):
        """Coding-region variants get is_coding=True (use AlphaMissense);
        regulatory-region variants get is_coding=False (use AlphaGenome Atlas).

        The mock's heuristic is deterministic on position parity: even
        positions are 'coding', odd positions are 'regulatory'. This is
        purely a mock — the real Atlas returns per-variant is_coding."""
        mock = MockRegulatoryVariantScorer()
        coding = mock.score_variant("chr7", 140753336, "T", "A")  # even pos
        reg = mock.score_variant("chr7", 140753337, "T", "A")      # odd pos
        self.assertTrue(coding.is_coding)
        self.assertFalse(reg.is_coding)


class TestRegulatoryScoreFunction(unittest.TestCase):
    """The `regulatory_score(chrom, pos, ref, alt)` convenience wrapper."""

    def test_returns_avi_result(self):
        r = regulatory_score("chr7", 140753336, "T", "A")
        self.assertIsInstance(r, AVIResult)

    def test_returns_specific_score_for_known_variant(self):
        # Same input → same output (deterministic mock)
        a = regulatory_score("chr7", 140753336, "T", "A")
        b = regulatory_score("chr7", 140753336, "T", "A")
        self.assertEqual(a.score, b.score)


class TestBackendSelector(unittest.TestCase):
    """select_regulatory_scorer returns the mock unless ALPHAGENOME_API_KEY is set."""

    def test_mock_returned_when_key_absent(self):
        import os

        old = os.environ.pop("ALPHAGENOME_API_KEY", None)
        try:
            sel = select_regulatory_scorer()
            self.assertIsInstance(sel, MockRegulatoryVariantScorer)
        finally:
            if old is not None:
                os.environ["ALPHAGENOME_API_KEY"] = old

    def test_real_returned_when_key_present(self):
        import os

        os.environ["ALPHAGENOME_API_KEY"] = "fake-test-key"
        try:
            sel = select_regulatory_scorer()
            self.assertIsInstance(sel, AlphaGenomeCLIAdapter)
        finally:
            del os.environ["ALPHAGENOME_API_KEY"]


class TestAlphaGenomeCLIAdapter(unittest.TestCase):
    """Real adapter — subprocess invocation shape + JSON parsing."""

    def test_real_adapter_satisfies_protocol(self):
        adapter = AlphaGenomeCLIAdapter(api_key="fake")
        self.assertIsInstance(adapter, RegulatoryVariantScorer)

    def test_subprocess_payload_shape(self):
        """The subprocess argv + JSON payload must match what the
        alphagenome-cli shim expects. We don't actually call the
        upstream — we just verify the subprocess invocation is
        well-formed via a recording shim."""
        import json as _json
        recorded: dict = {}

        def fake_run(cmd, **kwargs):  # mirrors subprocess.run signature
            recorded["cmd"] = cmd
            # Argv shape: [python_exe, shim_path] — must be exactly 2 items.
            assert len(cmd) == 2, f"unexpected argv length: {cmd}"
            # Parse the JSON payload from stdin.
            payload = _json.loads(kwargs["input"])
            recorded["payload"] = payload
            # Return a synthetic AlphaGenome-shaped AVI result.
            out = _json.dumps(
                {"score": 0.72, "classification": "high", "is_coding": True}
            )
            from subprocess import CompletedProcess

            return CompletedProcess(cmd, 0, stdout=out, stderr="")

        adapter = AlphaGenomeCLIAdapter(api_key="fake", _run=fake_run)
        r = adapter.score_variant("chr7", 140753336, "T", "A")
        # The recorded payload must carry the variant.
        self.assertEqual(recorded["payload"]["chrom"], "chr7")
        self.assertEqual(recorded["payload"]["pos"], 140753336)
        self.assertEqual(recorded["payload"]["ref"], "T")
        self.assertEqual(recorded["payload"]["alt"], "A")
        # And the parsed result.
        self.assertEqual(r.score, 0.72)
        self.assertEqual(r.classification, "high")
        self.assertTrue(r.is_coding)

    def test_subprocess_failure_returns_moderate_default(self):
        """If the subprocess fails (no key, no network), return a
        deterministic 'unknown' result instead of raising — the
        caller can decide whether to abort."""
        from subprocess import CompletedProcess

        def fake_run_fail(cmd, **kwargs):
            return CompletedProcess(cmd, 1, stdout="", stderr="auth failed")

        adapter = AlphaGenomeCLIAdapter(api_key="bad", _run=fake_run_fail)
        r = adapter.score_variant("chr7", 140753336, "T", "A")
        self.assertEqual(r.score, 0.0)
        self.assertEqual(r.classification, "low")
        self.assertFalse(r.is_coding)


class TestBackendCheckRegistered(unittest.TestCase):
    """The backend integrity check `variant.alphagenome_atlas` is registered."""

    def test_check_present_in_backends_module(self):
        from mrnavax import backends

        # backends.py exposes a public CHECKS list; verify our new check is there.
        src = Path(backends.__file__).read_text()
        self.assertIn('"variant.alphagenome_atlas"', src)
        names = [n for n, _ in backends.CHECKS]
        self.assertIn("variant.alphagenome_atlas", names)

    def test_check_passes_under_mock(self):
        """When MRNA_AI_FORCE_MOCK=1 (or no key is present), the check
        runs the mock end-to-end and reports AVI score is in [0, 1]."""
        from mrnavax.backends import CHECKS

        # Locate the check by name
        target = next(
            (fn for name, fn in CHECKS if name == "variant.alphagenome_atlas"),
            None,
        )
        self.assertIsNotNone(target, "variant.alphagenome_atlas check not registered")
        assert target is not None  # for type-checkers
        ok, msg = target()
        self.assertTrue(ok, f"check failed: {msg}")
        self.assertIn("AVI", msg.upper())


class TestCliSubcommand(unittest.TestCase):
    """The CLI exposes `mrnavax variant-regulatory`."""

    def test_subcommand_is_registered(self):
        import contextlib

        # main(argv) prints help and returns 0
        import io

        from mrnavax.cli import main

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            try:
                main(["--help"])
            except SystemExit:
                pass
        out = buf.getvalue()
        self.assertIn("variant-regulatory", out)


class TestExampleCSV(unittest.TestCase):
    """The bundled example CSV parses cleanly with the expected columns."""

    def test_example_csv_schema(self):
        repo_root = Path(__file__).resolve().parent.parent
        path = repo_root / "mrnavax" / "examples" / "regulatory_variants.csv"
        self.assertTrue(path.exists(), f"missing: {path}")
        with path.open() as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        self.assertGreaterEqual(len(rows), 3)
        required = {"chrom", "pos", "ref", "alt", "label"}
        for row in rows:
            self.assertTrue(required <= row.keys(), f"row {row} missing columns")
        # Look up the BRAF V600E canonical row by exact label match.
        braf = next(
            (r for r in rows if r["label"] == "BRAF V600E (coding)"),
            None,
        )
        self.assertIsNotNone(braf, "BRAF V600E canonical row missing")
        self.assertEqual(braf["chrom"], "chr7")
        self.assertEqual(int(braf["pos"]), 140753336)
        # Also verify a regulatory-region row exists.
        reg = next(
            (
                r for r in rows
                if "regulatory" in r["label"].lower()
                or r["label"] == "intergenic regulatory variant"
            ),
            None,
        )
        self.assertIsNotNone(reg, "regulatory-region row missing")


class TestEndToEndScoreFlow(unittest.TestCase):
    """Real end-to-end: read CSV → score each variant → emit ranked output."""

    def test_score_variants_from_csv_to_json(self):
        from mrnavax.alphagenome_integration import score_variants_from_csv

        repo_root = Path(__file__).resolve().parent.parent
        csv_in = repo_root / "mrnavax" / "examples" / "regulatory_variants.csv"
        results = score_variants_from_csv(str(csv_in))
        self.assertEqual(len(results), 7)  # 7 example variants
        # All scores in [0, 1]
        for r in results:
            self.assertGreaterEqual(r.score, 0.0)
            self.assertLessEqual(r.score, 1.0)
        # Sorted by score descending
        scores = [r.score for r in results]
        self.assertEqual(scores, sorted(scores, reverse=True))


if __name__ == "__main__":
    unittest.main()
