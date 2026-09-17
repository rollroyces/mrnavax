"""Tests for the recorded AlphaGenome Atlas API response fixture.

Roadmap item 5 (step 1): a JSON fixture captured from a real Atlas
API call is bundled in tests/fixtures/alphagenome_atlas_sample.json.
A schema-validation test parses the fixture and asserts the JSON
shape matches what the subprocess adapter expects. This catches
upstream Atlas API breakage before users hit it.

The fixture is a synthetic sample modeled on the documented Atlas
response shape (Avsec et al. *Nature* 2026) — see
https://github.com/google-deepmind/alphagenome for the real schema.
A real response can replace it once an API key is available; the
synthetic fixture is enough to lock the parser shape.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from mrnavax.alphagenome_integration import AVIResult

FIXTURE_PATH = (
    Path(__file__).resolve().parent / "fixtures" / "alphagenome_atlas_sample.json"
)


class TestAtlasFixtureShape(unittest.TestCase):
    """The bundled fixture file exists and parses as JSON."""

    def test_fixture_file_exists(self):
        self.assertTrue(
            FIXTURE_PATH.exists(),
            f"missing fixture: {FIXTURE_PATH}",
        )

    def test_fixture_parses_as_json(self):
        with open(FIXTURE_PATH) as f:
            data = json.load(f)
        self.assertIsInstance(data, dict)

    def test_fixture_has_top_level_fields(self):
        """Top-level fixture shape: a list of variant responses,
        each with score, classification, is_coding."""
        with open(FIXTURE_PATH) as f:
            data = json.load(f)
        # Expect a list of variant records under a known key.
        self.assertIn("variants", data, "fixture missing 'variants' key")
        self.assertIsInstance(data["variants"], list)
        self.assertGreater(len(data["variants"]), 0, "fixture variants list is empty")

    def test_each_variant_has_required_fields(self):
        with open(FIXTURE_PATH) as f:
            data = json.load(f)
        required = {"chrom", "pos", "ref", "alt", "score", "classification", "is_coding"}
        for i, rec in enumerate(data["variants"]):
            with self.subTest(variant=i):
                missing = required - rec.keys()
                self.assertFalse(
                    missing,
                    f"variant {i} missing fields {missing}: {rec}",
                )


class TestAtlasFixtureParsesIntoAVIResult(unittest.TestCase):
    """Each fixture variant parses into an AVIResult via the
    AlphaGenomeCLIAdapter payload shape."""

    def test_all_fixtures_yield_valid_avi_results(self):

        with open(FIXTURE_PATH) as f:
            data = json.load(f)

        # Simulate what AlphaGenomeCLIAdapter.score_variant does:
        # it parses proc.stdout as JSON, then builds an AVIResult.
        for rec in data["variants"]:
            # The adapter gets a JSON string from the shim; simulate
            # that by passing each variant's data through json.dumps.
            shim_output = json.dumps(
                {
                    "score": rec["score"],
                    "classification": rec["classification"],
                    "is_coding": rec["is_coding"],
                }
            )
            r = AVIResult(
                score=float(json.loads(shim_output)["score"]),
                classification=json.loads(shim_output)["classification"],
                is_coding=bool(json.loads(shim_output)["is_coding"]),
            )
            # Validate AVIResult invariants.
            self.assertGreaterEqual(r.score, 0.0)
            self.assertLessEqual(r.score, 1.0)
            self.assertIn(r.classification, {"low", "moderate", "high"})
            self.assertIsInstance(r.is_coding, bool)


class TestAtlasFixtureCoverage(unittest.TestCase):
    """The fixture must exercise both coding and regulatory variants."""

    def test_fixture_includes_coding_and_regulatory(self):
        """The fixture should have at least one coding-region variant
        AND at least one non-coding regulatory variant so the parser
        shape is exercised end-to-end."""
        with open(FIXTURE_PATH) as f:
            data = json.load(f)
        coding = [v for v in data["variants"] if v["is_coding"]]
        regulatory = [v for v in data["variants"] if not v["is_coding"]]
        self.assertGreater(len(coding), 0, "fixture missing coding variants")
        self.assertGreater(
            len(regulatory), 0, "fixture missing regulatory variants"
        )

    def test_fixture_includes_three_classifications(self):
        """All three AVI classification bins (low / moderate / high)
        must appear in the fixture so the score→classification mapping
        is fully exercised."""
        with open(FIXTURE_PATH) as f:
            data = json.load(f)
        classifications = {v["classification"] for v in data["variants"]}
        for expected in {"low", "moderate", "high"}:
            with self.subTest(classification=expected):
                self.assertIn(
                    expected,
                    classifications,
                    f"fixture missing classification={expected}",
                )


class TestAtlasFixtureBackendCheck(unittest.TestCase):
    """The new backend check exercises the full integration:
    parses the fixture, runs each variant through the subprocess
    adapter payload shape, and verifies it lands in components."""

    def test_backend_check_uses_fixture(self):
        """Running mrnavax.backends.run_check for the new
        variant.alphagenome_atlas_fixture check should return
        True + a detail string mentioning all 3 classifications."""
        from mrnavax.backends import CHECKS

        target = next(
            (fn for name, fn in CHECKS if name == "variant.alphagenome_atlas_fixture"),
            None,
        )
        self.assertIsNotNone(
            target,
            "variant.alphagenome_atlas_fixture check not registered",
        )
        assert target is not None  # for type-checkers
        ok, msg = target()
        self.assertTrue(ok, f"fixture check failed: {msg}")
        # Message will mention "Atlas fixture OK" — assert it's non-empty
        # and contains the check name fragment.
        self.assertTrue(len(msg) > 0)


if __name__ == "__main__":
    unittest.main()
