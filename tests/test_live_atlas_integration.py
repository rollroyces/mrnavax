"""Tests for the scheduled live Atlas integration harness.

Roadmap item 5b: a weekly GitHub Actions cron job hits the real
AlphaGenome Atlas API with a stable test variant, records the AVI
score, and compares against a baseline tolerance. Catches upstream
breakage in real time (the v0.18.0 fixture catches it in 7-day lag).

The harness has 3 modes:

* **mock** — runs against the stdlib mock. Always passes; no
  network. Used in unit tests and offline CI.
* **record** — calls the real Atlas, writes the response to
  tests/fixtures/alphagenome_atlas_live_<timestamp>.json. Triggered
  manually when refreshing the fixture.
* **regression** — calls the real Atlas, compares against a baseline
  score with ±5% tolerance. Triggered weekly.

Run modes::

    python -m mrnavax.live_atlas_integration --mode mock
    python -m mrnavax.live_atlas_integration --mode record \
        --output tests/fixtures/alphagenome_atlas_live.json
    python -m mrnavax.live_atlas_integration --mode regression \
        --baseline 0.72 --tolerance 0.05
"""

from __future__ import annotations

import json
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path

from mrnavax.live_atlas_integration import (
    LiveAtlasReport,
    compare_to_baseline,
    run_live_atlas_check,
)

# A canonical test variant: BRAF V600E coding-region. Stable position
# across releases (Avsec et al. cite it as a benchmark).
TEST_VARIANT = {
    "chrom": "chr7",
    "pos": 140753336,
    "ref": "T",
    "alt": "A",
}


def _make_record(score: float, classification: str = "high", is_coding: bool = True) -> dict:
    """Build a synthetic Atlas-like response record."""
    return {
        "chrom": TEST_VARIANT["chrom"],
        "pos": TEST_VARIANT["pos"],
        "ref": TEST_VARIANT["ref"],
        "alt": TEST_VARIANT["alt"],
        "score": score,
        "classification": classification,
        "is_coding": is_coding,
        "_captured_at": "2026-09-17T00:00:00Z",
    }


class TestLiveAtlasReport(unittest.TestCase):
    """LiveAtlasReport dataclass carries the captured + verified info."""

    def test_report_fields(self):
        r = LiveAtlasReport(
            chrom="chr7", pos=140753336, ref="T", alt="A",
            score=0.72, classification="high", is_coding=True,
            captured_at="2026-09-17",
            source="mock",
        )
        self.assertEqual(r.score, 0.72)
        self.assertEqual(r.source, "mock")


class TestCompareToBaseline(unittest.TestCase):
    """Baseline comparison with tolerance."""

    def test_within_tolerance_passes(self):
        result = compare_to_baseline(score=0.72, baseline=0.70, tolerance=0.05)
        self.assertTrue(result["within_tolerance"])
        self.assertAlmostEqual(result["delta"], 0.02, places=6)

    def test_outside_tolerance_fails(self):
        result = compare_to_baseline(score=0.50, baseline=0.70, tolerance=0.05)
        self.assertFalse(result["within_tolerance"])
        self.assertAlmostEqual(result["delta"], -0.20, places=6)

    def test_exact_match_passes(self):
        result = compare_to_baseline(score=0.70, baseline=0.70, tolerance=0.05)
        self.assertTrue(result["within_tolerance"])
        self.assertEqual(result["delta"], 0.0)

    def test_negative_delta_within_tolerance(self):
        """Score dropped 3%, well within 5% tolerance — passes."""
        result = compare_to_baseline(score=0.679, baseline=0.70, tolerance=0.05)
        self.assertTrue(result["within_tolerance"])


class TestRunLiveAtlasCheckMock(unittest.TestCase):
    """The harness runs in mock mode without network access."""

    def test_mock_mode_returns_valid_report(self):
        report = run_live_atlas_check(
            mode="mock",
            variant=TEST_VARIANT,
            api_key="dummy",  # ignored in mock mode
        )
        self.assertIsInstance(report, LiveAtlasReport)
        self.assertEqual(report.source, "mock")
        self.assertEqual(report.chrom, "chr7")
        self.assertEqual(report.pos, 140753336)
        self.assertGreaterEqual(report.score, 0.0)
        self.assertLessEqual(report.score, 1.0)
        self.assertIn(report.classification, {"low", "moderate", "high"})


class TestRunLiveAtlasCheckRegression(unittest.TestCase):
    """Regression mode compares against a baseline."""

    def test_regression_within_tolerance_returns_pass(self):
        baseline = 0.70
        report = run_live_atlas_check(
            mode="regression",
            variant=TEST_VARIANT,
            api_key="dummy",
            baseline=baseline,
            tolerance=0.05,
            # Use a stub mock_score_fn so we don't actually call the API.
            # Returns a dict (mirrors the real score_fn shape).
            score_fn=lambda: {
                "score": 0.72,  # within tolerance
                "classification": "high",
                "is_coding": True,
            },
        )
        self.assertTrue(report.within_tolerance)

    def test_regression_outside_tolerance_returns_fail(self):
        report = run_live_atlas_check(
            mode="regression",
            variant=TEST_VARIANT,
            api_key="dummy",
            baseline=0.70,
            tolerance=0.05,
            score_fn=lambda: {
                "score": 0.40,  # way outside tolerance
                "classification": "moderate",
                "is_coding": True,
            },
        )
        self.assertFalse(report.within_tolerance)


class TestRunLiveAtlasCheckRecord(unittest.TestCase):
    """Record mode writes the response to a JSON file."""

    def test_record_mode_writes_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "captured.json"

            def record_fn() -> dict:
                return _make_record(0.72)

            run_live_atlas_check(
                mode="record",
                variant=TEST_VARIANT,
                api_key="dummy",
                output_path=out_path,
                score_fn=record_fn,
            )
            self.assertTrue(out_path.exists())
            data = json.loads(out_path.read_text())
            self.assertEqual(data["score"], 0.72)
            self.assertEqual(data["chrom"], "chr7")


class TestCLI(unittest.TestCase):
    """CLI wrapper exposes the harness as a module."""

    def test_module_imports(self):
        from mrnavax import live_atlas_integration

        self.assertTrue(hasattr(live_atlas_integration, "run_live_atlas_check"))
        self.assertTrue(hasattr(live_atlas_integration, "compare_to_baseline"))


# Convenience wrapper around tempfile.TemporaryDirectory that yields a Path
@contextmanager
def _tempfile_TemporaryDirectory():
    with tempfile.TemporaryDirectory() as tmp:
        yield Path(tmp)


if __name__ == "__main__":
    unittest.main()
