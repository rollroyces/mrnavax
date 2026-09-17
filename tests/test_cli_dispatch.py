"""Tests for the unified CLI dispatcher.

Verifies that every documented tool subcommand is registered and
executable, and that error paths produce non-zero exit codes.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class TestCLISubcommandRegistration(unittest.TestCase):
    """The CLI dispatcher must register every tool we ship."""

    def test_all_seven_subcommands_listed_in_help(self) -> None:
        """`mrnavax --help` must list all 7 tools."""
        result = subprocess.run(
            [sys.executable, "-m", "mrnavax.cli", "--help"],
            capture_output=True,
            text=True,
            env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        for tool in (
            "codon", "neoantigen", "trial", "lnp",
            "scrna", "manufacture", "spatial",
        ):
            self.assertIn(tool, result.stdout, msg=f"missing tool: {tool}")

    def test_help_text_includes_spatial_subcommand(self) -> None:
        """Regression: the spatial subcommand is wired up.

        Bug history: prior to v0.13.0, README and docs mentioned
        `mrnavax spatial` but the CLI dispatcher did NOT register it.
        Verified by checking that `mrnavax --help` lists `spatial`
        AND that the spatial sub-parser's args are reachable (the
        dispatcher routes to `_spatial_run` which has the right args).
        """
        # 1. The --help output mentions spatial
        result = subprocess.run(
            [sys.executable, "-m", "mrnavax.cli", "--help"],
            capture_output=True,
            text=True,
            env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("spatial", result.stdout)

        # 2. The dispatcher routes `spatial` to _spatial_run. If the
        # spatial CLI is invoked with valid args, _spatial_run is
        # called and a real result is produced. (See TestSpatialCLIEndToEnd
        # below for the actual run.)
        # 3. If the spatial CLI is invoked with INVALID args
        # (missing --count-file), the dispatcher should call
        # _spatial_run which then errors — NOT silently route to a
        # different subcommand.
        result = subprocess.run(
            [
                sys.executable, "-m", "mrnavax.cli", "spatial",
                # Missing --count-file and --locations-file
            ],
            capture_output=True,
            text=True,
            env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
        )
        # Argparse in _spatial_run exits with code 2 (argparse error)
        self.assertNotEqual(result.returncode, 0)
        # Should be an argparse-style "required" error mentioning
        # --count-file
        self.assertIn("--count-file", result.stderr)


class TestSpatialCLIEndToEnd(unittest.TestCase):
    """End-to-end test of `mrnavax spatial` CLI invocation."""

    def _make_inputs(self, tmp: Path) -> tuple[Path, Path]:
        count = tmp / "counts.tsv"
        loc = tmp / "locs.tsv"
        count.write_text(
            "\t".join([""] + [f"G{j}" for j in range(5)]) + "\n"
            + "\n".join(
                f"loc_{i+1}\t1\t0\t1\t1\t0" for i in range(6)
            )
            + "\n"
        )
        loc.write_text(
            "\t".join(["", "x", "y"]) + "\n"
            + "\n".join(
                f"loc_{i+1}\t{17.0 + i*0.5}\t{5.0 + (i % 4)}" for i in range(6)
            )
            + "\n"
        )
        return count, loc

    def test_spatial_cli_runs_with_mock_backend(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            count, loc = self._make_inputs(tmp_path)
            result = subprocess.run(
                [
                    sys.executable, "-m", "mrnavax.cli", "spatial",
                    "--count-file", str(count),
                    "--locations-file", str(loc),
                    "--platform", "ST",
                    "--num-modules", "5",
                    "--backend", "mock",
                ],
                capture_output=True,
                text=True,
                env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            import json
            payload = json.loads(result.stdout)
            self.assertEqual(payload["platform"], "ST")
            self.assertEqual(len(payload["modules"]), 5)
            self.assertEqual(payload["backend"], "mock")
            self.assertEqual(payload["n_spots"], 6)

    def test_spatial_cli_missing_count_file_returns_error(self) -> None:
        """Missing --count-file must exit non-zero with a clear error."""
        result = subprocess.run(
            [
                sys.executable, "-m", "mrnavax.cli", "spatial",
                "--count-file", "/nonexistent.tsv",
                "--locations-file", "/nonexistent.tsv",
                "--platform", "ST",
            ],
            capture_output=True,
            text=True,
            env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
        )
        self.assertNotEqual(result.returncode, 0)
        # Error message should mention the missing file
        self.assertIn("/nonexistent", result.stderr)

    def test_spatial_cli_writes_to_out_file(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            count, loc = self._make_inputs(tmp_path)
            out_file = tmp_path / "out.json"
            result = subprocess.run(
                [
                    sys.executable, "-m", "mrnavax.cli", "spatial",
                    "--count-file", str(count),
                    "--locations-file", str(loc),
                    "--platform", "ST",
                    "--num-modules", "3",
                    "--backend", "mock",
                    "--out", str(out_file),
                ],
                capture_output=True,
                text=True,
                env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            self.assertTrue(out_file.exists())
            import json
            payload = json.loads(out_file.read_text())
            self.assertEqual(payload["platform"], "ST")
            self.assertEqual(len(payload["modules"]), 3)


class TestManufactureCLI(unittest.TestCase):
    """The manufacture CLI should still work after the spatial addition."""

    def test_manufacture_cli_runs(self) -> None:
        """A short CDS path produces valid JSON overall_score."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as f:
            f.write("ATGGCATCGATCGATCGCATCATCGCGCGCGCGCCGCATAA")
            tmp_path = f.name
        try:
            result = subprocess.run(
                [
                    sys.executable, "-m", "mrnavax.cli", "manufacture",
                    "--cds", tmp_path,
                ],
                capture_output=True,
                text=True,
                env={**os.environ, "MRNA_AI_FORCE_MOCK": "1"},
            )
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            import json
            payload = json.loads(result.stdout)
            self.assertIn("overall_score", payload)
        finally:
            os.unlink(tmp_path)

    def test_score_manufacturability_handles_none_utrs(self) -> None:
        """Regression: score_manufacturability(utr5=None, utr3=None)
        must not crash with AttributeError.

        Bug history: prior to v0.13.1, score_manufacturability called
        `utr5.upper()` on None, raising AttributeError. The CLI
        `_manufacture_run` passed None for empty --utr5/--utr3 flags,
        so every `mrnavax manufacture --cds ...` invocation crashed.
        """
        from mrnavax.manufacturability import score_manufacturability

        # Must not raise
        report = score_manufacturability(
            "ATGGCATCGATCGATCGCATCATCGCGCGCGCGCCGCATAA",
            utr5=None,  # was: crashed
            utr3=None,  # was: crashed
        )
        self.assertGreaterEqual(report.overall_score, 0.0)


if __name__ == "__main__":
    unittest.main()
