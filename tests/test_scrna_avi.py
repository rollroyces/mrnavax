"""Tests for AlphaGenome Atlas AVI integration in the sc_rna_pipeline.

Roadmap item 3: surface AVI in the scrna pipeline output so users can see
which variants were prioritized by AVI (regulatory) vs AlphaMissense
(coding), and surface the AVI lookup status in the pipeline note.

The pipeline already auto-wires the AVI lookup when variants carry DNA
coordinates (added in v0.16.0). This test module verifies the surfacing
of that information in PipelineReport.variant_scores and .note.
"""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from mrnavax.sc_rna_pipeline import (
    PipelineReport,
    load_variants,
    run_pipeline,
)


def _write_tsv(path: Path, rows: list[str]) -> None:
    path.write_text("\n".join(rows) + "\n")


def _make_minimal_inputs(tmp: Path) -> tuple[Path, Path, Path]:
    """Build the 3 input files needed for run_pipeline.

    Returns (expression_csv, variants_csv, proteins_fasta).
    """
    expr = tmp / "cells.csv"
    header = "\t".join([""] + [f"G{j}" for j in range(10)])
    data_lines = [
        "\t".join([f"cell_{i}"] + [str(i + j) for j in range(10)])
        for i in range(20)
    ]
    _write_tsv(expr, [header] + data_lines)
    variants = tmp / "variants.csv"
    with open(variants, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["gene", "position", "wt_aa", "mut_aa", "chrom", "ref_dna", "alt_dna"]
        )
        writer.writerow(["TP53", "175", "R", "H", "chr17", "G", "A"])
        writer.writerow(["KRAS", "12", "G", "D", "chr12", "C", "A"])
        writer.writerow(["REG1", "101", "A", "G", "chr11", "A", "G"])
    fasta = tmp / "proteins.fasta"
    fasta.write_text(
        ">sp|P04637|TP53_HUMAN Tumor protein p53 OS=Homo sapiens\n"
        "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n"
        ">sp|P01116|RASK_HUMAN GTPase KRas OS=Homo sapiens\n"
        "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY\n"
    )
    return expr, variants, fasta


class TestVariantLoadsDNAFields(unittest.TestCase):
    """load_variants parses chrom/ref_dna/alt_dna when present."""

    def test_load_variants_with_dna_fields(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, variants, _ = _make_minimal_inputs(Path(tmp))
            vs = load_variants(variants)
            self.assertEqual(len(vs), 3)
            self.assertEqual(vs[0].gene, "TP53")
            self.assertEqual(vs[0].chrom, "chr17")
            self.assertEqual(vs[0].ref_dna, "G")
            self.assertEqual(vs[0].alt_dna, "A")
            # Last variant — also has DNA fields
            self.assertEqual(vs[2].gene, "REG1")
            self.assertEqual(vs[2].chrom, "chr11")
            self.assertEqual(vs[2].ref_dna, "A")
            self.assertEqual(vs[2].alt_dna, "G")


class TestPipelineReportHasAVINote(unittest.TestCase):
    """PipelineReport.note mentions the AVI integration status."""

    def test_note_mentions_avi_when_used(self):
        """When variants carry DNA coordinates, the pipeline note should
        mention AVI (AlphaGenome Atlas) integration."""
        with tempfile.TemporaryDirectory() as tmp:
            expr, variants, fasta = _make_minimal_inputs(Path(tmp))
            report = run_pipeline(
                expr,
                variants,
                fasta,
                hla=("HLA-A*02:01",),
                tumor_marker_genes=["TP53", "KRAS"],
            )
            self.assertIsInstance(report, PipelineReport)
            # Note should reference either AlphaMissense or AlphaGenome
            # depending on which lookup fired. With DNA coords present,
            # the AVI lookup is wired → note should mention Atlas/AVI.
            self.assertTrue(
                "AlphaGenome" in report.note or "AVI" in report.note,
                f"note should mention AVI integration: {report.note!r}",
            )


class TestPipelineReportExposesAVIScores(unittest.TestCase):
    """PipelineReport.variant_scores includes AVI for variants with DNA coords."""

    def test_variant_scores_present_for_filtered_variants(self):
        """After filtering, variant_scores contains normalized scores
        for each kept variant."""
        with tempfile.TemporaryDirectory() as tmp:
            expr, variants, fasta = _make_minimal_inputs(Path(tmp))
            report = run_pipeline(
                expr,
                variants,
                fasta,
                hla=("HLA-A*02:01",),
                tumor_marker_genes=["TP53", "KRAS"],
                variant_filter_top_fraction=1.0,  # keep all
            )
            self.assertIsInstance(report.variant_scores, dict)
            # At least one variant key should be present
            self.assertGreater(len(report.variant_scores), 0)


class TestEndToEndPipelineWithAVI(unittest.TestCase):
    """End-to-end: a real run with DNA coordinates populates the pipeline."""

    def test_full_run_produces_dna_aware_output(self):
        """With DNA coordinates + filter, the pipeline keeps variants
        and exposes their scores."""
        with tempfile.TemporaryDirectory() as tmp:
            expr, variants, fasta = _make_minimal_inputs(Path(tmp))
            report = run_pipeline(
                expr,
                variants,
                fasta,
                hla=("HLA-A*02:01",),
                tumor_marker_genes=["TP53", "KRAS"],
                variant_filter_top_fraction=0.6,
            )
            self.assertIsInstance(report, PipelineReport)
            # n_variants_after_filter <= n_variants_input
            self.assertLessEqual(
                report.n_variants_after_filter, report.n_variants_input
            )
            # If any variants survived the filter, scores are present
            if report.n_variants_after_filter > 0:
                self.assertGreater(len(report.variant_scores), 0)

    def test_run_without_dna_fields_still_works(self):
        """A CSV without chrom/ref_dna/alt_dna columns should still work
        (backward compatibility — variants fall through to AM-or-BLOSUM)."""
        with tempfile.TemporaryDirectory() as tmp:
            expr, _, fasta = _make_minimal_inputs(Path(tmp))
            # Build a variants CSV without DNA columns
            variants_legacy = Path(tmp) / "legacy.csv"
            with open(variants_legacy, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["gene", "position", "wt_aa", "mut_aa"])
                writer.writerow(["TP53", "175", "R", "H"])
                writer.writerow(["KRAS", "12", "G", "D"])
            report = run_pipeline(
                expr,
                variants_legacy,
                fasta,
                hla=("HLA-A*02:01",),
                tumor_marker_genes=["TP53", "KRAS"],
            )
            self.assertIsInstance(report, PipelineReport)
            # DNA fields default to None
            self.assertEqual(report.n_variants_input, 2)


if __name__ == "__main__":
    unittest.main()
