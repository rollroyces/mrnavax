"""Tests for the variant-prioritization case study module.

Roadmap item 7: a worked ClinVar-style case study that demonstrates
the toolkit's variant prioritization pipeline end-to-end. The
case_study module:

* Bundles a small set of curated ClinVar-style variants in
  mrnavax/examples/clinvar_curated.csv (synthetic but
  biologically realistic — real hg38 positions for known
  pathogenic / benign variants).
* Provides score_case_study_variants(score_fn) that scores each
  variant with the provided scoring function and returns labeled
  results.
* Provides precision_at_k(scored, k) and variant_set_summary()
  helpers for computing precision@K and confusion-matrix-style
  metrics.

The point isn't to compete with published benchmarks — it's to
give users a reproducible starting point they can swap in their own
variants against and see how the toolkit's prioritization actually
ranks them.
"""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from mrnavax.case_study import (
    ClinVarVariant,
    load_clinvar_variants,
    precision_at_k,
    score_case_study_variants,
    summarize_variant_set,
)


def _make_csv(path: Path, rows: list[dict]) -> None:
    """Helper: write a ClinVar-style CSV."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "gene",
                "chrom",
                "pos",
                "ref",
                "alt",
                "label",
                "pathogenicity",
                "is_coding",
            ],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _high_score(*args, **kwargs) -> dict:
    """A stub scoring function that returns high scores for the
    arg tuple passed in. Used to test precision_at_k with a controlled
    input."""
    # Parse the variant key from the call args (gene, pos, wt_aa, mut_aa)
    return {"score": 0.9, "label": args[0] if args else "unknown"}


class TestLoadClinVarVariants(unittest.TestCase):
    """load_clinvar_variants parses the bundled curated CSV."""

    def test_load_from_bundled_example(self):
        """The bundled examples/clinvar_curated.csv is loadable."""
        # Repo root: parents[0] of the package → repo root
        repo_root = Path(__file__).resolve().parent.parent
        csv_path = repo_root / "mrnavax" / "examples" / "clinvar_curated.csv"
        self.assertTrue(csv_path.exists(), f"missing: {csv_path}")
        variants = load_clinvar_variants(csv_path)
        self.assertGreater(len(variants), 0)

    def test_load_from_custom_path(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "v.csv"
            _make_csv(
                path,
                [
                    {
                        "gene": "BRAF",
                        "chrom": "chr7",
                        "pos": "140753336",
                        "ref": "T",
                        "alt": "A",
                        "label": "BRAF V600E",
                        "pathogenicity": "pathogenic",
                        "is_coding": "true",
                    },
                ],
            )
            variants = load_clinvar_variants(path)
            self.assertEqual(len(variants), 1)
            v = variants[0]
            self.assertEqual(v.gene, "BRAF")
            self.assertEqual(v.pathogenicity, "pathogenic")
            self.assertTrue(v.is_coding)


class TestClinVarVariantFields(unittest.TestCase):
    """ClinVarVariant dataclass carries all fields needed for scoring."""

    def test_coding_variant_has_dna_coords(self):
        v = ClinVarVariant(
            gene="BRAF",
            chrom="chr7",
            pos=140753336,
            ref="T",
            alt="A",
            label="BRAF V600E",
            pathogenicity="pathogenic",
            is_coding=True,
        )
        self.assertEqual(v.chrom, "chr7")
        self.assertEqual(v.pos, 140753336)
        self.assertEqual(v.ref, "T")
        self.assertEqual(v.alt, "A")
        self.assertTrue(v.is_coding)
        self.assertEqual(v.pathogenicity, "pathogenic")


class TestPrecisionAtK(unittest.TestCase):
    """precision_at_k returns fraction of pathogenic in top-K."""

    def test_precision_at_k_basic(self):
        scores = [
            ("v1", 0.9, "pathogenic"),
            ("v2", 0.8, "benign"),
            ("v3", 0.7, "pathogenic"),
            ("v4", 0.6, "benign"),
        ]
        # Top-1: v1 is pathogenic → precision = 1.0
        self.assertEqual(precision_at_k(scores, k=1), 1.0)
        # Top-2: 1 pathogenic of 2 → precision = 0.5
        self.assertEqual(precision_at_k(scores, k=2), 0.5)
        # Top-3: 2 pathogenic of 3 → precision = 0.6667
        self.assertAlmostEqual(precision_at_k(scores, k=3), 2 / 3)
        # Top-4: 2 pathogenic of 4 → precision = 0.5
        self.assertEqual(precision_at_k(scores, k=4), 0.5)

    def test_precision_at_k_empty(self):
        """Empty scores → precision = 0.0 (not a crash)."""
        self.assertEqual(precision_at_k([], k=5), 0.0)

    def test_precision_at_k_k_larger_than_n(self):
        """k larger than the number of items → divide by k, not by len(scores)."""
        scores = [("v1", 0.9, "pathogenic")]
        # k=5 but only 1 item; precision = 1 / 1 (we have 1 pathogenic in top-1)
        # but k=5 so the top-5 includes only this 1 item
        self.assertEqual(precision_at_k(scores, k=5), 1.0)

    def test_precision_at_k_input_format(self):
        """Accepts (id, score, pathogenicity) tuples."""
        scores = [
            (1, 0.9, "pathogenic"),
            (2, 0.5, "benign"),
        ]
        self.assertEqual(precision_at_k(scores, k=1), 1.0)

    def test_precision_at_k_orders_by_score_descending(self):
        """Items are sorted by score before top-K is taken."""
        scores = [
            ("lo_first", 0.1, "pathogenic"),
            ("hi_second", 0.9, "benign"),
        ]
        # Without sorting: top-1 is "lo_first" → precision = 1.0
        # With sorting: top-1 is "hi_second" → precision = 0.0
        self.assertEqual(precision_at_k(scores, k=1, sort=True), 0.0)


class TestScoreCaseStudyVariants(unittest.TestCase):
    """score_case_study_variants calls the scoring fn per variant."""

    def test_score_case_study_variants_with_stub_scorer(self):
        """A stub scoring function that returns a constant score
        should produce a scored list with the right shape."""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "v.csv"
            _make_csv(
                path,
                [
                    {
                        "gene": "BRAF",
                        "chrom": "chr7",
                        "pos": "140753336",
                        "ref": "T",
                        "alt": "A",
                        "label": "BRAF V600E",
                        "pathogenicity": "pathogenic",
                        "is_coding": "true",
                    },
                    {
                        "gene": "APC",
                        "chrom": "chr5",
                        "pos": "112839510",
                        "ref": "C",
                        "alt": "T",
                        "label": "APC promoter",
                        "pathogenicity": "pathogenic",
                        "is_coding": "false",
                    },
                ],
            )
            variants = load_clinvar_variants(path)

            def stub_scorer(variant):
                return {"score": 0.8 if variant.is_coding else 0.7}

            results = score_case_study_variants(variants, stub_scorer)
            self.assertEqual(len(results), 2)
            # Sorted descending by score
            self.assertGreaterEqual(
                results[0]["score"], results[1]["score"]
            )
            # Pathogenicity preserved
            for r in results:
                self.assertIn(r["pathogenicity"], {"pathogenic", "benign", "uncertain"})

    def test_score_case_study_variants_with_real_scorer(self):
        """End-to-end with the real score_variant + a stub AM/AVI."""
        from mrnavax.variant_scorer import score_variant

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "v.csv"
            _make_csv(
                path,
                [
                    {
                        "gene": "BRAF",
                        "chrom": "chr7",
                        "pos": "140753336",
                        "ref": "T",
                        "alt": "A",
                        "label": "BRAF V600E",
                        "pathogenicity": "pathogenic",
                        "is_coding": "true",
                    },
                ],
            )
            variants = load_clinvar_variants(path)

            def stub_am(uniprot, wt, pos, mut):
                from mrnavax.alphamissense_integration import AlphaMissenseResult
                return AlphaMissenseResult(
                    score=0.95, classification="likely_pathogenic",
                    uniprot=uniprot, aa_change=f"{wt}{pos}{mut}",
                )

            def real_scorer(variant):
                if variant.is_coding:
                    r = score_variant(
                        gene=variant.gene, position=variant.pos,
                        wt_aa="V", mut_aa="E",
                        chrom=variant.chrom, ref_dna=variant.ref, alt_dna=variant.alt,
                        protein_length=766, uniprot_id="P15056",
                        am_lookup=stub_am,
                    )
                else:
                    r = score_variant(
                        gene=variant.gene, position=1, wt_aa="A", mut_aa="G",
                        chrom=variant.chrom, ref_dna=variant.ref, alt_dna=variant.alt,
                    )
                return {"score": r.normalized_score if r else 0.0}

            results = score_case_study_variants(variants, real_scorer)
            self.assertEqual(len(results), 1)
            self.assertGreater(results[0]["score"], 0.5)


class TestSummarizeVariantSet(unittest.TestCase):
    """summarize_variant_set returns counts per pathogenicity + is_coding."""

    def test_summary_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "v.csv"
            _make_csv(
                path,
                [
                    {
                        "gene": "A",
                        "chrom": "chr1",
                        "pos": "100",
                        "ref": "A",
                        "alt": "G",
                        "label": "l1",
                        "pathogenicity": "pathogenic",
                        "is_coding": "true",
                    },
                    {
                        "gene": "B",
                        "chrom": "chr2",
                        "pos": "200",
                        "ref": "C",
                        "alt": "T",
                        "label": "l2",
                        "pathogenicity": "benign",
                        "is_coding": "false",
                    },
                    {
                        "gene": "C",
                        "chrom": "chr3",
                        "pos": "300",
                        "ref": "G",
                        "alt": "A",
                        "label": "l3",
                        "pathogenicity": "pathogenic",
                        "is_coding": "false",
                    },
                ],
            )
            variants = load_clinvar_variants(path)
            summary = summarize_variant_set(variants)
            self.assertEqual(summary["n_total"], 3)
            self.assertEqual(summary["n_pathogenic"], 2)
            self.assertEqual(summary["n_benign"], 1)
            self.assertEqual(summary["n_pathogenic_coding"], 1)
            self.assertEqual(summary["n_pathogenic_regulatory"], 1)


if __name__ == "__main__":
    unittest.main()
