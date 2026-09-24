"""Tests for the v0.29.0 public-API refinements.

Three changes are tested here:
  1. ``VariantInputs`` dataclass + ``score_variant_from_inputs()`` —
     the new structured-input API for ``score_variant()``.
  2. ``_safe_selector`` helper — the deduplicated try-import-then-fall-back
     pattern used by ``_try_am_lookup``, ``_try_avi_lookup``, and
     ``_try_conservation_lookup``.
  3. ``filter_variants_with_lookups(score_only=True)`` — explicit
     score-only mode that removes the previous inference ambiguity.
"""

from __future__ import annotations

import unittest

from mrnavax._adapter_selectors import _safe_selector
from mrnavax._scrna_filter import filter_variants_with_lookups
from mrnavax.sc_rna_pipeline import Variant
from mrnavax.variant_scorer import (
    VariantInputs,
    score_variant,
    score_variant_from_inputs,
)


class TestVariantInputs(unittest.TestCase):
    """The new structured-input API for score_variant()."""

    def test_from_kwargs_minimal(self):
        inputs = VariantInputs.from_kwargs("BRAF", 600, "V", "E")
        self.assertEqual(inputs.gene, "BRAF")
        self.assertEqual(inputs.position, 600)
        self.assertEqual(inputs.wt_aa, "V")
        self.assertEqual(inputs.mut_aa, "E")
        # All optional fields default to None / False.
        self.assertIsNone(inputs.chrom)
        self.assertIsNone(inputs.ref_dna)
        self.assertIsNone(inputs.alt_dna)
        self.assertIsNone(inputs.protein_length)
        self.assertIsNone(inputs.protein_sequence)
        self.assertIsNone(inputs.uniprot_id)
        self.assertIsNone(inputs.am_lookup)
        self.assertIsNone(inputs.avi_lookup)
        self.assertIsNone(inputs.conservation_lookup)
        self.assertIsNone(inputs.driver_genes)
        self.assertFalse(inputs.strict)
        self.assertIsNone(inputs.utr5)
        self.assertIsNone(inputs.utr3)

    def test_from_kwargs_with_options(self):
        inputs = VariantInputs.from_kwargs(
            "BRAF", 600, "V", "E",
            chrom="chr7",
            pos=140753336,
            protein_length=766,
            utr5="GGGCGACGCGGTGGCGGCCACCAAT",
        )
        self.assertEqual(inputs.chrom, "chr7")
        self.assertEqual(inputs.pos, 140753336)
        self.assertEqual(inputs.protein_length, 766)
        self.assertIsNotNone(inputs.utr5)


class TestScoreVariantFromInputs(unittest.TestCase):
    """score_variant_from_inputs() produces the same result as score_variant()."""

    def test_equivalent_to_score_variant(self):
        r_kwargs = score_variant("BRAF", 600, "V", "E", protein_length=766)
        r_inputs = score_variant_from_inputs(
            VariantInputs.from_kwargs("BRAF", 600, "V", "E", protein_length=766)
        )
        self.assertIsNotNone(r_kwargs)
        self.assertIsNotNone(r_inputs)
        # Same logical result (normalized_score must match exactly).
        self.assertEqual(r_inputs.normalized_score, r_kwargs.normalized_score)
        self.assertEqual(r_inputs.gene, r_kwargs.gene)
        self.assertEqual(r_inputs.position, r_kwargs.position)

    def test_with_utr(self):
        r_inputs = score_variant_from_inputs(
            VariantInputs.from_kwargs(
                "BRAF", 600, "V", "E",
                protein_length=766,
                utr5="GGGCGACGCGGTGGCGGCCACCAAT",
                utr3="A" * 200,
            )
        )
        self.assertIsNotNone(r_inputs)
        self.assertIn("utr_context_score", r_inputs.components)


class TestSafeSelector(unittest.TestCase):
    """_safe_selector deduplicates the try-import-then-fall-back pattern."""

    def test_returns_real_when_module_imports(self):
        # ``os.getcwd()`` is callable and takes no required args.
        result, is_real = _safe_selector("os", "getcwd")
        self.assertTrue(is_real)
        self.assertIsNotNone(result)
        self.assertIsInstance(result, str)

    def test_returns_none_when_module_missing(self):
        result, is_real = _safe_selector(
            "definitely_not_a_real_module_xyz_123", "anything"
        )
        self.assertFalse(is_real)
        self.assertIsNone(result)

    def test_returns_none_when_factory_missing(self):
        result, is_real = _safe_selector("os", "not_a_real_attribute_xyz_456")
        self.assertFalse(is_real)
        self.assertIsNone(result)

    def test_passes_factory_args_and_kwargs(self):
        # ``dict`` is callable and accepts args/kwargs.
        result, is_real = _safe_selector(
            "builtins", "dict", factory_args=(), factory_kwargs={"a": 1}
        )
        self.assertTrue(is_real)
        self.assertEqual(result, {"a": 1})

    def test_extra_setup_runs_on_success(self):
        calls = []

        def setup(_result):
            calls.append("ran")

        _safe_selector("os", "getcwd", extra_setup=setup)
        self.assertEqual(calls, ["ran"])

    def test_extra_setup_skipped_on_failure(self):
        calls = []

        def setup(_result):
            calls.append("ran")

        _safe_selector(
            "definitely_not_a_real_module", "x", extra_setup=setup
        )
        self.assertEqual(calls, [])


class TestFilterVariantsScoreOnlyMode(unittest.TestCase):
    """score_only=True gives the score-each behavior explicitly."""

    def test_score_only_no_dna_coords(self):
        variants = [
            Variant(gene="BRAF", position=600, wt_aa="V", mut_aa="E"),
        ]
        kept, scores, _ = filter_variants_with_lookups(
            variants,
            top_fraction=1.0,
            min_score=0.0,
            proteins={},
            uniprot_ids={},
            score_only=True,
        )
        # No DNA coords → empty scores even with score_only=True
        self.assertEqual(kept, variants)
        self.assertEqual(scores, {})

    def test_score_only_with_dna_coords(self):
        variants = [
            Variant(
                gene="BRAF", position=600, wt_aa="V", mut_aa="E",
                chrom="chr7", ref_dna="T", alt_dna="A",
            ),
        ]
        kept, scores, _ = filter_variants_with_lookups(
            variants,
            top_fraction=1.0,
            min_score=0.0,
            proteins={"BRAF": "M" * 766},
            uniprot_ids={"BRAF": "P15056"},
            score_only=True,
        )
        # DNA coords + score_only → per-variant scores populated
        self.assertEqual(kept, variants)
        self.assertEqual(len(scores), 1)
        self.assertTrue(any("BRAF" in k for k in scores))
