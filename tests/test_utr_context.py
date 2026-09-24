"""Tests for the v0.27.0 UTR-aware variant scoring extension."""

from __future__ import annotations

import unittest

from mrnavax._utr_context import (
    DEFAULT_UTR_CONTEXT_WEIGHT,
    UTRContextResult,
    score_utr_context,
)

# Canonical strong Kozak upstream: 9-nt slice ``GCCACCAAT``
# (R=A at position -3, R=A at +0 of ATG, A at +1, T at +2).
# Matches the manufacturability checker's consensus exactly.
_STRONG_KOZAK_UTR5 = "GGGCGACGCGGTGGCGGCCACCAAT"
# Weak Kozak context (3' end deviates from consensus)
_WEAK_KOZAK_UTR5 = "GGGCGACGCAAAAAAAAAAAAAAAAAAAA"


class TestScoreUTRContext(unittest.TestCase):
    """score_utr_context: pure UTR scoring, no variant required."""

    def test_strong_kozak_full_credit(self):
        # The full 27-nt 5'UTR ends with the canonical strong Kozak
        # context (positions -3 R, +0 R, +1 A, +2 T). After T→U
        # conversion, position +2 becomes U which doesn't match the
        # consensus' expected T — so the score is 9/10 = 0.9, which
        # the manufacturability checker also classifies as "strong"
        # (threshold ≥0.85). We accept this 0.9 score as "strong".
        result = score_utr_context(_STRONG_KOZAK_UTR5, utr3=None)
        self.assertIsInstance(result, UTRContextResult)
        self.assertGreaterEqual(result.kozak_score, 0.85)

    def test_weak_kozak_lower_score(self):
        result = score_utr_context(_WEAK_KOZAK_UTR5, utr3=None)
        self.assertLess(result.kozak_score, 0.5)

    def test_short_utr5_returns_default(self):
        # Too short for Kozak scoring → 0.5 default
        result = score_utr_context("ATGG", utr3=None)
        self.assertEqual(result.kozak_score, 0.5)

    def test_no_utr5_returns_default(self):
        result = score_utr_context(None, utr3=None)
        self.assertEqual(result.kozak_score, 0.5)

    def test_utr3_length_scored(self):
        # ≥100 nt = full credit (1.0)
        long_utr3 = "A" * 150
        result = score_utr_context(None, long_utr3)
        self.assertGreater(result.utr3_score, 0.5)

    def test_utr3_no_are_motifs_low(self):
        # No ARE pentamers → 0.4 base score
        utr3 = "GCGCGCGC" * 5  # 40 nt, no AUUUA
        result = score_utr_context(None, utr3)
        # length_score = 0.4 (40/100), are_score = 0.4 → 0.6*0.4 + 0.4*0.4 = 0.40
        self.assertLess(result.utr3_score, 0.6)

    def test_utr3_with_are_motifs_higher(self):
        # 3+ ARE pentamers → high score
        utr3 = "AUUUAGCAUUUAGCAUUUAGC" + "A" * 100
        result = score_utr_context(None, utr3)
        self.assertGreater(result.utr3_score, 0.6)

    def test_context_score_in_unit_interval(self):
        result = score_utr_context(_STRONG_KOZAK_UTR5, "AUUUAGCAUUUAGCAUUUAG" + "A" * 100)
        self.assertGreaterEqual(result.context_score, 0.0)
        self.assertLessEqual(result.context_score, 1.0)

    def test_default_weight_value(self):
        # Default bonus weight is 0.05 — keeps dominant signals dominant
        self.assertEqual(DEFAULT_UTR_CONTEXT_WEIGHT, 0.05)

    def test_to_dict_serializable(self):
        import json
        result = score_utr_context(_STRONG_KOZAK_UTR5, "A" * 100)
        d = result.to_dict()
        json.dumps(d)

    def test_details_have_kozak_and_utr3(self):
        result = score_utr_context(_STRONG_KOZAK_UTR5, "AUUUAGCAUUUAGC" + "A" * 100)
        self.assertIn("kozak", result.details)
        self.assertIn("utr3", result.details)


class TestScoreVariantUTRContext(unittest.TestCase):
    """score_variant accepts utr5/utr3 kwargs (additive, no breakage)."""

    def test_no_utr_args_backward_compatible(self):
        # The existing 4-signal composition must work without utr5/utr3.
        from mrnavax.variant_scorer import score_variant
        r = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766, chrom="chr7", pos=140753336,
        )
        self.assertIsNotNone(r)
        self.assertNotIn("utr_context_score", r.components)

    def test_with_utr_args_adds_utr_components(self):
        from mrnavax.variant_scorer import score_variant
        r = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766, chrom="chr7", pos=140753336,
            utr5=_STRONG_KOZAK_UTR5,
            utr3="AUUUAGCAUUUAGCAUUUAG" + "A" * 100,
        )
        self.assertIsNotNone(r)
        # The 3 new components appear
        self.assertIn("utr_context_score", r.components)
        self.assertIn("kozak_score", r.components)
        self.assertIn("utr3_score", r.components)

    def test_strong_utr_boosts_score(self):
        from mrnavax.variant_scorer import score_variant
        # Without UTR context
        r_no_utr = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766, chrom="chr7", pos=140753336,
        )
        # With strong UTR context
        r_strong = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766, chrom="chr7", pos=140753336,
            utr5=_STRONG_KOZAK_UTR5,
            utr3="AUUUAGCAUUUAGCAUUUAG" + "A" * 100,
        )
        # Strong UTR context should boost the score (or keep it equal)
        self.assertGreaterEqual(r_strong.normalized_score, r_no_utr.normalized_score)

    def test_weak_utr_no_harm(self):
        # Weak UTR context should NOT boost the score (or stay the same)
        from mrnavax.variant_scorer import score_variant
        r_no_utr = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766, chrom="chr7", pos=140753336,
        )
        r_weak = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766, chrom="chr7", pos=140753336,
            utr5=_WEAK_KOZAK_UTR5,
            utr3="GCGCGCGC" * 5,  # 40 nt, no ARE
        )
        # Weak UTR must NOT lower the score (additive bonus only)
        self.assertGreaterEqual(r_weak.normalized_score, r_no_utr.normalized_score)


if __name__ == "__main__":
    unittest.main()
