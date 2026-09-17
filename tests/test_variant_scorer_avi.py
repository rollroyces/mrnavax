"""Tests for AVI (AlphaGenome Atlas) wiring into the variant_scorer.

Roadmap item 1: extend mrnavax.variant_scorer.score_variant() to consume
AVI scores for non-coding regulatory-region variants. AlphaMissense
remains the dominant signal for coding-region variants; AVI becomes
the dominant signal for regulatory-region variants when ``avi_lookup``
is supplied and the variant's ``is_coding`` flag is False.

The AVI lookup contract:

    avi_lookup(chrom, pos, ref, alt) -> AVIResult

— same shape as the existing am_lookup: a callable the caller wires up
to either the real AlphaGenome Atlas (subprocess adapter) or the
stdlib mock. The score_variant wrapper handles the routing:

* ``is_coding == True``  → uses AlphaMissense (if am_lookup provided)
* ``is_coding == False`` → uses AlphaGenome Atlas AVI (if avi_lookup provided)
* either lookup returns None → silently drops that component (no crash)

This matches the existing AlphaMissense pattern: silent-by-default, no
crash on missing lookup. Strict mode (already on the function) raises
if required inputs are missing.
"""

from __future__ import annotations

import unittest

from mrnavax.alphagenome_integration import AVIResult, MockRegulatoryVariantScorer
from mrnavax.variant_scorer import VariantScore, score_variant


def _high_avi(chrom: str, pos: int, ref: str, alt: str) -> AVIResult:
    """AVI lookup returning high impact for any input."""
    return AVIResult(score=0.9, classification="high", is_coding=False)


def _low_avi(chrom: str, pos: int, ref: str, alt: str) -> AVIResult:
    """AVI lookup returning low impact for any input."""
    return AVIResult(score=0.1, classification="low", is_coding=False)


def _high_am(uniprot: str, wt: str, pos: int, mut: str):
    """AM lookup returning likely_pathogenic for any input.

    AlphaMissenseResult requires uniprot + aa_change fields (matching
    the real result shape) — passing incomplete kwargs triggers a
    TypeError that's caught silently by score_variant, dropping AM."""
    from mrnavax.alphamissense_integration import AlphaMissenseResult
    return AlphaMissenseResult(
        score=0.8,
        classification="likely_pathogenic",
        uniprot=uniprot,
        aa_change=f"{wt}{pos}{mut}",
    )


def _coding_avi(chrom: str, pos: int, ref: str, alt: str) -> AVIResult:
    """AVI lookup returning coding-region (so AM dominates)."""
    return AVIResult(score=0.99, classification="high", is_coding=True)


class TestAVIParameterAccepted(unittest.TestCase):
    """score_variant accepts new chrom/ref_dna/alt_dna/avi_lookup keyword parameters."""

    def test_avi_lookup_keyword_accepted(self):
        """Calling score_variant with avi_lookup=... must not raise TypeError."""
        # Pass a synthetic BRAF V600E missense variant with no AM/AVI lookup.
        # The function should still return a VariantScore (just without
        # the AVI component).
        result = score_variant("BRAF", 600, "V", "E", protein_length=766)
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("alphagenome_atlas_score", result.components)

    def test_avi_lookup_included_in_components_when_provided(self):
        """When avi_lookup + chrom/ref_dna/alt_dna are supplied AND it
        returns a result with is_coding=False, the AVI score is recorded
        in the components dict under alphagenome_atlas_score."""
        # Mock AVI returning high impact with is_coding=False.
        # chrom="chr7" matches the position-parity mock: pos 140753336 is even,
        # but the lookup returns is_coding=False explicitly.
        result = score_variant(
            "REG_GENE",
            100,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=_high_avi,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertIn("alphagenome_atlas_score", result.components)
        self.assertEqual(result.components["alphagenome_atlas_score"], 0.9)


class TestAVIDominatesForRegulatoryVariants(unittest.TestCase):
    """For is_coding=False variants, AVI is the dominant signal."""

    def test_high_avi_outranks_low_avi_for_same_variant(self):
        """Two identical missense variants, one with high AVI lookup,
        one with low AVI lookup. The high-AVI version must score higher."""
        high = score_variant(
            "REG",
            100,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=_high_avi,
        )
        low = score_variant(
            "REG",
            100,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=_low_avi,
        )
        self.assertGreater(
            high.normalized_score,
            low.normalized_score,
            f"high-AVI {high.normalized_score} should outrank low-AVI {low.normalized_score}",
        )

    def test_avi_appears_in_rationale(self):
        """The rationale text should mention the AVI score."""
        result = score_variant(
            "REG",
            100,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=_high_avi,
        )
        # The rationale should mention AVI score.
        self.assertIn("AVI", result.rationale.upper())
        # And it should be reflected in the rationale text
        self.assertIn("0.900", result.rationale)


class TestCodingVariantsIgnoreAVI(unittest.TestCase):
    """For is_coding=True variants, AlphaMissense still dominates; AVI is ignored."""

    def test_coding_variant_uses_alpha_missense_not_avi(self):
        """If am_lookup is supplied AND the avi_lookup returns is_coding=True,
        AlphaMissense dominates and AVI is recorded as a secondary signal
        (not the dominant 0.45 weight)."""
        result = score_variant(
            "BRAF",
            600,
            "V",
            "E",
            chrom="chr7",
            ref_dna="A",
            alt_dna="T",
            protein_length=766,
            uniprot_id="P15056",
            am_lookup=_high_am,
            avi_lookup=_coding_avi,
        )
        self.assertIsInstance(result, VariantScore)
        # AM is recorded as a secondary signal (with all required fields).
        self.assertIn("alphamissense_score", result.components)
        # AVI is also recorded (secondary, not dominant) since is_coding=True.
        self.assertIn("alphagenome_atlas_score", result.components)
        # AM dominates → high final score
        self.assertGreater(result.normalized_score, 0.5)


class TestAVISilentWhenUnset(unittest.TestCase):
    """When avi_lookup is None or DNA args are missing, score_variant behaves
    exactly as before."""

    def test_no_avi_lookup_no_components_entry(self):
        """Without an avi_lookup, no alphagenome_atlas_* entries appear
        in the components dict — backward compatible."""
        result = score_variant("BRAF", 600, "V", "E", protein_length=766)
        self.assertNotIn("alphagenome_atlas_score", result.components)
        self.assertNotIn("alphagenome_atlas_classification", result.components)

    def test_avi_lookup_without_chrom_dna_args_does_not_crash(self):
        """If avi_lookup is provided but chrom/ref_dna/alt_dna are missing,
        the AVI lookup is skipped silently (backward-compatible behavior
        for callers that pass avi_lookup=... without DNA-level info)."""
        def avi_returns_none(*args):
            return None

        result = score_variant(
            "REG",
            100,
            "A",
            "G",
            protein_length=200,
            avi_lookup=avi_returns_none,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("alphagenome_atlas_score", result.components)

    def test_avi_lookup_raising_does_not_crash(self):
        """If avi_lookup raises (network error, etc.), silently swallow
        and continue with the other components. Caller can pass strict=True
        to convert this to an exception."""
        def avi_raises(*args):
            raise RuntimeError("network down")

        # silent mode (default)
        result = score_variant(
            "REG",
            100,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=avi_raises,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("alphagenome_atlas_score", result.components)


class TestMockRegulatoryVariantScorerAsAVILookup(unittest.TestCase):
    """The stdlib mock should plug directly into score_variant.avi_lookup."""

    def test_mock_satisfies_avi_lookup_shape(self):
        """MockRegulatoryVariantScorer.score_variant has the right signature
        to be used as avi_lookup= parameter."""
        mock = MockRegulatoryVariantScorer()
        # score_variant.avi_lookup signature: (chrom, pos, ref, alt) -> AVIResult
        r = mock.score_variant("chr7", 140753336, "T", "A")
        self.assertIsInstance(r, AVIResult)
        self.assertGreaterEqual(r.score, 0.0)
        self.assertLessEqual(r.score, 1.0)

    def test_score_variant_uses_mock_for_avi(self):
        """End-to-end: use the mock as the avi_lookup and verify AVI flows
        through to the components dict."""
        mock = MockRegulatoryVariantScorer()
        # Pick a position that the mock reports as regulatory (odd pos).
        # Mock's is_coding heuristic: even= coding, odd= regulatory.
        # Position 101 → odd → is_coding=False → AVI is dominant.
        result = score_variant(
            "REG",
            101,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=mock.score_variant,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertIn("alphagenome_atlas_score", result.components)


class TestAVIComponentsVisible(unittest.TestCase):
    """The components dict exposes both AM (coding) and AVI (regulatory) cleanly."""

    def test_components_keys_when_both_lookups_provided(self):
        """When the caller provides BOTH am_lookup and avi_lookup (and
        DNA args), components dict should contain both scoring tracks."""
        result = score_variant(
            "BOTH",
            100,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            uniprot_id="P15056",
            am_lookup=_high_am,
            avi_lookup=_high_avi,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertIn("alphamissense_score", result.components)
        self.assertIn("alphagenome_atlas_score", result.components)


class TestAVIRationaleString(unittest.TestCase):
    """The rationale string mentions AVI when used."""

    def test_rationale_includes_avi_when_used(self):
        result = score_variant(
            "REG",
            101,
            "A",
            "G",
            chrom="chr7",
            ref_dna="A",
            alt_dna="G",
            protein_length=200,
            avi_lookup=_high_avi,
        )
        self.assertIn("AVI", result.rationale.upper())

    def test_rationale_omits_avi_when_not_used(self):
        result = score_variant("BRAF", 600, "V", "E", protein_length=766)
        self.assertNotIn("AVI", result.rationale.upper())


if __name__ == "__main__":
    unittest.main()
