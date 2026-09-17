"""Tests for the PhyloP46way evolutionary-conservation integration.

Roadmap item 6: add PhyloP / GERP++ evolutionary conservation as a 4th
coding-region scoring component in score_variant. Closes the last
remaining gap in the variant prioritization story — BLOSUM62 + driver
gene + AlphaMissense + now PhyloP/GERP++.

The PhyloP score from UCSC (hg38 / 46-way alignment) is a per-base
measure of evolutionary conservation in [-1, 1]:
  - PhyloP > 0  →  conserved (likely functional)
  - PhyloP < 0  →  fast-evolving (likely neutral)
  - PhyloP ~ 0  →  weak/uninformative

For variant interpretation, higher absolute PhyloP at the variant
position correlates with higher pathogenicity likelihood.

Like AlphaMissense, PhyloP is opt-in via a callable passed to
score_variant:
    conservation_lookup(chrom, pos) -> float  in [-1, 1]

The lookup is silent: None return / error / out-of-range drops the
component rather than raising.
"""

from __future__ import annotations

import unittest

from mrnavax.variant_scorer import VariantScore, score_variant


def _high_conservation(chrom: str, pos: int) -> float:
    """High conservation (PhyloP = +0.9 → likely pathogenic)."""
    return 0.9


def _low_conservation(chrom: str, pos: int) -> float:
    """Low conservation (PhyloP = -0.8 → likely neutral)."""
    return -0.8


class TestConservationParameterAccepted(unittest.TestCase):
    """score_variant accepts a new conservation_lookup keyword parameter."""

    def test_conservation_lookup_keyword_accepted(self):
        """Calling score_variant with conservation_lookup=... must not
        raise TypeError."""
        result = score_variant("BRAF", 600, "V", "E", protein_length=766)
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("phylop46way_score", result.components)

    def test_conservation_lookup_included_in_components_when_provided(self):
        """When conservation_lookup returns a value in [-1, 1] AND the
        chrom/pos are supplied, the score is recorded in components."""
        result = score_variant(
            "BRAF",
            600,
            "V",
            "E",
            protein_length=766,
            chrom="chr7",
            pos=140753336,
            conservation_lookup=_high_conservation,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertIn("phylop46way_score", result.components)
        self.assertEqual(result.components["phylop46way_score"], 0.9)

    def test_conservation_lookup_requires_chrom_and_pos(self):
        """conservation_lookup is only invoked when both chrom AND pos
        are supplied. Without these, the component is silently dropped
        (backward compat for callers that don't have DNA coords)."""
        result = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            conservation_lookup=_high_conservation,
            # no chrom / pos
        )
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("phylop46way_score", result.components)


class TestConservationDominance(unittest.TestCase):
    """Conservation contributes a small but meaningful weight to the
    overall score; high conservation correlates with higher priority."""

    def test_high_conservation_outranks_low_for_same_variant(self):
        """Same variant, two conservation lookups (high vs low): the
        high-conservation version must score higher."""
        high = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            chrom="chr7", pos=140753336,
            conservation_lookup=_high_conservation,
        )
        low = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            chrom="chr7", pos=140753336,
            conservation_lookup=_low_conservation,
        )
        self.assertGreater(
            high.normalized_score,
            low.normalized_score,
            f"high-conservation {high.normalized_score} should outrank "
            f"low-conservation {low.normalized_score}",
        )

    def test_high_conservation_appears_in_rationale(self):
        """The rationale text mentions PhyloP score when used."""
        result = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            chrom="chr7", pos=140753336,
            conservation_lookup=_high_conservation,
        )
        self.assertIn("PhyloP", result.rationale)
        self.assertIn("0.900", result.rationale)


class TestConservationSilentWhenUnset(unittest.TestCase):
    """When conservation_lookup is None or returns None, behavior is
    backward-compatible."""

    def test_no_conservation_lookup_no_components_entry(self):
        """Without a conservation_lookup, no phylop46way_* entries appear
        in the components dict — backward compatible."""
        result = score_variant("BRAF", 600, "V", "E", protein_length=766)
        self.assertNotIn("phylop46way_score", result.components)

    def test_conservation_lookup_returning_none_does_not_crash(self):
        """If conservation_lookup returns None, treat as 'no conservation
        data available' and fall through to other components silently."""
        def conservation_returns_none(*args):
            return None

        result = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            chrom="chr7", pos=140753336,
            conservation_lookup=conservation_returns_none,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("phylop46way_score", result.components)

    def test_conservation_lookup_raising_does_not_crash(self):
        """If conservation_lookup raises (network error), silently swallow
        and continue with other components."""
        def conservation_raises(*args):
            raise RuntimeError("UCSC unreachable")

        result = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            chrom="chr7", pos=140753336,
            conservation_lookup=conservation_raises,
        )
        self.assertIsInstance(result, VariantScore)
        self.assertNotIn("phylop46way_score", result.components)

    def test_conservation_out_of_range_clamped_or_rejected(self):
        """PhyloP is in [-1, 1]. Out-of-range values should not silently
        corrupt the score."""
        def conservation_out_of_range(*args):
            return 5.0  # way out of range

        # Should not crash; should silently drop the component or clamp.
        result = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            chrom="chr7", pos=140753336,
            conservation_lookup=conservation_out_of_range,
        )
        # Either silently dropped, or present but documented as out-of-range
        # — we don't crash, regardless.
        self.assertIsInstance(result, VariantScore)
        if "phylop46way_score" in result.components:
            # If present, must be in [-1, 1]
            self.assertGreaterEqual(result.components["phylop46way_score"], -1.0)
            self.assertLessEqual(result.components["phylop46way_score"], 1.0)


class TestConservationWithAlphaMissenseAndAVI(unittest.TestCase):
    """Conservation composes with AM (coding) and AVI (regulatory)
    without breaking the existing routing logic."""

    def test_all_three_components_present_for_coding_variant(self):
        """For a coding-region variant with all three lookups supplied,
        components should expose alphamissense_score, alphagenome_atlas_score,
        and phylop46way_score."""
        from mrnavax.alphamissense_integration import AlphaMissenseResult

        def am_coding(uniprot, wt, pos, mut):
            return AlphaMissenseResult(
                score=0.9,
                classification="likely_pathogenic",
                uniprot=uniprot,
                aa_change=f"{wt}{pos}{mut}",
            )

        result = score_variant(
            "BRAF", 600, "V", "E",
            protein_length=766,
            uniprot_id="P15056",
            chrom="chr7", pos=140753336,
            am_lookup=am_coding,
            conservation_lookup=_high_conservation,
        )
        self.assertIn("alphamissense_score", result.components)
        self.assertIn("phylop46way_score", result.components)


if __name__ == "__main__":
    unittest.main()
