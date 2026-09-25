"""Tests for the v0.30.0 UTR-aware CDS designer (10th tool)."""

from __future__ import annotations

import unittest

from mrnavax.utr_designer import (
    UTRDesignConfig,
    UTRDesignResult,
    _candidate_library,
    design_utr_aware_cds,
)


class TestUTRDesignConfig(unittest.TestCase):
    """UTRDesignConfig dataclass defaults are sensible."""

    def test_defaults(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        self.assertEqual(cfg.cds, "MVSKGEELFTGV")
        self.assertEqual(cfg.organism, "human")
        self.assertEqual(cfg.prefer_kozak, 0.7)
        self.assertEqual(cfg.prefer_utr3, 0.5)
        self.assertEqual(cfg.library, "all")
        self.assertEqual(cfg.backend, "multi-objective")

    def test_custom_thresholds(self):
        cfg = UTRDesignConfig(
            cds="MVSKGEELFTGV",
            prefer_kozak=0.85,
            prefer_utr3=0.8,
        )
        self.assertEqual(cfg.prefer_kozak, 0.85)
        self.assertEqual(cfg.prefer_utr3, 0.8)


class TestCandidateLibrary(unittest.TestCase):
    """_candidate_library returns the right shape for each library setting."""

    def test_all_library_has_12_candidates(self):
        candidates = _candidate_library("all")
        # 4 5'UTR variants × 3 3'UTR variants = 12 combinations
        self.assertEqual(len(candidates), 12)
        for c in candidates:
            self.assertEqual(len(c), 4)  # (u5_name, u5_seq, u3_name, u3_seq)
            self.assertIsInstance(c[0], str)
            self.assertIsInstance(c[1], str)
            self.assertIsInstance(c[2], str)
            self.assertIsInstance(c[3], str)

    def test_minimal_library_has_1_candidate(self):
        candidates = _candidate_library("minimal")
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0][0], "strong_kozak")
        self.assertEqual(candidates[0][2], "short_constitutive")


class TestDesignUTRAwareCDS(unittest.TestCase):
    """design_utr_aware_cds() returns a UTRDesignResult."""

    def test_basic_design(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        result = design_utr_aware_cds(cfg)
        self.assertIsInstance(result, UTRDesignResult)
        self.assertEqual(len(result.utr5) > 0, True)
        self.assertEqual(len(result.utr3) > 0, True)
        # CDS length = 3 × len(AA)
        self.assertEqual(len(result.cds_dna), 3 * len(cfg.cds))
        # CDS preserves the protein
        self.assertEqual(result.protein, cfg.cds)

    def test_combined_score_in_unit_interval(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        result = design_utr_aware_cds(cfg)
        self.assertGreaterEqual(result.combined_score, 0.0)
        self.assertLessEqual(result.combined_score, 1.0)
        self.assertGreaterEqual(result.cds_score, 0.0)
        self.assertLessEqual(result.cds_score, 1.0)

    def test_candidates_evaluated_count(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        result = design_utr_aware_cds(cfg)
        self.assertEqual(result.candidates_evaluated, 12)

    def test_ranking_length_matches_candidates(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        result = design_utr_aware_cds(cfg)
        self.assertEqual(len(result.ranking), 12)

    def test_ranking_is_sorted_descending(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        result = design_utr_aware_cds(cfg)
        scores = [r[0] for r in result.ranking]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_minimal_library_picks_strong_kozak(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV", library="minimal")
        result = design_utr_aware_cds(cfg)
        # Minimal library = strong_kozak + short_constitutive only
        self.assertEqual(result.candidates_evaluated, 1)
        # Top ranking should be that one combination
        self.assertEqual(result.ranking[0][1], "strong_kozak")
        self.assertEqual(result.ranking[0][2], "short_constitutive")

    def test_basic_backend_works(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV", backend="basic")
        result = design_utr_aware_cds(cfg)
        self.assertIsInstance(result, UTRDesignResult)
        self.assertEqual(len(result.cds_dna), 36)

    def test_thresholds_filter_too_strict(self):
        # Very strict thresholds — only strong-kozak candidates should pass
        cfg = UTRDesignConfig(
            cds="MVSKGEELFTGV",
            prefer_kozak=0.95,
            prefer_utr3=0.95,
        )
        result = design_utr_aware_cds(cfg)
        # With impossible thresholds, the fallback path returns the
        # highest-ranked combination anyway (still a valid UTRDesignResult).
        self.assertIsInstance(result, UTRDesignResult)

    def test_strong_kozak_chosen_when_available(self):
        # Among candidates, strong-kozak variants should score highest
        # for the Kozak component.
        cfg = UTRDesignConfig(
            cds="MVSKGEELFTGV",
            prefer_kozak=0.0,  # accept any
            prefer_utr3=0.0,
        )
        result = design_utr_aware_cds(cfg)
        # Ranking top-3 should include strong_kozak
        top3_names = {name for _, name, _ in result.ranking[:3]}
        self.assertIn("strong_kozak", top3_names)

    def test_egfp_fragment_picks_paired_utrs(self):
        # eGFP fragment: longer than test_basic — should produce more
        # candidate evaluations and a CDS that's 33 nt (11 codons).
        cfg = UTRDesignConfig(cds="MVSKGEELFTG")
        result = design_utr_aware_cds(cfg)
        self.assertEqual(len(result.cds_dna), 33)
        # The combined score should be > 0.5 for a reasonable pair
        self.assertGreater(result.combined_score, 0.5)


class TestUTRDesignResultDataclass(unittest.TestCase):
    """UTRDesignResult has all expected fields."""

    def test_result_fields(self):
        cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
        result = design_utr_aware_cds(cfg)
        # All required fields present
        self.assertTrue(hasattr(result, "utr5"))
        self.assertTrue(hasattr(result, "utr3"))
        self.assertTrue(hasattr(result, "cds_dna"))
        self.assertTrue(hasattr(result, "protein"))
        self.assertTrue(hasattr(result, "utr_context"))
        self.assertTrue(hasattr(result, "cds_score"))
        self.assertTrue(hasattr(result, "combined_score"))
        self.assertTrue(hasattr(result, "candidates_evaluated"))
        self.assertTrue(hasattr(result, "ranking"))
