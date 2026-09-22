"""Tests for mrnavax/codon_multi_objective.py — SOTA multi-objective codon optimizer."""

from __future__ import annotations

import unittest

from mrnavax.codon_multi_objective import (
    ComponentBreakdown,
    MultiObjectiveConfig,
    MultiObjectiveResult,
    _component_breakdown,
    _overall_score,
    multi_objective_optimize,
    multi_objective_score,
)

# A short, real CDS (eGFP fragment, ~30 AA) — small enough for unit tests,
# long enough to exercise per-component breakdown.
_EGFP_FRAGMENT = (
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGT"
    "CGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGG"
    "CGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCAC"
    "CGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGG"
)


class TestComponentBreakdown(unittest.TestCase):
    """Per-component scorers return values in [0, 1]."""

    def test_all_components_in_unit_interval(self):
        cfg = MultiObjectiveConfig()
        b = _component_breakdown(_EGFP_FRAGMENT, cfg)
        self.assertGreaterEqual(b.cai, 0.0)
        self.assertLessEqual(b.cai, 1.0)
        self.assertGreaterEqual(b.gc_score, 0.0)
        self.assertLessEqual(b.gc_score, 1.0)
        self.assertGreaterEqual(b.cpg_score, 0.0)
        self.assertLessEqual(b.cpg_score, 1.0)
        self.assertGreaterEqual(b.rare_run_penalty, 0.0)
        self.assertLessEqual(b.rare_run_penalty, 1.0)
        self.assertGreaterEqual(b.structure_proxy, 0.0)
        self.assertLessEqual(b.structure_proxy, 1.0)

    def test_to_dict_round_trip(self):
        cfg = MultiObjectiveConfig()
        b = _component_breakdown(_EGFP_FRAGMENT, cfg)
        d = b.to_dict()
        self.assertEqual(set(d.keys()), {
            "cai", "gc_score", "cpg_score", "rare_run_penalty", "structure_proxy"
        })


class TestOverallScore(unittest.TestCase):
    """Weighted-sum is clamped to [0, 1]."""

    def test_all_weights_zero_returns_zero(self):
        cfg = MultiObjectiveConfig(
            cai_weight=0, gc_weight=0, cpg_weight=0,
            rare_run_weight=0, structure_weight=0,
        )
        b = ComponentBreakdown(cai=0.5, gc_score=0.5, cpg_score=0.5,
                                rare_run_penalty=0.5, structure_proxy=0.5)
        self.assertEqual(_overall_score(b, cfg), 0.0)

    def test_high_components_high_score(self):
        cfg = MultiObjectiveConfig()
        # With all components at 1.0 except rare_run_penalty (which is
        # subtracted) at 0.0, the weighted sum is 0.40+0.20+0.15+0.10
        # = 0.85 (the rare_run_weight cancels because penalty=0).
        b = ComponentBreakdown(cai=1.0, gc_score=1.0, cpg_score=1.0,
                                rare_run_penalty=0.0, structure_proxy=1.0)
        self.assertEqual(_overall_score(b, cfg), 0.85)

    def test_low_components_low_score(self):
        cfg = MultiObjectiveConfig()
        b = ComponentBreakdown(cai=0.0, gc_score=0.0, cpg_score=0.0,
                                rare_run_penalty=1.0, structure_proxy=0.0)
        # With default weights, the rare_run penalty subtracts 0.15,
        # so the result is 0 (clamped, not negative).
        score = _overall_score(b, cfg)
        self.assertEqual(score, 0.0)

    def test_score_clamped_to_unit_interval(self):
        cfg = MultiObjectiveConfig(rare_run_weight=100.0)
        b = ComponentBreakdown(cai=1.0, gc_score=1.0, cpg_score=1.0,
                                rare_run_penalty=1.0, structure_proxy=1.0)
        score = _overall_score(b, cfg)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)


class TestMultiObjectiveScore(unittest.TestCase):
    """multi_objective_score() — score-only entry point."""

    def test_valid_cds_returns_float_in_unit_interval(self):
        score = multi_objective_score(_EGFP_FRAGMENT)
        self.assertIsInstance(score, float)
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    def test_short_cds_returns_zero(self):
        # Less than one codon after cleaning
        self.assertEqual(multi_objective_score("AT"), 0.0)
        self.assertEqual(multi_objective_score(""), 0.0)

    def test_custom_config_changes_score(self):
        cfg = MultiObjectiveConfig(cai_weight=1.0, gc_weight=0,
                                    cpg_weight=0, rare_run_weight=0,
                                    structure_weight=0)
        # CAI-only: should give a non-trivial score on eGFP fragment
        score_cai_only = multi_objective_score(_EGFP_FRAGMENT, cfg)
        self.assertGreater(score_cai_only, 0.5)

    def test_u_to_t_normalization(self):
        # Same sequence in RNA form (U instead of T) should score the same
        rna_form = _EGFP_FRAGMENT.replace("T", "U")
        score_dna = multi_objective_score(_EGFP_FRAGMENT)
        score_rna = multi_objective_score(rna_form)
        self.assertEqual(score_dna, score_rna)


class TestMultiObjectiveOptimize(unittest.TestCase):
    """multi_objective_optimize() — the main entry point."""

    def test_optimization_preserves_amino_acid_sequence(self):
        from mrnavax.codon_optimizer import CODON_TO_AA
        result = multi_objective_optimize(_EGFP_FRAGMENT)
        input_cds = result.input_cds  # use the cleaned CDS (input may have leftover after trimming)
        input_aas = [
            CODON_TO_AA[input_cds[i:i+3]]
            for i in range(0, len(input_cds), 3)
        ]
        output_aas = [
            CODON_TO_AA[result.optimized_cds[i:i+3]]
            for i in range(0, len(result.optimized_cds), 3)
        ]
        self.assertEqual(input_aas, output_aas,
                         "fidelity: AA sequence must be preserved")

    def test_optimized_cds_has_same_length_as_input(self):
        result = multi_objective_optimize(_EGFP_FRAGMENT)
        self.assertEqual(len(result.optimized_cds), len(result.input_cds))

    def test_overall_after_gte_overall_before(self):
        # Optimization should not regress the standard: either improve
        # or stay equal (if the input is already locally optimal).
        result = multi_objective_optimize(_EGFP_FRAGMENT)
        self.assertGreaterEqual(result.overall_after, result.overall_before)

    def test_n_changes_matches_actual_diff(self):
        result = multi_objective_optimize(_EGFP_FRAGMENT)
        actual_diff = sum(
            1 for i in range(0, len(result.input_cds), 3)
            if result.input_cds[i:i+3] != result.optimized_cds[i:i+3]
        )
        self.assertEqual(result.n_changes, actual_diff)

    def test_result_has_correct_types(self):
        result = multi_objective_optimize(_EGFP_FRAGMENT)
        self.assertIsInstance(result, MultiObjectiveResult)
        self.assertIsInstance(result.before, ComponentBreakdown)
        self.assertIsInstance(result.after, ComponentBreakdown)
        self.assertIsInstance(result.config, MultiObjectiveConfig)

    def test_to_dict_serializable(self):
        import json
        result = multi_objective_optimize(_EGFP_FRAGMENT)
        d = result.to_dict()
        # Must be JSON-serializable
        json.dumps(d)

    def test_short_cds_raises(self):
        with self.assertRaises(ValueError):
            multi_objective_optimize("AT")

    def test_cai_only_config_increases_cai(self):
        # With CAI-only config, the post-optimization CAI should be >= pre.
        cfg = MultiObjectiveConfig(
            cai_weight=1.0, gc_weight=0, cpg_weight=0,
            rare_run_weight=0, structure_weight=0,
        )
        result = multi_objective_optimize(_EGFP_FRAGMENT, cfg)
        self.assertGreaterEqual(result.after.cai, result.before.cai)

    def test_rare_run_only_config_reduces_penalty(self):
        # With rare_run_weight dominant, the post-opt penalty should be
        # <= pre-opt penalty (fewer long rare-codon runs).
        cfg = MultiObjectiveConfig(
            cai_weight=0, gc_weight=0, cpg_weight=0,
            rare_run_weight=1.0, structure_weight=0,
        )
        result = multi_objective_optimize(_EGFP_FRAGMENT, cfg)
        self.assertLessEqual(result.after.rare_run_penalty,
                             result.before.rare_run_penalty)


class TestMultiObjectiveConfig(unittest.TestCase):
    """Configuration dataclass has sensible defaults."""

    def test_default_weights_sum_to_one(self):
        cfg = MultiObjectiveConfig()
        total = (cfg.cai_weight + cfg.gc_weight + cfg.cpg_weight
                 + cfg.rare_run_weight + cfg.structure_weight)
        self.assertAlmostEqual(total, 1.0, places=2)

    def test_default_gc_band_is_45_to_60(self):
        cfg = MultiObjectiveConfig()
        self.assertEqual(cfg.target_gc_min, 45.0)
        self.assertEqual(cfg.target_gc_max, 60.0)


if __name__ == "__main__":
    unittest.main()
