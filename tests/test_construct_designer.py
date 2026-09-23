"""Tests for the mRNA construct designer (mrnavax/construct_designer.py)."""

from __future__ import annotations

import json
import unittest

from mrnavax.construct_designer import (
    _POLYA_SIGNAL,
    ConstructConfig,
    ConstructResult,
    _reverse_translate,
    design_construct,
)

# A small test protein (eGFP N-terminus fragment, 12 AA).
_TEST_PROTEIN = "MVSKGEELFTGV"


class TestReverseTranslate(unittest.TestCase):
    """Default reverse translation uses the most-frequent human codon."""

    def test_known_codons(self):
        # M = ATG (single), V = GTG (highest freq in human), K = AAG,
        # G = GGC, E = GAG, L = CTG, F = TTC, T = ACC
        result = _reverse_translate("MVKGE")
        self.assertEqual(result, "ATGGTGAAGGGCGAG")

    def test_all_20_amino_acids(self):
        # Each AA should reverse-translate to a known codon
        for aa in "ACDEFGHIKLMNPQRSTVWY":
            cds = _reverse_translate(aa)
            self.assertEqual(len(cds), 3, f"AA {aa}: expected 3 nt, got {cds!r}")
            # Should be uppercase ACGT only
            self.assertTrue(set(cds) <= set("ACGT"))

    def test_empty_protein_raises(self):
        with self.assertRaises(ValueError):
            _reverse_translate("")

    def test_unknown_amino_acid_raises(self):
        with self.assertRaises(ValueError):
            _reverse_translate("MVAU")  # U is not an AA


class TestConstructConfig(unittest.TestCase):
    """ConstructConfig defaults are sensible."""

    def test_default_species_is_human(self):
        self.assertEqual(ConstructConfig().species, "human")

    def test_default_poly_a_length(self):
        self.assertEqual(ConstructConfig().poly_a_length, 120)

    def test_default_backend_is_multi_objective(self):
        # The v0.25.0 multi-objective backend is the SOTA default.
        self.assertEqual(ConstructConfig().optimize_backend, "multi-objective")

    def test_default_utrs_are_bundled(self):
        cfg = ConstructConfig()
        self.assertGreater(len(cfg.utr5), 10)
        self.assertGreater(len(cfg.utr3), 50)


class TestDesignConstruct(unittest.TestCase):
    """The main entry point."""

    def test_basic_design_produces_full_construct(self):
        result = design_construct(_TEST_PROTEIN)
        self.assertIsInstance(result, ConstructResult)
        # Construct length = utr5 + cds + utr3 + polyA_signal + polyA
        expected_len = (
            len(result.utr5)
            + result.cds_length
            + len(result.utr3)
            + len(_POLYA_SIGNAL)
            + len(result.poly_a_tail)
        )
        self.assertEqual(len(result.construct_dna), expected_len)
        # Construct starts with 5'UTR and ends with poly-A tail
        self.assertTrue(result.construct_dna.startswith(result.utr5))
        self.assertTrue(result.construct_dna.endswith(result.poly_a_tail))

    def test_cds_preserves_amino_acid_sequence(self):
        from mrnavax.codon_optimizer import CODON_TO_AA
        result = design_construct(_TEST_PROTEIN)
        cds_aas = [
            CODON_TO_AA[result.cds[i : i + 3]]
            for i in range(0, len(result.cds), 3)
        ]
        self.assertEqual("".join(cds_aas), _TEST_PROTEIN)

    def test_cds_length_multiple_of_three(self):
        result = design_construct(_TEST_PROTEIN)
        self.assertEqual(result.cds_length % 3, 0)
        self.assertEqual(result.n_codons * 3, result.cds_length)

    def test_n_codons_matches_protein_length(self):
        result = design_construct(_TEST_PROTEIN)
        self.assertEqual(result.n_codons, len(_TEST_PROTEIN))

    def test_poly_a_tail_length_matches_config(self):
        result = design_construct(_TEST_PROTEIN, ConstructConfig(poly_a_length=150))
        self.assertEqual(len(result.poly_a_tail), 150)
        self.assertTrue(result.poly_a_tail == "A" * 150)
        self.assertTrue(result.construct_dna.endswith("A" * 150))

    def test_backend_basic_works(self):
        cfg = ConstructConfig(optimize_backend="basic")
        result = design_construct(_TEST_PROTEIN, cfg)
        self.assertEqual(result.cds_length, len(_TEST_PROTEIN) * 3)

    def test_backend_multi_objective_adds_note(self):
        cfg = ConstructConfig(optimize_backend="multi-objective")
        result = design_construct(_TEST_PROTEIN, cfg)
        self.assertGreater(len(result.optimization_notes), 0)
        self.assertIn("multi-objective", result.optimization_notes[0])

    def test_unknown_backend_raises(self):
        cfg = ConstructConfig(optimize_backend="not-a-real-backend")
        with self.assertRaises(ValueError):
            design_construct(_TEST_PROTEIN, cfg)

    def test_to_dict_serializable(self):
        result = design_construct(_TEST_PROTEIN)
        d = result.to_dict()
        json.dumps(d)  # must be JSON-serializable

    def test_cds_analysis_populated(self):
        result = design_construct(_TEST_PROTEIN)
        self.assertIn("cai", result.cds_analysis)
        self.assertIn("gc_percent", result.cds_analysis)
        self.assertIn("n_codons", result.cds_analysis)
        # CAI must be in (0, 1]
        self.assertGreater(result.cds_analysis["cai"], 0.0)
        self.assertLessEqual(result.cds_analysis["cai"], 1.0)

    def test_gc_percent_in_range(self):
        result = design_construct(_TEST_PROTEIN)
        self.assertGreaterEqual(result.gc_percent, 0.0)
        self.assertLessEqual(result.gc_percent, 100.0)

    def test_empty_protein_raises(self):
        with self.assertRaises(ValueError):
            design_construct("")

    def test_invalid_protein_raises(self):
        with self.assertRaises(ValueError):
            design_construct("MVAU")


class TestLargerConstructs(unittest.TestCase):
    """Larger inputs (full-length eGFP) work end-to-end."""

    def test_full_length_egfp(self):
        # Full-length eGFP, 239 AA
        egfp = (
            "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
            "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
            "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHY"
            "QQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
        )
        result = design_construct(egfp)
        self.assertEqual(result.n_codons, 239)
        self.assertEqual(result.cds_length, 239 * 3)
        # Construct should be substantial
        self.assertGreater(result.construct_length, 1000)


if __name__ == "__main__":
    unittest.main()
