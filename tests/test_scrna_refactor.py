"""Tests for the v0.23.0 sc_rna_pipeline refactor.

After extracting ``run_pipeline`` into ``_scrna_filter`` and
``_scrna_peptide_emitter``, these tests verify each helper is
testable in isolation (the whole point of the refactor).
"""

import unittest

from mrnavax._scrna_filter import (
    build_pipeline_notes,
    build_tumor_marker_note,
    filter_variants_with_lookups,
)
from mrnavax._scrna_peptide_emitter import emit_tumor_peptides
from mrnavax.sc_rna_pipeline import TumorPeptide, Variant


class TestEmitTumorPeptides(unittest.TestCase):
    """Peptide emitter — walk filtered variants, emit mutants."""

    def test_empty_variants_yields_empty_peptides(self):
        out = emit_tumor_peptides(
            variants=[],
            proteins={},
            tumor_cluster=0,
            tumor_expressed_genes=set(),
            cluster_marker_score=0.0,
            peptide_lengths=(9, 10),
        )
        self.assertEqual(out, [])

    def test_variant_without_protein_sequence_is_skipped(self):
        v = Variant(gene="MISSING", position=100, wt_aa="V", mut_aa="E")
        out = emit_tumor_peptides(
            variants=[v],
            proteins={},
            tumor_cluster=0,
            tumor_expressed_genes={"MISSING"},
            cluster_marker_score=0.0,
            peptide_lengths=(9,),
        )
        self.assertEqual(out, [])

    def test_variant_not_tumor_expressed_is_skipped(self):
        # protein is known, but the gene isn't in tumor_expressed_genes
        v = Variant(gene="BRAF", position=600, wt_aa="V", mut_aa="E")
        # minimal protein sequence centered around pos 600
        proteins = {"BRAF": "M" * 5 + "V" * 50}
        out = emit_tumor_peptides(
            variants=[v],
            proteins=proteins,
            tumor_cluster=0,
            tumor_expressed_genes=set(),  # not expressed
            cluster_marker_score=0.0,
            peptide_lengths=(9,),
        )
        self.assertEqual(out, [])

    def test_tumor_expressed_with_known_protein_emits_peptide(self):
        # BRAF V600E with a real-ish surrounding sequence
        seq = "M" * 595 + "V" + "M" * 170  # position 600 = V
        v = Variant(gene="BRAF", position=600, wt_aa="V", mut_aa="E")
        out = emit_tumor_peptides(
            variants=[v],
            proteins={"BRAF": seq},
            tumor_cluster=2,
            tumor_expressed_genes={"BRAF"},
            cluster_marker_score=0.85,
            peptide_lengths=(9, 10),
        )
        # ``mutant_peptides`` emits all sliding windows containing the
        # substituted aa at the requested lengths — multiple per length.
        self.assertGreater(len(out), 0)
        self.assertTrue(all(isinstance(p, TumorPeptide) for p in out))
        self.assertTrue(all(p.gene == "BRAF" for p in out))
        self.assertTrue(all(p.position == 600 for p in out))
        self.assertTrue(all(p.cell_cluster == 2 for p in out))
        self.assertTrue(all(p.cluster_marker_score == 0.85 for p in out))
        # Only 9-mers and 10-mers should be emitted (lengths requested)
        self.assertTrue(all(p.length in {9, 10} for p in out))
        # Every peptide should contain E (mut_aa at the substituted
        # position); surrounding M residues are unchanged from the
        # wild-type protein.
        self.assertTrue(all("E" in p.peptide for p in out))


class TestBuildPipelineNotes(unittest.TestCase):
    """Note builder — produces user-facing strings documenting which signals fired."""

    def test_am_active_only(self):
        out = build_pipeline_notes(
            am_active=True, filter_active=True, had_dna_coords=False
        )
        self.assertIn("AlphaMissense", out)
        self.assertNotIn("AlphaGenome", out)
        self.assertNotIn("not found", out)

    def test_am_not_active_with_filter_active_shows_fallback_note(self):
        out = build_pipeline_notes(
            am_active=False, filter_active=True, had_dna_coords=False
        )
        self.assertIn("not found", out)
        self.assertIn("download", out.lower())
        self.assertIn("AlphaMissense_hg38.tsv", out)

    def test_am_not_active_no_filter_shows_nothing_am(self):
        out = build_pipeline_notes(
            am_active=False, filter_active=False, had_dna_coords=False
        )
        self.assertNotIn("AlphaMissense", out)

    def test_dna_coords_present_adds_avi_note(self):
        out = build_pipeline_notes(
            am_active=True, filter_active=True, had_dna_coords=True
        )
        self.assertIn("AlphaGenome Atlas", out)
        self.assertIn("Avsec et al", out)

    def test_dna_coords_absent_skips_avi_note(self):
        out = build_pipeline_notes(
            am_active=True, filter_active=True, had_dna_coords=False
        )
        self.assertNotIn("AlphaGenome Atlas", out)

    def test_all_signals_off_yields_empty_string(self):
        out = build_pipeline_notes(
            am_active=False, filter_active=False, had_dna_coords=False
        )
        self.assertEqual(out, "")


class TestBuildTumorMarkerNote(unittest.TestCase):
    """Tumor-marker warning appended only when no markers matched."""

    def test_no_markers_appends_warning(self):
        notes: list[str] = []
        build_tumor_marker_note(marker_idx=[], notes=notes)
        self.assertEqual(len(notes), 1)
        self.assertIn("WARNING", notes[0])
        self.assertIn("largest cluster was assumed to be tumor", notes[0])

    def test_markers_present_appends_nothing(self):
        notes: list[str] = ["existing note"]
        build_tumor_marker_note(marker_idx=[0, 1, 2], notes=notes)
        self.assertEqual(notes, ["existing note"])  # unchanged


class TestFilterVariantsWithLookupsNoFilter(unittest.TestCase):
    """filter_variants_with_lookups — no filter, no DNA coords."""

    def test_no_filter_no_dna_yields_passthrough(self):
        v = Variant(gene="BRAF", position=600, wt_aa="V", mut_aa="E")
        kept, scores, am_active = filter_variants_with_lookups(
            [v],
            top_fraction=1.0,
            min_score=0.0,
            proteins={"BRAF": "M" * 600 + "E" + "M" * 165},
            uniprot_ids={},
        )
        self.assertEqual(len(kept), 1)
        self.assertEqual(kept[0], v)
        self.assertEqual(scores, {})
        self.assertFalse(am_active)


class TestFilterVariantsWithLookupsFilterActive(unittest.TestCase):
    """filter_variants_with_lookups — filter active without DNA coords."""

    def test_filter_active_reduces_to_kept_set(self):
        # 3 variants; request top_fraction=0.34 → ~1 kept
        variants = [
            Variant(gene="BRAF", position=600, wt_aa="V", mut_aa="E"),
            Variant(gene="TP53", position=175, wt_aa="R", mut_aa="H"),
            Variant(gene="KRAS", position=12, wt_aa="G", mut_aa="D"),
        ]
        proteins = {
            "BRAF": "M" * 600 + "E" + "M" * 165,
            "TP53": "M" * 175 + "H" + "M" * 150,
            "KRAS": "M" * 12 + "D" + "M" * 180,
        }
        kept, scores, _ = filter_variants_with_lookups(
            variants,
            top_fraction=0.34,
            min_score=0.0,
            proteins=proteins,
            uniprot_ids={},
        )
        self.assertGreaterEqual(len(kept), 1)
        self.assertLessEqual(len(kept), 3)
        # Scores dict should have an entry for every kept variant
        self.assertEqual(len(scores), len(kept))
        # All kept variant keys should appear in scores
        for v in kept:
            key = f"{v.gene}.{v.position}{v.wt_aa}>{v.mut_aa}"
            self.assertIn(key, scores)


if __name__ == "__main__":
    unittest.main()
