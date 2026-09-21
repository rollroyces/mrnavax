"""Backend checks for the ``scrna.*`` family.

Extracted from ``mrnavax/backends.py`` in v0.24.0. Each check is
registered with the global ``CHECKS`` registry on import.
"""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

from ._backends_registry import _example_path, register


@register("scrna.pipeline")
def _check_scrna() -> tuple[bool, str]:
    from .sc_rna_pipeline import run_pipeline

    r = run_pipeline(
        _example_path("cells.csv"),
        _example_path("variants_coding.csv"),
        _example_path("proteins.fasta"),
        hla=("HLA-A*02:01",),
        tumor_marker_genes=["TP53", "KRAS", "BRAF"],
    )
    ok = r.n_cells > 0 and r.n_candidate_peptides > 0 and r.tumor_cluster >= 0
    return ok, (
        f"cells={r.n_cells} clusters={len(set(r.cluster_labels))} "
        f"tumor={r.tumor_cluster} peptides={r.n_candidate_peptides}"
    )


@register("scrna.variant_filter")
def _check_scrna_variant_filter() -> tuple[bool, str]:
    from .sc_rna_pipeline import run_pipeline

    full = run_pipeline(
        _example_path("cells.csv"),
        _example_path("variants_coding.csv"),
        _example_path("proteins.fasta"),
        hla=("HLA-A*02:01",),
        tumor_marker_genes=["TP53", "KRAS", "BRAF"],
    )
    filt = run_pipeline(
        _example_path("cells.csv"),
        _example_path("variants_coding.csv"),
        _example_path("proteins.fasta"),
        hla=("HLA-A*02:01",),
        tumor_marker_genes=["TP53", "KRAS", "BRAF"],
        variant_filter_top_fraction=0.30,
    )
    ok = (
        filt.n_variants_after_filter < full.n_variants_after_filter
        and filt.n_candidate_peptides < full.n_candidate_peptides
    )
    return ok, (
        f"top-30% filter: {full.n_variants_input} variants → "
        f"{filt.n_variants_after_filter} kept → "
        f"{filt.n_candidate_peptides} peptides (was {full.n_candidate_peptides})"
    )


@register("scrna.pipeline_with_avi")
def _check_scrna_pipeline_with_avi() -> tuple[bool, str]:
    """End-to-end: scrna pipeline auto-wires AlphaGenome Atlas AVI for
    variants carrying DNA coordinates.

    Validates:
      1. When the variants CSV has chrom/ref_dna/alt_dna columns, the
         pipeline auto-wires the AVI lookup (select_regulatory_scorer).
      2. variant_scores dict is populated even when no filter is applied
         (so users see per-variant scores for regulatory + coding variants).
      3. The pipeline note mentions AlphaGenome Atlas AVI integration.
      4. Without DNA coordinates, the pipeline still works (backward
         compat — falls back to AM/BLOSUM62 only).
    """
    # Build a minimal variants CSV with DNA coordinates for some variants.
    # We synthesize the file using tempfile because the bundled example
    # CSV (variants_coding.csv) doesn't have DNA columns.

    from .sc_rna_pipeline import run_pipeline

    with tempfile.TemporaryDirectory() as tmp:
        dna_path = Path(tmp) / "variants_dna.csv"
        with open(dna_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["gene", "position", "wt_aa", "mut_aa", "chrom", "ref_dna", "alt_dna"])
            # Mix: one coding-region variant (chr17 TP53 R175H), one regulatory (odd pos)
            w.writerow(["TP53", "175", "R", "H", "chr17", "G", "A"])
            w.writerow(["REG1", "101", "A", "G", "chr11", "A", "G"])

        # 1 + 2 + 3: full pipeline with DNA coords; variant_scores should be populated
        # and note should mention AVI.
        report = run_pipeline(
            _example_path("cells.csv"),
            dna_path,
            _example_path("proteins.fasta"),
            hla=("HLA-A*02:01",),
            tumor_marker_genes=["TP53"],
            variant_filter_top_fraction=1.0,  # no filter, so we test
            # the no-filter-but-with-DNA path that fills variant_scores
        )
        bad: list[str] = []
        if not report.variant_scores:
            bad.append("variant_scores empty when DNA coords present")
        if "AlphaGenome" not in report.note and "AVI" not in report.note:
            bad.append(f"note missing AVI mention: {report.note!r}")

        # 4: without DNA coords, pipeline still works
        # Build a CSV without chrom/ref_dna/alt_dna.
        legacy_path = Path(tmp) / "legacy.csv"
        with open(legacy_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["gene", "position", "wt_aa", "mut_aa"])
            w.writerow(["TP53", "175", "R", "H"])
        legacy_report = run_pipeline(
            _example_path("cells.csv"),
            legacy_path,
            _example_path("proteins.fasta"),
            hla=("HLA-A*02:01",),
            tumor_marker_genes=["TP53"],
        )
        if "AlphaGenome" in legacy_report.note or "AVI" in legacy_report.note:
            bad.append(
                f"note wrongly includes AVI for legacy CSV without DNA coords: "
                f"{legacy_report.note!r}"
            )

    if bad:
        return False, "scrna+AVI integration issues: " + "; ".join(bad)
    return True, (
        f"scrna+AVI OK: DNA-coord variants → variant_scores populated "
        f"({len(report.variant_scores)} entries), note mentions "
        f"AlphaGenome Atlas AVI; legacy CSV (no DNA coords) → note "
        f"correctly omits AVI mention"
    )


@register("scrna.variant_scorer_alphamissense")
def _check_variant_scorer() -> tuple[bool, str]:
    from .variant_scorer import score_variant

    # Known driver mutations should rank highest
    braf = score_variant("BRAF", 600, "V", "E", protein_length=766)
    kras = score_variant("KRAS", 12, "G", "V", protein_length=189)
    benign = score_variant("MYC", 100, "A", "A", protein_length=439)  # silent
    ok = braf.normalized_score > kras.normalized_score > 0
    return ok, (
        f"BRAF.V600E={braf.normalized_score} > "
        f"KRAS.G12V={kras.normalized_score} > 0 (silent={benign.normalized_score})"
    )


@register("scrna.scgpt_integration")
def _check_scgpt_integration() -> tuple[bool, str]:
    """scGPT integration loads real weights and produces meaningful embeddings.

    Skip when:
      - MRNA_AI_FORCE_MOCK=1 (CI without scGPT weights)
      - MRNA_AI_SKIP_SCGPT_CHECK=1 (CI escape hatch)
      - scGPT weights are not on disk (first-run before download)

    When run:
      1. Load weights + vocab
      2. Embed 30 synthetic cells with named gene IDs (tumor vs normal markers)
      3. Verify output shape (n_cells, 512) and variance > 0
    """

    if os.environ.get("MRNA_AI_FORCE_MOCK"):
        return True, "skipped (MRNA_AI_FORCE_MOCK=1)"
    if os.environ.get("MRNA_AI_SKIP_SCGPT_CHECK"):
        return True, "skipped (MRNA_AI_SKIP_SCGPT_CHECK=1)"

    from .scgpt_integration import (
        ScGPTConfig,
        embed_with_scgpt,
        scgpt_available,
    )

    if not scgpt_available():
        return True, (
            "skipped (scGPT weights not present at ~/.cache/mrnavax/). "
            "Download from https://huggingface.co/perturblab/scgpt-human to enable."
        )

    cfg = ScGPTConfig.from_json()
    gene_set = [
        "TP53",
        "MYC",
        "KRAS",
        "EGFR",
        "EPCAM",
        "CD8A",
        "CD4",
        "CD3E",
        "PTPRC",
        "ACTB",
    ]
    import random

    random.seed(42)
    matrix = []
    for i in range(30):
        is_tumor = i < 15
        row = []
        for g in gene_set:
            if is_tumor and g in {"TP53", "MYC", "KRAS", "EGFR", "EPCAM"}:
                row.append(random.gauss(5.0, 0.5))
            elif (not is_tumor) and g in {"CD8A", "CD4", "CD3E", "PTPRC"}:
                row.append(random.gauss(5.0, 0.5))
            else:
                row.append(random.gauss(2.0, 0.3))
        matrix.append([max(0.0, v) for v in row])
    emb = embed_with_scgpt(matrix, gene_names=gene_set, max_cells=30)
    assert len(emb) == 30, f"got {len(emb)} embeddings"
    assert len(emb[0]) == cfg.d_hid, f"got dim {len(emb[0])}, expected {cfg.d_hid}"
    var = sum(sum(x * x for x in row) for row in emb) / (len(emb) * len(emb[0]))
    assert var > 1e-6, f"all-zero embeddings (var={var})"
    return True, (
        f"scGPT OK: 30 cells -> 30x{cfg.d_hid} embeddings (var={var:.3f}, "
        f"nlayers={cfg.nlayers}, vocab={cfg.ntoken})"
    )


@register("scrna.alphamissense_integration")
def _check_alphamissense_integration() -> tuple[bool, str]:
    """Confirm the AlphaMissense plug-in is wired end-to-end.

    Always runs (no network, no model download) by exercising the
    synthetic test index. Verifies:
      1. build_test_index() returns the expected number of entries.
      2. lookup() returns an AlphaMissenseResult with the right fields.
      3. variant_scorer integrates the score and reports it in rationale.
    """
    from .alphamissense_integration import build_test_index, lookup
    from .variant_scorer import score_variant

    idx = build_test_index()
    assert len(idx) == 5, f"expected 5 test entries, got {len(idx)}"
    r = lookup("P01116", "G", 12, "D", index=idx)
    assert r is not None, "KRAS.G12D lookup returned None"
    assert r.classification == "likely_pathogenic"
    assert 0.0 <= r.score <= 1.0
    # End-to-end: BRAF.V600E with synthetic AlphaMissense
    r2 = score_variant(
        "BRAF",
        600,
        "V",
        "E",
        protein_length=766,
        uniprot_id="P15056",
        am_lookup=lambda u, w, p, m: lookup(u, w, p, m, index=idx),
    )
    assert "alphamissense_score" in r2.components
    assert r2.components["alphamissense_score"] > 0.5
    return True, (
        f"AlphaMissense lookup + integration OK: "
        f"5-entry test index, BRAF.V600E norm={r2.normalized_score}, "
        f"AM={r2.components['alphamissense_score']:.3f}"
    )


@register("scrna.structural_disruption_chou_fasman")
def _check_structural_disruption() -> tuple[bool, str]:
    """L→P in a helix context must score higher than L→P in a coil context.

    Proline is a helix breaker; in a real alpha-helix, L→P scores higher
    because of the additional structural-disruption component.
    """
    from .variant_scorer import score_variant

    helix_seq = "LAELAEKLAEEK"
    coil_seq = "GGPGGPPPGGPG"
    v_helix = score_variant(
        "FAKE",
        2,
        "L",
        "P",
        protein_length=12,
        protein_sequence=helix_seq,
    )
    v_coil = score_variant(
        "FAKE",
        2,
        "L",
        "P",
        protein_length=12,
        protein_sequence=coil_seq,
    )
    ok = v_helix.normalized_score > v_coil.normalized_score
    return ok, (
        f"L→P in helix={v_helix.normalized_score:.3f} (struct={v_helix.components['structural_disruption']}) "
        f"> L→P in coil={v_coil.normalized_score:.3f} (struct={v_coil.components['structural_disruption']})"
    )


@register("scrna.embedding_tfidf_svd")
def _check_embedding() -> tuple[bool, str]:
    import math

    from .foundation_embedder import embed_cells
    from .sc_rna_pipeline import _synthetic

    _, _, matrix = _synthetic()
    emb = embed_cells(matrix, model="tfidf-svd", n_components=8)
    # Cluster 0 (cells 0-19) vs cluster 2 (cells 35-49)
    intra = sum(
        math.sqrt(sum((emb[i][k] - emb[j][k]) ** 2 for k in range(8)))
        for i in range(19)
        for j in range(i + 1, 20)
    ) / (19 * 18 / 2)
    inter = (
        sum(math.sqrt(sum((emb[0][k] - emb[j][k]) ** 2 for k in range(8))) for j in range(35, 50))
        / 15
    )
    sep = inter / intra if intra else 0
    ok = sep > 1.5  # require at least 1.5x separation
    return ok, f"shape={len(emb)}x{len(emb[0])} inter/intra separation={sep:.2f}x"


__all__ = []
