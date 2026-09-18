"""Internal peptide-candidate emitter for ``sc_rna_pipeline.run_pipeline``.

Extracted from ``mrnavax/sc_rna_pipeline.py`` in v0.23.0. Walks the
filtered variants and emits mutant-peptide candidates for those that
are both (a) tumor-cluster-expressed and (b) have a known protein
sequence. Stdlib-only.
"""

from __future__ import annotations

from .sc_rna_pipeline import (
    TumorPeptide,
    Variant,
    mutant_peptides,
)


def emit_tumor_peptides(
    variants: list[Variant],
    proteins: dict[str, str],
    tumor_cluster: int,
    tumor_expressed_genes: set[str],
    cluster_marker_score: float,
    peptide_lengths: tuple[int, ...],
) -> list[TumorPeptide]:
    """Walk the filtered variants and emit mutant-peptide candidates.

    For each variant that:
      * has a known protein sequence in ``proteins`` (else skipped —
        peptide enumeration requires the wild-type sequence), AND
      * is expressed in the tumor cluster (else skipped — the
        downstream peptide list is for tumor-present mutations only)
    …emit one ``TumorPeptide`` per mutant peptide at each requested
    ``peptide_lengths`` size.
    """
    peptides: list[TumorPeptide] = []
    for v in variants:
        if v.gene not in proteins:
            continue
        if v.gene not in tumor_expressed_genes:
            continue
        for pep in mutant_peptides(
            proteins[v.gene], v.position, v.mut_aa, lengths=peptide_lengths
        ):
            peptides.append(
                TumorPeptide(
                    cell_cluster=tumor_cluster,
                    gene=v.gene,
                    position=v.position,
                    wt_aa=v.wt_aa,
                    mut_aa=v.mut_aa,
                    peptide=pep,
                    length=len(pep),
                    cluster_marker_score=cluster_marker_score,
                )
            )
    return peptides


__all__ = ["emit_tumor_peptides"]
