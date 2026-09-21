"""Internal variant-filter helper for ``sc_rna_pipeline.run_pipeline``.

Extracted from ``mrnavax/sc_rna_pipeline.py`` in v0.23.0 to keep
``run_pipeline`` focused on orchestration rather than detail.
Private API; downstream code should call
``mrnavax.sc_rna_pipeline.run_pipeline`` directly.

This module owns:
  * Wiring the optional AlphaMissense lookup (real index if
    cached; mock otherwise).
  * Wiring the optional AlphaGenome Atlas AVI lookup.
  * Wiring the optional PhyloP46way conservation lookup.
  * Calling ``filter_variants`` when filtering is enabled.
  * Calling ``score_variant`` per-variant when filtering is disabled
    but DNA coords are present (so the report still has per-variant
    scores for regulatory variants).
  * Building the user-facing ``notes`` strings that document what
    signals fired during the run.

Stdlib-only except for the opt-in upstream-model adapters, which
are themselves wrapped with silent failure at the boundary.
"""

from __future__ import annotations

from typing import Callable

from .sc_rna_pipeline import Variant


def _try_am_lookup() -> tuple[Callable | None, bool]:
    """Build the AlphaMissense lookup, preferring the real cached index.

    Returns (lookup_fn, is_real_index). lookup_fn is None if no real
    index is available; is_real_index is True only if a real
    AlphaMissense TSV cache was loaded.
    """
    try:
        from .alphamissense_integration import (
            CACHE_FILE,
            load_index,
            lookup,
        )

        if CACHE_FILE.exists():
            load_index()
            return lookup, True
    except Exception:
        pass
    return None, False


def _try_avi_lookup() -> Callable | None:
    """Build the AlphaGenome Atlas AVI lookup, preferring the real Atlas
    when ``ALPHAGENOME_API_KEY`` is set + ``[variant-alphagenome]``
    extra installed; mock otherwise.
    """
    try:
        from .alphagenome_integration import select_regulatory_scorer

        return select_regulatory_scorer().score_variant
    except Exception:
        return None


def _try_conservation_lookup() -> Callable | None:
    """Build the PhyloP46way conservation lookup. Mock by default."""
    try:
        from .conservation import select_conservation_lookup

        return select_conservation_lookup().lookup
    except Exception:
        return None


def _build_v_dicts(variants: list[Variant]) -> list[dict]:
    """Convert Variant objects to dicts compatible with score_variant /
    filter_variants."""
    return [
        {
            "gene": v.gene,
            "position": v.position,
            "wt_aa": v.wt_aa,
            "mut_aa": v.mut_aa,
            "chrom": v.chrom,
            "ref_dna": v.ref_dna,
            "alt_dna": v.alt_dna,
        }
        for v in variants
    ]


def _score_one_each(
    variants: list[Variant],
    proteins: dict[str, str],
    uniprot_ids: dict[str, str],
    am_lookup_fn: Callable | None,
    avi_lookup_fn: Callable | None,
    conservation_lookup_fn: Callable | None,
) -> dict[str, float]:
    """Score each variant individually (no filtering).

    Used when ``variant_filter_top_fraction == 1.0`` AND any variant
    carries DNA coordinates. The result populates
    ``filter_scores`` with per-variant scores so the report carries
    AVI for regulatory variants even when no filtering is applied.
    """
    from .variant_scorer import score_variant

    out: dict[str, float] = {}
    for v in variants:
        try:
            r = score_variant(
                gene=v.gene,
                position=v.position,
                wt_aa=v.wt_aa,
                mut_aa=v.mut_aa,
                chrom=v.chrom,
                ref_dna=v.ref_dna,
                alt_dna=v.alt_dna,
                protein_length=(
                    len(proteins[v.gene]) if v.gene in proteins else None
                ),
                uniprot_id=uniprot_ids.get(v.gene),
                am_lookup=am_lookup_fn,
                avi_lookup=avi_lookup_fn,
                conservation_lookup=conservation_lookup_fn,
            )
        except Exception:
            r = None
        if r is not None:
            out[f"{r.gene}.{r.position}{r.wt_aa}>{r.mut_aa}"] = r.normalized_score
    return out


def filter_variants_with_lookups(
    variants: list[Variant],
    *,
    top_fraction: float,
    min_score: float,
    proteins: dict[str, str],
    uniprot_ids: dict[str, str],
) -> tuple[list[Variant], dict[str, float], bool]:
    """Filter variants with all lookups wired, returning (kept, scores,
    am_active).

    ``kept`` is the variants list reduced to the kept set when
    filtering is enabled; otherwise it's the original list.
    ``scores`` is the per-variant normalized-score dict (empty when
    filtering is enabled AND no DNA coords present).
    ``am_active`` is True if a real AlphaMissense index was loaded.
    """
    from .variant_scorer import filter_variants

    am_lookup_fn, am_active = _try_am_lookup()

    filter_active = (
        top_fraction < 1.0 or min_score > 0.0
    )
    has_dna_coords = any(v.chrom is not None for v in variants)

    if filter_active:
        # Wire the optional lookups (mock by default; real Atlas /
        # UCSC when secrets are set + extras installed).
        avi_lookup_fn = _try_avi_lookup() if has_dna_coords else None
        conservation_lookup_fn = (
            _try_conservation_lookup() if has_dna_coords else None
        )

        scored = filter_variants(
            _build_v_dicts(variants),
            top_fraction=top_fraction,
            min_score=min_score,
            protein_lengths={g: len(p) for g, p in proteins.items()},
            protein_sequences=proteins,
            uniprot_ids=uniprot_ids,
            am_lookup=am_lookup_fn,
            avi_lookup=avi_lookup_fn,
            conservation_lookup=conservation_lookup_fn,
            strict=False,
        )
        keep_keys = {(s.gene, s.position, s.wt_aa, s.mut_aa) for s in scored}
        kept = [
            v for v in variants
            if (v.gene, v.position, v.wt_aa, v.mut_aa) in keep_keys
        ]
        scores = {
            f"{s.gene}.{s.position}{s.wt_aa}>{s.mut_aa}": s.normalized_score
            for s in scored
        }
        return kept, scores, am_active

    if has_dna_coords:
        # No filter requested but DNA coordinates present — score
        # each variant individually so the report carries per-variant
        # AVI scores for regulatory + coding variants.
        avi_lookup_fn = _try_avi_lookup()
        conservation_lookup_fn = _try_conservation_lookup()
        scores = _score_one_each(
            variants,
            proteins=proteins,
            uniprot_ids=uniprot_ids,
            am_lookup_fn=am_lookup_fn,
            avi_lookup_fn=avi_lookup_fn,
            conservation_lookup_fn=conservation_lookup_fn,
        )
        return variants, scores, am_active

    # No filter, no DNA coords → no scoring happens; report just
    # carries n_variants_input and the peptides list.
    return variants, {}, am_active


def build_pipeline_notes(
    *,
    am_active: bool,
    filter_active: bool,
    had_dna_coords: bool,
) -> str:
    """Build the user-facing ``note`` string for ``PipelineReport``.

    Documents which signals fired:
      * AlphaMissense used vs not (with download instructions)
      * AlphaGenome Atlas used (when DNA coords present)
      * tumor-marker warning (when no tumor-marker genes matched
        the expression matrix)
    """
    notes: list[str] = []
    if am_active:
        notes.append("AlphaMissense pathogenicity scores used in variant filtering.")
    elif filter_active:
        notes.append(
            "AlphaMissense predictions not found — heuristic BLOSUM62 + driver "
            "gene + structural score used. To enable AlphaMissense, download "
            "the predictions TSV from "
            "https://storage.googleapis.com/dm_alphamissense/ "
            "to ~/.cache/mrnavax/AlphaMissense_hg38.tsv "
            "(licensed CC BY-NC-SA 4.0, non-commercial)."
        )
    if had_dna_coords:
        notes.append(
            "AlphaGenome Atlas AVI scores used for non-coding regulatory "
            "variants (Avsec et al. Nature 2026). Coding-region variants "
            "are still scored by AlphaMissense (when available)."
        )
    return " | ".join(notes)


def build_tumor_marker_note(marker_idx: list[int], notes: list[str]) -> None:
    """Append the tumor-marker warning to ``notes`` if no markers found."""
    if not marker_idx:
        notes.append(
            "WARNING: no tumor-marker genes supplied; the largest cluster was "
            "assumed to be tumor. Pass --tumor-markers GENE1,GENE2 for "
            "accurate selection. The downstream peptide list may contain "
            "many false positives."
        )


__all__ = [
    "filter_variants_with_lookups",
    "build_pipeline_notes",
    "build_tumor_marker_note",
]
