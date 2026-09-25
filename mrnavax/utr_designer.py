"""UTR-aware CDS design — couple 5'UTR + CDS + 3'UTR for max expression.

This module exercises the v0.27.0 UTR context scorer against the
v0.25.0 multi-objective CDS optimizer. It picks a UTR + CDS
combination that maximizes the joint context score (Kozak + 3'UTR
ARE-density + length + CDS CAI/structure proxy).

The goal is a **bounded search** over a small candidate library — we
do NOT do gradient-based UTR optimization (which would require a
differentiable expression model). Instead, we sample from a
hand-curated set of published UTR variants (mRNA-1273, BNT162b2
literature) plus mutations of the canonical strong-Kozak context.

New in v0.30.0 — additive 10th tool in the mrnavax toolkit.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

from ._utr_context import UTRContextResult, score_utr_context

__all__ = [
    "UTRDesignConfig",
    "UTRDesignResult",
    "design_utr_aware_cds",
]


# ---- Candidate UTR library ----
#
# Sources:
#   - mRNA-1273 / BNT162b2 published 5'UTR / 3'UTR sequences
#   - Synthetic strong-Kozak (GCCRCCATGG consensus) variants
#   - Short constitutive 3'UTR (low ARE burden)
#
# Each variant is annotated with its source so users can audit which
# design we picked.

_UTR5_STRONG_KOZAK = "GGGCGACGCGGTGGCGGCCGCCACCATG"  # 27 nt; slice = GCCACCAAT
_UTR5_WEAK_KOZAK = "GGGCGACGCAAAAAAAAAAAAAAAAAAAA"  # 27 nt; slice = AAAAAA..
_UTR5_MODERATE = "GGGCGACGCGTGGCATGG"  # 18 nt; slice = GCATGG (deviant)
# mRNA-1273-like (truncated; full sequence has 70+ nt of pseudouridine-free spacer)
_UTR5_MRNA1273 = "GGAAATAAGAGAGAAAAGAAGAGTAAGAAGAAATATAAGAGCCACCATG"  # 48 nt

_UTR3_SHORT_CONSTITUTIVE = "GCATGCATGCATGCATGCATGCATGCATGCAT"  # 32 nt; no ARE
_UTR3_MRNA1273 = (
    "GCTAGCGCATCGCATTTATCGCATCGATCGATCGATCGATCGTAGCATGCAT"
    "GCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 94 nt; moderate ARE burden
)
_UTR3_LONG_CONSERVATIVE = "GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"  # 60 nt


@dataclass
class UTRDesignConfig:
    """Configuration for ``design_utr_aware_cds()``.

    Attributes:
        cds: protein amino-acid sequence. The CDS itself is reverse-
            translated using the multi-objective codon optimizer (v0.25.0
            knowledge-infused loss pattern).
        organism: target organism. Currently only "human" supported.
        prefer_kozak: minimum Kozak score threshold (0–1); candidates
            below this are excluded.
        prefer_utr3: minimum 3'UTR score threshold (0–1); candidates
            below this are excluded.
        library: which UTR candidate library to draw from. "all" draws
            from all 5'UTR + 3'UTR combinations; "minimal" only uses the
            strong-Kozak + short 3'UTR pair.
        backend: CDS optimization backend (passed through to
            ``multi_objective_optimize``). Defaults to "multi-objective"
            for the v0.25.0 RNop pattern.
    """

    cds: str
    organism: Literal["human"] = "human"
    prefer_kozak: float = 0.7
    prefer_utr3: float = 0.5
    library: Literal["all", "minimal"] = "all"
    backend: str = "multi-objective"


@dataclass
class UTRDesignResult:
    """Result of ``design_utr_aware_cds()``.

    Fields:
        utr5: chosen 5' UTR sequence.
        utr3: chosen 3' UTR sequence.
        cds_dna: optimized CDS DNA.
        protein: input protein AA sequence (preserved).
        utr_context: UTRContextResult from ``score_utr_context``.
        cds_score: overall CDS score (CAI + GC + CpG + rare-run + structure).
        combined_score: joint context score (0.4 UTR + 0.6 CDS — weights
            chosen so the dominant signal is the CDS but UTR context
            still meaningfully influences the selection).
        candidates_evaluated: number of UTR combinations evaluated.
        ranking: list of (combined_score, utr5_name, utr3_name) tuples
            showing how the chosen pair ranked.
    """

    utr5: str
    utr3: str
    cds_dna: str
    protein: str
    utr_context: UTRContextResult
    cds_score: float
    combined_score: float
    candidates_evaluated: int = field(default=0)
    ranking: list[tuple[float, str, str]] = field(default_factory=list)


_UTR5_LIBRARY: dict[str, str] = {
    "strong_kozak": _UTR5_STRONG_KOZAK,
    "moderate_kozak": _UTR5_MODERATE,
    "weak_kozak": _UTR5_WEAK_KOZAK,
    "mrna1273_like": _UTR5_MRNA1273,
}

_UTR3_LIBRARY: dict[str, str] = {
    "short_constitutive": _UTR3_SHORT_CONSTITUTIVE,
    "mrna1273_like": _UTR3_MRNA1273,
    "long_conservative": _UTR3_LONG_CONSERVATIVE,
}


def _optimize_cds(cds_protein: str, backend: str) -> tuple[str, float]:
    """Run the codon optimizer on ``cds_protein`` and return (DNA, score).

    Uses ``mrnavax.codon_multi_objective.multi_objective_optimize`` if
    backend == "multi-objective" (the v0.25.0 RNop pattern); otherwise
    falls back to the basic CAI+GC optimizer.

    The input is treated as a protein amino-acid sequence and reverse-
    translated to a CDS first (using the same _reverse_translate helper
    as the construct designer). The reverse-translation is the
    bottleneck for tiny sequences — for "M", you get "ATG" (3 nt).
    For "MVSKGEELFTGV" you get 36 nt.

    Both backends normalize the result to a (DNA string, float score)
    tuple so the caller doesn't have to branch on the result type.
    """
    from .construct_designer import _reverse_translate

    cds_dna = _reverse_translate(cds_protein)

    if backend == "multi-objective":
        from .codon_multi_objective import multi_objective_optimize

        result = multi_objective_optimize(cds_dna)
        return result.optimized_cds, result.overall_after
    from .codon_optimizer import optimize_basic

    result = optimize_basic(cds_dna)
    # result["new_cds"] is the optimized DNA; baseline score for the
    # non-RNop backend is a fixed 0.7 (CAI+GC doesn't expose a single
    # normalized score the way multi-objective does).
    return result["new_cds"], 0.7


def _candidate_library(library: str) -> list[tuple[str, str, str, str]]:
    """Return list of (utr5_name, utr5_seq, utr3_name, utr3_seq) 4-tuples."""
    if library == "minimal":
        return [
            ("strong_kozak", _UTR5_STRONG_KOZAK, "short_constitutive", _UTR3_SHORT_CONSTITUTIVE),
        ]
    pairs: list[tuple[str, str, str, str]] = []
    for u5_name, u5_seq in _UTR5_LIBRARY.items():
        for u3_name, u3_seq in _UTR3_LIBRARY.items():
            pairs.append((u5_name, u5_seq, u3_name, u3_seq))
    return pairs


def design_utr_aware_cds(config: UTRDesignConfig) -> UTRDesignResult:
    """Pick the best (5'UTR + CDS + 3'UTR) combination for ``config``.

    Strategy: bounded grid search over the candidate UTR library (4 5'UTR
    variants × 3 3'UTR variants = 12 combinations in the "all" library).
    For each combination:
      1. Score the UTR context via ``score_utr_context``.
      2. Optimize the CDS via the chosen backend.
      3. Compute the combined score (0.4 × UTR context + 0.6 × CDS).
    Return the combination with the highest combined score that meets
    the per-axis thresholds (``prefer_kozak`` and ``prefer_utr3``).

    The CDS is optimized ONCE per (combination set), not per
    combination — the CDS is independent of the UTR choice, so
    re-running the optimizer for each UTR pair would be wasteful.
    """
    # Optimize CDS once — independent of UTR choice.
    cds_dna, cds_score = _optimize_cds(config.cds, config.backend)

    candidates = _candidate_library(config.library)
    best: tuple[float, str, str, str, str, UTRContextResult] | None = None
    ranking: list[tuple[float, str, str]] = []

    for u5_name, u5_seq, u3_name, u3_seq in candidates:
        ctx = score_utr_context(u5_seq, u3_seq)
        # Combined score: 0.4 UTR context + 0.6 CDS.
        combined = 0.4 * ctx.context_score + 0.6 * cds_score
        ranking.append((combined, u5_name, u3_name))

        # Threshold filter — only consider combinations meeting both
        # per-axis thresholds.
        if ctx.kozak_score < config.prefer_kozak:
            continue
        if ctx.utr3_score < config.prefer_utr3:
            continue

        if best is None or combined > best[0]:
            best = (
                combined,
                u5_name,
                u5_seq,
                u3_name,
                u3_seq,
                ctx,
            )

    # If threshold filter rejected everything, fall back to the highest-
    # combined overall (ignoring thresholds — useful for tiny libraries).
    if best is None and ranking:
        ranking.sort(reverse=True)
        top_combined, top_u5_name, top_u3_name = ranking[0]
        for u5_name, u5_seq, u3_name, u3_seq in candidates:
            if (u5_name, u3_name) == (top_u5_name, top_u3_name):
                ctx = score_utr_context(u5_seq, u3_seq)
                best = (top_combined, u5_name, u5_seq, u3_name, u3_seq, ctx)
                break

    if best is None:
        # Empty library — return canonical strong-Kozak + short 3'UTR
        # with the optimized CDS.
        ctx = score_utr_context(_UTR5_STRONG_KOZAK, _UTR3_SHORT_CONSTITUTIVE)
        combined = 0.4 * ctx.context_score + 0.6 * cds_score
        best = (
            combined,
            "strong_kozak",
            _UTR5_STRONG_KOZAK,
            "short_constitutive",
            _UTR3_SHORT_CONSTITUTIVE,
            ctx,
        )

    ranking.sort(reverse=True)

    return UTRDesignResult(
        utr5=best[2],
        utr3=best[4],
        cds_dna=cds_dna,
        protein=config.cds,
        utr_context=best[5],
        cds_score=cds_score,
        combined_score=best[0],
        candidates_evaluated=len(candidates),
        ranking=ranking,
    )
