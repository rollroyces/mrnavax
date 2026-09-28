"""RiboDecode-style codon backend.

A data-driven codon optimizer that improves on the classic per-codon greedy
swap (``optimize_basic``) by modeling codon-pair interactions and ribosome
stalling. Falls back gracefully when no ribosome-profiling data is provided.

The original RiboDecode (Fu et al., *Nat Commun* 2025) trains a deep model on
Ribo-seq data. This module captures the same operational idea — *context-aware,
data-driven codon choice* — without requiring GPU inference:

1. **Codon-pair penalty.** Adjacent codon pairs that are both rare in human
   usage stall ribosomes. Penalize (rare, rare) neighbors.
2. **Rare-run penalty.** Consecutive runs of rare codons (>3 in a row) are
   strongly disfavored.
3. **Local GC smoothness.** Ribosomes stall at sharp GC transitions; the
   optimizer softens the GC% window profile.
4. **Optional empirical weights.** If the caller passes a ``ribo_weights``
   dict mapping ``{codon: relative_translation_rate}``, those override the
   CAI table. This is the integration point for real Ribo-seq data.

References
----------
Fu et al., *Deep generative optimization of mRNA codon sequences for enhanced
mRNA translation and therapeutic efficacy*, Nat Commun 16, 9957 (2025).
"""

from __future__ import annotations

import math
import random  # noqa: F401  — kept for downstream code that may import from this module
from dataclasses import dataclass, field

from .codon_optimizer import (
    _AA_MAX_FREQ,
    CODON_TO_AA,
    HUMAN_CODON_FREQ,
    _gc_window_stddev,
    analyze_cds,
)

# Default penalty weights (tuned heuristically — the values don't need to be
# perfect, they just need to bias the optimizer away from pathological
# patterns the basic greedy swap allows).
DEFAULT_WEIGHTS: dict[str, float] = {
    "rare_run_threshold": 0.30,  # codons below this usage count as "rare"
    "rare_run_max_length": 3,  # > this many rare codons in a row = penalty
    "rare_run_penalty": 0.15,  # additive penalty per codon over the threshold
    "rare_pair_penalty": 0.05,  # penalty for each (rare, rare) adjacent pair
    "gc_window_stddev_max": 6.0,  # penalty kicks in above this stddev
    "gc_window_stddev_penalty": 0.02,  # per unit of stddev over the threshold
}


@dataclass
class OptimizationResult:
    new_cds: str
    before: dict
    after: dict
    changes: int
    rare_run_count: int
    rare_pair_count: int
    gc_window_stddev_before: float
    gc_window_stddev_after: float
    weights: dict[str, float] = field(default_factory=dict)
    note: str = ""

    def to_dict(self) -> dict:
        return {
            "new_cds": self.new_cds,
            "before": self.before,
            "after": self.after,
            "changes": self.changes,
            "rare_run_count": self.rare_run_count,
            "rare_pair_count": self.rare_pair_count,
            "gc_window_stddev_before": self.gc_window_stddev_before,
            "gc_window_stddev_after": self.gc_window_stddev_after,
            "weights": self.weights,
            "note": self.note,
        }


def _is_rare(codon: str, threshold: float, weights: dict[str, float]) -> bool:
    aa = CODON_TO_AA.get(codon)
    if aa is None or aa == "*":
        return False
    # empirical weights override table
    ribo_rate = weights.get("ribo_rates", {}).get(codon)
    if ribo_rate is not None:
        return ribo_rate < threshold
    return HUMAN_CODON_FREQ[aa].get(codon, 0.0) < threshold


def _codon_score(codon: str, weights: dict[str, float]) -> float:
    """Higher score = better codon."""
    aa = CODON_TO_AA.get(codon)
    if aa is None or aa == "*":
        return 0.0
    ribo_rate = weights.get("ribo_rates", {}).get(codon)
    if ribo_rate is not None:
        return float(ribo_rate)
    freq = HUMAN_CODON_FREQ[aa].get(codon, 0.0)
    max_freq = _AA_MAX_FREQ[aa]
    if max_freq <= 0 or freq <= 0:
        return 0.0
    # log-scaled relative frequency — more stable than raw CAI
    return math.log(freq / max_freq + 1e-9) + 2.0  # shift to positive


def _count_rare_runs(codons: list[str], threshold: float, weights: dict[str, float]) -> int:
    """Total number of codons in rare-runs (consecutive rare codons)."""
    in_run = 0
    total = 0
    max_run = int(weights.get("rare_run_max_length", DEFAULT_WEIGHTS["rare_run_max_length"]))
    for c in codons:
        if _is_rare(c, threshold, weights):
            in_run += 1
        else:
            if in_run > max_run:
                total += in_run
            in_run = 0
    if in_run > max_run:
        total += in_run
    return total


def _count_rare_pairs(codons: list[str], threshold: float, weights: dict[str, float]) -> int:
    """Count of adjacent (rare, rare) codon pairs."""
    n = 0
    for i in range(len(codons) - 1):
        if _is_rare(codons[i], threshold, weights) and _is_rare(codons[i + 1], threshold, weights):
            n += 1
    return n


def _optimize_with_score(
    cds: str,
    weights: dict[str, float],
    *,
    n_passes: int = 3,
) -> str:
    """Hill-climb codon optimization under a custom scoring function.

    The scorer is the sum of:
      - per-codon translation score (CAI / ribo-rate)
      - rare-run penalty (per extra rare codon in a long run)
      - rare-pair penalty (per adjacent rare pair)
      - GC window stddev penalty (smooths out GC spikes)

    Deterministic: ties are broken by preferring the codon that appears
    first in the codon-usage table. For a stochastic variant, swap in a
    different optimizer — this function is the deterministic baseline.
    """
    cds = cds.upper().replace("U", "T")
    if cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
        cds = cds[:-3]
    # Reject internal stop codons — they'd truncate the protein. Caller's job
    # to clean up, but we surface a clear error here.
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    for i, c in enumerate(codons):
        aa = CODON_TO_AA.get(c)
        if aa == "*":
            raise ValueError(
                f"internal stop codon at position {i + 1} (codon {c!r}); "
                "refusing to optimize. Remove the stop or split the CDS."
            )
    threshold = weights.get("rare_run_threshold", DEFAULT_WEIGHTS["rare_run_threshold"])
    rare_run_pen = weights.get("rare_run_penalty", DEFAULT_WEIGHTS["rare_run_penalty"])
    rare_pair_pen = weights.get("rare_pair_penalty", DEFAULT_WEIGHTS["rare_pair_penalty"])
    gc_std_max = weights.get("gc_window_stddev_max", DEFAULT_WEIGHTS["gc_window_stddev_max"])
    gc_std_pen = weights.get(
        "gc_window_stddev_penalty", DEFAULT_WEIGHTS["gc_window_stddev_penalty"]
    )

    def score(codons: list[str]) -> float:
        # sum of per-codon scores
        s = sum(_codon_score(c, weights) for c in codons)
        # rare-run penalty
        s -= rare_run_pen * _count_rare_runs(codons, threshold, weights)
        # rare-pair penalty
        s -= rare_pair_pen * _count_rare_pairs(codons, threshold, weights)
        # GC window smoothness
        gc_std = _gc_window_stddev("".join(codons))
        if gc_std > gc_std_max:
            s -= gc_std_pen * (gc_std - gc_std_max)
        return s

    best_codons = list(codons)
    best_score = score(best_codons)

    for _pass in range(n_passes):
        improved = False
        for i in range(len(best_codons)):
            codon = best_codons[i]
            aa = CODON_TO_AA.get(codon)
            if aa is None or aa == "*" or aa == "M" or aa == "W":
                continue
            # Try every synonymous codon and keep the best
            current_local = best_codons[:]
            best_local = current_local[:]
            best_local_score = score(best_local)
            for cand in HUMAN_CODON_FREQ[aa]:
                if cand == codon:
                    continue
                candidate = best_codons[:]
                candidate[i] = cand
                cs = score(candidate)
                if cs > best_local_score + 1e-9:
                    best_local = candidate
                    best_local_score = cs
            if best_local_score > best_score + 1e-9:
                best_codons = best_local
                best_score = best_local_score
                improved = True
        if not improved:
            break

    return "".join(best_codons)


def optimize_ribodecode(
    cds: str,
    *,
    ribo_weights: dict[str, float] | None = None,
    custom_weights: dict[str, float] | None = None,
) -> OptimizationResult:
    """RiboDecode-style codon optimization with context awareness.

    Deterministic hill-climb. Parameters
    ----------
    cds : str
        Coding sequence (uppercase, ACGT). May include U (converted to T).
        Must not contain internal stop codons.
    ribo_weights : dict, optional
        Per-codon relative translation rates (e.g., from a Ribo-seq experiment).
        Keys are codons (uppercase, T-form), values are positive floats. Higher
        is better. If omitted, falls back to the human codon-usage table.
    custom_weights : dict, optional
        Override the penalty weights. Defaults to ``DEFAULT_WEIGHTS``.
    """
    weights = dict(DEFAULT_WEIGHTS)
    if custom_weights:
        weights.update(custom_weights)
    if ribo_weights is not None:
        weights["ribo_rates"] = dict(ribo_weights)

    before = analyze_cds(cds).to_dict()

    # Strip trailing stop if present
    working = cds.upper().replace("U", "T")
    if working[-3:] in CODON_TO_AA and CODON_TO_AA[working[-3:]] == "*":
        working = working[:-3]

    new_cds = _optimize_with_score(working, weights)
    after = analyze_cds(new_cds).to_dict()

    codons_before = [working[i : i + 3] for i in range(0, len(working), 3)]
    codons_after = [new_cds[i : i + 3] for i in range(0, len(new_cds), 3)]
    threshold = weights["rare_run_threshold"]
    rare_runs_after = _count_rare_runs(codons_after, threshold, weights)
    rare_pairs_after = _count_rare_pairs(codons_after, threshold, weights)
    gc_std_before = _gc_window_stddev(working)
    gc_std_after = _gc_window_stddev(new_cds)

    changes = sum(1 for a, b in zip(codons_before, codons_after) if a != b)

    note = ""
    if ribo_weights is not None:
        note = "Used provided ribosome-profiling weights; treat as data-driven optimization."
    else:
        note = "Used human codon-usage table (no Ribo-seq data provided). Pass ribo_weights={} for data-driven mode."

    return OptimizationResult(
        new_cds=new_cds,
        before=before,
        after=after,
        changes=changes,
        rare_run_count=rare_runs_after,
        rare_pair_count=rare_pairs_after,
        gc_window_stddev_before=round(gc_std_before, 2),
        gc_window_stddev_after=round(gc_std_after, 2),
        weights=weights,
        note=note,
    )


__all__ = [
    # Default penalty weights.
    "DEFAULT_WEIGHTS",
    # Result dataclass.
    "OptimizationResult",
    # Public API.
    "optimize_ribodecode",
]

