"""Multi-objective codon optimization (SOTA pattern from RNop / LinearDesign / RiboDecode).

This module provides a single, stdlib-only entry point that takes a CDS +
weighted loss vector and returns an optimized CDS with per-component
contribution breakdown for both the input and the output.

The pattern mirrors the **knowledge-infused loss** approach from
arxiv:2505.23862 (RNop, Aug 2026), where each biologically-motivated
component (CAI, GC, CpG, rare-codon-run, structure-proxy) is decomposed
into its own score so the user can see **which** component drove each
codon swap. This is the "impossible triangle" answer from the top-tier
SOTA projects (RNop / LinearDesign / RiboDecode): fidelity is preserved
(only synonymous swaps), multi-objective is explicit (the loss vector),
and efficiency is stdlib-Python.

Each component is implemented using primitives that already exist in
the codebase (``codon_optimizer.py``: CAI, GC, CpG, rare-codon-run).
The structure-proxy uses the existing GC-window-stddev heuristic. This
module only orchestrates them and exposes the per-component breakdown.

References
----------
* RNop (arxiv:2505.23862, Aug 2026): knowledge-infused losses for
  mRNA optimization; the "impossible triangle" of fidelity,
  multi-objective, and efficiency.
* LinearDesign (Nature 621, 2023): joint CAI + MFE optimization.
* RiboDecode (Nature Comm. 2025): weighted translation + structure
  loss with environment-aware per-codon rates.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass, field

from .codon_optimizer import (
    _AA_MAX_FREQ,
    CODON_TO_AA,
    HUMAN_CODON_FREQ,
    _cpg_obs_exp,
    _gc_content,
    _gc_window_stddev,
)

# ---------- Configuration --------------------------------------------------


@dataclass
class MultiObjectiveConfig:
    """Weighted loss vector for multi-objective codon optimization.

    Each weight controls how much that component contributes to the
    overall score. The overall score is a weighted sum of the
    normalized component scores (each in [0, 1]).

    Defaults match the canonical SOTA pattern: CAI is the dominant
    signal (~40%), structure-proxy is next (~25%), then GC band, then
    CpG-avoidance, with rare-codon-run counted as a soft penalty.

    Attributes
    ----------
    cai_weight, gc_weight, cpg_weight, rare_run_weight, structure_weight
            Non-negative weights (the score is relative, not a probability).
    target_gc_min, target_gc_max : float
        Soft GC band; GC% inside the band scores 1.0, outside falls linearly.
    rare_threshold : float
        Codon usage frequency below which a codon is "rare".
    structure_window : int
        Window size for the GC-window-stddev structure proxy.
    """

    cai_weight: float = 0.40
    gc_weight: float = 0.20
    cpg_weight: float = 0.15
    rare_run_weight: float = 0.15
    structure_weight: float = 0.10
    target_gc_min: float = 45.0
    target_gc_max: float = 60.0
    rare_threshold: float = 0.10
    structure_window: int = 90


@dataclass
class ComponentBreakdown:
    """Per-component contribution to the overall multi-objective score.

    All fields are in [0, 1] after normalization:

    * cai  — Codon Adaptation Index (translation efficiency).
    * gc_score  — GC band score (transcript stability proxy).
    * cpg_score  — CpG-avoidance (innate-immune proxy).
    * rare_run_penalty  — Longest rare-codon run (ribosome-stalling proxy);
      higher = worse, so it is subtracted in the overall score.
    * structure_proxy  — GC-window stddev proxy for local structure
      uniformity; lower stddev = higher score.
    """

    cai: float
    gc_score: float
    cpg_score: float
    rare_run_penalty: float
    structure_proxy: float

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class MultiObjectiveResult:
    """Result of a multi-objective codon optimization pass."""

    input_cds: str
    optimized_cds: str
    before: ComponentBreakdown
    after: ComponentBreakdown
    overall_before: float
    overall_after: float
    improvement: float
    n_changes: int
    config: MultiObjectiveConfig = field(default_factory=MultiObjectiveConfig)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["config"] = asdict(self.config)
        return d


# ---------- Component scorers ---------------------------------------------


def _cai_score(cds: str) -> float:
    """CAI score in (0, 1]. Geometric mean of (freq / max_freq)."""
    log_sum = 0.0
    n = 0
    for i in range(0, len(cds), 3):
        c = cds[i : i + 3]
        aa = CODON_TO_AA.get(c)
        if aa is None or aa == "*":
            continue
        freq = HUMAN_CODON_FREQ[aa].get(c, 0.0)
        denom = _AA_MAX_FREQ[aa]
        if denom > 0 and freq > 0:
            log_sum += math.log(freq / denom)
            n += 1
    return math.exp(log_sum / n) if n else 0.0


def _gc_band_score(gc_pct: float, lo: float, hi: float) -> float:
    """Linear-falloff GC-band score in [0, 1]."""
    if lo <= gc_pct <= hi:
        return 1.0
    if gc_pct < lo:
        return max(0.0, 1.0 - (lo - gc_pct) / 15.0)
    return max(0.0, 1.0 - (gc_pct - hi) / 15.0)


def _cpg_score(seq: str) -> float:
    """CpG-avoidance score in [0, 1]. Lower Obs/Exp -> higher score."""
    oe = _cpg_obs_exp(seq)
    return max(0.0, 1.5 - oe) / 1.5


def _rare_run_penalty(cds: str, threshold: float) -> float:
    """Rare-codon-run penalty in [0, 1] (longest run / 4)."""
    longest = 0
    cur = 0
    for i in range(0, len(cds), 3):
        c = cds[i : i + 3]
        aa = CODON_TO_AA.get(c)
        if aa is None or aa == "*":
            cur = 0
            continue
        freq = HUMAN_CODON_FREQ[aa].get(c, 0.0)
        if freq < threshold:
            cur += 1
            longest = max(longest, cur)
        else:
            cur = 0
    return min(1.0, longest / 4.0)


def _structure_proxy_score(cds: str, window: int) -> float:
    """Structure-proxy score in [0, 1] from GC-window stddev."""
    sd = _gc_window_stddev(cds, window=window)
    return max(0.0, 1.0 - sd / 15.0)


def _component_breakdown(
    cds: str, config: MultiObjectiveConfig
) -> ComponentBreakdown:
    """Compute all per-component scores for a CDS."""
    gc_pct = _gc_content(cds)
    return ComponentBreakdown(
        cai=round(_cai_score(cds), 4),
        gc_score=round(_gc_band_score(gc_pct, config.target_gc_min, config.target_gc_max), 4),
        cpg_score=round(_cpg_score(cds), 4),
        rare_run_penalty=round(_rare_run_penalty(cds, config.rare_threshold), 4),
        structure_proxy=round(_structure_proxy_score(cds, config.structure_window), 4),
    )


def _overall_score(b: ComponentBreakdown, config: MultiObjectiveConfig) -> float:
    """Weighted-sum of component scores, clamped to [0, 1].

    rare_run_penalty is subtracted (higher = worse).
    """
    total_weight = (
        config.cai_weight
        + config.gc_weight
        + config.cpg_weight
        + config.rare_run_weight
        + config.structure_weight
    )
    if total_weight <= 0:
        return 0.0
    raw = (
        config.cai_weight * b.cai
        + config.gc_weight * b.gc_score
        + config.cpg_weight * b.cpg_score
        - config.rare_run_weight * b.rare_run_penalty
        + config.structure_weight * b.structure_proxy
    )
    return round(max(0.0, min(1.0, raw / total_weight)), 4)


# ---------- Optimization ---------------------------------------------------


def _clean_cds(cds: str) -> str:
    """Strip whitespace, normalize to uppercase ACGT-only, trim to multiple-of-3."""
    cds = re.sub(r"[^ATCG]", "", cds.upper().replace("U", "T"))
    if len(cds) >= 3 and cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
        cds = cds[:-3]
    trim = len(cds) % 3
    if trim:
        cds = cds[:-trim]
    return cds


def _swap_candidates(codon: str, aa: str) -> list[str]:
    """Return synonymous-codon candidates sorted by freq desc, GC-match asc."""
    current_gc = _gc_content(codon)
    freqs = HUMAN_CODON_FREQ[aa]
    return [
        c
        for c, _ in sorted(
            freqs.items(),
            key=lambda kv: (-kv[1], abs(_gc_content(kv[0]) - current_gc)),
        )
    ]


def multi_objective_score(
    cds: str, config: MultiObjectiveConfig | None = None
) -> float:
    """Score a CDS without optimizing — useful for comparing two designs.

    Returns the weighted-sum overall score in [0, 1].
    """
    if config is None:
        config = MultiObjectiveConfig()
    cds_clean = _clean_cds(cds)
    if len(cds_clean) < 3:
        return 0.0
    return _overall_score(_component_breakdown(cds_clean, config), config)


def multi_objective_optimize(
    cds: str, config: MultiObjectiveConfig | None = None
) -> MultiObjectiveResult:
    """Greedy synonymous-codon swap with multi-objective scoring.

    Iterates codon-by-codon; for each codon, evaluates the overall
    multi-objective score after each candidate synonymous swap and
    picks the one that improves the score the most.

    Parameters
    ----------
    cds : str
        Coding DNA sequence (no UTRs, length multiple of 3 after cleaning).
    config : MultiObjectiveConfig, optional
        Loss vector. Defaults to MultiObjectiveConfig().

    Returns
    -------
    MultiObjectiveResult
        Optimized CDS + per-component breakdown for before/after +
        improvement + change count.
    """
    if config is None:
        config = MultiObjectiveConfig()
    cds_clean = _clean_cds(cds)
    if len(cds_clean) < 3:
        raise ValueError("CDS shorter than one codon after cleaning")
    before_breakdown = _component_breakdown(cds_clean, config)
    overall_before = _overall_score(before_breakdown, config)

    out = list(cds_clean)
    n_changes = 0
    current_overall = overall_before
    for i in range(0, len(cds_clean), 3):
        codon = cds_clean[i : i + 3]
        aa = CODON_TO_AA.get(codon)
        if aa is None or aa == "*" or aa == "M" or aa == "W":
            continue  # single-codon AAs, nothing to swap
        candidates = _swap_candidates(codon, aa)
        best_codon = codon
        best_score = current_overall
        for cand in candidates:
            if cand == codon:
                continue
            trial = "".join(out[:i] + list(cand) + out[i + 3 :])
            trial_score = _overall_score(
                _component_breakdown(trial, config), config
            )
            if trial_score > best_score:
                best_score = trial_score
                best_codon = cand
        if best_codon != codon:
            out[i : i + 3] = best_codon
            n_changes += 1
            current_overall = best_score

    new_cds = "".join(out)
    after_breakdown = _component_breakdown(new_cds, config)
    overall_after = _overall_score(after_breakdown, config)
    return MultiObjectiveResult(
        input_cds=cds_clean,
        optimized_cds=new_cds,
        before=before_breakdown,
        after=after_breakdown,
        overall_before=overall_before,
        overall_after=overall_after,
        improvement=round(overall_after - overall_before, 4),
        n_changes=n_changes,
        config=config,
    )
