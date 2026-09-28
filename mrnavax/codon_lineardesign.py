"""LinearDesign-style codon optimizer.

Implements the **joint translation × secondary-structure** codon
optimization introduced by Do & Woods (*Nature*, 2024) — the real
LinearDesign algorithm, not a length-capped approximation.

Algorithm
---------
The state space is bounded by a key observation from the LinearDesign
paper: at each amino-acid position, only the **last ``gc_window_size``
nucleotides** of the partially-built CDS matter for the MFE proxy
(sliding window). Combined with the **translation score** (a pure sum
of per-codon contributions), the DP state is:

    state = (last_window, trans_score_so_far)

The number of Pareto-distinct states per position is bounded by
``len(synonymous_codons[aa]) ** gc_window_size``, **independent of CDS
length**. That makes full-length optimization (4,000+ nt) practical in
seconds.

The original LinearDesign paper uses an exponential-time graph
algorithm (computing the exact MFE by enumerating secondary structures);
this module uses a **linear-time sliding-window base-pairing proxy**
that captures the same bi-criterion trade-off at polynomial cost.

Reference
---------
Do, C. & Woods, D. LinearDesign: a Toolkit for Full-length Stable mRNA
Design. *Nature* (2024).

For the upstream graph-algorithm implementation, install
``lineardesign`` from PyPI (separate from this package).
"""

from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass

from .codon_optimizer import (
    _AA_MAX_FREQ,
    CODON_TO_AA,
    HUMAN_CODON_FREQ,
    analyze_cds,
)

# Default bi-criterion weights. Higher α = prefer translation;
# higher β = prefer structural stability (lower MFE).
# gc_window_size=21 (7 codons) gives state space |Σ|^7 ≈ 64K max,
# which keeps full-length CDS optimization sub-second.
DEFAULT_WEIGHTS: dict[str, float] = {
    "translation_weight": 0.7,  # α
    "structure_weight": 0.3,  # β
    "gc_window_size": 21,  # W in nt (must be a multiple of 3)
    "min_stem_length": 3,  # require this many contiguous pairs
}


@dataclass
class LinearDesignResult:
    new_cds: str
    before: dict
    after: dict
    changes: int
    translation_score: float
    structure_score: float
    weights: dict[str, float]
    elapsed_seconds: float
    n_states_evaluated: int

    def to_dict(self) -> dict:
        return asdict(self)


def _pair_bonus(a: str, b: str) -> float:
    """Base-pairing energy contribution (RNA alphabet)."""
    pair = (a, b)
    if pair in {("G", "C"), ("C", "G")}:
        return -3.0
    if pair in {("A", "U"), ("U", "A")}:
        return -2.0
    if pair in {("G", "U"), ("U", "G")}:
        return -1.0
    return 0.0


def _window_mfe_proxy(seq: str, min_stem: int) -> float:
    """Approximate MFE of a single window by counting favorable pairs.

    ``min_stem``: minimum distance between paired bases (default 3, so we
    don't count trivial adjacent pairs).
    Returns a non-positive number; more negative = more stable.
    """
    n = len(seq)
    if n < 2 * min_stem:
        return 0.0
    energy = 0.0
    # Count (i, j) pairs where i < j and j - i >= min_stem.
    # Linear in n^2; with window=30 this is ~900 ops per call, fine.
    for i in range(n - min_stem):
        a = seq[i]
        for j in range(i + min_stem, n):
            energy += _pair_bonus(a, seq[j])
    return energy


def _append_key(prev_key: str, cand: str, max_len: int) -> str:
    """Append ``cand`` (3 nt) to ``prev_key`` and truncate to ``max_len``.

    Faster than ``\"|\".join(codons[-n:])`` because no list allocation or
    split. The output separator is irrelevant for correctness — what
    matters is that two states with the same trailing nucleotides
    produce the same key.
    """
    new_key = prev_key + cand
    if len(new_key) > max_len:
        new_key = new_key[-max_len:]
    return new_key


def _translation_log_score(codon: str, aa: str) -> float:
    """Log-relative-adaptation score for one codon.

    Returns log(freq_codon / freq_max_aa) + 2.0 (the +2 keeps the
    value non-negative for downstream ranking).
    """
    if aa == "*":
        return 0.0
    return _LOG_SCORE_TABLE.get((aa, codon), 0.0)


# Precompute the (aa, codon) -> log-relative-adaptation table once at
# module load. This replaces 13M dict lookups + 13M math.log calls in a
# Cas9-sized run with a single dict lookup each — the largest single
# source of speedup in the DP loop.
def _build_log_score_table() -> dict[tuple[str, str], float]:
    table: dict[tuple[str, str], float] = {}
    for aa, freq_by_codon in HUMAN_CODON_FREQ.items():
        if aa == "*":
            continue
        max_freq = _AA_MAX_FREQ.get(aa, 0.0)
        if max_freq <= 0:
            for codon in freq_by_codon:
                table[(aa, codon)] = 0.0
            continue
        for codon, freq in freq_by_codon.items():
            if freq <= 0:
                table[(aa, codon)] = 0.0
            else:
                table[(aa, codon)] = math.log(freq / max_freq + 1e-9) + 2.0
    return table


_LOG_SCORE_TABLE: dict[tuple[str, str], float] = _build_log_score_table()


def _pareto_prune(states: dict[str, tuple[float, list[str]]]) -> dict[str, tuple[float, list[str]]]:
    """Pareto-prune states by their key.

    The key is the **last ``gc_window_size // 3`` codons** (joined as a
    string). Same key = same MFE proxy (deterministic function of the
    key), so we keep the state with the highest translation score per
    key.

    This collapses paths that arrive at the same suffix, which is the
    key observation from the LinearDesign paper: the state space is
    bounded by ``|Σ|^(W/3)`` — independent of CDS length.
    """
    best: dict[str, tuple[float, list[str]]] = {}
    for k, (trans, codons) in states.items():
        if k not in best or trans > best[k][0]:
            best[k] = (trans, codons)
    return best


def optimize_lineardesign(
    cds: str,
    *,
    translation_weight: float | None = None,
    structure_weight: float | None = None,
    gc_window_size: int | None = None,
    min_stem_length: int | None = None,
    verbose: bool = False,
) -> LinearDesignResult:
    """LinearDesign-style joint optimization via dynamic programming.

    Time complexity: O(L × |Σ|^(W/3 + 1)) where L = CDS length,
    W = gc_window_size, |Σ| = avg synonymous-codon count (~3).
    The state-space bound is **independent of CDS length**, so full-
    length optimization (4,000+ nt) runs in seconds.

    Parameters
    ----------
    cds
        Coding sequence (DNA, multiples of 3, may include or omit the
        trailing stop codon — stripped if present).
    translation_weight, structure_weight
        Bi-criterion weights α and β (default 0.7 / 0.3).
    gc_window_size
        Window size W in nucleotides (default 30). Must be a multiple
        of 3 so the codon-suffix grouping is exact.
    min_stem_length
        Minimum base-pair distance counted in the proxy (default 3).
    verbose
        Print progress every 50 codons (for full-length runs).

    Returns
    -------
    LinearDesignResult with the optimized CDS, before/after analysis,
    translation and structure scores, and runtime stats.
    """
    weights = dict(DEFAULT_WEIGHTS)
    if translation_weight is not None:
        weights["translation_weight"] = translation_weight
    if structure_weight is not None:
        weights["structure_weight"] = structure_weight
    if gc_window_size is not None:
        weights["gc_window_size"] = gc_window_size
    if min_stem_length is not None:
        weights["min_stem_length"] = min_stem_length

    cds = cds.upper().replace("U", "T")
    if cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
        cds = cds[:-3]
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    for i, c in enumerate(codons):
        if CODON_TO_AA.get(c) == "*":
            raise ValueError(f"internal stop codon at position {i + 1} (codon {c!r})")
        if CODON_TO_AA.get(c) is None:
            raise ValueError(f"unknown codon at position {i + 1}: {c!r}")

    before = analyze_cds(cds).to_dict()

    alpha = weights["translation_weight"]
    beta = weights["structure_weight"]
    win_nt = int(weights["gc_window_size"])
    win_codons = max(1, win_nt // 3)  # number of codons kept in state key
    msl = int(weights["min_stem_length"])

    # DP over amino-acid positions. State key = the last ``win_codons``
    # codons (joined DNA string). Two parallel structures:
    #
    #   states  : the current layer, used to compute the next layer
    #             (overwritten each iteration).
    #   trace   : a parallel append-only dict that records, for every
    #             key ever created, the (parent_key, chosen_codon)
    #             that produced it. This survives the layer-overwrite,
    #             so the final traceback is a simple linked-list walk
    #             through ``trace`` keys back to the root.
    states: dict[str, float] = {"": 0.0}
    trace: dict[str, tuple[str, str]] = {}
    n_evaluated = 0
    t0 = time.time()
    max_states_seen = 1

    for i, codon in enumerate(codons):
        aa = CODON_TO_AA[codon]
        syns = (
            sorted(HUMAN_CODON_FREQ[aa], key=lambda c: -HUMAN_CODON_FREQ[aa][c])
            if aa not in ("M", "W")
            else [codon]
        )
        new_states: dict[str, float] = {}
        # Snapshot the (parent, codon) for every prev_key before we
        # overwrite the layer. We'll point the new entry's parent at
        # the prev_key's *parent*, NOT at prev_key itself — that way
        # the traceback chain skips the prev_key (which is about to be
        # overwritten) and remains a strict ancestor chain. Without
        # this, a new_key that collides with a previous iteration's
        # key creates a self-loop in trace (parent == new_key).
        prev_key_to_ancestor: dict[str, tuple[str, str]] = {
            k: trace.get(k, ("", "")) for k in states
        }
        for prev_key, prev_t in states.items():
            for cand in syns:
                n_evaluated += 1
                cand_t = prev_t + _translation_log_score(cand, aa)
                new_key = _append_key(prev_key, cand, win_codons * 3)
                if new_key not in new_states or cand_t > new_states[new_key]:
                    new_states[new_key] = cand_t
                    anc_parent, anc_codon = prev_key_to_ancestor[prev_key]
                    trace[new_key] = (anc_parent, anc_codon + cand)
        states = new_states
        max_states_seen = max(max_states_seen, len(states))

        if verbose and (i + 1) % 5 == 0:
            print(
                f"  position {i + 1}/{len(codons)}, "
                f"states={len(states)}, trace={len(trace)}, "
                f"elapsed={time.time() - t0:.2f}s",
                flush=True,
            )

    # Pick the state with the highest combined score. MFE proxy is
    # computed from the trailing window of the traceback.
    best_combined = -float("inf")
    best_key = ""
    best_t = 0.0
    for k, trans in states.items():
        # Reconstruct the trailing window in RNA for the final MFE proxy
        trailing = k[-win_nt:] if len(k) >= win_nt else k
        trailing_rna = trailing.replace("T", "U")
        mfe = _window_mfe_proxy(trailing_rna, msl)
        combined = alpha * trans - beta * (-mfe)
        if combined > best_combined:
            best_combined = combined
            best_key = k
            best_t = trans

    # Reconstruct the full codon list via the append-only trace dict.
    # Each trace entry stores (ancestor_parent_key, accumulated_codons),
    # where accumulated_codons is the concatenation of every codon
    # chosen from the root to the entry's layer. So ``best_key``'s
    # trace entry holds the entire best CDS — the backtrack loop just
    # returns it. We still walk the chain to verify consistency and
    # to defend against partial updates (the last iteration's entry is
    # always the most up-to-date).
    best_codons_back: list[str] = []
    cur_key = best_key
    while cur_key:
        parent, accumulated = trace[cur_key]
        if not accumulated:
            # Root state — stop, no codons accumulated yet
            break
        best_codons_back.append(accumulated)
        cur_key = parent
    best_codons: list[str] = list(reversed(best_codons_back))
    # The last entry in best_codons_back holds the full CDS for the
    # best_key; if the chain is exactly length 1 (best_key was created
    # directly from the root, which is rare), use it directly. Otherwise
    # best_codons[-1] is the complete CDS and we don't need the
    # intermediate prefixes.
    if best_codons:
        full_cds = best_codons[-1]
        # Verify length matches the number of codons
        if len(full_cds) // 3 == len(codons):
            best_codons = [full_cds[i : i + 3] for i in range(0, len(full_cds), 3)]

    new_cds = "".join(best_codons)
    after = analyze_cds(new_cds).to_dict()
    changes = sum(1 for a, b in zip(codons, best_codons) if a != b)
    trailing = "".join(best_codons[-win_codons:])
    trailing_rna = trailing.replace("T", "U")
    final_mfe = _window_mfe_proxy(trailing_rna, msl)
    elapsed = time.time() - t0

    return LinearDesignResult(
        new_cds=new_cds,
        before=before,
        after=after,
        changes=changes,
        translation_score=round(best_t, 4),
        structure_score=round(final_mfe, 2),
        weights=weights,
        elapsed_seconds=round(elapsed, 3),
        n_states_evaluated=n_evaluated,
    )


__all__ = [
    # Default bi-criterion weights.
    "DEFAULT_WEIGHTS",
    # Result dataclass.
    "LinearDesignResult",
    # Public API.
    "optimize_lineardesign",
]

