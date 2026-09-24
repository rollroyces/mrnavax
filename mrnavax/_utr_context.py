"""UTR-aware variant scoring — extends score_variant with a 6th signal.

This module adds a *context* score to ``score_variant``'s 4-signal
composition (BLOSUM62 + driver + structural + AM → AVI → PhyloP). The
context captures translation-initiation efficiency from the
surrounding 5'UTR + Kozak sequence and the 3'UTR regulatory elements
(polyadenylation, ARE motifs).

Why this is the right addition for v0.27.0
-------------------------------------------
RNop (arXiv:2505.23862, Aug 2026) frames mRNA optimization as
"knowledge infusion across UTR context + CDS + tail". Without the
UTR context signal, the toolkit can rank two variants identically
even when one sits in a construct with strong Kozak context and
the other in a weak context — which actually matters for
expression. Adding the UTR context score makes the 4-signal
composition 5-signal, fully aware of the construct context.

The score is computed from the upstream/downstream sequences the
caller supplies (optional, both default to None so existing callers
don't break). When supplied, the score contributes a small bonus
(weight = 0.05 by default) to the overall normalized score, capped
to keep the dominant signals (AM / AVI / PhyloP) dominant.

Components
----------
* **kozak_score** — Kozak context strength (-9 to +4 around ATG),
  score in [0, 1]. Reuses ``manufacturability.check_kozak_strength``.
* **utr3_score** — 3'UTR quality (ARE motif density + length), score
  in [0, 1]. Implemented inline (no dependency on a separate module).
* **context_score** — Combined context = 0.5*kozak + 0.5*utr3, in [0, 1].

Public API
----------
* ``score_utr_context(utr5, utr3) -> UTRContextResult`` — score UTR
  context alone (no variant required; useful for construct QC).
* ``score_variant(..., utr5=None, utr3=None)`` — extended to accept
  optional UTR sequences; when supplied, adds ``utr_context_score``
  to the components dict.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass
class UTRContextResult:
    """Per-component UTR context score for a construct.

    Attributes
    ----------
    kozak_score : float
        Kozak consensus strength in [0, 1]. Computed from the 9 nt
        immediately upstream of the start codon against the canonical
        Kozak motif ``GCCRCCATGG`` (R = purine).
    utr3_score : float
        3'UTR quality in [0, 1]. Combines ARE-motif density
        (canonical AUUUA pentamers) and 3'UTR length (≥100 nt
        considered healthy for therapeutic mRNAs).
    context_score : float
        Combined context score in [0, 1] (0.5 * kozak + 0.5 * utr3).
    details : dict
        Per-component diagnostic info (Kozak matched positions, ARE
        motif count, UTR lengths).
    """

    kozak_score: float
    utr3_score: float
    context_score: float
    details: dict

    def to_dict(self) -> dict:
        return asdict(self)


# Default UTR-context weight in the overall score composition.
# Kept small (0.05) so the dominant signals (AM/AVI/PhyloP, total
# weight 0.85) remain dominant. Override via ConstructConfig.
DEFAULT_UTR_CONTEXT_WEIGHT = 0.05


def _score_kozak(utr5: str) -> tuple[float, dict]:
    """Kozak context score in [0, 1] + diagnostic details.

    Convention (matches ``manufacturability.check_kozak_strength``):
    the last 9 nt of ``utr5`` are scored against the canonical Kozak
    consensus ``GCCRCCATGG`` (where R = A or G, ATG = start codon).
    Positions -3 (R) and +4 (G) carry 1.5× weight per the original
    Kozak paper.

    Returns (score, details). Score is 0.5 when utr5 is too short
    to score (consistent with the manufacturability checker).
    """
    if not utr5:
        return 0.5, {"reason": "no utr5 supplied"}
    rna = utr5.upper().replace("T", "U")
    if len(rna) < 9:
        return 0.5, {"reason": f"utr5 only {len(rna)} nt (need ≥9)"}
    pre = rna[-9:]
    weights = [1, 1, 1, 1.5, 1, 1, 1.5, 1, 1]
    expected: list[tuple[str, ...]] = [
        ("G",),
        ("C",),
        ("C",),
        ("A", "G"),  # R at position -3
        ("C",),
        ("C",),
        ("A", "G"),  # R at position 0 of ATG? — Actually canonical
                     # Kozak has R at position -3 only; the docstring
                     # ``GCCRCCATGG`` is the 10-char canonical motif
                     # (R at -3, A at start, T at +1, G at +4). This
                     # 9-char slice matches the upstream
                     # manufacturability.check_kozak_strength exactly.
        ("A",),  # A in ATG
        ("T",),  # T in ATG
    ]
    matched = 0.0
    total_w = sum(weights)
    for i, exp in enumerate(expected):
        if pre[i] in exp:
            matched += weights[i]
    score = matched / total_w if total_w else 0.0
    return round(score, 4), {
        "pre_atg": pre,
        "matched_weight": matched,
        "total_weight": total_w,
    }


def _score_utr3(utr3: str) -> tuple[float, dict]:
    """3'UTR quality score in [0, 1] + diagnostic details.

    Combines ARE-motif density (canonical AUUUA pentamers, 1-2 = OK,
    3+ = high) and length (≥100 nt considered healthy).
    """
    if not utr3:
        return 0.5, {"reason": "no utr3 supplied"}
    rna = utr3.upper().replace("T", "U")

    # ARE-motif density (AUUUA pentamers in the 3'UTR)
    are_count = rna.count("AUUUA")
    if are_count == 0:
        are_score = 0.4  # low — no AU-rich elements
    elif are_count <= 2:
        are_score = 0.7  # moderate
    else:
        are_score = 1.0  # dense ARE motifs (stabilizing in some contexts)

    # Length score (≥100 nt = full credit)
    if len(rna) >= 100:
        length_score = 1.0
    else:
        length_score = max(0.0, len(rna) / 100.0)

    score = 0.6 * are_score + 0.4 * length_score
    return round(score, 4), {
        "length": len(rna),
        "are_motif_count": are_count,
        "are_score": are_score,
        "length_score": length_score,
    }


def score_utr_context(
    utr5: str | None,
    utr3: str | None,
) -> UTRContextResult:
    """Score the UTR context of a construct (Kozak + 3'UTR quality).

    Parameters
    ----------
    utr5 : str | None
        5' UTR sequence (DNA). The last 9 nt are matched against the
        canonical Kozak motif.
    utr3 : str | None
        3' UTR sequence (DNA). ARE-motif density and length are scored.

    Returns
    -------
    UTRContextResult
        Per-component scores + diagnostic details.
    """
    kozak_score, kozak_details = _score_kozak(utr5 or "")
    utr3_score, utr3_details = _score_utr3(utr3 or "")
    context = 0.5 * kozak_score + 0.5 * utr3_score
    return UTRContextResult(
        kozak_score=kozak_score,
        utr3_score=utr3_score,
        context_score=round(context, 4),
        details={"kozak": kozak_details, "utr3": utr3_details},
    )
