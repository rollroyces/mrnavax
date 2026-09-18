"""Internal upstream-model lookup helpers for ``variant_scorer.score_variant``.

Extracted from ``variant_scorer.py`` in v0.22.0 to keep the public
``score_variant`` function focused on orchestration rather than
detail. Private API; downstream code should call
``mrnavax.variant_scorer.score_variant`` directly.

Three upstream-model integrations live here, each with the same
"silent failure" contract:

  * Network error / exception → component dropped, weight = 0
  * Returns None → component dropped, weight = 0
  * Returns invalid range → component dropped (caught by the
    individual validators)

Each function returns a small dataclass with the score and any
classification metadata, plus a weight that ``score_variant`` uses
to decide which signal dominates the weighted combination.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

# ---------------------------------------------------------------------------
# AlphaMissense
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AlphaMissenseComponent:
    """Result of an AlphaMissense lookup, normalized for scoring.

    score : float
        Pathogenicity probability in [0, 1].
    classification : str
        "likely_benign" / "ambiguous" / "likely_pathogenic".
    weight : float
        0.45 when AM is available, 0 otherwise.
    """

    score: float
    classification: str
    weight: float


def compute_am_component(
    uniprot_id: str | None,
    wt_aa: str,
    position: int,
    mut_aa: str,
    am_lookup: Callable | None,
) -> AlphaMissenseComponent:
    """Resolve AlphaMissense pathogenicity for a coding-region variant.

    Silent failure: returns weight=0 if either ``uniprot_id`` or
    ``am_lookup`` is missing, or if the lookup raises / returns None.
    """
    if uniprot_id is None or am_lookup is None:
        return AlphaMissenseComponent(score=0.0, classification="", weight=0.0)
    try:
        result = am_lookup(uniprot_id, wt_aa, position, mut_aa)
    except Exception:
        return AlphaMissenseComponent(score=0.0, classification="", weight=0.0)
    if result is None:
        return AlphaMissenseComponent(score=0.0, classification="", weight=0.0)
    # AlphaMissense is the highest-fidelity signal here when present.
    return AlphaMissenseComponent(
        score=result.score,
        classification=result.classification,
        weight=0.45,
    )


# ---------------------------------------------------------------------------
# AlphaGenome Atlas AVI
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AVIComponent:
    """Result of an AlphaGenome Atlas lookup, normalized for scoring.

    score : float
        AlphaGenome Variant Impact score in [0, 1].
    classification : str
        "low" / "moderate" / "high".
    is_coding : bool
        True if Atlas reports the variant as protein-coding.
    weight : float
        0.45 when is_coding=False (regulatory variant dominates);
        0 otherwise (AlphaMissense wins for coding variants).
    """

    score: float
    classification: str
    is_coding: bool
    weight: float


def compute_avi_component(
    chrom: str | None,
    position: int,
    ref_dna: str | None,
    alt_dna: str | None,
    avi_lookup: Callable | None,
) -> AVIComponent | None:
    """Resolve AlphaGenome Atlas AVI for a DNA-level variant.

    Returns ``None`` if any required argument is missing, or if the
    lookup fails silently. The returned component's ``weight`` is 0.45
    only when the variant is non-coding regulatory — for coding-region
    variants the weight is 0 (AVI is a secondary observation but does
    not compete with AlphaMissense for the dominant-signal slot).
    """
    if (
        avi_lookup is None
        or chrom is None
        or ref_dna is None
        or alt_dna is None
    ):
        return None
    try:
        result = avi_lookup(chrom, position, ref_dna, alt_dna)
    except Exception:
        return None
    if result is None:
        return None
    return AVIComponent(
        score=result.score,
        classification=result.classification,
        is_coding=result.is_coding,
        weight=0.45 if not result.is_coding else 0.0,
    )


# ---------------------------------------------------------------------------
# PhyloP46way conservation
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ConservationComponent:
    """Result of a PhyloP46way lookup, normalized for scoring.

    score : float
        PhyloP score in [-1, 1]. Out-of-range values are dropped
        (see ``compute_conservation_component``).
    """

    score: float


def compute_conservation_component(
    chrom: str | None,
    position: int,
    conservation_lookup: Callable | None,
) -> ConservationComponent | None:
    """Resolve PhyloP46way evolutionary-conservation at a DNA position.

    Returns ``None`` if ``chrom`` or ``conservation_lookup`` is
    missing, if the lookup raises, or if the returned value is
    outside [-1, 1].
    """
    if chrom is None or conservation_lookup is None:
        return None
    try:
        raw = conservation_lookup(chrom, position)
    except Exception:
        return None
    if raw is None:
        return None
    try:
        val = float(raw)
    except (TypeError, ValueError):
        return None
    if not -1.0 <= val <= 1.0:
        return None
    return ConservationComponent(score=val)


__all__ = [
    "AlphaMissenseComponent",
    "AVIComponent",
    "ConservationComponent",
    "compute_am_component",
    "compute_avi_component",
    "compute_conservation_component",
]
