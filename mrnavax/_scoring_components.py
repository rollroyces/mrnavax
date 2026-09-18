"""Internal scoring component helpers for ``variant_scorer.score_variant``.

Extracted from ``variant_scorer.py`` in v0.22.0 to keep the public
``score_variant`` function focused on orchestration rather than
detail. Private API; downstream code should call
``mrnavax.variant_scorer.score_variant`` directly.

These helpers encapsulate the four "local" scoring components
(BLOSUM62, driver-gene boost, hydrophobicity change, structural
disruption via Chou-Fasman) plus the position-penalty logic. They
are stdlib-only, deterministic, and have no upstream-model
dependencies — the upstream-model integration (AlphaMissense,
AlphaGenome Atlas AVI, PhyloP46way) lives in ``_scoring_lookups.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

from .variant_scorer import BLOSUM62, HYDROPHOBICITY

# ---------------------------------------------------------------------------
# Local (no-upstream) components
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class LocalComponents:
    """Per-variant local scoring components (no upstream model).

    Attributes
    ----------
    blosum_score : int
        Raw BLOSUM62 score for (wt_aa, mut_aa) substitution.
    blosum_norm : float
        BLOSUM62 score mapped to [0, 1] (high = disruptive).
    driver_gene : bool
        True if the variant's gene is a known cancer driver.
    hydro_delta : float
        |Δhydrophobicity| between wt_aa and mut_aa (in [-0, ~9]).
    hydro_norm : float
        hydro_delta mapped to [0, 1] (high = large change).
    position_penalty : float
        Negative for variants near N/C termini (likely cleaved); 0
        otherwise.
    structural_disruption : float
        Chou-Fasman structural-disruption score in [0, 1] (high =
        likely to break the local secondary structure).
    """

    blosum_score: int
    blosum_norm: float
    driver_gene: bool
    hydro_delta: float
    hydro_norm: float
    position_penalty: float
    structural_disruption: float


def _structural_disruption_for_test(*args, **kwargs):
    """Lazy import of structural_disruption_penalty (defined later in
    variant_scorer) so this module has no circular dependency.

    Resolved at call time, not at import time, via the public name
    on the variant_scorer module."""
    from .variant_scorer import structural_disruption_penalty

    return structural_disruption_penalty(*args, **kwargs)


def compute_local_components(
    gene: str,
    wt_aa: str,
    mut_aa: str,
    position: int,
    *,
    protein_length: int | None = None,
    protein_sequence: str | None = None,
    driver_genes: set[str] | None = None,
) -> LocalComponents:
    """Compute the local (no-upstream) scoring components.

    See :class:`LocalComponents` for field descriptions. The
    Chou-Fasman structural-disruption component is computed only
    when ``protein_sequence`` is supplied.
    """
    from .variant_scorer import DRIVER_GENES

    driver = driver_genes if driver_genes is not None else DRIVER_GENES

    # 1. Substitution severity (BLOSUM62). Low BLOSUM = high priority.
    blosum = BLOSUM62.get((wt_aa, mut_aa), -4)
    # BLOSUM ranges roughly -4 (disruptive) to +11 (identity).
    # Map to [0, 1] where low BLOSUM = high score.
    blosum_norm = 1.0 - (blosum + 4) / 15.0
    blosum_norm = max(0.0, min(1.0, blosum_norm))

    # 2. Driver-gene boost (binary on/off).
    driver_gene = gene in driver

    # 3. Position penalty (N/C termini are often cleaved during
    # antigen processing).
    if protein_length is not None:
        n_term = position <= 10
        c_term = position >= protein_length - 10
        position_penalty = -0.15 if (n_term or c_term) else 0.0
    else:
        position_penalty = 0.0

    # 4. Hydrophobicity change.
    h_wt = HYDROPHOBICITY.get(wt_aa, 0.0)
    h_mut = HYDROPHOBICITY.get(mut_aa, 0.0)
    hydro_delta = abs(h_wt - h_mut)
    hydro_norm = min(1.0, hydro_delta / 9.0)

    # 5. Structural disruption (Chou-Fasman).
    struct = 0.0
    if protein_sequence is not None and len(protein_sequence) >= position:
        struct = _structural_disruption_for_test(
            protein_sequence, position, wt_aa, mut_aa
        )

    return LocalComponents(
        blosum_score=blosum,
        blosum_norm=blosum_norm,
        driver_gene=driver_gene,
        hydro_delta=hydro_delta,
        hydro_norm=hydro_norm,
        position_penalty=position_penalty,
        structural_disruption=struct,
    )


# ---------------------------------------------------------------------------
# Rationale parts for the local components
# ---------------------------------------------------------------------------

def local_rationale(local: LocalComponents, gene: str, wt_aa: str, mut_aa: str) -> list[str]:
    """Build the rationale string fragments for the local components.

    Returns lines like::

        BLOSUM62 V→E = -2
        Δhydrophobicity = 7.7
        BRAF is a known driver gene (+boost)
        position 600 near terminus (penalty)
    """
    parts = [
        f"BLOSUM62 {wt_aa}→{mut_aa} = {local.blosum_score}",
        f"Δhydrophobicity = {local.hydro_delta:.1f}",
    ]
    if local.driver_gene:
        parts.append(f"{gene} is a known driver gene (+boost)")
    if local.position_penalty < 0:
        parts.append(
            "position near terminus (penalty)"
        )
    return parts


__all__ = [
    "LocalComponents",
    "compute_local_components",
    "local_rationale",
]
