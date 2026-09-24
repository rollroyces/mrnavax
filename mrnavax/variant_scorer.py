"""AlphaMissense-style variant pathogenicity scorer.

Lightweight, deterministic, stdlib-only implementation that captures the
operational idea behind AlphaMissense (Cheng et al. *Science* 2023): score
each missense variant by how likely it is to be functionally impactful, then
filter to the top-N candidates before downstream peptide enumeration.

Scoring components (all stdlib, no model downloads):

1. **Position in protein** — N- and C-terminal residues get a mild penalty
   (peptides from the termini are often cleaved during antigen processing).
2. **Substitution severity** — BLOSUM62-style substitution matrix embedded
   below; lower score = more disruptive substitution = higher priority.
3. **Driver-gene boost** — known oncogenes / tumor suppressors get a flat
   +0.2 priority boost.
4. **Hydrophobicity change** — dramatic Δhydrophobicity suggests a
   conformational effect.

The output is a per-variant score in [0, 1] that the `scrna` tool uses to
filter the candidate list before peptide enumeration.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Callable

# ---------- BLOSUM62-style substitution matrix (subset, stdlib-only) -------
# Rows = wildtype, cols = mutation. Standard biology-textbook values.
# Negative = disruptive, positive = conservative.
BLOSUM62: dict[tuple[str, str], int] = {
    ("A", "A"): 4,
    ("A", "R"): -1,
    ("A", "N"): -2,
    ("A", "D"): -2,
    ("A", "C"): 0,
    ("A", "Q"): -1,
    ("A", "E"): -1,
    ("A", "G"): 0,
    ("A", "H"): -2,
    ("A", "I"): -1,
    ("A", "L"): -1,
    ("A", "K"): -1,
    ("A", "M"): -1,
    ("A", "F"): -2,
    ("A", "P"): -1,
    ("A", "S"): 1,
    ("A", "T"): 0,
    ("A", "W"): -3,
    ("A", "Y"): -2,
    ("A", "V"): 0,
    ("R", "A"): -1,
    ("R", "R"): 5,
    ("R", "N"): 0,
    ("R", "D"): -2,
    ("R", "C"): -3,
    ("R", "Q"): 1,
    ("R", "E"): 0,
    ("R", "G"): -2,
    ("R", "H"): 0,
    ("R", "I"): -3,
    ("R", "L"): -2,
    ("R", "K"): 2,
    ("R", "M"): -1,
    ("R", "F"): -3,
    ("R", "P"): -2,
    ("R", "S"): -1,
    ("R", "T"): -1,
    ("R", "W"): -3,
    ("R", "Y"): -2,
    ("R", "V"): -3,
    ("N", "A"): -2,
    ("N", "R"): 0,
    ("N", "N"): 6,
    ("N", "D"): 1,
    ("N", "C"): -3,
    ("N", "Q"): 0,
    ("N", "E"): 0,
    ("N", "G"): 0,
    ("N", "H"): 1,
    ("N", "I"): -3,
    ("N", "L"): -3,
    ("N", "K"): 0,
    ("N", "M"): -2,
    ("N", "F"): -3,
    ("N", "P"): -2,
    ("N", "S"): 1,
    ("N", "T"): 0,
    ("N", "W"): -4,
    ("N", "Y"): -2,
    ("N", "V"): -3,
    ("D", "A"): -2,
    ("D", "R"): -2,
    ("D", "N"): 1,
    ("D", "D"): 6,
    ("D", "C"): -3,
    ("D", "Q"): 0,
    ("D", "E"): 2,
    ("D", "G"): -1,
    ("D", "H"): -1,
    ("D", "I"): -3,
    ("D", "L"): -4,
    ("D", "K"): -1,
    ("D", "M"): -3,
    ("D", "F"): -3,
    ("D", "P"): -1,
    ("D", "S"): 0,
    ("D", "T"): -1,
    ("D", "W"): -4,
    ("D", "Y"): -3,
    ("D", "V"): -3,
    ("C", "A"): 0,
    ("C", "R"): -3,
    ("C", "N"): -3,
    ("C", "D"): -3,
    ("C", "C"): 9,
    ("C", "Q"): -3,
    ("C", "E"): -4,
    ("C", "G"): -3,
    ("C", "H"): -3,
    ("C", "I"): -1,
    ("C", "L"): -1,
    ("C", "K"): -3,
    ("C", "M"): -1,
    ("C", "F"): -2,
    ("C", "P"): -3,
    ("C", "S"): -1,
    ("C", "T"): -1,
    ("C", "W"): -2,
    ("C", "Y"): -2,
    ("C", "V"): -1,
    ("Q", "A"): -1,
    ("Q", "R"): 1,
    ("Q", "N"): 0,
    ("Q", "D"): 0,
    ("Q", "C"): -3,
    ("Q", "Q"): 5,
    ("Q", "E"): 2,
    ("Q", "G"): -2,
    ("Q", "H"): 0,
    ("Q", "I"): -3,
    ("Q", "L"): -2,
    ("Q", "K"): 1,
    ("Q", "M"): 0,
    ("Q", "F"): -3,
    ("Q", "P"): -1,
    ("Q", "S"): 0,
    ("Q", "T"): -1,
    ("Q", "W"): -2,
    ("Q", "Y"): -1,
    ("Q", "V"): -2,
    ("E", "A"): -1,
    ("E", "R"): 0,
    ("E", "N"): 0,
    ("E", "D"): 2,
    ("E", "C"): -4,
    ("E", "Q"): 2,
    ("E", "E"): 5,
    ("E", "G"): -2,
    ("E", "H"): 0,
    ("E", "I"): -3,
    ("E", "L"): -3,
    ("E", "K"): 1,
    ("E", "M"): -2,
    ("E", "F"): -3,
    ("E", "P"): -1,
    ("E", "S"): 0,
    ("E", "T"): -1,
    ("E", "W"): -3,
    ("E", "Y"): -2,
    ("E", "V"): -2,
    ("G", "A"): 0,
    ("G", "R"): -2,
    ("G", "N"): 0,
    ("G", "D"): -1,
    ("G", "C"): -3,
    ("G", "Q"): -2,
    ("G", "E"): -2,
    ("G", "G"): 6,
    ("G", "H"): -2,
    ("G", "I"): -4,
    ("G", "L"): -4,
    ("G", "K"): -2,
    ("G", "M"): -3,
    ("G", "F"): -3,
    ("G", "P"): -2,
    ("G", "S"): 0,
    ("G", "T"): -2,
    ("G", "W"): -2,
    ("G", "Y"): -3,
    ("G", "V"): -3,
    ("H", "A"): -2,
    ("H", "R"): 0,
    ("H", "N"): 1,
    ("H", "D"): -1,
    ("H", "C"): -3,
    ("H", "Q"): 0,
    ("H", "E"): 0,
    ("H", "G"): -2,
    ("H", "H"): 8,
    ("H", "I"): -3,
    ("H", "L"): -3,
    ("H", "K"): -1,
    ("H", "M"): -2,
    ("H", "F"): -1,
    ("H", "P"): -2,
    ("H", "S"): -1,
    ("H", "T"): -2,
    ("H", "W"): -2,
    ("H", "Y"): 2,
    ("H", "V"): -3,
    ("I", "A"): -1,
    ("I", "R"): -3,
    ("I", "N"): -3,
    ("I", "D"): -3,
    ("I", "C"): -1,
    ("I", "Q"): -3,
    ("I", "E"): -3,
    ("I", "G"): -4,
    ("I", "H"): -3,
    ("I", "I"): 4,
    ("I", "L"): 2,
    ("I", "K"): -3,
    ("I", "M"): 1,
    ("I", "F"): 0,
    ("I", "P"): -3,
    ("I", "S"): -2,
    ("I", "T"): -1,
    ("I", "W"): -3,
    ("I", "Y"): -1,
    ("I", "V"): 3,
    ("L", "A"): -1,
    ("L", "R"): -2,
    ("L", "N"): -3,
    ("L", "D"): -4,
    ("L", "C"): -1,
    ("L", "Q"): -2,
    ("L", "E"): -3,
    ("L", "G"): -4,
    ("L", "H"): -3,
    ("L", "I"): 2,
    ("L", "L"): 4,
    ("L", "K"): -2,
    ("L", "M"): 2,
    ("L", "F"): 0,
    ("L", "P"): -3,
    ("L", "S"): -2,
    ("L", "T"): -1,
    ("L", "W"): -2,
    ("L", "Y"): -1,
    ("L", "V"): 1,
    ("K", "A"): -1,
    ("K", "R"): 2,
    ("K", "N"): 0,
    ("K", "D"): -1,
    ("K", "C"): -3,
    ("K", "Q"): 1,
    ("K", "E"): 1,
    ("K", "G"): -2,
    ("K", "H"): -1,
    ("K", "I"): -3,
    ("K", "L"): -2,
    ("K", "K"): 5,
    ("K", "M"): -1,
    ("K", "F"): -3,
    ("K", "P"): -1,
    ("K", "S"): 0,
    ("K", "T"): -1,
    ("K", "W"): -3,
    ("K", "Y"): -2,
    ("K", "V"): -2,
    ("M", "A"): -1,
    ("M", "R"): -1,
    ("M", "N"): -2,
    ("M", "D"): -3,
    ("M", "C"): -1,
    ("M", "Q"): 0,
    ("M", "E"): -2,
    ("M", "G"): -3,
    ("M", "H"): -2,
    ("M", "I"): 1,
    ("M", "L"): 2,
    ("M", "K"): -1,
    ("M", "M"): 5,
    ("M", "F"): 0,
    ("M", "P"): -2,
    ("M", "S"): -1,
    ("M", "T"): -1,
    ("M", "W"): -1,
    ("M", "Y"): -1,
    ("M", "V"): 1,
    ("F", "A"): -2,
    ("F", "R"): -3,
    ("F", "N"): -3,
    ("F", "D"): -3,
    ("F", "C"): -2,
    ("F", "Q"): -3,
    ("F", "E"): -3,
    ("F", "G"): -3,
    ("F", "H"): -1,
    ("F", "I"): 0,
    ("F", "L"): 0,
    ("F", "K"): -3,
    ("F", "M"): 0,
    ("F", "F"): 6,
    ("F", "P"): -4,
    ("F", "S"): -2,
    ("F", "T"): -2,
    ("F", "W"): 1,
    ("F", "Y"): 3,
    ("F", "V"): -1,
    ("P", "A"): -1,
    ("P", "R"): -2,
    ("P", "N"): -2,
    ("P", "D"): -1,
    ("P", "C"): -3,
    ("P", "Q"): -1,
    ("P", "E"): -1,
    ("P", "G"): -2,
    ("P", "H"): -2,
    ("P", "I"): -3,
    ("P", "L"): -3,
    ("P", "K"): -1,
    ("P", "M"): -2,
    ("P", "F"): -4,
    ("P", "P"): 7,
    ("P", "S"): -1,
    ("P", "T"): -1,
    ("P", "W"): -4,
    ("P", "Y"): -3,
    ("P", "V"): -2,
    ("S", "A"): 1,
    ("S", "R"): -1,
    ("S", "N"): 1,
    ("S", "D"): 0,
    ("S", "C"): -1,
    ("S", "Q"): 0,
    ("S", "E"): 0,
    ("S", "G"): 0,
    ("S", "H"): -1,
    ("S", "I"): -2,
    ("S", "L"): -2,
    ("S", "K"): 0,
    ("S", "M"): -1,
    ("S", "F"): -2,
    ("S", "P"): -1,
    ("S", "S"): 4,
    ("S", "T"): 1,
    ("S", "W"): -3,
    ("S", "Y"): -2,
    ("S", "V"): -2,
    ("T", "A"): 0,
    ("T", "R"): -1,
    ("T", "N"): 0,
    ("T", "D"): -1,
    ("T", "C"): -1,
    ("T", "Q"): -1,
    ("T", "E"): -1,
    ("T", "G"): -2,
    ("T", "H"): -2,
    ("T", "I"): -1,
    ("T", "L"): -1,
    ("T", "K"): -1,
    ("T", "M"): -1,
    ("T", "F"): -2,
    ("T", "P"): -1,
    ("T", "S"): 1,
    ("T", "T"): 4,
    ("T", "W"): -2,
    ("T", "Y"): -2,
    ("T", "V"): 0,
    ("W", "A"): -3,
    ("W", "R"): -3,
    ("W", "N"): -4,
    ("W", "D"): -4,
    ("W", "C"): -2,
    ("W", "Q"): -2,
    ("W", "E"): -3,
    ("W", "G"): -2,
    ("W", "H"): -2,
    ("W", "I"): -3,
    ("W", "L"): -2,
    ("W", "K"): -3,
    ("W", "M"): -1,
    ("W", "F"): 1,
    ("W", "P"): -4,
    ("W", "S"): -3,
    ("W", "T"): -2,
    ("W", "W"): 11,
    ("W", "Y"): 2,
    ("W", "V"): -3,
    ("Y", "A"): -2,
    ("Y", "R"): -2,
    ("Y", "N"): -2,
    ("Y", "D"): -3,
    ("Y", "C"): -2,
    ("Y", "Q"): -1,
    ("Y", "E"): -2,
    ("Y", "G"): -3,
    ("Y", "H"): 2,
    ("Y", "I"): -1,
    ("Y", "L"): -1,
    ("Y", "K"): -2,
    ("Y", "M"): -1,
    ("Y", "F"): 3,
    ("Y", "P"): -3,
    ("Y", "S"): -2,
    ("Y", "T"): -2,
    ("Y", "W"): 2,
    ("Y", "Y"): 7,
    ("Y", "V"): -1,
    ("V", "A"): 0,
    ("V", "R"): -3,
    ("V", "N"): -3,
    ("V", "D"): -3,
    ("V", "C"): -1,
    ("V", "Q"): -2,
    ("V", "E"): -2,
    ("V", "G"): -3,
    ("V", "H"): -3,
    ("V", "I"): 3,
    ("V", "L"): 1,
    ("V", "K"): -2,
    ("V", "M"): 1,
    ("V", "F"): -1,
    ("V", "P"): -2,
    ("V", "S"): -2,
    ("V", "T"): 0,
    ("V", "W"): -3,
    ("V", "Y"): -1,
    ("V", "V"): 4,
}

# Known oncogenes / tumor suppressors. Boost priority for variants in these.
DRIVER_GENES: set[str] = {
    "TP53",
    "KRAS",
    "BRAF",
    "EGFR",
    "MYC",
    "PTEN",
    "PIK3CA",
    "APC",
    "BRCA1",
    "BRCA2",
    "CDKN2A",
    "RB1",
    "NRAS",
    "HRAS",
    "ALK",
    "ROS1",
    "MET",
    "RET",
    "HER2",
    "ERBB2",
    "IDH1",
    "IDH2",
    "FLT3",
    "JAK2",
    "KIT",
    "PDGFRA",
    "CTNNB1",
    "SMAD4",
    "STK11",
    "KEAP1",
    "ARID1A",
}

# Kyte-Doolittle hydrophobicity (subset)
HYDROPHOBICITY: dict[str, float] = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}

# Chou-Fasman helix/strand propensity (P_alpha, P_beta). Higher = more
# likely to be in that secondary structure. From Chou & Fasman (1978),
# normalized to ~[0.4, 1.7].
HELIX_PROPENSITY: dict[str, float] = {
    "A": 1.42,
    "R": 0.98,
    "N": 0.67,
    "D": 1.01,
    "C": 0.70,
    "Q": 1.11,
    "E": 1.51,
    "G": 0.57,
    "H": 1.00,
    "I": 1.08,
    "L": 1.21,
    "K": 1.16,
    "M": 1.45,
    "F": 1.13,
    "P": 0.57,
    "S": 0.77,
    "T": 0.83,
    "W": 1.08,
    "Y": 0.69,
    "V": 1.06,
}

STRAND_PROPENSITY: dict[str, float] = {
    "A": 0.83,
    "R": 0.93,
    "N": 0.89,
    "D": 0.54,
    "C": 1.19,
    "Q": 1.10,
    "E": 0.37,
    "G": 0.75,
    "H": 0.87,
    "I": 1.60,
    "L": 1.30,
    "K": 0.74,
    "M": 1.05,
    "F": 1.38,
    "P": 0.55,
    "S": 0.75,
    "T": 1.19,
    "W": 1.37,
    "Y": 1.47,
    "V": 1.70,
}


@dataclass
class VariantScore:
    gene: str
    position: int
    wt_aa: str
    mut_aa: str
    raw_score: float
    normalized_score: float  # 0–1
    components: dict[str, float]  # sub-scores
    rationale: str


@dataclass
class VariantInputs:
    """Structured input for ``score_variant()`` and ``score_variant_from_inputs()``.

    The four positional fields of ``score_variant()`` (gene, position,
    wt_aa, mut_aa) are required; all other fields are optional. Pass
    this dataclass as the single argument to
    ``score_variant_from_inputs(inputs)`` to avoid the 14-kwarg
    signature of the legacy ``score_variant()`` function.

    New in v0.29.0 — the dataclass is additive; existing kwargs callers
    of ``score_variant()`` are unaffected.
    """
    gene: str
    position: int
    wt_aa: str
    mut_aa: str
    chrom: str | None = None
    ref_dna: str | None = None
    alt_dna: str | None = None
    pos: int | None = None
    protein_length: int | None = None
    protein_sequence: str | None = None
    uniprot_id: str | None = None
    am_lookup: Callable | None = None
    avi_lookup: Callable | None = None
    conservation_lookup: Callable | None = None
    driver_genes: set[str] | None = None
    strict: bool = False
    utr5: str | None = None
    utr3: str | None = None

    @classmethod
    def from_kwargs(
        cls,
        gene: str,
        position: int,
        wt_aa: str,
        mut_aa: str,
        **kwargs,
    ) -> "VariantInputs":
        """Build from the same kwargs shape as ``score_variant()``."""
        return cls(
            gene=gene,
            position=position,
            wt_aa=wt_aa,
            mut_aa=mut_aa,
            **kwargs,
        )


def score_variant_from_inputs(inputs: VariantInputs) -> VariantScore | None:
    """Score a single missense variant from a structured ``VariantInputs``.

    Functionally identical to ``score_variant()`` but accepts the full
    input as a single dataclass argument. New in v0.29.0.
    """
    return score_variant(
        inputs.gene,
        inputs.position,
        inputs.wt_aa,
        inputs.mut_aa,
        chrom=inputs.chrom,
        ref_dna=inputs.ref_dna,
        alt_dna=inputs.alt_dna,
        pos=inputs.pos,
        protein_length=inputs.protein_length,
        protein_sequence=inputs.protein_sequence,
        uniprot_id=inputs.uniprot_id,
        am_lookup=inputs.am_lookup,
        avi_lookup=inputs.avi_lookup,
        conservation_lookup=inputs.conservation_lookup,
        driver_genes=inputs.driver_genes,
        strict=inputs.strict,
        utr5=inputs.utr5,
        utr3=inputs.utr3,
    )


def score_variant(
    gene: str,
    position: int,
    wt_aa: str,
    mut_aa: str,
    *,
    chrom: str | None = None,
    ref_dna: str | None = None,
    alt_dna: str | None = None,
    pos: int | None = None,
    protein_length: int | None = None,
    protein_sequence: str | None = None,
    uniprot_id: str | None = None,
    am_lookup: Callable | None = None,
    avi_lookup: Callable | None = None,
    conservation_lookup: Callable | None = None,
    driver_genes: set[str] | None = None,
    strict: bool = False,
    utr5: str | None = None,
    utr3: str | None = None,
) -> VariantScore | None:
    """Score a single missense variant.

    Higher ``normalized_score`` = higher priority for downstream analysis.
    Returns ``None`` for unknown amino acids (silent-mode default) or raises
    ``ValueError`` when ``strict=True``.

    The 5-signal composition (v0.27.0 adds UTR context as a 5th signal):
      * AlphaMissense (Cheng 2023) — dominant for coding variants (0.45)
      * AlphaGenome Atlas AVI (Avsec 2026) — dominant for non-coding
        regulatory variants (0.45)
      * BLOSUM62 + driver-gene + hydrophobicity + structural
        disruption — always computed, split into the residual weight
      * PhyloP46way conservation (Pollard 2010) — small bonus
        when conservation > 0
      * UTR context (Kozak + ARE-motif density) — optional bonus,
        default weight 0.05, only when ``utr5``/``utr3`` supplied

    All upstream-model lookups are silent: a network error / None
    return / out-of-range value drops the component without raising.
    See ``mrnavax/_scoring_components.py`` and
    ``mrnavax/_scoring_lookups.py`` for the implementation.

    Call ``avi_lookup`` as ``avi_lookup(chrom, position, ref_dna,
    alt_dna)``. Call ``conservation_lookup`` as
    ``conservation_lookup(chrom, pos)`` (pos defaults to ``position``
    if not supplied separately). If any of the DNA-level arguments
    are missing, the corresponding lookup is skipped silently.
    """
    from ._scoring_components import (
        compute_local_components,
        local_rationale,
    )
    from ._scoring_lookups import (
        compute_am_component,
        compute_avi_component,
        compute_conservation_component,
    )

    wt_aa = wt_aa.upper()
    mut_aa = mut_aa.upper()

    # Validate amino acids — silent default returns None rather than
    # silently maxing the BLOSUM penalty, which would inflate the score.
    if wt_aa not in HYDROPHOBICITY or mut_aa not in HYDROPHOBICITY:
        if strict:
            raise ValueError(
                f"unknown amino acid(s): wt={wt_aa!r}, mut={mut_aa!r}; "
                "expected standard 20-letter amino acid codes"
            )
        return None

    # ---- Local (no-upstream) components ----
    local = compute_local_components(
        gene=gene,
        wt_aa=wt_aa,
        mut_aa=mut_aa,
        position=position,
        protein_length=protein_length,
        protein_sequence=protein_sequence,
        driver_genes=driver_genes,
    )

    # ---- Upstream-model components (silent-failure) ----
    am = compute_am_component(uniprot_id, wt_aa, position, mut_aa, am_lookup)
    avi = compute_avi_component(chrom, position, ref_dna, alt_dna, avi_lookup)
    cons = compute_conservation_component(
        chrom, pos if pos is not None else position, conservation_lookup
    )

    # ---- Weighted combination ----
    # Priority of dominant signal: AVI > AM > neither. Conservation
    # adds a small bonus when present.
    utr_ctx = None
    if utr5 is not None or utr3 is not None:
        from ._utr_context import score_utr_context
        utr_ctx = score_utr_context(utr5, utr3)
    raw, components = _combine_components(local, am, avi, cons, utr_ctx)

    # ---- Rationale ----
    rationale_parts = local_rationale(local, gene, wt_aa, mut_aa)
    if am.weight > 0 or am.classification:
        rationale_parts.append(
            f"AlphaMissense pathogenicity={am.score:.3f} ({am.classification})"
        )
    if avi is not None:
        tag = " (dominant)" if avi.weight > 0 else ""
        rationale_parts.append(
            f"AlphaGenome AVI={avi.score:.3f} ({avi.classification}){tag}"
        )
    if cons is not None:
        rationale_parts.append(f"PhyloP46way conservation={cons.score:+.3f}")
    if utr_ctx is not None:
        rationale_parts.append(
            f"UTR context={utr_ctx.context_score:.3f} "
            f"(Kozak={utr_ctx.kozak_score:.3f}, 3'UTR={utr_ctx.utr3_score:.3f})"
        )

    normalized = max(0.0, min(1.0, raw))
    return VariantScore(
        gene=gene,
        position=position,
        wt_aa=wt_aa,
        mut_aa=mut_aa,
        raw_score=round(raw, 4),
        normalized_score=round(normalized, 4),
        components=components,
        rationale="; ".join(rationale_parts),
    )


def _combine_components(
    local,
    am,
    avi,
    cons,
    utr_ctx=None,
) -> tuple[float, dict]:
    """Combine the 4 (or 5, when ``utr_ctx`` is supplied) scoring
    components into a single raw score + components dict.

    Returns (raw, components). Components dict matches the previous
    contract: BLOSUM62 + driver + hydrophobicity + structural
    disruption always present; AM / AVI / PhyloP entries added when
    available; ``utr_context_score`` added when ``utr_ctx`` is supplied.

    The UTR context contributes a small bonus (weight 0.05) when
    supplied, capped to keep the dominant signals (AM/AVI) dominant.
    """
    components: dict = {
        "blosum62_score": local.blosum_score,
        "blosum62_norm": round(local.blosum_norm, 3),
        "driver_gene": local.driver_gene,
        "hydrophobicity_delta": round(local.hydro_delta, 2),
        "position_penalty": local.position_penalty,
        "structural_disruption": local.structural_disruption,
    }
    if am.weight > 0 or am.classification:
        components["alphamissense_score"] = am.score
        components["alphamissense_classification"] = am.classification

    if avi is not None:
        components["alphagenome_atlas_score"] = avi.score
        components["alphagenome_atlas_classification"] = avi.classification
        if avi.is_coding is not None:
            components["alphagenome_atlas_is_coding"] = avi.is_coding

    if cons is not None:
        components["phylop46way_score"] = cons.score

    if utr_ctx is not None:
        components["utr_context_score"] = utr_ctx.context_score
        components["kozak_score"] = utr_ctx.kozak_score
        components["utr3_score"] = utr_ctx.utr3_score

    # Weighted combination. The dominant signal slot takes 0.45;
    # the residual 0.55 is split among the local components plus a
    # small conservation bonus.
    #
    # PhyloP contribution: gradient ``0.10 * cons.score`` when the
    # PhyloP score is positive (highly-conserved positions add bonus),
    # zero otherwise. Used identically across all three branches so
    # the conservation signal contributes consistently regardless of
    # which dominant signal fires.
    if cons is not None and cons.score > 0:
        phylo_boost = 0.10 * cons.score
    else:
        phylo_boost = 0.0

    # UTR context contribution: small additive bonus when supplied.
    # Weight 0.05 keeps the dominant signals dominant. Caps the
    # total raw score at 1.0 via the max() in score_variant().
    utr_boost = 0.0
    if utr_ctx is not None:
        utr_boost = 0.05 * utr_ctx.context_score

    local_weighted = (
        0.35 * local.blosum_norm
        + 0.20 * (0.2 if local.driver_gene else 0.0)
        + 0.15 * local.hydro_norm
        + 0.20 * local.structural_disruption
        + 0.10 * (1.0 if local.position_penalty == 0 else 0.0)
        + local.position_penalty
    )

    if avi is not None and avi.weight > 0:
        # AVI dominates (regulatory-region variant).
        raw = avi.weight * avi.score + (1.0 - avi.weight) * local_weighted + phylo_boost + utr_boost
    elif am.weight > 0:
        # AlphaMissense dominates (coding-region variant).
        raw = am.weight * am.score + (1.0 - am.weight) * local_weighted + phylo_boost + utr_boost
    else:
        # Neither AVI nor AM — local-only weighted combination + PhyloP + UTR.
        raw = local_weighted + phylo_boost + utr_boost

    return raw, components



def predict_secondary_structure(protein: str, *, window: int = 6) -> dict[int, str]:
    """Chou-Fasman secondary-structure prediction.

    Returns ``{position_1based: 'H'|'E'|'C'}`` where H = helix, E = strand,
    C = coil. Per-residue prediction averaged over a sliding window.

    Reference: Chou PY, Fasman GD. "Prediction of the secondary structure of
    proteins from their amino acid sequence." Adv Enzymol Relat Areas Mol
    Biol 47, 45–148 (1978).

    Memoized via functools.lru_cache: filtering N variants on the same
    protein recomputes the structure prediction N times otherwise
    (O(N × L × W) total). The cache key is the (protein, window) tuple;
    bound at 1024 entries to cap memory.
    """
    if not protein:
        return {}
    return _predict_secondary_structure_inner(protein, window)


@lru_cache(maxsize=1024)
def _predict_secondary_structure_inner(protein: str, window: int) -> dict[int, str]:
    pred: dict[int, str] = {}
    n = len(protein)
    half = window // 2
    for i in range(n):
        lo = max(0, i - half)
        hi = min(n, i + half + 1)
        sub = protein[lo:hi]
        helix_score = sum(HELIX_PROPENSITY.get(a, 1.0) for a in sub) / max(1, len(sub))
        strand_score = sum(STRAND_PROPENSITY.get(a, 1.0) for a in sub) / max(1, len(sub))
        if helix_score >= 1.03 and helix_score > strand_score:
            pred[i + 1] = "H"
        elif strand_score >= 1.05 and strand_score > helix_score:
            pred[i + 1] = "E"
        else:
            pred[i + 1] = "C"
    return pred


def structural_disruption_penalty(
    protein: str,
    position: int,
    wt_aa: str,
    mut_aa: str,
) -> float:
    """Score how much a substitution would disrupt local secondary structure.

    Returns a value in [0, 1]. 0 = predicted loop region. 1 = strong
    helix-breaker or beta-sheet disruptor in a structured region.
    """
    if not protein or position < 1 or position > len(protein):
        return 0.0
    pred = predict_secondary_structure(protein, window=6)
    local = pred.get(position, "C")
    if local == "C":
        return 0.0
    if local == "H":
        delta = HELIX_PROPENSITY.get(wt_aa, 1.0) - HELIX_PROPENSITY.get(mut_aa, 1.0)
    else:
        delta = STRAND_PROPENSITY.get(wt_aa, 1.0) - STRAND_PROPENSITY.get(mut_aa, 1.0)
    return round(max(0.0, min(1.0, delta / 0.8)), 4)


def filter_variants(
    variants: list[dict],
    *,
    top_fraction: float = 0.20,
    min_score: float = 0.4,
    protein_lengths: dict[str, int] | None = None,
    protein_sequences: dict[str, str] | None = None,
    uniprot_ids: dict[str, str] | None = None,
    am_lookup: Callable | None = None,
    avi_lookup: Callable | None = None,
    conservation_lookup: Callable | None = None,
    strict: bool = False,
) -> list[VariantScore]:
    """Score and filter a list of variants to the top ``top_fraction``.

    Returns a list of ``VariantScore`` sorted by ``normalized_score`` descending.
    Variants with unknown amino acids are skipped (or raise if ``strict=True``).

    If ``protein_sequences`` is provided, the scorer also computes a
    Chou-Fasman structural-disruption component.
    If ``uniprot_ids`` and ``am_lookup`` are both provided, AlphaMissense
    pathogenicity scores are included as the dominant signal.
    If ``avi_lookup`` is provided AND the variant dict has
    ``chrom``/``ref_dna``/``alt_dna`` keys, the AlphaGenome Atlas AVI
    score is used as the dominant signal for non-coding regulatory
    variants (Avsec et al. *Nature* 2026). For coding-region variants
    where AVI returns ``is_coding=True``, AVI is recorded as a secondary
    signal and AlphaMissense still dominates.
    If ``conservation_lookup`` is provided AND the variant dict has
    ``chrom`` (and optionally ``pos``), the PhyloP46way conservation
    score is included as a 4th signal (Pollard et al. 2010).
    """
    scored: list[VariantScore] = []
    for v in variants:
        prot_len = None
        prot_seq = None
        uniprot_id = None
        gene = v.get("gene")
        if gene:
            if protein_lengths and gene in protein_lengths:
                prot_len = protein_lengths[gene]
            if protein_sequences and gene in protein_sequences:
                prot_seq = protein_sequences[gene]
            if uniprot_ids and gene in uniprot_ids:
                uniprot_id = uniprot_ids[gene]
        result = score_variant(
            gene=gene,
            position=v["position"],
            wt_aa=v["wt_aa"],
            mut_aa=v["mut_aa"],
            # Pass DNA-level coordinates through to score_variant; it
            # skips the AVI lookup silently if any are missing.
            chrom=v.get("chrom"),
            ref_dna=v.get("ref_dna"),
            alt_dna=v.get("alt_dna"),
            # Pass conservation lookup through; falls back to protein
            # position if pos isn't supplied separately in the dict.
            pos=v.get("pos"),
            protein_length=prot_len,
            protein_sequence=prot_seq,
            uniprot_id=uniprot_id,
            am_lookup=am_lookup,
            avi_lookup=avi_lookup,
            conservation_lookup=conservation_lookup,
            strict=strict,
        )
        if result is not None:
            scored.append(result)
    scored.sort(key=lambda s: -s.normalized_score)

    n_keep = max(1, int(len(scored) * top_fraction))
    threshold_idx = min(n_keep, len(scored))
    threshold_score = scored[threshold_idx - 1].normalized_score if scored else 1.0
    keep = [
        s
        for s in scored
        if s.normalized_score >= threshold_score and s.normalized_score >= min_score
    ]
    return keep
