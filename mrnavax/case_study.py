"""Variant-prioritization case study — worked example grounded in
ClinVar.

This module ships a small curated set of ClinVar-style variants and
provides scoring + precision@K helpers for running a reproducible
benchmark. The point isn't to compete with published benchmarks
(those use thousands of variants and per-variant gold-standard
labels); the point is to give users a starting point they can
swap in their own variants against and see how the toolkit's
prioritization actually ranks them.

The bundled CSV at ``mrnavax/examples/clinvar_curated.csv`` contains
a mix of:

* Known driver-coding variants (BRAF V600E, KRAS G12D, TP53 R175H)
  — coding pathogenic, should rank high via AlphaMissense
* Regulatory-region pathogenic variants (chr5 APC promoter,
  chr11 intergenic) — non-coding pathogenic, should rank high
  via AlphaGenome Atlas AVI
* Benign / likely-benign variants — should rank low

Each row carries:

* ``gene`` / ``chrom`` / ``pos`` / ``ref`` / ``alt`` — DNA-level
  coordinates for AlphaGenome Atlas lookup
* ``label`` — human-readable identifier
* ``pathogenicity`` — one of "pathogenic", "benign", "uncertain"
* ``is_coding`` — "true" / "false" hint; the real Atlas API
  reports per-variant is_coding
* The CSV is intentionally small (a dozen rows) so the example
  runs in seconds without any network access.

Usage
-----
    from mrnavax.case_study import (
        load_clinvar_variants,
        score_case_study_variants,
        precision_at_k,
    )
    variants = load_clinvar_variants()  # default path
    scored = score_case_study_variants(variants, my_scorer)
    p_at_5 = precision_at_k([(s["label"], s["score"], s["pathogenicity"])
                            for s in scored], k=5)

Reference
---------
The ClinVar variant set is curated from real hg38 positions reported
in Landrum et al., "ClinGen's Variant Curation Interface,"
*Genetics in Medicine* (2024) and the AlphaMissense paper
(Cheng et al., *Science* 2023). Pathogenicity labels follow ClinVar's
3-tier review status.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ClinVarVariant:
    """One curated ClinVar-style variant for the case study.

    Attributes
    ----------
    gene : str
        HGNC gene symbol (e.g. "BRAF").
    chrom : str
        Chromosome with "chr" prefix (e.g. "chr7").
    pos : int
        1-based hg38 position.
    ref : str
        Reference allele (single nucleotide).
    alt : str
        Alternate allele (single nucleotide).
    label : str
        Human-readable identifier (e.g. "BRAF V600E").
    pathogenicity : str
        One of "pathogenic", "benign", "uncertain".
    is_coding : bool
        Whether the variant sits in a protein-coding region.
        For the case study this is a curated hint; the real Atlas
        API reports per-variant is_coding.
    """

    gene: str
    chrom: str
    pos: int
    ref: str
    alt: str
    label: str
    pathogenicity: str
    is_coding: bool


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

DEFAULT_CLINVAR_CSV = (
    Path(__file__).resolve().parent / "examples" / "clinvar_curated.csv"
)


def load_clinvar_variants(path: str | Path | None = None) -> list[ClinVarVariant]:
    """Load a ClinVar-style CSV into a list of ``ClinVarVariant``.

    The CSV must have these columns: gene, chrom, pos, ref, alt,
    label, pathogenicity, is_coding. See the bundled
    ``mrnavax/examples/clinvar_curated.csv`` for the canonical
    shape.

    Parameters
    ----------
    path : str | Path | None
        Path to a CSV. Defaults to the bundled example.
    """
    if path is None:
        path = DEFAULT_CLINVAR_CSV
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"ClinVar CSV not found: {path}. "
            f"Pass a path to a custom CSV or use the bundled example."
        )
    out: list[ClinVarVariant] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            out.append(
                ClinVarVariant(
                    gene=row["gene"],
                    chrom=row["chrom"],
                    pos=int(row["pos"]),
                    ref=row["ref"],
                    alt=row["alt"],
                    label=row["label"],
                    pathogenicity=row["pathogenicity"],
                    is_coding=row["is_coding"].strip().lower() in ("true", "1", "yes"),
                )
            )
    return out


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

ScoringResult = dict  # type alias for clarity: {label, score, pathogenicity, is_coding, gene}


def score_case_study_variants(
    variants: list[ClinVarVariant],
    score_fn: Callable[[ClinVarVariant], dict],
) -> list[ScoringResult]:
    """Run ``score_fn`` on each variant and return a sorted result list.

    The score_fn callable should accept a ClinVarVariant and return
    a dict with at least ``{"score": float}``. Any extra keys in
    the returned dict (e.g. ``"rationale"``) are preserved in the
    output for downstream reporting.

    Results are sorted by ``score`` descending (most-impactful first).
    """
    out: list[ScoringResult] = []
    for v in variants:
        result = score_fn(v)
        result.setdefault("label", v.label)
        result.setdefault("pathogenicity", v.pathogenicity)
        result.setdefault("gene", v.gene)
        result.setdefault("is_coding", v.is_coding)
        result.setdefault("chrom", v.chrom)
        result.setdefault("pos", v.pos)
        out.append(result)
    out.sort(key=lambda r: r.get("score", 0.0), reverse=True)
    return out


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def precision_at_k(
    scored: Iterable,
    k: int,
    *,
    sort: bool = True,
) -> float:
    """Precision@K over pathogenicity-labelled scored variants.

    Parameters
    ----------
    scored : iterable of (id, score, pathogenicity) tuples
        Each element is a 3-tuple. ``id`` is an arbitrary label,
        ``score`` is a float, ``pathogenicity`` is a string with
        values "pathogenic" / "benign" / "uncertain".
    k : int
        Top-K. If k > len(scored), uses min(k, len(scored)).
    sort : bool
        If True (default), sort by score descending before taking
        top-K. If False, assume the iterable is already sorted.

    Returns
    -------
    float
        Fraction of pathogenic variants in top-K. Returns 0.0 for
        empty input.
    """
    items = list(scored)
    if not items:
        return 0.0
    if sort:
        items.sort(key=lambda x: x[1], reverse=True)
    n = min(k, len(items))
    if n == 0:
        return 0.0
    pathogenic_count = sum(1 for _, _, patho in items[:n] if patho == "pathogenic")
    return pathogenic_count / n


def summarize_variant_set(variants: list[ClinVarVariant]) -> dict:
    """Return counts of variants by pathogenicity and is_coding bucket."""
    out = {
        "n_total": len(variants),
        "n_pathogenic": sum(1 for v in variants if v.pathogenicity == "pathogenic"),
        "n_benign": sum(1 for v in variants if v.pathogenicity == "benign"),
        "n_uncertain": sum(1 for v in variants if v.pathogenicity == "uncertain"),
        "n_pathogenic_coding": sum(
            1
            for v in variants
            if v.pathogenicity == "pathogenic" and v.is_coding
        ),
        "n_pathogenic_regulatory": sum(
            1
            for v in variants
            if v.pathogenicity == "pathogenic" and not v.is_coding
        ),
        "n_benign_coding": sum(
            1 for v in variants if v.pathogenicity == "benign" and v.is_coding
        ),
        "n_benign_regulatory": sum(
            1 for v in variants if v.pathogenicity == "benign" and not v.is_coding
        ),
    }
    return out


__all__ = [
    "ClinVarVariant",
    "DEFAULT_CLINVAR_CSV",
    "ScoringResult",
    "load_clinvar_variants",
    "score_case_study_variants",
    "precision_at_k",
    "summarize_variant_set",
]
