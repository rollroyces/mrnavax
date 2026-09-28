"""Neoantigen screener — peptide→MHC binding + immunogenicity.

Uses the LLM as the predictor (mimicking how a practitioner would integrate
TrambaHLApan / DeepHLApan / DeepNeo in production without a GPU). The LLM is
given a small, evidence-anchored prompt with the variant's 9-11-mer peptides
and the patient's HLA alleles, and asked to score each.

A simple **heuristic fallback** is included (HLA-A*02:01 anchor matrix) so the
tool runs end-to-end without an API key.

References
----------
TrambaHLApan: https://link.springer.com/article/10.1007/s12539-025-00777-5
DeepNeo:      https://pmc.ncbi.nlm.nih.gov/articles/PMC10320182/
DeepHLApan:   Wu et al. Front Immunol 2019
NetMHCpan:    Reynisson et al. NAR 2020
"""

from __future__ import annotations

import csv
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Iterable

# ---------- heuristic predictor (offline fallback) -------------------------

# Classical HLA-A*02:01 anchor matrix: P2 = L/M, P9 = V/L/I.
# Position-specific weights (rough log-affinity) from IEDB / SYFPEITHI motifs.
A0201_ANCHORS = {
    2: {"L": 1.0, "M": 0.9, "I": 0.5, "V": 0.4},
    9: {"V": 1.0, "L": 0.9, "I": 0.8, "A": 0.4, "M": 0.4},
}


# Hydrophobicity at secondary anchors (P1, P3) helps. Crude penalty for charged.
def _heuristic_a0201(peptide: str) -> tuple[float, bool, str]:
    pep = peptide.upper()
    if len(pep) != 9:
        return (9999.0, False, "non-9mer skipped")
    score = 0.0
    # primary anchors
    for pos, weights in A0201_ANCHORS.items():
        aa = pep[pos - 1]  # 1-indexed
        score += weights.get(aa, 0.0)
    # penalize charged residues outside anchors
    for i, aa in enumerate(pep, start=1):
        if i in A0201_ANCHORS:
            continue
        if aa in "DE":
            score -= 0.3
        if aa in "KR":
            score -= 0.2
    # map to nM affinity — heuristic, not calibrated
    affinity = max(5.0, 5000.0 * 2 ** (-3.0 * score))
    binder = affinity < 500
    return (round(affinity, 1), binder, "heuristic HLA-A*02:01 anchor matrix")


# ---------- LLM prompt ------------------------------------------------------

PROMPT_TEMPLATE = """You are an immunogenicity-prediction model. Score the
neoantigen candidates below for HLA binding and T-cell immunogenicity.

Patient HLA alleles: {hla}

For EACH candidate, return a JSON object with:
  peptide (str):       the peptide sequence
  hla (str):           the HLA allele
  binding_affinity_nM (float): predicted MHC binding affinity in nM (<50=strong)
  binder (bool):       true if predicted to bind
  immunogenicity_score (float): 0.0–1.0, probability of T-cell recognition
  rationale (str):     one sentence grounded in anchor motifs / TCR contact
                       preferences

Candidates (CSV):
{csv}

Return a JSON array, one object per candidate, in input order. No prose."""


@dataclass
class NeoantigenCall:
    peptide: str
    hla: str
    binding_affinity_nM: float
    binder: bool
    immunogenicity_score: float
    rationale: str
    source: str = "heuristic"  # or "llm"


@dataclass
class ScreenReport:
    hla: list[str]
    candidates: list[NeoantigenCall] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "hla": self.hla,
            "candidates": [asdict(c) for c in self.candidates],
        }


def screen_peptide_llm(
    peptide: str,
    hla: str,
    *,
    backend: str | None = None,
) -> NeoantigenCall:
    """Screen one peptide × one HLA via LLM (or fallback)."""
    peptide = peptide.upper().strip()
    hla = hla.strip()

    # Backend resolution: explicit > MRNA_AI_LLM_BACKEND > auto
    if backend is None:
        backend = os.environ.get("MRNA_AI_LLM_BACKEND", "auto").lower()
        if backend == "auto":
            if _mhcflurry_available():
                backend = "mhcflurry"
            elif _openai_available():
                backend = "openai"
            else:
                backend = "mock"

    if backend == "mhcflurry":
        return _screen_peptide_mhcflurry(peptide, hla)

    if backend == "mock" or backend is None:
        # heuristic
        if hla.startswith("HLA-A*02"):
            aff, binder, rat = _heuristic_a0201(peptide)
        else:
            aff, binder, rat = (9999.0, False, "no heuristic for non-A2 allele")
        return NeoantigenCall(
            peptide=peptide,
            hla=hla,
            binding_affinity_nM=aff,
            binder=binder,
            immunogenicity_score=0.5 if binder else 0.05,
            rationale=rat,
            source="heuristic",
        )

    if backend == "openai":
        from .llm import llm_json

        prompt = (
            f"peptide={peptide} hla={hla}\n\n"
            'Return JSON: {"peptide":..,"hla":..,"binding_affinity_nM":..,'
            '"binder":..,"immunogenicity_score":..,"rationale":..}'
        )
        data = llm_json(prompt, backend="openai")
        return NeoantigenCall(
            peptide=peptide,
            hla=hla,
            binding_affinity_nM=float(data.get("binding_affinity_nM", 9999)),
            binder=bool(data.get("binder", False)),
            immunogenicity_score=float(data.get("immunogenicity_score", 0.0)),
            rationale=str(data.get("rationale", "")),
            source="openai",
        )

    raise ValueError(f"unknown backend: {backend!r}")


def _openai_available() -> bool:
    import os

    return bool(os.environ.get("OPENAI_API_KEY"))


# ---------- mhcflurry backend (optional) -----------------------------------

_MHCFLURRY_PREDICTOR = None
_MHCFLURRY_AVAILABLE: bool | None = None


def _mhcflurry_available() -> bool:
    """Return True iff ``mhcflurry`` is importable (downloads models on first use)."""
    global _MHCFLURRY_AVAILABLE
    if _MHCFLURRY_AVAILABLE is None:
        try:
            import mhcflurry  # noqa: F401

            _MHCFLURRY_AVAILABLE = True
        except Exception:
            _MHCFLURRY_AVAILABLE = False
    return _MHCFLURRY_AVAILABLE


def _get_mhcflurry_predictor():
    global _MHCFLURRY_PREDICTOR
    if _MHCFLURRY_PREDICTOR is None:
        from mhcflurry import Class1PresentationPredictor

        _MHCFLURRY_PREDICTOR = Class1PresentationPredictor.load()
    return _MHCFLURRY_PREDICTOR


def _screen_peptide_mhcflurry(peptide: str, hla: str) -> NeoantigenCall:
    """Real binding-affinity prediction via mhcflurry Class1PresentationPredictor."""
    if not _mhcflurry_available():
        raise RuntimeError("mhcflurry not installed")
    predictor = _get_mhcflurry_predictor()
    # Class1PresentationPredictor returns a DataFrame with `affinity` and `presentation_score`.
    df = predictor.predict(
        peptides=[peptide],
        alleles=[hla],
        verbose=0,
    )
    row = df.iloc[0]
    affinity = float(row["affinity"])
    presentation = float(row.get("presentation_score", 0.0))
    binder = affinity < 500
    return NeoantigenCall(
        peptide=peptide,
        hla=hla,
        binding_affinity_nM=round(affinity, 1),
        binder=binder,
        immunogenicity_score=round(min(1.0, presentation), 3),
        rationale=(
            f"mhcflurry predicted presentation_score={presentation:.3f}, affinity={affinity:.0f} nM"
        ),
        source="mhcflurry",
    )


def screen_csv(
    csv_path: str | Path,
    *,
    hla_alleles: Iterable[str],
    backend: str | None = None,
) -> ScreenReport:
    """Screen all peptides in a CSV.

    CSV must have a column ``peptide``. The ``variant`` column (if present) is
    preserved in the rationale.
    """
    hla_list = list(hla_alleles)
    report = ScreenReport(hla=hla_list)
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pep = row["peptide"].upper().strip()
            variant = row.get("variant", "")
            for hla in hla_list:
                call = screen_peptide_llm(pep, hla, backend=backend)
                if variant and variant not in call.rationale:
                    call.rationale = f"{variant}: {call.rationale}"
                report.candidates.append(call)
    return report


# ---------- Protein-LM immunogenicity scoring ------------------------------


def lm_immunogenicity_score(
    peptide: str,
    *,
    embedder: object | None = None,
    model_id: str = "facebook/esm2_t12_35M_UR50D",
) -> dict:
    """Score peptide immunogenicity via frozen protein-LM embeddings.

    Implements the Wong et al. 2025 Applm pattern
    (arXiv 2508.10541): use a frozen protein LM to embed the
    candidate peptide, then a lightweight downstream classifier
    scores it for immunogenicity.

    Default backend is the stdlib mock embedder (no torch dep).
    Pass ``embedder=ESM2Embedder(model_id=...)`` to use the real
    ESM2 model.

    Parameters
    ----------
    peptide
        Amino-acid sequence (9-11 mer typical for MHC-I epitopes).
    embedder
        Optional :class:`ProteinLMEmbedder`. When None, loads via
        :func:`select_protein_lm_embedder` which picks mock when
        transformers is unavailable.
    model_id
        HuggingFace model identifier. Ignored when ``embedder`` is
        supplied.

    Returns
    -------
    dict with keys:
        - ``score``: float in [0, 1], higher = more immunogenic
        - ``dim``: embedding dimension used
        - ``model_id``: model that produced the embedding
        - ``backend``: ``"transformers"`` or ``"mock"``
        - ``elapsed_seconds``: wall-clock time
        - ``notes``: tuple of debug notes
    """
    from .protein_lm_adapter import (
        ApplmStyleClassifier,
        select_protein_lm_embedder,
    )
    from .protein_lm_protocols import EmbeddingRequest

    if embedder is None:
        embedder = select_protein_lm_embedder(model_id=model_id)
    req = EmbeddingRequest(sequences=(peptide,), model_id=model_id)
    emb = embedder.embed(req)
    clf = ApplmStyleClassifier(embedder=embedder)
    score = clf.score(emb.embeddings[0])
    return {
        "score": score,
        "dim": emb.dim,
        "model_id": emb.model_id,
        "backend": emb.backend,
        "elapsed_seconds": emb.elapsed_seconds,
        "notes": emb.notes,
    }


# ---------- CLI ------------------------------------------------------------


def _run_cli(argv: list[str]) -> int:
    import argparse
    import json as _json

    p = argparse.ArgumentParser(prog="mrnavax neoantigen")
    p.add_argument("--variants", required=True, help="CSV with a 'peptide' column")
    p.add_argument("--hla", action="append", required=True, help="repeatable, e.g. HLA-A*02:01")
    p.add_argument("--backend", choices=["auto", "mock", "openai", "mhcflurry"], default="auto")
    p.add_argument("--out")
    args = p.parse_args(argv)

    # In CI / sandboxed environments, mhcflurry is unavailable and we want
    # the deterministic mock path. Honor MRNA_AI_FORCE_MOCK=1 to override
    # the auto-detect preference for mhcflurry.
    backend = None if args.backend == "auto" else args.backend
    if backend is None and os.environ.get("MRNA_AI_FORCE_MOCK") == "1":
        backend = "mock"

    rep = screen_csv(args.variants, hla_alleles=args.hla, backend=backend)
    out_text = _json.dumps(rep.to_dict(), indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_run_cli(sys.argv[1:]))


__all__ = [
    # Reference data (heuristic HLA-A*02:01 anchor matrix).
    "A0201_ANCHORS",
    # Dataclasses.
    "NeoantigenCall",
    "ScreenReport",
    # Public API.
    "screen_peptide_llm",
    "screen_csv",
    "lm_immunogenicity_score",
]

