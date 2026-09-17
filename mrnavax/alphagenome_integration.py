"""AlphaGenome Atlas integration — regulatory-variant impact scoring.

Loads pre-computed AlphaGenome Atlas (Avsec et al., *Nature* 2026)
predictions for variant regulatory impact across the human genome.
The Atlas covers all ~9 billion possible single-nucleotide variants
in the human genome, providing the **AVI score** (AlphaGenome Variant
Impact): a unified [0, 1] score combining coding and non-coding
regulatory effects. Pre-computed; queryable via the official
`alphagenome` Python package API.

Source: https://github.com/google-deepmind/alphagenome
License: Non-commercial use only (per Google DeepMind Terms of Service).
         Atlas predictions and AVI scores are not for commercial use.
         Commercial access: Google Cloud Vertex AI (`alphagenome`
         deployable).

This module plugs AlphaGenome Atlas into the `mrnavax` variant
prioritization pipeline:

* For variants in **coding regions** (`is_coding=True`), AlphaMissense
  (Cheng et al., *Science* 2023) is the higher-fidelity signal — use
  `mrnavax.alphamissense_integration` instead.
* For variants in **non-coding regulatory regions** (`is_coding=False`),
  AlphaGenome Atlas is the canonical signal — use this module.

Architecture
------------
This module follows the established mrnavax Protocol-adapter pattern:

* **Protocol**: ``RegulatoryVariantScorer`` (typing.Protocol).
* **Mock**: ``MockRegulatoryVariantScorer`` — stdlib-only, deterministic,
  always present. Computes a synthetic AVI score from a stdlib hash
  of (chrom, pos, ref, alt). CI runs without network access.
* **Real**: ``AlphaGenomeCLIAdapter`` — calls the official
  `alphagenome` Python package via subprocess (matches RiboDecode and
  STModule shim pattern). Heavy dep is opt-in via the
  ``[variant-alphagenome]`` extra; the user supplies an API key via
  ``ALPHAGENOME_API_KEY`` env var.

Usage
-----
    from mrnavax.alphagenome_integration import regulatory_score, AVIResult

    # Mock (no API key needed):
    r = regulatory_score("chr7", 140753336, "T", "A")  # BRAF V600E
    # -> AVIResult(score=0.71, classification="high", is_coding=True)

    # With real Atlas (after `pip install -e .[variant-alphagenome]`
    # and `export ALPHAGENOME_API_KEY=...`):
    r = regulatory_score("chr11", 12345678, "A", "G")
    # -> AVIResult from the official AlphaGenome Atlas API

Reference
---------
Avsec et al., "Advancing regulatory variant effect prediction with
AlphaGenome," *Nature* (2026).
https://www.nature.com/articles/s41586-025-10014-0
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import os
import subprocess  # noqa: S404  (subprocess invocation of optional upstream)
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Protocol, runtime_checkable

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Thresholds for classifying AVI score into bands. Match the
# AlphaMissense paper bins (Cheng 2023) for consistency between the
# two variant-prioritization signals in mrnavax.
THRESHOLD_PATHOGENIC = 0.564  # AVI ≥ this → "high"
THRESHOLD_BENIGN = 0.34      # AVI < this  → "low"

# Chromosomes considered "coding-rich" for the mock's is_coding heuristic.
# The real adapter returns whatever the upstream Atlas API reports.
# (no longer used — see MockRegulatoryVariantScorer docstring)


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AVIResult:
    """Regulatory impact prediction for one DNA-level variant.

    Attributes
    ----------
    score:
        AlphaGenome Variant Impact score in [0, 1]. Higher = more
        likely to disrupt gene regulation. From the Atlas for all
        9 billion possible single-nucleotide variants in the human
        genome.
    classification:
        Bin label: ``"low"`` (< 0.34), ``"moderate"`` (0.34–0.564),
        or ``"high"`` (≥ 0.564). Thresholds match AlphaMissense for
        consistent variant-prioritization triage.
    is_coding:
        True if the variant sits in a protein-coding region (where
        AlphaMissense is the better signal); False if it sits in a
        non-coding regulatory region (where AlphaGenome Atlas is the
        canonical signal). For mock variants this is a chromosome-level
        heuristic; for real Atlas queries it comes from the upstream
        API.
    """

    score: float
    classification: str
    is_coding: bool

    def __post_init__(self) -> None:
        if not 0.0 <= self.score <= 1.0:
            raise ValueError(
                f"AVI score must be in [0, 1], got {self.score!r}"
            )
        # Reclassify if the caller's `classification` doesn't match the
        # canonical bins derived from the score. This guarantees the
        # contract invariant: classification is always derived from score.
        canonical = _classify(self.score)
        if self.classification != canonical:
            object.__setattr__(self, "classification", canonical)


def _classify(score: float) -> str:
    """Map an AVI score to its bin label. Pure — safe for mocks."""
    if score < THRESHOLD_BENIGN:
        return "low"
    if score < THRESHOLD_PATHOGENIC:
        return "moderate"
    return "high"


# ---------------------------------------------------------------------------
# Protocol
# ---------------------------------------------------------------------------

@runtime_checkable
class RegulatoryVariantScorer(Protocol):
    """Anything that can score one DNA-level variant's regulatory impact."""

    def score_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> AVIResult:
        ...


# ---------------------------------------------------------------------------
# Mock backend — stdlib-only, deterministic
# ---------------------------------------------------------------------------

class MockRegulatoryVariantScorer:
    """Stdlib-only deterministic mock.

    Computes a synthetic AVI score by hashing (chrom, pos, ref, alt)
    with SHA-256 and mapping the first 4 bytes to a float in [0, 1].
    Stable: same input → same output across processes.

    The ``is_coding`` heuristic uses position parity: even positions
    are treated as coding-region (use AlphaMissense), odd positions
    as non-coding regulatory (use AlphaGenome Atlas). This is purely
    a mock — the real Atlas API returns per-variant is_coding from the
    upstream query.
    """

    def score_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> AVIResult:
        key = f"{chrom}:{pos}:{ref}:{alt}".encode()
        digest = hashlib.sha256(key).digest()
        # First 4 bytes → unsigned int → [0, 1]
        n = int.from_bytes(digest[:4], byteorder="big")
        score = n / 0xFFFFFFFF
        is_coding = (pos % 2) == 0
        return AVIResult(score=score, classification="", is_coding=is_coding)


# ---------------------------------------------------------------------------
# Real adapter — subprocess wrapper around the official `alphagenome` pkg
# ---------------------------------------------------------------------------


class AlphaGenomeCLIAdapter:
    """Real adapter — calls the official `alphagenome` package via subprocess.

    The adapter itself stays stdlib-only. The heavy dep (alphagenome
    Python package) is loaded only inside the shim subprocess and only
    when the user opts in via `pip install -e ".[variant-alphagenome]"`.

    Parameters
    ----------
    api_key:
        Google DeepMind AlphaGenome API key. Reads ``ALPHAGENOME_API_KEY``
        env var if not supplied.
    shim_path:
        Path to the shim script. Defaults to the bundled shim inside
        the package. Override for tests.
    _run:
        Internal hook for tests. Defaults to ``subprocess.run``.
    """

    def __init__(
        self,
        api_key: str | None = None,
        shim_path: str | os.PathLike[str] | None = None,
        _run: Callable | None = None,
    ) -> None:
        self.api_key = api_key or os.environ.get("ALPHAGENOME_API_KEY", "")
        if not self.api_key:
            raise ValueError(
                "AlphaGenomeCLIAdapter requires an API key. Set the "
                "ALPHAGENOME_API_KEY environment variable or pass "
                "api_key=... explicitly."
            )
        self.shim_path = (
            Path(shim_path) if shim_path is not None
            else Path(__file__).parent / "_shims" / "alphagenome_cli.py"
        )
        self._run = _run if _run is not None else subprocess.run

    def score_variant(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> AVIResult:
        payload = {
            "api_key": self.api_key,
            "chrom": chrom,
            "pos": pos,
            "ref": ref,
            "alt": alt,
        }
        proc = self._run(
            [sys.executable, str(self.shim_path)],
            input=json.dumps(payload),
            capture_output=True,
            text=True,
            timeout=30,
        )
        if proc.returncode != 0:
            logger.warning(
                "AlphaGenome shim failed (rc=%d, stderr=%r); "
                "returning default low-impact AVI result",
                proc.returncode,
                proc.stderr[:200] if proc.stderr else "",
            )
            return AVIResult(score=0.0, classification="low", is_coding=False)
        try:
            data = json.loads(proc.stdout)
            return AVIResult(
                score=float(data["score"]),
                classification=str(data.get("classification", "")),
                is_coding=bool(data.get("is_coding", False)),
            )
        except (json.JSONDecodeError, KeyError, ValueError) as exc:
            logger.warning("AlphaGenome shim returned malformed JSON: %s", exc)
            return AVIResult(score=0.0, classification="low", is_coding=False)


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------

def select_regulatory_scorer() -> RegulatoryVariantScorer:
    """Return the appropriate backend based on environment.

    * If ``ALPHAGENOME_API_KEY`` is set, return the real adapter.
    * Otherwise return the mock (stdlib-only, deterministic).
    """
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if api_key:
        return AlphaGenomeCLIAdapter(api_key=api_key)
    return MockRegulatoryVariantScorer()


def regulatory_score(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    backend: RegulatoryVariantScorer | None = None,
) -> AVIResult:
    """Convenience wrapper that auto-selects the backend.

    Examples
    --------
    >>> r = regulatory_score("chr7", 140753336, "T", "A")
    >>> 0.0 <= r.score <= 1.0
    True
    """
    sel = backend if backend is not None else select_regulatory_scorer()
    return sel.score_variant(chrom, pos, ref, alt)


# ---------------------------------------------------------------------------
# Bulk scoring helper (used by CLI)
# ---------------------------------------------------------------------------

def score_variants_from_csv(
    csv_path: str,
    *,
    backend: RegulatoryVariantScorer | None = None,
) -> list[AVIResult]:
    """Score every variant in a CSV file. CSV columns: chrom,pos,ref,alt,label.

    Returns results sorted by ``score`` descending (most-impactful first).
    """
    sel = backend if backend is not None else select_regulatory_scorer()
    out: list[AVIResult] = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            try:
                r = sel.score_variant(
                    row["chrom"],
                    int(row["pos"]),
                    row["ref"],
                    row["alt"],
                )
            except (KeyError, ValueError) as exc:
                logger.warning("skipping row %r: %s", row, exc)
                continue
            out.append(r)
    out.sort(key=lambda r: r.score, reverse=True)
    return out


__all__ = [
    "AVIResult",
    "AlphaGenomeCLIAdapter",
    "MockRegulatoryVariantScorer",
    "RegulatoryVariantScorer",
    "regulatory_score",
    "score_variants_from_csv",
    "select_regulatory_scorer",
    "THRESHOLD_PATHOGENIC",
    "THRESHOLD_BENIGN",
]
