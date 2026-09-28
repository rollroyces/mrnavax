"""Protocol interfaces for codon optimization backends.

This module defines the typed contracts that any codon-optimization
backend must satisfy to plug into the toolkit. The actual implementations
live elsewhere (e.g. ``codon_ribodecode_adapter.py`` for the
CLI-shell-out adapter to the published RiboDecode package, and
``codon_optimizer.py`` for the stdlib LinearDesign + heuristic
implementations).

Two protocols are defined:

- :class:`TranslationPredictor`: predicts the translation level of a
  given coding sequence, optionally conditioned on a cellular
  environment. Maps to the published ``TranslationModel`` package
  (Li, Wang, Yang et al., *Nat Commun* 16, 9957 (2025)) when installed.

- :class:`CodonOptimizer`: jointly optimizes a CDS for translation and
  secondary-structure stability. Maps to the published ``ribo-decode``
  CLI when installed.

Both protocols are deliberately narrow — they expose only what the
toolkit needs to wire a backend into the ``codon`` CLI and into the
``manufacture`` post-optimization check pipeline. Heavier features
(batch optimization, custom env files, multi-epoch sampling) are
supported via constructor args on the implementing class.

Why a Protocol rather than ABC? The whole point is that we ship the
toolkit without ViennaRNA / torch / CUDA installed, but want users who
*do* install the official RiboDecode package to plug it in with a
two-line subclass — no inheritance contract to satisfy, just duck
typing that the type checker can verify.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RiboDecodeRequest:
    """Input contract for joint translation × MFE codon optimization.

    Mirrors the CLI flags of the published ``ribo-decode`` tool so the
    adapter is a thin pass-through. The string ``name`` becomes the
    output directory name (``results_<name>/...``) in the upstream
    package; the toolkit ignores it and only consumes the optimized
    CDS string.
    """

    cds: str
    env: str = "HEK293T"  # "HEK293T" | "A549" | "HeLa" | "custom"
    mfe_weight: float = 0.0  # 0=translation-only, 1=MFE-only
    optim_epoch: int = 10
    alpha: float = 100.0
    beta: float = 100.0
    custom_env_csv: Path | None = None
    name: str = "example"

    def __post_init__(self) -> None:
        # Defensive validation — these are caught at the public-API
        # boundary so a downstream adapter never has to worry about
        # them. Mirrors the upstream CLI's documented constraints.
        if not self.cds:
            raise ValueError("cds must be non-empty")
        if len(self.cds) % 3 != 0:
            raise ValueError(f"cds length {len(self.cds)} is not a multiple of 3")
        if len(self.cds) > 4500:
            raise ValueError(
                f"cds length {len(self.cds)} exceeds 4500 nt cap (upstream RiboDecode limit)"
            )
        if not 0.0 <= self.mfe_weight <= 1.0:
            raise ValueError(f"mfe_weight {self.mfe_weight} not in [0, 1]")
        if self.env == "custom" and self.custom_env_csv is None:
            raise ValueError(
                "env='custom' requires custom_env_csv (a path to a "
                "CSV with gene-ID, RPKM columns per upstream docs)"
            )

    def to_dict(self) -> dict:
        d = asdict(self)
        if d["custom_env_csv"] is not None:
            d["custom_env_csv"] = str(d["custom_env_csv"])
        return d


@dataclass(frozen=True)
class RiboDecodeResult:
    """Output of a single codon-optimization run."""

    optimized_cds: str
    predicted_translation: float
    predicted_mfe: float | None  # None when mfe_weight == 0
    epoch: int
    elapsed_seconds: float
    backend: str  # "ribodecode" | "mock" | "heuristic" | "lineardesign"
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "optimized_cds": self.optimized_cds,
            "predicted_translation": self.predicted_translation,
            "predicted_mfe": self.predicted_mfe,
            "epoch": self.epoch,
            "elapsed_seconds": round(self.elapsed_seconds, 3),
            "backend": self.backend,
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class TranslationPrediction:
    """Output of a translation-level prediction."""

    cds: str
    translation_level: float
    env: str
    backend: str  # "translationmodel" | "mock"
    elapsed_seconds: float
    notes: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return {
            "cds": self.cds,
            "translation_level": round(self.translation_level, 4),
            "env": self.env,
            "backend": self.backend,
            "elapsed_seconds": round(self.elapsed_seconds, 3),
            "notes": list(self.notes),
        }


# ---------------------------------------------------------------------------
# Protocol interfaces
# ---------------------------------------------------------------------------


@runtime_checkable
class TranslationPredictor(Protocol):
    """Predict translation level for a CDS, optionally per cellular env.

    Implementations:
        - ``MockTranslationPredictor``: deterministic stdlib stub for
          tests and offline runs. Returns a simple CAI-derived score.
        - ``TranslationModelCLIAdapter``: shells out to the upstream
          ``pred-translation`` CLI when the ``ribodecode`` wheel is
          installed.

    Contract: ``predict()`` must never raise on input that's a valid
    DNA string (length a multiple of 3). On a malformed CDS (unknown
    bases, ambiguous IUPAC, internal stop), implementations should
    return a result with ``notes`` explaining what was rejected
    rather than raising — the toolkit's codon CLI aggregates over
    many candidates and one bad input should not abort the whole
    batch.
    """

    def predict(
        self,
        cds: str,
        env: str = "HEK293T",
        custom_env_csv: Path | None = None,
    ) -> TranslationPrediction: ...


@runtime_checkable
class CodonOptimizer(Protocol):
    """Joint translation × secondary-structure codon optimizer.

    Implementations:
        - ``MockCodonOptimizer``: deterministic stdlib stub. Returns
          the LinearDesign result with a synthetic translation score
          derived from CAI.
        - ``RiboDecodeCLIAdapter``: shells out to the upstream
          ``ribo-decode`` CLI when the ``ribodecode`` wheel is
          installed. Heavy: requires ViennaRNA + CUDA.

    Contract: ``optimize()`` must preserve the protein sequence. The
    output CDS, when translated with the standard genetic code,
    must equal the translation of the input CDS (excluding any
    trailing stop codon).
    """

    def optimize(self, req: RiboDecodeRequest) -> RiboDecodeResult: ...


__all__ = [
    # Input/output dataclasses.
    "RiboDecodeRequest",
    "RiboDecodeResult",
    "TranslationPrediction",
    # Protocol interfaces.
    "TranslationPredictor",
    "CodonOptimizer",
]

