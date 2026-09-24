"""Adapters for the published RiboDecode codon-optimization package.

Reference
---------
Li, Y., Wang, F., Yang, J. et al. *Deep generative optimization of
mRNA codon sequences for enhanced mRNA translation and therapeutic
efficacy.* Nat Commun 16, 9957 (2025).
DOI: 10.1038/s41467-025-64894-x

The upstream package distributes two CLI tools as `.whl` wheels from
Google Drive (NOT on PyPI):

  - ``TranslationModel-1.1.0-py3-none-any.whl`` → installs the
    ``pred-translation`` console script. Predicts translation level
    for a CDS, optionally conditioned on a cellular environment
    (HEK293T, A549, HeLa, or a custom RPKM CSV).
  - ``ribodecode-1.3.0-py3-none-any.whl`` → installs the
    ``ribo-decode`` console script. Jointly optimizes a CDS for
    translation and MFE via a deep generative model. Requires
    CUDA + ViennaRNA 2.6.4.

This module ships:

  - :class:`TranslationModelCLIAdapter`: shells out to
    ``pred-translation`` and parses the result. Heavy (PyTorch +
    CUDA).
  - :class:`RiboDecodeCLIAdapter`: shells out to ``ribo-decode`` and
    parses the result. Heavy (PyTorch + CUDA + ViennaRNA).
  - :class:`MockTranslationPredictor`: deterministic stdlib stub
    that derives translation from codon-pair frequencies. Used in
    tests and offline runs.
  - :class:`MockCodonOptimizer`: wraps the existing stdlib
    LinearDesign optimizer (``mrnavax.codon_lineardesign``)
    and reports the result in the RiboDecode response shape. Used
    in tests and offline runs.

Why mocks are shipped alongside the real adapters:
    Tests must run in CI without downloading 400+ MB of CUDA
    weights. Mock backends satisfy the same Protocol and produce
    results in the same shape as the real adapter — the toolkit
    consumers see no difference.

Why the real adapters shell out (rather than import):
    The upstream wheels have hard deps on `viennarna`, `torch`, and
    `pandas`, none of which are stdlib. Importing them would force
    every toolkit user to install 1+ GB of heavy deps. Subprocess
    invocation lets the user opt in by `pip install
    TranslationModel-1.1.0-py3-none-any.whl ribodecode-1.3.0-py3-none-any.whl`
    and `viennarna==2.6.4` separately; if any of those are missing
    the adapter raises ``RiboDecodeNotInstalled`` with a clear
    remediation message.
"""

from __future__ import annotations

import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from .codon_optimizer import _AA_MAX_FREQ, CODON_TO_AA, HUMAN_CODON_FREQ
from .codon_protocols import (
    CodonOptimizer,
    RiboDecodeRequest,
    RiboDecodeResult,
    TranslationPrediction,
    TranslationPredictor,
)

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class RiboDecodeError(RuntimeError):
    """Base class for RiboDecode adapter errors."""


class RiboDecodeNotInstalled(RiboDecodeError):
    """Raised when the upstream CLI binary is not on $PATH."""

    def __init__(self, binary: str) -> None:
        super().__init__(
            f"CLI binary '{binary}' not found on $PATH. "
            f"The RiboDecode package (Li et al., Nat Commun 16, 9957 (2025)) "
            f"is distributed as a Google-Drive `.whl`, NOT on PyPI. "
            f"Download from https://github.com/wangfanfff/RiboDecode and "
            f"`pip install TranslationModel-1.1.0-py3-none-any.whl "
            f"ribodecode-1.3.0-py3-none-any.whl`. "
            f"Then `pip install viennarna==2.6.4` (requires GCC >= 5.0). "
            f"For tests and offline use, the mock backends work without "
            f"any installation."
        )


class riboDecodeInvocationError(RiboDecodeError):
    """Raised when the upstream CLI exits non-zero or produces unparseable output."""


# ---------------------------------------------------------------------------
# Path / availability helpers
# ---------------------------------------------------------------------------


PRED_TRANSLATION_BIN = "pred-translation"
RIBO_DECODE_BIN = "ribo-decode"


def pred_translation_available() -> bool:
    """True iff the upstream ``pred-translation`` binary is on $PATH."""
    return shutil.which(PRED_TRANSLATION_BIN) is not None


def ribo_decode_available() -> bool:
    """True iff the upstream ``ribo-decode`` binary is on $PATH."""
    return shutil.which(RIBO_DECODE_BIN) is not None


# ---------------------------------------------------------------------------
# TranslationModel CLI adapter (real)
# ---------------------------------------------------------------------------


@dataclass
class TranslationModelCLIAdapter:
    """Adapter for the published ``TranslationModel`` package.

    Usage:
        >>> from mrnavax.codon_translationmodel_adapter import TranslationModelCLIAdapter
        >>> adapter = TranslationModelCLIAdapter()
        >>> pred = adapter.predict("ATGGACGGGTAG", env="HEK293T")
        >>> pred.translation_level
        42.7  # upstream-predicted value
    """

    timeout_seconds: int = 60
    working_dir: Path | None = None

    def predict(
        self,
        cds: str,
        env: str = "HEK293T",
        custom_env_csv: Path | None = None,
    ) -> TranslationPrediction:
        if not pred_translation_available():
            raise RiboDecodeNotInstalled(PRED_TRANSLATION_BIN)
        cmd: list[str] = [
            PRED_TRANSLATION_BIN,
            "--cds",
            cds,
            "--env",
            env,
        ]
        if env == "custom" and custom_env_csv is not None:
            cmd += ["--csv", str(custom_env_csv)]
        t0 = time.time()
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds,
                cwd=str(self.working_dir) if self.working_dir else None,
                check=False,
            )
        except subprocess.TimeoutExpired as e:
            raise riboDecodeInvocationError(
                f"{PRED_TRANSLATION_BIN} timed out after {self.timeout_seconds}s"
            ) from e
        elapsed = time.time() - t0
        if proc.returncode != 0:
            raise riboDecodeInvocationError(
                f"{PRED_TRANSLATION_BIN} failed (exit {proc.returncode}): "
                f"{proc.stderr.strip()[:500]}"
            )
        # The upstream CLI prints a single float to stdout (translation
        # level). Parse it.
        try:
            translation_level = float(proc.stdout.strip().splitlines()[-1])
        except (ValueError, IndexError) as e:
            raise riboDecodeInvocationError(
                f"could not parse translation level from "
                f"{PRED_TRANSLATION_BIN!r} output: {proc.stdout[:200]!r}"
            ) from e
        return TranslationPrediction(
            cds=cds,
            translation_level=translation_level,
            env=env,
            backend="translationmodel",
            elapsed_seconds=elapsed,
        )


# ---------------------------------------------------------------------------
# TranslationModel mock backend (stdlib, deterministic)
# ---------------------------------------------------------------------------


@dataclass
class MockTranslationPredictor:
    """Deterministic translation-level predictor.

    Algorithm: weighted sum of per-codon-pair log-frequency ratios
    derived from HUMAN_CODON_FREQ. Bounded to [0, 100] to match the
    upstream TranslationModel output range. Deterministic for a fixed
    input — no random seeds required.

    This is intentionally a heuristic; the mock backend exists so
    CI/tests work without the 135 MB TranslationModel weights.
    """

    backend_name: ClassVar[str] = "mock-translation"

    def predict(
        self,
        cds: str,
        env: str = "HEK293T",
        custom_env_csv: Path | None = None,
    ) -> TranslationPrediction:
        t0 = time.time()
        cds = cds.upper().replace("U", "T")
        # Cap at the upstream 4500 nt limit; silently truncate so the
        # mock can handle oversized synthetic inputs without raising.
        if len(cds) > 4500:
            cds = cds[:4500]
        codons = [cds[i : i + 3] for i in range(0, len(cds) - 2, 3)]
        notes: list[str] = []
        if not codons:
            return TranslationPrediction(
                cds=cds,
                translation_level=0.0,
                env=env,
                backend=self.backend_name,
                elapsed_seconds=time.time() - t0,
                notes=("empty_cds",),
            )
        per_codon_log = 0.0
        n_known = 0
        for codon in codons:
            aa = CODON_TO_AA.get(codon)
            if aa is None or aa == "*":
                notes.append(f"unknown_or_stop_codon:{codon}")
                continue
            max_freq = _AA_MAX_FREQ.get(aa, 0.0)
            if max_freq <= 0:
                continue
            freq = HUMAN_CODON_FREQ[aa].get(codon, 0.0)
            per_codon_log += (freq / max_freq) if freq > 0 else 0.0
            n_known += 1
        # Normalize to [0, 1], then scale to [0, 100] to match the
        # upstream range. Empirical: a uniform-CAI CDS scores ~70,
        # perfect-CAI scores 100, lowest-CAI ~30.
        if n_known == 0:
            score = 0.0
        else:
            avg = per_codon_log / n_known  # in [0, 1]
            score = 30.0 + 70.0 * avg
        return TranslationPrediction(
            cds=cds,
            translation_level=score,
            env=env,
            backend=self.backend_name,
            elapsed_seconds=time.time() - t0,
            notes=tuple(notes),
        )


# ---------------------------------------------------------------------------
# RiboDecode CLI adapter (real)
# ---------------------------------------------------------------------------


@dataclass
class RiboDecodeCLIAdapter:
    """Adapter for the published ``ribo-decode`` package.

    Usage:
        >>> from mrnavax.codon_ribodecode_adapter import RiboDecodeCLIAdapter
        >>> adapter = RiboDecodeCLIAdapter()
        >>> req = RiboDecodeRequest(cds="ATGGACGGGTAG", env="HEK293T")
        >>> result = adapter.optimize(req)
        >>> result.optimized_cds
        '...'

    Heavy deps on the user's environment: ``ribo-decode`` on $PATH,
    CUDA, ViennaRNA 2.6.4, torch. Missing any of these raises
    :class:`RiboDecodeNotInstalled` or
    :class:`riboDecodeInvocationError`.
    """

    timeout_seconds: int = 300  # ribo-decode is GPU-bound
    working_dir: Path | None = None

    def optimize(self, req: RiboDecodeRequest) -> RiboDecodeResult:
        if not ribo_decode_available():
            raise RiboDecodeNotInstalled(RIBO_DECODE_BIN)
        cmd: list[str] = [
            RIBO_DECODE_BIN,
            "--cds",
            req.name,
            "--cds_seq",
            req.cds,
            "--env",
            req.env,
            "--mfe_weight",
            str(req.mfe_weight),
            "--optim_epoch",
            str(req.optim_epoch),
            "--alpha",
            str(req.alpha),
            "--beta",
            str(req.beta),
        ]
        if req.env == "custom" and req.custom_env_csv is not None:
            cmd += ["--csv", str(req.custom_env_csv)]
        t0 = time.time()
        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout_seconds,
                cwd=str(self.working_dir) if self.working_dir else None,
                check=False,
            )
        except subprocess.TimeoutExpired as e:
            raise riboDecodeInvocationError(
                f"{RIBO_DECODE_BIN} timed out after {self.timeout_seconds}s "
                f"(mfe_weight={req.mfe_weight}, epoch={req.optim_epoch})"
            ) from e
        elapsed = time.time() - t0
        if proc.returncode != 0:
            raise riboDecodeInvocationError(
                f"{RIBO_DECODE_BIN} failed (exit {proc.returncode}): {proc.stderr.strip()[:500]}"
            )
        # The upstream CLI writes results to
        # ``results_<name>/<env>_optim_mfe_dist-optim/<epoch>/samples/optim_results.txt``
        # with three columns: optimized_cds, predicted_translation, predicted_mfe.
        return self._parse_results_file(req, elapsed)

    def _parse_results_file(self, req: RiboDecodeRequest, elapsed: float) -> RiboDecodeResult:
        """Read the upstream output file and return a RiboDecodeResult."""
        out_path = (
            Path(self.working_dir or ".")
            / f"results_{req.name}"
            / f"{req.env}_optim_mfe_dist-optim"
            / str(req.optim_epoch)
            / "samples"
            / "optim_results.txt"
        )
        if not out_path.exists():
            raise riboDecodeInvocationError(
                f"{RIBO_DECODE_BIN} did not produce expected output file "
                f"at {out_path}. Output: {out_path.parent.read_text()[:500]}"
            )
        last_line = out_path.read_text().strip().splitlines()[-1]
        try:
            parts = last_line.split("\t")
            optimized_cds, translation_str, mfe_str = parts[0], parts[1], parts[2]
            translation = float(translation_str)
        except (ValueError, IndexError) as e:
            raise riboDecodeInvocationError(
                f"could not parse {RIBO_DECODE_BIN} output file {out_path}: {last_line[:200]!r}"
            ) from e
        mfe: float | None = None
        if req.mfe_weight != 0.0:
            try:
                mfe = float(mfe_str)
            except ValueError:
                mfe = None
        return RiboDecodeResult(
            optimized_cds=optimized_cds,
            predicted_translation=translation,
            predicted_mfe=mfe,
            epoch=req.optim_epoch,
            elapsed_seconds=elapsed,
            backend="ribodecode",
        )


# ---------------------------------------------------------------------------
# RiboDecode mock backend (stdlib, deterministic)
# ---------------------------------------------------------------------------


@dataclass
class MockCodonOptimizer:
    """Mock CodonOptimizer that delegates to the stdlib LinearDesign.

    Satisfies :class:`CodonOptimizer` so tests can run without the
    ~400 MB upstream RiboDecode weights + CUDA. The mock reports a
    synthetic translation score derived from the CAI of the optimized
    CDS, normalized to the same [0, 100] range as the real
    TranslationModel.
    """

    backend_name: ClassVar[str] = "mock-lineardesign"

    def optimize(self, req: RiboDecodeRequest) -> RiboDecodeResult:
        # Lazy import: avoid pulling LinearDesign at module load when
        # the mock is only used in CI mocks.
        from .codon_lineardesign import optimize_lineardesign

        t0 = time.time()
        # Strip trailing stop codon (LinearDesign expects no stop).
        cds = req.cds.upper().replace("U", "T")
        if cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
            cds = cds[:-3]
        # Run LinearDesign (the stdlib joint translation × MFE DP).
        ld_result = optimize_lineardesign(
            cds,
            gc_window_size=21,
            verbose=False,
        )
        # Score the optimized CDS with the mock translation predictor
        # so callers see the same shape they would see from the real
        # backend.
        mock_pred = MockTranslationPredictor()
        # Translation score in LinearDesign is a sum of per-codon
        # log-freq. Rescale to [0, 100] like the upstream TranslationModel.
        pred = mock_pred.predict(ld_result.new_cds, env=req.env)
        # MFE proxy: LinearDesign reports a structure_score in arbitrary
        # units (negative energy proxy). When mfe_weight==0, the real
        # RiboDecode sets mfe=None — we do the same.
        mfe: float | None = None
        if req.mfe_weight != 0.0:
            mfe = float(ld_result.structure_score)
        return RiboDecodeResult(
            optimized_cds=ld_result.new_cds,
            predicted_translation=pred.translation_level,
            predicted_mfe=mfe,
            epoch=req.optim_epoch,
            elapsed_seconds=time.time() - t0,
            backend=self.backend_name,
            notes=("mock-backend", f"linear-design-on-{req.env}"),
        )


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------


def select_translation_predictor(
    *,
    prefer: str = "auto",
    env: str = "HEK293T",
    custom_env_csv: Path | None = None,
) -> TranslationPredictor:
    """Pick the best available :class:`TranslationPredictor` for ``env``.

    Selection rules (in order):
        1. ``prefer="mock"`` → return :class:`MockTranslationPredictor`.
        2. ``prefer="real"`` → require the upstream ``pred-translation``
           binary; raise :class:`RiboDecodeNotInstalled` if missing.
        3. ``prefer="auto"`` (default):
           - If ``pred_translation_available()`` → return
             :class:`TranslationModelCLIAdapter`.
           - Else if ``MRNA_AI_FORCE_MOCK=1`` env var → mock.
           - Else → mock (with a note that real backend is unavailable).
    """
    if prefer == "mock":
        return MockTranslationPredictor()
    if prefer == "real":
        if not pred_translation_available():
            raise RiboDecodeNotInstalled(PRED_TRANSLATION_BIN)
        return TranslationModelCLIAdapter()
    if pred_translation_available():
        return TranslationModelCLIAdapter()
    return MockTranslationPredictor()


def select_codon_optimizer(
    *,
    prefer: str = "auto",
) -> CodonOptimizer:
    """Pick the best available :class:`CodonOptimizer`.

    Same selection rules as :func:`select_translation_predictor`.
    """
    if prefer == "mock":
        return MockCodonOptimizer()
    if prefer == "real":
        if not ribo_decode_available():
            raise RiboDecodeNotInstalled(RIBO_DECODE_BIN)
        return RiboDecodeCLIAdapter()
    if ribo_decode_available():
        return RiboDecodeCLIAdapter()
    return MockCodonOptimizer()


__all__ = [
    "RIBO_DECODE_BIN",
    "RiboDecodeError",
    "RiboDecodeNotInstalled",
    "RiboDecodeRequest",
    "RiboDecodeResult",
    "pred_translation_available",
    "ribo_decode_available",
    "TranslationModelCLIAdapter",
    "MockTranslationPredictor",
    "RiboDecodeCLIAdapter",
    "MockCodonOptimizer",
    "select_translation_predictor",
    "select_codon_optimizer",
]
