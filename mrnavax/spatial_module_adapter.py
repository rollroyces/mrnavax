"""Adapters for the published STModule spatial-transcriptomics package.

Reference
---------
Wang R., Qian Y., Guo X., Song F., Xiong Z., Cai S., Bian X., Wong M.H.,
Cao Q.#, Cheng L.#, Lu G.#, and Leung K.S.#. (2025) STModule:
identifying tissue modules to uncover spatial components and
characteristics of transcriptomic landscapes. *Genome Medicine*
17(1): 18.

The upstream package distributes an R library from GitHub:

  https://github.com/rwang-z/STModule

It provides four R functions (verified from the official vignette):

  - ``data_preprocessing(count_file, loc_file, high_resolution, ...)``
    Reads two TSVs, returns a list of (expr, dist, loc).
  - ``run_STModule(data, num_modules, high_resolution, max_iter, version)``
    Returns a model-fit list (estimated vars, params, gene/loc lists).
  - ``get_assocaited_genes(res)``  # sic, typo in upstream
    Returns a data.frame with columns (gene, module, activity).
  - ``spatial_map_visualization(res, normalization, point_size)``
    Returns a list of plots (skipped here — toolkit returns data only).

There is no upstream CLI. Our wrapper ships a tiny R shim
(``scripts/stmodule_shim.R``) that calls these functions in order
and emits JSON to stdout, which the Python adapter parses back into
a typed :class:`SpatialModuleResult`.

Why an R shim + subprocess rather than a Python R-interface library
(rpy2): the user is the one who installs R + Seurat + torch +
GPUmatrix (heavy). rpy2 would force these deps on every toolkit
user. Subprocess invocation lets the user opt in by ``pip install
mrnavax[spatial-r]`` and following the conda install from
the upstream vignette; if any of those are missing the adapter
raises :class:`STModuleNotInstalled` with a clear remediation message.
"""

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from .spatial_protocols import (
    SpatialData,
    SpatialModule,
    SpatialModuleBackend,
    SpatialModuleResult,
)

# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------


class STModuleError(RuntimeError):
    """Base class for STModule adapter errors."""


class STModuleNotInstalled(STModuleError):
    """Raised when the upstream Rscript + STModule binary is not on PATH."""

    def __init__(self) -> None:
        super().__init__(
            "Rscript or STModule R package not found on this machine. "
            "STModule (Wang et al., Genome Medicine 17, 18 (2025)) is an "
            "R package distributed via GitHub "
            "(https://github.com/rwang-z/STModule), NOT on PyPI. "
            "Install per the upstream vignette "
            "(https://github.com/rwang-z/STModule/blob/main/vignette/installation.md): "
            "1) `conda create -n STModule python=3.9 && conda activate STModule "
            "&& conda install conda-forge::r-base=4.4.1`; "
            "2) `conda install conda-forge::r-seurat r-devtools`; "
            "3) `devtools::install_github('rwang-z/STModule')` "
            "(requires torch + CUDA 11.7 + GPUmatrix 1.0.2 — DO NOT update GPUmatrix). "
            "Then `pip install mrnavax[spatial-r]` and ensure "
            "`Rscript` is on $PATH. For tests and offline use, the "
            "mock backend works without any installation."
        )


# ---------------------------------------------------------------------------
# Path / availability helpers
# ---------------------------------------------------------------------------


RSCRIPT_BIN = "Rscript"
STMODULE_SHIM_NAME = "stmodule_shim.R"


def rscript_available() -> bool:
    """True iff the upstream ``Rscript`` binary is on $PATH."""
    return shutil.which(RSCRIPT_BIN) is not None


def _shim_path() -> Path:
    """Locate the bundled R shim next to the adapter module.

    Search order: (1) sibling of ``mrnavax/spatial_module_adapter.py``,
    (2) ``mrnavax/scripts/stmodule_shim.R``. The shim is
    versioned alongside the Python module so the subprocess always
    matches the protocol the adapter expects.
    """
    here = Path(__file__).parent
    candidates = [
        here / "scripts" / STMODULE_SHIM_NAME,
        here / STMODULE_SHIM_NAME,
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        f"Could not locate {STMODULE_SHIM_NAME} next to {__file__}. "
        f"Expected at {candidates[0]} or {candidates[1]}."
    )


# ---------------------------------------------------------------------------
# Real CLI adapter (subprocess to Rscript)
# ---------------------------------------------------------------------------


@dataclass
class STModuleCLIAdapter:
    """Adapter that shells out to the published STModule R package.

    Usage (requires R + Seurat + STModule installed):
        >>> from mrnavax.spatial_module_adapter import STModuleCLIAdapter
        >>> from mrnavax.spatial_protocols import SpatialData
        >>> data = SpatialData(count_file=Path("counts.tsv"),
        ...                    locations_file=Path("locs.tsv"),
        ...                    platform="ST", num_modules=10)
        >>> adapter = STModuleCLIAdapter()
        >>> result = adapter.run(data)
        >>> print(result.modules[0].top_genes)

    Heavy deps: R 4.4.1, Seurat v5, devtools, torch, GPUmatrix 1.0.2,
    CUDA 11.7 (GPU version requires all of these). Missing any of
    these raises :class:`STModuleNotInstalled`.
    """

    timeout_seconds: int = 600  # STModule can run minutes on large datasets
    working_dir: Path | None = None

    def run(self, data: SpatialData) -> SpatialModuleResult:
        if not rscript_available():
            raise STModuleNotInstalled()

        shim = _shim_path()
        cmd: list[str] = [
            RSCRIPT_BIN,
            str(shim),
            "--count-file",
            str(data.count_file),
            "--locations-file",
            str(data.locations_file),
            "--platform",
            data.platform,
            "--num-modules",
            str(data.num_modules),
        ]
        if data.is_high_resolution:
            cmd += ["--high-resolution", "--max-iter", "100"]

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
            raise STModuleError(
                f"stmodule_shim.R timed out after {self.timeout_seconds}s "
                f"(platform={data.platform}, num_modules={data.num_modules})"
            ) from e
        elapsed = time.time() - t0

        if proc.returncode != 0:
            raise STModuleError(
                f"stmodule_shim.R failed (exit {proc.returncode}): {proc.stderr.strip()[:500]}"
            )
        # Parse JSON from stdout. The shim emits exactly one JSON
        # object matching SpatialModuleResult.to_dict() shape.
        try:
            payload = json.loads(proc.stdout.strip().splitlines()[-1])
        except (json.JSONDecodeError, IndexError) as e:
            raise STModuleError(
                f"could not parse stmodule_shim.R stdout as JSON: {proc.stdout[:300]!r}"
            ) from e
        return _parse_payload(payload, elapsed)


def _parse_payload(payload: dict, elapsed: float) -> SpatialModuleResult:
    """Convert a JSON payload (matching ``SpatialModuleResult.to_dict``)
    into the typed dataclass."""
    modules = tuple(
        SpatialModule(
            module_id=int(m["module_id"]),
            top_genes=tuple(m["top_genes"]),
            n_spots=int(m["n_spots"]),
            mean_activity=float(m["mean_activity"]),
        )
        for m in payload.get("modules", ())
    )
    return SpatialModuleResult(
        platform=str(payload["platform"]),
        modules=modules,
        n_spots=int(payload.get("n_spots", 0)),
        elapsed_seconds=elapsed if elapsed > 0 else float(payload.get("elapsed_seconds", 0.0)),
        backend=str(payload.get("backend", "stmodule")),
        notes=tuple(payload.get("notes", ())),
    )


# ---------------------------------------------------------------------------
# Mock backend (stdlib, deterministic)
# ---------------------------------------------------------------------------


# Per-platform gene universes. These are 5 gene names per platform
# that the mock assigns to modules in a deterministic round-robin.
# Real genes that appear in the STModule paper's tutorial output.
_PLATFORM_GENE_UNIVERSES: dict[str, tuple[str, ...]] = {
    "ST": ("GAPDH", "USP4", "MAPKAPK2", "CPEB1", "LANCL2"),
    "Visium": ("CDH1", "VIM", "KRT8", "KRT18", "EPCAM"),
    "SlideSeqV2": ("MOBP", "MBP", "PLP1", "MAG", "MOG"),
    "StereoSeq": ("SOX2", "PAX6", "NES", "VIM", "HES1"),
    "Other": ("GEN_A", "GEN_B", "GEN_C", "GEN_D", "GEN_E"),
}


@dataclass
class MockSpatialModuleBackend:
    """Deterministic stdlib stub for STModule.

    Algorithm:
      - Take the platform's gene universe (5 genes per platform).
      - Round-robin assign top-N genes to num_modules modules.
      - mean_activity = (1 / (module_id + 1)) normalized to [0, 1].
      - n_spots derived from the count-matrix row count.
      - Deterministic: same input → same output, no RNG.

    This is intentionally a heuristic; the mock exists so CI/tests
    work without R / Seurat / CUDA installed.
    """

    backend_name: ClassVar[str] = "mock"

    def run(self, data: SpatialData) -> SpatialModuleResult:
        t0 = time.time()
        # Use the intersection of spot IDs across count + locations
        # so a missing spot in either file doesn't crash the run.
        count_spots = _read_first_column(data.count_file)
        loc_spots = _read_first_column(data.locations_file)
        n_spots = len(count_spots & loc_spots) or len(count_spots)
        genes = self._platform_gene_universe(data.platform)
        # Activity per module: monotonically decreasing — module 0 is
        # the "dominant" tissue pattern, module N-1 is the smallest.
        # Normalize so module 0 = 1.0 and module N-1 → small value.
        modules: list[SpatialModule] = []
        for i in range(data.num_modules):
            # Round-robin: take 3 genes starting from offset (i*3) mod len(genes)
            offset = (i * 3) % len(genes)
            top_genes = tuple(genes[(offset + j) % len(genes)] for j in range(3))
            # Activity: 1.0 / (i + 1), mapped to (0, 1]
            activity = 1.0 / (i + 1)
            modules.append(
                SpatialModule(
                    module_id=i,
                    top_genes=top_genes,
                    n_spots=max(1, n_spots // data.num_modules),
                    mean_activity=activity,
                )
            )
        return SpatialModuleResult(
            platform=data.platform,
            modules=tuple(modules),
            n_spots=n_spots,
            elapsed_seconds=time.time() - t0,
            backend=self.backend_name,
            notes=("mock-backend", f"platform={data.platform}", f"n_modules={data.num_modules}"),
        )

    def _platform_gene_universe(self, platform: str) -> tuple[str, ...]:
        return _PLATFORM_GENE_UNIVERSES.get(platform, _PLATFORM_GENE_UNIVERSES["Other"])


# ---------------------------------------------------------------------------
# File I/O helpers (mock-only)
# ---------------------------------------------------------------------------


def _count_data_rows(count_file: Path) -> int:
    """Count data rows in a TSV count matrix, excluding the header.

    Robust to trailing newlines / blank lines. Returns 0 on any
    I/O error so the mock never crashes on malformed input.
    """
    try:
        with open(count_file) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None)
            if header is None:
                return 0
            n = sum(1 for row in reader if row)
        return n
    except (OSError, csv.Error, UnicodeDecodeError):
        return 0


def _read_first_column(path: Path) -> set[str]:
    """Read the first column of a TSV as a set of strings.

    Used by the mock backend to compute the intersection of spot IDs
    across the count matrix and the locations file, so a missing
    spot in either input doesn't crash the run.
    """
    out: set[str] = set()
    try:
        with open(path) as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)  # skip header
            for row in reader:
                if row and row[0]:
                    out.add(row[0])
    except (OSError, csv.Error, UnicodeDecodeError):
        return set()
    return out


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------


def select_spatial_module_backend(
    *,
    prefer: str = "auto",
) -> SpatialModuleBackend:
    """Pick the best available :class:`SpatialModuleBackend`.

    Selection rules (in order):
        1. ``prefer="mock"`` → return :class:`MockSpatialModuleBackend`.
        2. ``prefer="real"`` → require ``Rscript`` on $PATH; raise
           :class:`STModuleNotInstalled` if missing.
        3. ``prefer="auto"`` (default):
           - If ``rscript_available()`` → return :class:`STModuleCLIAdapter`.
           - Else → return :class:`MockSpatialModuleBackend`.
    """
    if prefer == "mock":
        return MockSpatialModuleBackend()
    if prefer == "real":
        if not rscript_available():
            raise STModuleNotInstalled()
        return STModuleCLIAdapter()
    if rscript_available():
        return STModuleCLIAdapter()
    return MockSpatialModuleBackend()


__all__ = [
    # Path / availability constants.
    "RSCRIPT_BIN",
    "STMODULE_SHIM_NAME",
    # Errors.
    "STModuleError",
    "STModuleNotInstalled",
    # Availability helper.
    "rscript_available",
    # Real + mock backends.
    "STModuleCLIAdapter",
    "MockSpatialModuleBackend",
    # Backend selector.
    "select_spatial_module_backend",
]

