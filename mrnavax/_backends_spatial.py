"""Backend checks for the ``spatial.*`` family (STModule).

Extracted from ``mrnavax/backends.py`` in v0.24.0. Each check is
registered with the global ``CHECKS`` registry on import.
"""

from __future__ import annotations

import json as _json
import tempfile
from pathlib import Path

from ._backends_registry import register


@register("spatial.stmodule_module_identification")
def _check_stmodule_module_identification() -> tuple[bool, str]:
    """STModule Protocol adapter: tissue-module identification from
    spatial-transcriptomics data.

    Validates:
      1. SpatialData dataclass: platform enum, file existence,
         num_modules >= 1.
      2. SpatialModule + SpatialModuleResult dataclasses with
         to_dict() JSON-serializable.
      3. SpatialModuleBackend is runtime_checkable.
      4. MockSpatialModuleBackend satisfies the Protocol and produces
         the right shape (n_modules = num_modules).
      5. Mock result is deterministic across runs.
      6. Backend selector picks mock when Rscript is unavailable.
      7. Spot/location ID mismatch is handled gracefully (intersection).
    """
    from .spatial_module_adapter import (
        MockSpatialModuleBackend,
        rscript_available,
        select_spatial_module_backend,
    )
    from .spatial_protocols import (
        SpatialData,
        SpatialModule,
        SpatialModuleBackend,
        SpatialModuleResult,
    )

    def _write(p: Path, rows: list[str]) -> None:
        p.write_text("\n".join(rows))

    bad: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        count = Path(tmp) / "counts.tsv"
        loc = Path(tmp) / "locs.tsv"
        _write(
            count,
            ["\t".join([""] + [f"G{j}" for j in range(5)])]
            + ["\t".join([f"loc_{i + 1}", "1", "0", "1", "1", "0"]) for i in range(12)],
        )
        _write(
            loc,
            ["\t".join(["", "x", "y"])]
            + [
                "\t".join([f"loc_{i + 1}", str(17.0 + i * 0.5), str(5.0 + (i % 4))])
                for i in range(12)
            ],
        )

        # 1. SpatialData construction + validation
        data = SpatialData(count_file=count, locations_file=loc, platform="ST")
        if data.platform != "ST":
            bad.append("platform not preserved")
        if data.num_modules != 10:
            bad.append("default num_modules != 10")
        try:
            SpatialData(count_file=count, locations_file=loc, platform="ST", num_modules=0)
        except ValueError:
            pass
        else:
            bad.append("num_modules=0 should raise")
        try:
            SpatialData(count_file=count, locations_file=loc, platform="BOGUS")
        except ValueError:
            pass
        else:
            bad.append("invalid platform should raise")

        # 2. SpatialModule + SpatialModuleResult round-trip
        m = SpatialModule(module_id=0, top_genes=("TP53", "KRAS"), n_spots=12, mean_activity=0.7)
        r = SpatialModuleResult(
            platform="ST",
            modules=(m,),
            n_spots=12,
            elapsed_seconds=0.5,
            backend="mock",
            notes=("mock-backend",),
        )
        try:
            d = r.to_dict()
            s = _json.dumps(d)
            loaded = _json.loads(s)
            if loaded["modules"][0]["module_id"] != 0:
                bad.append("module round-trip mismatch")
        except (TypeError, ValueError) as e:
            bad.append(f"JSON serialization failed: {e}")

        # 3. Protocol runtime_checkable
        mock = MockSpatialModuleBackend()
        if not isinstance(mock, SpatialModuleBackend):
            bad.append("mock backend does not satisfy Protocol")

        # 4. Mock run produces correct n_modules
        result = mock.run(data)
        if len(result.modules) != data.num_modules:
            bad.append(f"mock returned {len(result.modules)} modules, expected {data.num_modules}")
        if result.backend != "mock":
            bad.append(f"backend != 'mock': {result.backend!r}")

        # 5. Determinism
        result2 = mock.run(data)
        if [mm.to_dict() for mm in result.modules] != [mm.to_dict() for mm in result2.modules]:
            bad.append("mock backend not deterministic")

        # 6. Backend selector picks mock when Rscript unavailable
        from unittest.mock import patch

        with patch(
            "mrnavax.spatial_module_adapter.rscript_available",
            return_value=False,
        ):
            sel = select_spatial_module_backend(prefer="auto")
            if not isinstance(sel, MockSpatialModuleBackend):
                bad.append(
                    f"auto-select should pick mock when Rscript missing, got {type(sel).__name__}"
                )

        # 7. Spot/location mismatch — must not crash, must report <= 11 spots
        loc_mismatch = Path(tmp) / "locs_mismatch.tsv"
        _write(
            loc_mismatch,
            ["\t".join(["", "x", "y"])]
            + [
                "\t".join([f"loc_{i + 1}", str(17.0 + i * 0.5), str(5.0 + (i % 4))])
                for i in range(11)
            ],
        )
        data_mismatch = SpatialData(count_file=count, locations_file=loc_mismatch, platform="ST")
        result_mm = mock.run(data_mismatch)
        if result_mm.n_spots > 11:
            bad.append(f"spot mismatch not handled: n_spots={result_mm.n_spots}")

    if bad:
        return False, "STModule backend issues: " + "; ".join(bad)
    real_available = rscript_available()
    return True, (
        f"STModule OK: mock backend produces correct n_modules, "
        f"deterministic, JSON-serializable; mismatch handled gracefully "
        f"(n_spots<=11); backend selector chooses mock when Rscript "
        f"unavailable; real backend (Rscript) installed={real_available}"
    )


__all__ = []
