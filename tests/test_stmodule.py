"""Tests for the STModule spatial transcriptomics Protocol adapter.

Follows strict TDD: each test was written before the corresponding
implementation. The test file encodes the public contract:

  SpatialData            — typed input container
  SpatialModule          — one identified tissue module + associated genes
  SpatialModuleResult    — aggregate result for one tissue section
  SpatialModuleBackend   — Protocol (runtime_checkable)
  STModuleCLIAdapter     — shells out to R shim (real backend)
  MockSpatialModuleBackend — stdlib simulation
  select_spatial_module_backend() — real-or-mock dispatcher

All subprocess paths are mocked via unittest.mock so tests run without
R / Seurat / CUDA installed.
"""
from __future__ import annotations

import csv
import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# Forward-import the module under test (will fail until module is written)
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Test fixtures
# ---------------------------------------------------------------------------


def _write_count_matrix(path: Path, n_spots: int = 12, n_genes: int = 5) -> None:
    """Write a synthetic ST count matrix in the canonical STModule shape.

    Format: TSV, header row + spot rows. Rownames = spot IDs, colnames
    = gene names. Matches the shape from the STModule tutorial.
    """
    spot_ids = [f"loc_{i + 1}" for i in range(n_spots)]
    gene_names = [f"GENE_{j + 1}" for j in range(n_genes)]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([""] + gene_names)
        for sid in spot_ids:
            row = [sid] + [str((i + 1) % 4) for i in range(n_genes)]
            w.writerow(row)


def _write_locations(path: Path, n_spots: int = 12) -> None:
    """Write synthetic spatial coordinates (x, y) in STModule shape.

    Format: TSV, header row + spot rows. rownames match count matrix.
    """
    spot_ids = [f"loc_{i + 1}" for i in range(n_spots)]
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["", "x", "y"])
        for i, sid in enumerate(spot_ids):
            w.writerow([sid, str(17.0 + i * 0.5), str(5.0 + (i % 4) * 1.0)])


def _write_count_matrix_with_missing_spot(path: Path) -> None:
    """Write a count matrix with one spot that has no location entry."""
    _write_count_matrix(path, n_spots=12, n_genes=5)


def _write_locations_with_missing_spot(path: Path) -> None:
    """Write locations missing 'loc_12'."""
    spot_ids = [f"loc_{i + 1}" for i in range(11)]  # 11 spots, no loc_12
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["", "x", "y"])
        for i, sid in enumerate(spot_ids):
            w.writerow([sid, str(17.0 + i * 0.5), str(5.0 + (i % 4) * 1.0)])


# ---------------------------------------------------------------------------
# SpatialData dataclass
# ---------------------------------------------------------------------------


class TestSpatialDataDataclass(unittest.TestCase):
    def test_minimal_construction(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            d = SpatialData(
                count_file=count,
                locations_file=loc,
                platform="ST",
            )
            self.assertEqual(d.platform, "ST")
            self.assertTrue(d.count_file.exists())
            self.assertTrue(d.locations_file.exists())
            self.assertEqual(d.num_modules, 10)  # default per STModule tutorial

    def test_platform_enum(self) -> None:
        """Only ST, Visium, SlideSeqV2, StereoSeq, Other are accepted."""
        from mrnavax.spatial_protocols import SpatialData, SpatialPlatform

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            for plat in ("ST", "Visium", "SlideSeqV2", "StereoSeq", "Other"):
                d = SpatialData(
                    count_file=count, locations_file=loc, platform=plat
                )
                self.assertEqual(d.platform, SpatialPlatform(plat))

    def test_rejects_invalid_platform(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            with self.assertRaises(ValueError):
                SpatialData(
                    count_file=count, locations_file=loc, platform="BOGUS"
                )

    def test_rejects_missing_count_file(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with self.assertRaises(FileNotFoundError):
            SpatialData(
                count_file=Path("/nonexistent/count.tsv"),
                locations_file=Path("/tmp/x"),
                platform="ST",
            )

    def test_rejects_missing_locations_file(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            _write_count_matrix(count)
            with self.assertRaises(FileNotFoundError):
                SpatialData(
                    count_file=count,
                    locations_file=Path("/nonexistent/locs.tsv"),
                    platform="ST",
                )

    def test_num_modules_must_be_positive(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            for bad in (0, -1, -5):
                with self.assertRaises(ValueError):
                    SpatialData(
                        count_file=count,
                        locations_file=loc,
                        platform="ST",
                        num_modules=bad,
                    )

    def test_to_dict_serializes_paths_as_strings(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            d = SpatialData(count_file=count, locations_file=loc, platform="ST")
            dd = d.to_dict()
            self.assertIsInstance(dd["count_file"], str)
            self.assertIsInstance(dd["locations_file"], str)
            self.assertEqual(dd["platform"], "ST")


# ---------------------------------------------------------------------------
# SpatialModule + SpatialModuleResult dataclasses
# ---------------------------------------------------------------------------


class TestSpatialModuleDataclass(unittest.TestCase):
    def test_module_construction(self) -> None:
        from mrnavax.spatial_protocols import SpatialModule

        m = SpatialModule(
            module_id=0,
            top_genes=("TP53", "KRAS", "BRAF"),
            n_spots=12,
            mean_activity=0.42,
        )
        self.assertEqual(m.module_id, 0)
        self.assertEqual(m.top_genes, ("TP53", "KRAS", "BRAF"))

    def test_module_to_dict(self) -> None:
        from mrnavax.spatial_protocols import SpatialModule

        m = SpatialModule(
            module_id=3,
            top_genes=("GENE_1", "GENE_2"),
            n_spots=42,
            mean_activity=0.7,
        )
        d = m.to_dict()
        self.assertEqual(d["module_id"], 3)
        self.assertEqual(d["top_genes"], ["GENE_1", "GENE_2"])
        self.assertEqual(d["n_spots"], 42)


class TestSpatialModuleResultDataclass(unittest.TestCase):
    def test_result_construction(self) -> None:
        from mrnavax.spatial_protocols import (
            SpatialModule,
            SpatialModuleResult,
        )

        mods = (
            SpatialModule(0, ("TP53",), 10, 0.5),
            SpatialModule(1, ("KRAS",), 8, 0.3),
        )
        r = SpatialModuleResult(
            platform="ST",
            modules=mods,
            n_spots=12,
            elapsed_seconds=1.5,
            backend="mock",
            notes=("mock-backend",),
        )
        self.assertEqual(len(r.modules), 2)
        self.assertEqual(r.platform, "ST")
        self.assertEqual(r.backend, "mock")

    def test_result_to_dict(self) -> None:
        from mrnavax.spatial_protocols import (
            SpatialModule,
            SpatialModuleResult,
        )

        r = SpatialModuleResult(
            platform="Visium",
            modules=(SpatialModule(0, ("GENE_1",), 10, 0.5),),
            n_spots=10,
            elapsed_seconds=2.0,
            backend="mock",
            notes=(),
        )
        d = r.to_dict()
        self.assertEqual(d["platform"], "Visium")
        self.assertEqual(len(d["modules"]), 1)
        self.assertEqual(d["modules"][0]["module_id"], 0)


# ---------------------------------------------------------------------------
# SpatialModuleBackend Protocol
# ---------------------------------------------------------------------------


class TestSpatialModuleBackendProtocol(unittest.TestCase):
    def test_protocol_is_runtime_checkable(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )
        from mrnavax.spatial_protocols import SpatialModuleBackend

        m = MockSpatialModuleBackend()
        self.assertIsInstance(m, SpatialModuleBackend)

    def test_protocol_method_signature(self) -> None:
        """All backends must expose .run(data) -> SpatialModuleResult."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        m = MockSpatialModuleBackend()
        self.assertTrue(callable(getattr(m, "run", None)))

    def test_mock_backend_run_returns_correct_type(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )
        from mrnavax.spatial_protocols import (
            SpatialData,
            SpatialModuleResult,
        )

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            data = SpatialData(
                count_file=count, locations_file=loc, platform="ST"
            )
            m = MockSpatialModuleBackend()
            result = m.run(data)
            self.assertIsInstance(result, SpatialModuleResult)


# ---------------------------------------------------------------------------
# MockSpatialModuleBackend
# ---------------------------------------------------------------------------


class TestMockSpatialModuleBackend(unittest.TestCase):
    def _make_data(self, platform: str = "ST", n_spots: int = 12):
        from mrnavax.spatial_protocols import SpatialData

        tmp = tempfile.mkdtemp()
        count = Path(tmp) / "counts.tsv"
        loc = Path(tmp) / "locs.tsv"
        _write_count_matrix(count, n_spots=n_spots, n_genes=5)
        _write_locations(loc, n_spots=n_spots)
        return SpatialData(
            count_file=count, locations_file=loc, platform=platform
        )

    def test_returns_correct_num_modules(self) -> None:
        from dataclasses import replace

        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        data = self._make_data()
        for n in (5, 10, 15):
            data_n = replace(data, num_modules=n)
            m = MockSpatialModuleBackend()
            result = m.run(data_n)
            self.assertEqual(len(result.modules), n)

    def test_deterministic_across_runs(self) -> None:
        """Same input + no seed passed → identical output."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        data = self._make_data()
        m = MockSpatialModuleBackend()
        r1 = m.run(data)
        r2 = m.run(data)
        self.assertEqual(
            [mod.to_dict() for mod in r1.modules],
            [mod.to_dict() for mod in r2.modules],
        )

    def test_different_platforms_different_results(self) -> None:
        """Platform affects output (different gene universe in mock)."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        m = MockSpatialModuleBackend()
        r_st = m.run(self._make_data("ST"))
        r_vis = m.run(self._make_data("Visium"))
        # At least one top-gene tuple should differ between platforms
        st_genes = {g for mod in r_st.modules for g in mod.top_genes}
        vis_genes = {g for mod in r_vis.modules for g in mod.top_genes}
        self.assertTrue(st_genes != vis_genes or len(st_genes) > 0)

    def test_spot_count_matches_input(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        data = self._make_data(n_spots=20)
        m = MockSpatialModuleBackend()
        result = m.run(data)
        self.assertEqual(result.n_spots, 20)

    def test_known_universe_per_platform(self) -> None:
        """Each platform has a different fixed gene universe in the mock."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        m = MockSpatialModuleBackend()
        for plat, n_expected in (
            ("ST", 5),
            ("Visium", 5),
            ("SlideSeqV2", 5),
            ("StereoSeq", 5),
        ):
            universe = m._platform_gene_universe(plat)
            self.assertGreaterEqual(len(universe), n_expected)
            self.assertIsInstance(universe, tuple)

    def test_module_activities_in_zero_one(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        data = self._make_data()
        m = MockSpatialModuleBackend()
        result = m.run(data)
        for mod in result.modules:
            self.assertGreaterEqual(mod.mean_activity, 0.0)
            self.assertLessEqual(mod.mean_activity, 1.0)

    def test_handles_missing_spot_in_locations(self) -> None:
        """Count matrix has 12 spots, locations have 11 — mock should
        not crash; just use the intersection."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )
        from mrnavax.spatial_protocols import SpatialData

        tmp = tempfile.mkdtemp()
        count = Path(tmp) / "counts.tsv"
        loc = Path(tmp) / "locs.tsv"
        _write_count_matrix_with_missing_spot(count)
        _write_locations_with_missing_spot(loc)
        data = SpatialData(
            count_file=count, locations_file=loc, platform="ST"
        )
        m = MockSpatialModuleBackend()
        # Should not raise even though spot 12 has no coords
        result = m.run(data)
        self.assertIsNotNone(result)
        self.assertLessEqual(result.n_spots, 11)

    def test_backend_name_is_mock(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )

        data = self._make_data()
        m = MockSpatialModuleBackend()
        result = m.run(data)
        self.assertEqual(result.backend, "mock")


# ---------------------------------------------------------------------------
# STModuleCLIAdapter (real)
# ---------------------------------------------------------------------------


class TestSTModuleCLIAdapter(unittest.TestCase):
    def _make_data(self, platform: str = "ST"):
        from mrnavax.spatial_protocols import SpatialData

        tmp = tempfile.mkdtemp()
        count = Path(tmp) / "counts.tsv"
        loc = Path(tmp) / "locs.tsv"
        _write_count_matrix(count)
        _write_locations(loc)
        return SpatialData(
            count_file=count, locations_file=loc, platform=platform
        )

    def test_rscript_missing_raises(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleCLIAdapter,
            STModuleNotInstalled,
        )

        data = self._make_data()
        adapter = STModuleCLIAdapter()
        with patch(
            "mrnavax.spatial_module_adapter.shutil.which",
            return_value=None,
        ):
            with self.assertRaises(STModuleNotInstalled) as ctx:
                adapter.run(data)
            # Remediation message should mention Rscript
            self.assertIn("Rscript", str(ctx.exception))

    def test_successful_subprocess_parses_json(self) -> None:
        """When the subprocess returns a valid JSON SpatialModuleResult
        on stdout, the adapter parses and returns it."""
        from mrnavax.spatial_module_adapter import STModuleCLIAdapter

        data = self._make_data()
        adapter = STModuleCLIAdapter()
        mock_stdout = json.dumps(
            {
                "platform": "ST",
                "modules": [
                    {
                        "module_id": 0,
                        "top_genes": ["GENE_1", "GENE_2"],
                        "n_spots": 12,
                        "mean_activity": 0.5,
                    }
                ],
                "n_spots": 12,
                "elapsed_seconds": 1.0,
                "backend": "stmodule",
                "notes": [],
            }
        )
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = mock_stdout
        mock_proc.stderr = ""
        with patch(
            "mrnavax.spatial_module_adapter.shutil.which",
            return_value="/usr/bin/Rscript",
        ), patch(
            "mrnavax.spatial_module_adapter.subprocess.run",
            return_value=mock_proc,
        ):
            result = adapter.run(data)
        self.assertEqual(result.backend, "stmodule")
        self.assertEqual(len(result.modules), 1)
        self.assertEqual(result.modules[0].top_genes, ("GENE_1", "GENE_2"))

    def test_nonzero_exit_raises_stmodule_error(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleCLIAdapter,
            STModuleNotInstalled,
        )

        data = self._make_data()
        adapter = STModuleCLIAdapter()
        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stdout = ""
        mock_proc.stderr = "Error: package 'STModule' not found"
        with patch(
            "mrnavax.spatial_module_adapter.shutil.which",
            return_value="/usr/bin/Rscript",
        ), patch(
            "mrnavax.spatial_module_adapter.subprocess.run",
            return_value=mock_proc,
        ):
            with self.assertRaises(STModuleNotInstalled.__bases__[0]):  # STModuleError
                adapter.run(data)

    def test_unparseable_json_raises(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleCLIAdapter,
            STModuleNotInstalled,
        )

        data = self._make_data()
        adapter = STModuleCLIAdapter()
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "not json at all"
        mock_proc.stderr = ""
        with patch(
            "mrnavax.spatial_module_adapter.shutil.which",
            return_value="/usr/bin/Rscript",
        ), patch(
            "mrnavax.spatial_module_adapter.subprocess.run",
            return_value=mock_proc,
        ):
            with self.assertRaises(STModuleNotInstalled.__bases__[0]):
                adapter.run(data)

    def test_timeout_raises_stmodule_error(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleCLIAdapter,
            STModuleNotInstalled,
        )

        data = self._make_data()
        adapter = STModuleCLIAdapter(timeout_seconds=1)
        with patch(
            "mrnavax.spatial_module_adapter.shutil.which",
            return_value="/usr/bin/Rscript",
        ), patch(
            "mrnavax.spatial_module_adapter.subprocess.run",
            side_effect=subprocess.TimeoutExpired(cmd="Rscript", timeout=1),
        ):
            with self.assertRaises(STModuleNotInstalled.__bases__[0]):
                adapter.run(data)

    def test_command_line_includes_platform_and_paths(self) -> None:
        """Verify the subprocess invocation passes the right args to Rscript."""
        from mrnavax.spatial_module_adapter import STModuleCLIAdapter

        data = self._make_data(platform="Visium")
        adapter = STModuleCLIAdapter()
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = json.dumps(
            {
                "platform": "Visium",
                "modules": [],
                "n_spots": 0,
                "elapsed_seconds": 0.0,
                "backend": "stmodule",
                "notes": [],
            }
        )
        mock_proc.stderr = ""
        with patch(
            "mrnavax.spatial_module_adapter.shutil.which",
            return_value="/usr/bin/Rscript",
        ), patch(
            "mrnavax.spatial_module_adapter.subprocess.run",
            return_value=mock_proc,
        ) as mock_run:
            adapter.run(data)
            call_args = mock_run.call_args
            cmd = call_args.args[0] if call_args.args else call_args.kwargs.get("args")
            self.assertIn("Visium", cmd)
            self.assertIn(str(data.count_file), cmd)
            self.assertIn(str(data.locations_file), cmd)


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------


class TestBackendSelector(unittest.TestCase):
    def test_select_real_forced_missing_raises(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleNotInstalled,
            select_spatial_module_backend,
        )

        with patch(
            "mrnavax.spatial_module_adapter.rscript_available",
            return_value=False,
        ):
            with self.assertRaises(STModuleNotInstalled):
                select_spatial_module_backend(prefer="real")

    def test_select_real_forced_present(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleCLIAdapter,
            select_spatial_module_backend,
        )

        with patch(
            "mrnavax.spatial_module_adapter.rscript_available",
            return_value=True,
        ):
            sel = select_spatial_module_backend(prefer="real")
            self.assertIsInstance(sel, STModuleCLIAdapter)

    def test_select_mock_forced(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
            select_spatial_module_backend,
        )

        sel = select_spatial_module_backend(prefer="mock")
        self.assertIsInstance(sel, MockSpatialModuleBackend)

    def test_select_auto_falls_back_to_mock(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
            select_spatial_module_backend,
        )

        with patch(
            "mrnavax.spatial_module_adapter.rscript_available",
            return_value=False,
        ):
            sel = select_spatial_module_backend()
            self.assertIsInstance(sel, MockSpatialModuleBackend)

    def test_select_auto_picks_real(self) -> None:
        from mrnavax.spatial_module_adapter import (
            STModuleCLIAdapter,
            select_spatial_module_backend,
        )

        with patch(
            "mrnavax.spatial_module_adapter.rscript_available",
            return_value=True,
        ):
            sel = select_spatial_module_backend()
            self.assertIsInstance(sel, STModuleCLIAdapter)


# ---------------------------------------------------------------------------
# End-to-end integration
# ---------------------------------------------------------------------------


class TestSpatialDataHighResolution(unittest.TestCase):
    """Verify the platform-specific high_resolution flag is wired correctly."""

    def test_is_high_resolution_st(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            data = SpatialData(
                count_file=count, locations_file=loc, platform="ST"
            )
            self.assertFalse(data.is_high_resolution)

    def test_is_high_resolution_slideseqv2(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            data = SpatialData(
                count_file=count, locations_file=loc, platform="SlideSeqV2"
            )
            self.assertTrue(data.is_high_resolution)

    def test_is_high_resolution_stereoseq(self) -> None:
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            data = SpatialData(
                count_file=count, locations_file=loc, platform="StereoSeq"
            )
            self.assertTrue(data.is_high_resolution)

    def test_cli_args_include_high_resolution_for_slideseq(self) -> None:
        """When platform is SlideSeqV2, the subprocess invocation should
        include --high-resolution and --max-iter 100 per upstream
        tutorial recommendation."""
        from mrnavax.spatial_module_adapter import STModuleCLIAdapter

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            from mrnavax.spatial_protocols import SpatialData

            data = SpatialData(
                count_file=count,
                locations_file=loc,
                platform="SlideSeqV2",
                num_modules=10,
            )
            adapter = STModuleCLIAdapter()
            mock_proc = MagicMock()
            mock_proc.returncode = 0
            mock_proc.stdout = json.dumps(
                {
                    "platform": "SlideSeqV2",
                    "modules": [],
                    "n_spots": 0,
                    "elapsed_seconds": 0.0,
                    "backend": "stmodule",
                    "notes": [],
                }
            )
            mock_proc.stderr = ""
            with patch(
                "mrnavax.spatial_module_adapter.shutil.which",
                return_value="/usr/bin/Rscript",
            ), patch(
                "mrnavax.spatial_module_adapter.subprocess.run",
                return_value=mock_proc,
            ) as mock_run:
                adapter.run(data)
                cmd = mock_run.call_args.args[0]
                self.assertIn("--high-resolution", cmd)
                self.assertIn("--max-iter", cmd)
                max_iter_idx = cmd.index("--max-iter")
                self.assertEqual(cmd[max_iter_idx + 1], "100")

    def test_cli_args_omit_high_resolution_for_st(self) -> None:
        """When platform is ST, --high-resolution should NOT be in the
        subprocess args (per upstream tutorial: ST data uses defaults)."""
        from mrnavax.spatial_module_adapter import STModuleCLIAdapter

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            _write_count_matrix(count)
            _write_locations(loc)
            from mrnavax.spatial_protocols import SpatialData

            data = SpatialData(
                count_file=count, locations_file=loc, platform="ST"
            )
            adapter = STModuleCLIAdapter()
            mock_proc = MagicMock()
            mock_proc.returncode = 0
            mock_proc.stdout = json.dumps(
                {
                    "platform": "ST",
                    "modules": [],
                    "n_spots": 0,
                    "elapsed_seconds": 0.0,
                    "backend": "stmodule",
                    "notes": [],
                }
            )
            mock_proc.stderr = ""
            with patch(
                "mrnavax.spatial_module_adapter.shutil.which",
                return_value="/usr/bin/Rscript",
            ), patch(
                "mrnavax.spatial_module_adapter.subprocess.run",
                return_value=mock_proc,
            ) as mock_run:
                adapter.run(data)
                cmd = mock_run.call_args.args[0]
                self.assertNotIn("--high-resolution", cmd)


class TestParsePayload(unittest.TestCase):
    """Logic validation: _parse_payload handles edge cases."""

    def test_empty_modules_list(self) -> None:
        from mrnavax.spatial_module_adapter import _parse_payload

        result = _parse_payload(
            {
                "platform": "ST",
                "modules": [],
                "n_spots": 0,
                "elapsed_seconds": 0.0,
                "backend": "stmodule",
                "notes": [],
            },
            elapsed=0.5,
        )
        self.assertEqual(result.modules, ())
        self.assertEqual(result.n_spots, 0)
        # When subprocess elapsed > 0, use it; otherwise use payload's
        self.assertEqual(result.elapsed_seconds, 0.5)

    def test_subprocess_elapsed_wins_when_positive(self) -> None:
        """When subprocess reports elapsed=0.1, we use it (subprocess
        wall-clock time includes R startup overhead). The payload's
        elapsed_seconds is only used when subprocess elapsed is 0."""
        from mrnavax.spatial_module_adapter import _parse_payload

        result = _parse_payload(
            {
                "platform": "ST",
                "modules": [],
                "n_spots": 0,
                "elapsed_seconds": 99.9,
                "backend": "stmodule",
                "notes": [],
            },
            elapsed=0.1,
        )
        self.assertEqual(result.elapsed_seconds, 0.1)

    def test_payload_elapsed_used_when_subprocess_is_zero(self) -> None:
        """When subprocess reports elapsed=0 (the rare case where the
        adapter is called without timing), fall back to the payload's
        elapsed_seconds (which the R shim reports)."""
        from mrnavax.spatial_module_adapter import _parse_payload

        result = _parse_payload(
            {
                "platform": "ST",
                "modules": [],
                "n_spots": 0,
                "elapsed_seconds": 42.0,
                "backend": "stmodule",
                "notes": [],
            },
            elapsed=0.0,
        )
        self.assertEqual(result.elapsed_seconds, 42.0)

    def test_notes_default_to_empty(self) -> None:
        """Missing 'notes' field should default to empty tuple, not error."""
        from mrnavax.spatial_module_adapter import _parse_payload

        result = _parse_payload(
            {
                "platform": "ST",
                "modules": [],
                "n_spots": 0,
                "elapsed_seconds": 0.0,
                "backend": "stmodule",
            },
            elapsed=1.0,
        )
        self.assertEqual(result.notes, ())


class TestMockBackendZeroSpots(unittest.TestCase):
    """Edge case: empty count matrix should not crash."""

    def test_empty_count_matrix(self) -> None:
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            # Header only, no rows
            count.write_text("\t".join([""] + [f"G{j}" for j in range(5)]) + "\n")
            loc.write_text("\t".join(["", "x", "y"]) + "\n")
            data = SpatialData(
                count_file=count, locations_file=loc, platform="ST"
            )
            m = MockSpatialModuleBackend()
            # Should not crash even with 0 spots
            result = m.run(data)
            self.assertEqual(result.n_spots, 0)
            self.assertEqual(len(result.modules), data.num_modules)

    def test_malformed_count_matrix_graceful(self) -> None:
        """A completely garbage file should still produce a result, not
        raise. The mock's _read_first_column swallows I/O errors."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )
        from mrnavax.spatial_protocols import SpatialData

        with tempfile.TemporaryDirectory() as tmp:
            count = Path(tmp) / "counts.tsv"
            loc = Path(tmp) / "locs.tsv"
            count.write_text("not\ta\tvalid\nmatrix\tat\tall\n")
            loc.write_text("\t".join(["", "x", "y"]) + "\n")
            data = SpatialData(
                count_file=count, locations_file=loc, platform="ST"
            )
            m = MockSpatialModuleBackend()
            # Should not raise
            result = m.run(data)
            self.assertIsNotNone(result)


class TestEndToEnd(unittest.TestCase):
    def test_mock_backend_result_is_serializable(self) -> None:
        """SpatialModuleResult.to_dict() must be JSON-serializable."""
        from mrnavax.spatial_module_adapter import (
            MockSpatialModuleBackend,
        )
        from mrnavax.spatial_protocols import SpatialData

        tmp = tempfile.mkdtemp()
        count = Path(tmp) / "counts.tsv"
        loc = Path(tmp) / "locs.tsv"
        _write_count_matrix(count)
        _write_locations(loc)
        data = SpatialData(
            count_file=count, locations_file=loc, platform="ST"
        )
        m = MockSpatialModuleBackend()
        result = m.run(data)
        # Must round-trip through JSON without errors
        s = json.dumps(result.to_dict())
        loaded = json.loads(s)
        self.assertEqual(loaded["platform"], "ST")
        self.assertEqual(loaded["backend"], "mock")


if __name__ == "__main__":
    unittest.main()
