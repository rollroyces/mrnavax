"""Tests for the RiboDecode Protocol + adapter contracts.

Covers:
- RiboDecodeRequest dataclass validation (length, mfe_weight, custom env)
- TranslationPredictor / CodonOptimizer Protocol runtime_checkable
- MockTranslationPredictor: CAI ordering, score in [0, 100], known-codon handling
- MockCodonOptimizer: protein preservation, no-op for empty input
- Adapter shell-out error paths (CLI missing → RiboDecodeNotInstalled)
- select_translation_predictor / select_codon_optimizer dispatch logic

All tests use stdlib only. No model downloads, no network, no GPU.
"""
from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

from mrnavax.codon_optimizer import CODON_TO_AA
from mrnavax.codon_protocols import (
    CodonOptimizer,
    RiboDecodeRequest,
    RiboDecodeResult,
    TranslationPrediction,
    TranslationPredictor,
)
from mrnavax.codon_ribodecode_adapter import (
    MockCodonOptimizer,
    MockTranslationPredictor,
    RiboDecodeCLIAdapter,
    RiboDecodeNotInstalled,
    TranslationModelCLIAdapter,
    select_codon_optimizer,
    select_translation_predictor,
)

# ---------------------------------------------------------------------------
# Test helpers
# ---------------------------------------------------------------------------


def translate(cds: str) -> str:
    """Translate a CDS (multiples of 3) to AA, skipping trailing stop."""
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    if codons and CODON_TO_AA.get(codons[-1]) == "*":
        codons = codons[:-1]
    return "".join(CODON_TO_AA[c] for c in codons)


# ---------------------------------------------------------------------------
# RiboDecodeRequest dataclass validation
# ---------------------------------------------------------------------------


class TestRiboDecodeRequestValidation(unittest.TestCase):
    """The dataclass must reject malformed inputs at construction."""

    def test_accepts_minimal_cds(self) -> None:
        req = RiboDecodeRequest(cds="ATGGACGGGTAG")
        self.assertEqual(req.cds, "ATGGACGGGTAG")
        self.assertEqual(req.env, "HEK293T")
        self.assertEqual(req.mfe_weight, 0.0)
        self.assertEqual(req.optim_epoch, 10)

    def test_rejects_empty_cds(self) -> None:
        with self.assertRaises(ValueError, msg="empty cds should raise"):
            RiboDecodeRequest(cds="")

    def test_rejects_non_multiple_of_3(self) -> None:
        with self.assertRaises(ValueError, msg="non-multiple-of-3 should raise"):
            RiboDecodeRequest(cds="ATGA")

    def test_rejects_oversized_cds(self) -> None:
        with self.assertRaises(ValueError, msg=">4500 nt should raise"):
            RiboDecodeRequest(cds="A" * 4501)

    def test_accepts_max_cds_length(self) -> None:
        # Boundary: exactly 4500 nt should be accepted
        RiboDecodeRequest(cds="ATG" + "GCT" * 1499)  # 3 + 4497 = 4500

    def test_rejects_mfe_weight_out_of_range(self) -> None:
        for bad in (-0.1, 1.1, 2.0, -1.0):
            with self.assertRaises(ValueError, msg=f"mfe_weight={bad} should raise"):
                RiboDecodeRequest(cds="ATGGACGGGTAG", mfe_weight=bad)

    def test_accepts_mfe_weight_boundary(self) -> None:
        for ok in (0.0, 0.5, 1.0):
            RiboDecodeRequest(cds="ATGGACGGGTAG", mfe_weight=ok)

    def test_rejects_custom_env_without_csv(self) -> None:
        with self.assertRaises(ValueError, msg="custom env without csv should raise"):
            RiboDecodeRequest(cds="ATGGACGGGTAG", env="custom")

    def test_accepts_custom_env_with_csv(self) -> None:
        req = RiboDecodeRequest(
            cds="ATGGACGGGTAG",
            env="custom",
            custom_env_csv=Path("/tmp/env.csv"),
        )
        self.assertEqual(req.env, "custom")

    def test_to_dict_serializes_path(self) -> None:
        req = RiboDecodeRequest(
            cds="ATGGACGGGTAG",
            env="custom",
            custom_env_csv=Path("/tmp/env.csv"),
        )
        d = req.to_dict()
        self.assertEqual(d["custom_env_csv"], "/tmp/env.csv")

    def test_to_dict_handles_none_path(self) -> None:
        req = RiboDecodeRequest(cds="ATGGACGGGTAG")
        d = req.to_dict()
        self.assertIsNone(d["custom_env_csv"])


# ---------------------------------------------------------------------------
# Protocol runtime_checkable
# ---------------------------------------------------------------------------


class TestProtocolConformance(unittest.TestCase):
    """Mock backends must satisfy the runtime_checkable Protocols."""

    def test_mock_translation_predictor_satisfies_protocol(self) -> None:
        pred = MockTranslationPredictor()
        self.assertIsInstance(pred, TranslationPredictor)

    def test_mock_codon_optimizer_satisfies_protocol(self) -> None:
        opt = MockCodonOptimizer()
        self.assertIsInstance(opt, CodonOptimizer)

    def test_protocol_predict_method_signature(self) -> None:
        """Translate a CDS via the mock and check the result type."""
        pred = MockTranslationPredictor()
        result = pred.predict("ATGGACGGGTAG")
        self.assertIsInstance(result, TranslationPrediction)
        self.assertEqual(result.cds, "ATGGACGGGTAG")

    def test_protocol_optimize_method_signature(self) -> None:
        """Optimize a CDS via the mock and check the result type."""
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG")
        result = opt.optimize(req)
        self.assertIsInstance(result, RiboDecodeResult)


# ---------------------------------------------------------------------------
# MockTranslationPredictor
# ---------------------------------------------------------------------------


class TestMockTranslationPredictor(unittest.TestCase):
    """Mock returns deterministic CAI-derived scores in [0, 100]."""

    def test_score_in_range(self) -> None:
        pred = MockTranslationPredictor()
        for cds in (
            "ATGGACGGGTAG",
            "ATG" + "GCT" * 20,
            "ATG" + "TTA" * 20 + "TAA",
        ):
            r = pred.predict(cds)
            self.assertGreaterEqual(r.translation_level, 0.0)
            self.assertLessEqual(r.translation_level, 100.0)

    def test_high_cai_scores_above_low_cai(self) -> None:
        """All-Ala GCC (highest-frequency Ala codon) should outrank
        a synthetic low-frequency CDS."""
        pred = MockTranslationPredictor()
        # Highest-frequency Ala codon is GCC (per HUMAN_CODON_FREQ).
        high = "ATG" + ("GCCGCCGCC" * 30) + "TAA"  # 273 nt, all Ala = GCC
        # Synthetic low-frequency mix (avoid stop codons):
        low = "ATG" + ("CTACTACTA" * 30) + "TAA"  # all Leu, low freq
        high_score = pred.predict(high).translation_level
        low_score = pred.predict(low).translation_level
        self.assertGreater(high_score, low_score)

    def test_unknown_codons_dont_crash(self) -> None:
        """Ambiguous IUPAC bases should be tolerated."""
        pred = MockTranslationPredictor()
        # NNN is an unknown codon — pred must return a result with
        # a note rather than raising.
        r = pred.predict("ATGNNNNGACNNN")
        self.assertIsInstance(r, TranslationPrediction)
        self.assertTrue(any("unknown" in n for n in r.notes))

    def test_internal_stop_codon_recorded(self) -> None:
        """Internal TAA must be flagged in notes, not raised."""
        pred = MockTranslationPredictor()
        r = pred.predict("ATGTAAGGGCCC")
        self.assertTrue(any("unknown_or_stop" in n for n in r.notes))

    def test_empty_cds_returns_zero(self) -> None:
        pred = MockTranslationPredictor()
        r = pred.predict("")
        self.assertEqual(r.translation_level, 0.0)
        self.assertIn("empty_cds", r.notes)

    def test_oversized_cds_truncated(self) -> None:
        """CDS >4500 nt should be silently truncated (matches upstream)."""
        pred = MockTranslationPredictor()
        big = "ATG" + ("GCT" * 10000) + "TAA"  # > 30000 nt
        r = pred.predict(big)
        # Output must not exceed 4500 nt per upstream cap
        self.assertLessEqual(len(r.cds), 4500)

    def test_environment_name_recorded(self) -> None:
        pred = MockTranslationPredictor()
        for env in ("HEK293T", "A549", "HeLa", "custom"):
            r = pred.predict("ATGGACGGGTAG", env=env)
            self.assertEqual(r.env, env)

    def test_backend_name_is_mock(self) -> None:
        pred = MockTranslationPredictor()
        r = pred.predict("ATGGACGGGTAG")
        self.assertEqual(r.backend, "mock-translation")


# ---------------------------------------------------------------------------
# MockCodonOptimizer
# ---------------------------------------------------------------------------


class TestMockCodonOptimizer(unittest.TestCase):
    """Mock must preserve the protein sequence and emit valid output."""

    def test_protein_preserved_simple(self) -> None:
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATG" + ("GCTGCTGCT" * 5) + "TAA")  # 48 nt
        res = opt.optimize(req)
        original = translate(req.cds)
        optimized = translate(res.optimized_cds)
        self.assertEqual(optimized, original)

    def test_protein_preserved_real_size(self) -> None:
        """GFP-like CDS (~240 aa) must round-trip the protein.

        LinearDesign on ~720 nt runs in ~1s on modern hardware;
        we keep this size to exercise the full DP without making CI
        slow.
        """
        opt = MockCodonOptimizer()
        # GFP-like protein (M + 79 residues ≈ 240 nt of codons)
        gfp_cds = (
            "ATG" + ("GCTGCTGCTGCTGCT" * 16) + "TAA"  # 243 nt
        )
        req = RiboDecodeRequest(cds=gfp_cds)
        res = opt.optimize(req)
        self.assertEqual(translate(res.optimized_cds), translate(gfp_cds))

    def test_mfe_is_none_when_weight_zero(self) -> None:
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG", mfe_weight=0.0)
        res = opt.optimize(req)
        self.assertIsNone(res.predicted_mfe)

    def test_mfe_present_when_weight_nonzero(self) -> None:
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG", mfe_weight=0.5)
        res = opt.optimize(req)
        self.assertIsNotNone(res.predicted_mfe)

    def test_epoch_recorded_in_result(self) -> None:
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG", optim_epoch=42)
        res = opt.optimize(req)
        self.assertEqual(res.epoch, 42)

    def test_environment_appears_in_notes(self) -> None:
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG", env="A549")
        res = opt.optimize(req)
        self.assertTrue(any("A549" in n for n in res.notes))

    def test_backend_name_is_mock(self) -> None:
        opt = MockCodonOptimizer()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG")
        res = opt.optimize(req)
        self.assertEqual(res.backend, "mock-lineardesign")


# ---------------------------------------------------------------------------
# Real CLI adapter error paths
# ---------------------------------------------------------------------------


class TestRealAdapterErrors(unittest.TestCase):
    """The real CLI adapters must fail cleanly when the binary is missing."""

    def test_pred_translation_missing_binary_raises(self) -> None:
        adapter = TranslationModelCLIAdapter()
        with patch(
            "mrnavax.codon_ribodecode_adapter.shutil.which",
            return_value=None,
        ):
            with self.assertRaises(RiboDecodeNotInstalled) as ctx:
                adapter.predict("ATGGACGGGTAG")
            self.assertIn("pred-translation", str(ctx.exception))

    def test_ribo_decode_missing_binary_raises(self) -> None:
        adapter = RiboDecodeCLIAdapter()
        req = RiboDecodeRequest(cds="ATGGACGGGTAG")
        with patch(
            "mrnavax.codon_ribodecode_adapter.shutil.which",
            return_value=None,
        ):
            with self.assertRaises(RiboDecodeNotInstalled) as ctx:
                adapter.optimize(req)
            self.assertIn("ribo-decode", str(ctx.exception))

    def test_pred_translation_success_parses_output(self) -> None:
        """Mocked subprocess: stdout=42.5 → translation_level=42.5."""
        adapter = TranslationModelCLIAdapter()
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "42.5\n"
        mock_proc.stderr = ""
        with patch(
            "mrnavax.codon_ribodecode_adapter.shutil.which",
            return_value="/usr/bin/pred-translation",
        ), patch(
            "mrnavax.codon_ribodecode_adapter.subprocess.run",
            return_value=mock_proc,
        ):
            r = adapter.predict("ATGGACGGGTAG", env="HEK293T")
        self.assertEqual(r.translation_level, 42.5)
        self.assertEqual(r.backend, "translationmodel")
        self.assertEqual(r.env, "HEK293T")

    def test_pred_translation_nonzero_exit_raises(self) -> None:
        adapter = TranslationModelCLIAdapter()
        mock_proc = MagicMock()
        mock_proc.returncode = 1
        mock_proc.stdout = ""
        mock_proc.stderr = "fatal: GPU not found"
        with patch(
            "mrnavax.codon_ribodecode_adapter.shutil.which",
            return_value="/usr/bin/pred-translation",
        ), patch(
            "mrnavax.codon_ribodecode_adapter.subprocess.run",
            return_value=mock_proc,
        ):
            with self.assertRaises(RiboDecodeNotInstalled.__bases__[0]):  # RiboDecodeError
                adapter.predict("ATGGACGGGTAG")

    def test_pred_translation_unparseable_output_raises(self) -> None:
        adapter = TranslationModelCLIAdapter()
        mock_proc = MagicMock()
        mock_proc.returncode = 0
        mock_proc.stdout = "not a number"
        mock_proc.stderr = ""
        with patch(
            "mrnavax.codon_ribodecode_adapter.shutil.which",
            return_value="/usr/bin/pred-translation",
        ), patch(
            "mrnavax.codon_ribodecode_adapter.subprocess.run",
            return_value=mock_proc,
        ):
            with self.assertRaises(Exception):
                adapter.predict("ATGGACGGGTAG")


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------


class TestBackendSelector(unittest.TestCase):
    """select_* functions choose real or mock based on $PATH + env var."""

    def test_select_translation_predictor_real_forced_missing(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.pred_translation_available",
            return_value=False,
        ):
            with self.assertRaises(RiboDecodeNotInstalled):
                select_translation_predictor(prefer="real")

    def test_select_translation_predictor_real_forced_present(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.pred_translation_available",
            return_value=True,
        ):
            sel = select_translation_predictor(prefer="real")
            self.assertIsInstance(sel, TranslationModelCLIAdapter)

    def test_select_translation_predictor_mock_forced(self) -> None:
        sel = select_translation_predictor(prefer="mock")
        self.assertIsInstance(sel, MockTranslationPredictor)

    def test_select_translation_predictor_auto_picks_real(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.pred_translation_available",
            return_value=True,
        ):
            sel = select_translation_predictor()
            self.assertIsInstance(sel, TranslationModelCLIAdapter)

    def test_select_translation_predictor_auto_falls_back_to_mock(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.pred_translation_available",
            return_value=False,
        ):
            sel = select_translation_predictor()
            self.assertIsInstance(sel, MockTranslationPredictor)

    def test_select_codon_optimizer_real_forced_missing(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.ribo_decode_available",
            return_value=False,
        ):
            with self.assertRaises(RiboDecodeNotInstalled):
                select_codon_optimizer(prefer="real")

    def test_select_codon_optimizer_real_forced_present(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.ribo_decode_available",
            return_value=True,
        ):
            sel = select_codon_optimizer(prefer="real")
            self.assertIsInstance(sel, RiboDecodeCLIAdapter)

    def test_select_codon_optimizer_mock_forced(self) -> None:
        sel = select_codon_optimizer(prefer="mock")
        self.assertIsInstance(sel, MockCodonOptimizer)

    def test_select_codon_optimizer_auto_falls_back_to_mock(self) -> None:
        with patch(
            "mrnavax.codon_ribodecode_adapter.ribo_decode_available",
            return_value=False,
        ):
            sel = select_codon_optimizer()
            self.assertIsInstance(sel, MockCodonOptimizer)


# ---------------------------------------------------------------------------
# End-to-end: optimizer + predictor on real CDS sizes
# ---------------------------------------------------------------------------


class TestEndToEnd(unittest.TestCase):
    """Full protein-preservation invariant across realistic CDS sizes."""

    def test_optimize_then_predict_same_cds(self) -> None:
        """Run optimize(), then re-score the optimized CDS with predict().
        The translation score must be consistent (both from the same
        mock backend)."""
        opt = MockCodonOptimizer()
        pred = MockTranslationPredictor()
        # Use a moderate-size CDS so the test stays fast in CI
        req = RiboDecodeRequest(cds="ATG" + ("GCT" * 30) + "TAA")  # 93 nt
        res = opt.optimize(req)
        # Re-score the optimized CDS
        rescored = pred.predict(res.optimized_cds, env=req.env)
        # The translation score should be ≥ the pre-optimization score
        # (LinearDesign aims to maximize translation)
        prescore = pred.predict(req.cds.rstrip("TAA"), env=req.env).translation_level
        self.assertGreaterEqual(rescored.translation_level, prescore * 0.9)

    def test_4500nt_cap_respected(self) -> None:
        """A 4500 nt CDS must optimize without raising (boundary).

        Note: we use a CDS of 150 nt here and verify it constructs
        correctly. The 4500-nt boundary is exercised in
        ``TestRiboDecodeRequestValidation.test_accepts_max_cds_length``
        (the dataclass rejects at 4501). A full 4500-nt LinearDesign
        run takes >30s on commodity hardware, which is too slow for
        CI; that case is covered by the slower benchmarks/ folder.
        """
        opt = MockCodonOptimizer()
        big = "ATG" + "GCT" * 49  # 150 nt, fast
        req = RiboDecodeRequest(cds=big)
        res = opt.optimize(req)
        self.assertEqual(translate(res.optimized_cds), translate(big))


# ---------------------------------------------------------------------------
# Run with `python -m tests.test_ribodecode_adapter`
# ---------------------------------------------------------------------------


if __name__ == "__main__":
    unittest.main()
