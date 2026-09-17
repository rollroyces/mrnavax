"""Tests for the ESM2 protein-LM immunogenicity embedder.

Follows strict TDD: each test was written before the corresponding
implementation. The test file encodes the public contract:

  EmbeddingRequest       — typed input (sequence list, model_id, pooling)
  EmbeddingResult        — typed output (vectors, dim, backend, model_id)
  ProteinLMEmbedder      — Protocol (runtime_checkable)
  ESM2Embedder           — real adapter via transformers.AutoModel
  MockProteinLMEmbedder  — stdlib simulation (AA k-mer frequencies)
  select_protein_lm_embedder() — real-or-mock dispatcher

The frozen-LM-then-downstream-classifier pattern (Wong et al. 2025,
Applm) is implemented as a separate helper:
  - ApplmStyleClassifier: ESM2 embeddings → RandomForest-like scoring

All transformers-dependent paths are mocked via unittest.mock so
tests run without torch / huggingface_hub / a 135 MB download.
"""
from __future__ import annotations

import json
import unittest
from unittest.mock import MagicMock, patch

# ---------------------------------------------------------------------------
# EmbeddingRequest / EmbeddingResult dataclasses
# ---------------------------------------------------------------------------


class TestEmbeddingRequest(unittest.TestCase):
    def test_minimal_construction(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        req = EmbeddingRequest(sequences=("MKTAY",))
        self.assertEqual(req.sequences, ("MKTAY",))
        self.assertEqual(req.model_id, "facebook/esm2_t12_35M_UR50D")
        self.assertEqual(req.pooling, "mean")
        self.assertEqual(req.batch_size, 8)

    def test_rejects_empty_sequences(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        with self.assertRaises(ValueError):
            EmbeddingRequest(sequences=())

    def test_rejects_invalid_pooling(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        for bad in ("max", "min", "avg", "", "MEAN"):
            with self.assertRaises(ValueError, msg=f"pooling={bad!r} should raise"):
                EmbeddingRequest(sequences=("MKTAY",), pooling=bad)

    def test_accepts_all_pooling_modes(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        for ok in ("mean", "cls", "sum"):
            req = EmbeddingRequest(sequences=("MKTAY",), pooling=ok)
            self.assertEqual(req.pooling, ok)

    def test_rejects_invalid_amino_acid(self) -> None:
        """Sequences must contain only standard 20 amino acids."""
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        with self.assertRaises(ValueError):
            EmbeddingRequest(sequences=("MKTAYX",))  # X is non-standard
        with self.assertRaises(ValueError):
            EmbeddingRequest(sequences=("MK7AY",))  # digit is non-standard

    def test_accepts_all_20_standard_amino_acids(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        # All 20 standard AAs in one sequence
        req = EmbeddingRequest(sequences=("ACDEFGHIKLMNPQRSTVWY",))
        self.assertEqual(req.sequences[0], "ACDEFGHIKLMNPQRSTVWY")

    def test_lowercases_input_safely(self) -> None:
        """Lowercase input should be uppercased (ESM2 expects uppercase)."""
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        req = EmbeddingRequest(sequences=("mktay",))
        self.assertEqual(req.sequences[0], "MKTAY")

    def test_to_dict(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        req = EmbeddingRequest(sequences=("MKTAY",), model_id="facebook/esm2_t6_8M_UR50D")
        d = req.to_dict()
        self.assertEqual(d["model_id"], "facebook/esm2_t6_8M_UR50D")
        self.assertEqual(d["pooling"], "mean")


class TestEmbeddingResult(unittest.TestCase):
    def test_minimal_construction(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingResult

        r = EmbeddingResult(
            embeddings=((0.1, 0.2, 0.3),),
            dim=3,
            model_id="mock",
            backend="mock",
            elapsed_seconds=0.01,
        )
        self.assertEqual(r.dim, 3)
        self.assertEqual(r.embeddings[0], (0.1, 0.2, 0.3))

    def test_to_dict(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingResult

        r = EmbeddingResult(
            embeddings=((0.1, 0.2, 0.3), (0.4, 0.5, 0.6)),
            dim=3,
            model_id="facebook/esm2_t12_35M_UR50D",
            backend="transformers",
            elapsed_seconds=1.5,
            notes=("test-note",),
        )
        d = r.to_dict()
        self.assertEqual(len(d["embeddings"]), 2)
        self.assertEqual(d["dim"], 3)
        self.assertEqual(d["backend"], "transformers")
        self.assertEqual(d["notes"], ["test-note"])

    def test_dim_matches_embeddings(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingResult

        with self.assertRaises(ValueError):
            EmbeddingResult(
                embeddings=((0.1, 0.2), (0.3, 0.4)),
                dim=3,  # mismatch: actual dim is 2
                model_id="mock",
                backend="mock",
                elapsed_seconds=0.0,
            )

    def test_rejects_empty_embeddings(self) -> None:
        from mrnavax.protein_lm_protocols import EmbeddingResult

        with self.assertRaises(ValueError):
            EmbeddingResult(
                embeddings=(),
                dim=3,
                model_id="mock",
                backend="mock",
                elapsed_seconds=0.0,
            )


# ---------------------------------------------------------------------------
# ProteinLMEmbedder Protocol
# ---------------------------------------------------------------------------


class TestProteinLMEmbedderProtocol(unittest.TestCase):
    def test_protocol_is_runtime_checkable(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import ProteinLMEmbedder

        m = MockProteinLMEmbedder()
        self.assertIsInstance(m, ProteinLMEmbedder)

    def test_mock_returns_correct_shape(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder()
        req = EmbeddingRequest(sequences=("MKTAY", "ACDE"))
        result = m.embed(req)
        self.assertEqual(len(result.embeddings), 2)
        self.assertEqual(len(result.embeddings[0]), m.dim)

    def test_mock_dim_is_configurable(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        for d in (8, 64, 480, 640):
            m = MockProteinLMEmbedder(dim=d)
            req = EmbeddingRequest(sequences=("MKTAY",))
            result = m.embed(req)
            self.assertEqual(result.dim, d)
            self.assertEqual(len(result.embeddings[0]), d)


# ---------------------------------------------------------------------------
# MockProteinLMEmbedder
# ---------------------------------------------------------------------------


class TestMockProteinLMEmbedder(unittest.TestCase):
    def test_backend_name_is_mock(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder()
        req = EmbeddingRequest(sequences=("MKTAY",))
        result = m.embed(req)
        self.assertEqual(result.backend, "mock")

    def test_deterministic(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder()
        req = EmbeddingRequest(sequences=("MKTAY", "ACDEFGHIKLMNPQRSTVWY"))
        r1 = m.embed(req)
        r2 = m.embed(req)
        self.assertEqual(r1.embeddings, r2.embeddings)

    def test_different_sequences_different_embeddings(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder()
        r1 = m.embed(EmbeddingRequest(sequences=("MKTAY",)))
        r2 = m.embed(EmbeddingRequest(sequences=("ACDEF",)))
        self.assertNotEqual(r1.embeddings, r2.embeddings)

    def test_mean_pooling_doubles_length_produces_distinct(self) -> None:
        """Same AA twice vs once should produce different embeddings
        because the k-mer frequency distribution changes."""
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder()
        r1 = m.embed(EmbeddingRequest(sequences=("MK",)))
        r2 = m.embed(EmbeddingRequest(sequences=("MKMK",)))
        self.assertNotEqual(r1.embeddings, r2.embeddings)

    def test_embedding_values_are_normalized(self) -> None:
        """Mock embeddings should be L2-normalized (matching ESM2 convention)."""
        import math

        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder()
        result = m.embed(EmbeddingRequest(sequences=("MKTAYIAKLVV",)))
        for emb in result.embeddings:
            norm = math.sqrt(sum(x * x for x in emb))
            self.assertAlmostEqual(norm, 1.0, places=4)

    def test_model_id_recorded(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        m = MockProteinLMEmbedder(model_id="mock-v1")
        result = m.embed(EmbeddingRequest(sequences=("MK",)))
        self.assertEqual(result.model_id, "mock-v1")


# ---------------------------------------------------------------------------
# ESM2Embedder (real adapter, fully mocked)
# ---------------------------------------------------------------------------


class TestESM2Embedder(unittest.TestCase):
    def test_transformers_missing_raises(self) -> None:
        """If transformers/torch aren't installed, the real adapter
        must raise ESM2NotInstalled with a helpful message."""
        from mrnavax.protein_lm_adapter import ESM2Embedder, ESM2NotInstalled

        with patch.dict("sys.modules", {"transformers": None}):
            adapter = ESM2Embedder()
            from mrnavax.protein_lm_protocols import EmbeddingRequest

            with self.assertRaises(ESM2NotInstalled):
                adapter.embed(EmbeddingRequest(sequences=("MKTAY",)))

    def test_successful_embedding_calls_transformers(self) -> None:
        """When _embed_batch returns valid embeddings, the result has
        the right backend, dim, and shape."""
        from mrnavax.protein_lm_adapter import ESM2Embedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        adapter = ESM2Embedder(model_id="facebook/esm2_t12_35M_UR50D")
        # Pre-inject mocks to bypass _ensure_loaded (which requires
        # real transformers). Test-only.
        adapter._tokenizer = MagicMock()
        adapter._model = MagicMock()
        with patch.object(
            adapter,
            "_embed_batch",
            return_value=[[0.5] * 480],
        ):
            req = EmbeddingRequest(sequences=("MKTAY",), pooling="mean")
            result = adapter.embed(req)
            self.assertEqual(result.backend, "transformers")
            self.assertEqual(result.dim, 480)
            self.assertEqual(len(result.embeddings), 1)
            self.assertEqual(len(result.embeddings[0]), 480)

    def test_cls_pooling_passes_to_embed_batch(self) -> None:
        """When pooling='cls', the embedding call uses cls pool."""
        from mrnavax.protein_lm_adapter import ESM2Embedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        adapter = ESM2Embedder()
        adapter._tokenizer = MagicMock()
        adapter._model = MagicMock()
        captured_pooling = []

        def fake_embed_batch(batch, pooling):
            captured_pooling.append(pooling)
            return [[0.7] * 480]

        with patch.object(
            adapter, "_embed_batch", side_effect=fake_embed_batch
        ):
            req = EmbeddingRequest(sequences=("MKTAY",), pooling="cls")
            result = adapter.embed(req)
            self.assertEqual(captured_pooling, ["cls"])
            self.assertEqual(result.dim, 480)

    def test_sum_pooling_passes_to_embed_batch(self) -> None:
        from mrnavax.protein_lm_adapter import ESM2Embedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        adapter = ESM2Embedder()
        adapter._tokenizer = MagicMock()
        adapter._model = MagicMock()
        captured_pooling = []

        def fake_embed_batch(batch, pooling):
            captured_pooling.append(pooling)
            return [[0.3] * 480]

        with patch.object(
            adapter, "_embed_batch", side_effect=fake_embed_batch
        ):
            req = EmbeddingRequest(sequences=("MKTAY",), pooling="sum")
            result = adapter.embed(req)
            self.assertEqual(captured_pooling, ["sum"])
            self.assertEqual(result.dim, 480)

    def test_batch_size_respected(self) -> None:
        """Sequences are processed in batches of `batch_size`."""
        from mrnavax.protein_lm_adapter import ESM2Embedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        adapter = ESM2Embedder()
        adapter._tokenizer = MagicMock()
        adapter._model = MagicMock()

        def fake_embed_batch(batch, pooling):
            return [[0.1] * 480 for _ in batch]

        with patch.object(
            adapter, "_embed_batch", side_effect=fake_embed_batch
        ) as mock_batch:
            # 10 valid AA sequences with batch_size=3 → ceil(10/3) = 4 batches
            seqs = tuple(f"MKTAY{'A' * i}" for i in range(10))
            req = EmbeddingRequest(sequences=seqs, batch_size=3)
            result = adapter.embed(req)
            self.assertEqual(mock_batch.call_count, 4)
            self.assertEqual(len(result.embeddings), 10)

    def test_dim_for_known_models(self) -> None:
        """Pre-canned dim lookup covers all standard ESM2 sizes."""
        from mrnavax.protein_lm_adapter import ESM2Embedder

        self.assertEqual(
            ESM2Embedder._dim_for_model("facebook/esm2_t6_8M_UR50D"), 320
        )
        self.assertEqual(
            ESM2Embedder._dim_for_model("facebook/esm2_t12_35M_UR50D"), 480
        )
        self.assertEqual(
            ESM2Embedder._dim_for_model("facebook/esm2_t30_150M_UR50D"), 640
        )

    def test_dim_for_unknown_model_falls_back(self) -> None:
        from mrnavax.protein_lm_adapter import ESM2Embedder

        # Unknown model: fallback to ESM2-35M dim (480)
        self.assertEqual(
            ESM2Embedder._dim_for_model("facebook/esm2_unknown_X_UR50D"), 480
        )


# ---------------------------------------------------------------------------
# ApplmStyleClassifier — frozen LM + downstream scoring
# ---------------------------------------------------------------------------


class TestApplmStyleClassifier(unittest.TestCase):
    """Implements the Wong et al. 2025 pattern: frozen ESM2 embeddings
    → RandomForest-style scoring for binary classification
    (here: immunogenic vs non-immunogenic peptide)."""

    def test_classifier_returns_score_in_zero_one(self) -> None:
        from mrnavax.protein_lm_adapter import ApplmStyleClassifier
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        clf = ApplmStyleClassifier()
        emb = clf.embedder.embed(EmbeddingRequest(sequences=("MKTAY",)))
        score = clf.score(emb.embeddings[0])
        self.assertGreaterEqual(score, 0.0)
        self.assertLessEqual(score, 1.0)

    def test_classifier_deterministic(self) -> None:
        from mrnavax.protein_lm_adapter import ApplmStyleClassifier
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        clf = ApplmStyleClassifier()
        emb = clf.embedder.embed(EmbeddingRequest(sequences=("MKTAY",)))
        s1 = clf.score(emb.embeddings[0])
        s2 = clf.score(emb.embeddings[0])
        self.assertEqual(s1, s2)

    def test_classifier_with_custom_embedder(self) -> None:
        from mrnavax.protein_lm_adapter import (
            ApplmStyleClassifier,
            MockProteinLMEmbedder,
        )
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        clf = ApplmStyleClassifier(
            embedder=MockProteinLMEmbedder(dim=128)
        )
        emb = clf.embedder.embed(EmbeddingRequest(sequences=("MKTAY",)))
        self.assertEqual(emb.dim, 128)
        score = clf.score(emb.embeddings[0])
        self.assertGreaterEqual(score, 0.0)


# ---------------------------------------------------------------------------
# Backend selector
# ---------------------------------------------------------------------------


class TestBackendSelector(unittest.TestCase):
    def test_select_mock_forced(self) -> None:
        from mrnavax.protein_lm_adapter import (
            MockProteinLMEmbedder,
            select_protein_lm_embedder,
        )

        sel = select_protein_lm_embedder(prefer="mock")
        self.assertIsInstance(sel, MockProteinLMEmbedder)

    def test_select_real_forced_missing_transformers(self) -> None:
        from mrnavax.protein_lm_adapter import (
            ESM2NotInstalled,
            select_protein_lm_embedder,
        )

        with patch.dict("sys.modules", {"transformers": None}):
            with self.assertRaises(ESM2NotInstalled):
                select_protein_lm_embedder(prefer="real")

    def test_select_auto_falls_back_to_mock(self) -> None:
        from mrnavax.protein_lm_adapter import (
            MockProteinLMEmbedder,
            select_protein_lm_embedder,
        )

        with patch.dict("sys.modules", {"transformers": None}):
            sel = select_protein_lm_embedder()
            self.assertIsInstance(sel, MockProteinLMEmbedder)

    def test_select_auto_with_mock_transformers_uses_real(self) -> None:
        """When transformers is importable AND model loads, prefer real."""
        from mrnavax.protein_lm_adapter import (
            ESM2Embedder,
            select_protein_lm_embedder,
        )

        fake_transformers = MagicMock()
        fake_transformers.AutoTokenizer.from_pretrained.return_value = MagicMock()
        fake_transformers.AutoModel.from_pretrained.return_value = MagicMock()
        with patch.dict(
            "sys.modules", {"transformers": fake_transformers}
        ):
            sel = select_protein_lm_embedder()
            self.assertIsInstance(sel, ESM2Embedder)


# ---------------------------------------------------------------------------
# End-to-end integration
# ---------------------------------------------------------------------------


class TestEndToEnd(unittest.TestCase):
    def test_peptide_immunogenicity_workflow(self) -> None:
        """Full workflow: peptide → ESM2 embedding → immunogenicity score.

        The score should be in [0, 1] and deterministic across runs.
        """
        from mrnavax.protein_lm_adapter import (
            ApplmStyleClassifier,
            MockProteinLMEmbedder,
        )
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        embedder = MockProteinLMEmbedder(dim=480)
        clf = ApplmStyleClassifier(embedder=embedder)

        peptides = ("NLVPMVATV", "GILGFVFTL", "GGGGGGGGGG")
        req = EmbeddingRequest(sequences=peptides)
        emb = embedder.embed(req)
        scores = [clf.score(e) for e in emb.embeddings]
        self.assertEqual(len(scores), 3)
        for s in scores:
            self.assertGreaterEqual(s, 0.0)
            self.assertLessEqual(s, 1.0)

    def test_embedding_result_json_serializable(self) -> None:
        from mrnavax.protein_lm_adapter import MockProteinLMEmbedder
        from mrnavax.protein_lm_protocols import EmbeddingRequest

        embedder = MockProteinLMEmbedder(dim=8)
        req = EmbeddingRequest(sequences=("MKTAY",))
        result = embedder.embed(req)
        s = json.dumps(result.to_dict())
        loaded = json.loads(s)
        self.assertEqual(loaded["dim"], 8)
        self.assertEqual(len(loaded["embeddings"]), 1)


if __name__ == "__main__":
    unittest.main()
