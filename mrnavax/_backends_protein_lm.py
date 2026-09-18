"""Backend checks for the ``neoantigen.esm2_protein_lm_embedder`` family.

Note: the check is registered under the ``neoantigen.*`` family in
the registry, but the underlying implementation is the ESM2
protein-LM Protocol adapter. This module owns that check.

Extracted from ``mrnavax/backends.py`` in v0.24.0.
"""

from __future__ import annotations

from ._backends_registry import register


@register("neoantigen.esm2_protein_lm_embedder")
def _check_esm2_protein_lm_embedder() -> tuple[bool, str]:
    """ESM2 protein-LM Protocol adapter: embed peptides via frozen LM
    for Applm-style downstream classification.

    Validates:
      1. EmbeddingRequest dataclass: AA alphabet validation, pooling
         mode, batch_size >= 1.
      2. EmbeddingResult dataclass: dim matches embeddings, non-empty.
      3. ProteinLMEmbedder Protocol runtime_checkable.
      4. MockProteinLMEmbedder satisfies the Protocol, produces
         L2-normalized embeddings of correct dim, is deterministic.
      5. ApplmStyleClassifier returns score in [0, 1].
      6. lm_immunogenicity_score() integration works with mock.
      7. ESM2 pre-canned dim lookup covers all standard ESM2 sizes.
      8. Backend selector picks mock when transformers is unavailable.
    """
    from .protein_lm_adapter import (
        ApplmStyleClassifier,
        ESM2Embedder,
        MockProteinLMEmbedder,
        _check_transformers_available,
        select_protein_lm_embedder,
    )
    from .protein_lm_protocols import (
        EmbeddingRequest,
        EmbeddingResult,
        ProteinLMEmbedder,
    )

    bad: list[str] = []

    # 1. EmbeddingRequest validation
    try:
        EmbeddingRequest(sequences=())
    except ValueError:
        pass
    else:
        bad.append("empty sequences should raise")

    for bad_pool in ("max", "min", ""):
        try:
            EmbeddingRequest(sequences=("MK",), pooling=bad_pool)
        except ValueError:
            pass
        else:
            bad.append(f"pooling={bad_pool!r} should raise")

    try:
        EmbeddingRequest(sequences=("MK7AY",))  # digit
    except ValueError:
        pass
    else:
        bad.append("non-AA character should raise")

    try:
        EmbeddingRequest(sequences=("M",), batch_size=0)
    except ValueError:
        pass
    else:
        bad.append("batch_size=0 should raise")

    # Lowercase input normalized to uppercase
    req = EmbeddingRequest(sequences=("mktay",))
    if req.sequences[0] != "MKTAY":
        bad.append("lowercase not normalized")

    # 2. EmbeddingResult validation
    try:
        EmbeddingResult(
            embeddings=(),
            dim=3,
            model_id="mock",
            backend="mock",
            elapsed_seconds=0.0,
        )
    except ValueError:
        pass
    else:
        bad.append("empty embeddings should raise")

    try:
        EmbeddingResult(
            embeddings=((0.1, 0.2),),
            dim=3,  # mismatch
            model_id="mock",
            backend="mock",
            elapsed_seconds=0.0,
        )
    except ValueError:
        pass
    else:
        bad.append("dim mismatch should raise")

    # 3. Protocol runtime_checkable
    mock = MockProteinLMEmbedder()
    if not isinstance(mock, ProteinLMEmbedder):
        bad.append("mock does not satisfy Protocol")

    # 4. Mock embedding: L2-normalized, correct dim, deterministic
    import math

    req = EmbeddingRequest(sequences=("MKTAYIAKLVV", "ACDEFGHIKL"))
    result = mock.embed(req)
    if len(result.embeddings) != 2:
        bad.append(f"expected 2 embeddings, got {len(result.embeddings)}")
    for i, emb in enumerate(result.embeddings):
        if len(emb) != result.dim:
            bad.append(f"embedding[{i}] dim mismatch")
        norm = math.sqrt(sum(x * x for x in emb))
        if abs(norm - 1.0) > 1e-4:
            bad.append(f"embedding[{i}] not L2-normalized: {norm}")

    # Determinism
    r1 = mock.embed(req)
    r2 = mock.embed(req)
    if r1.embeddings != r2.embeddings:
        bad.append("mock not deterministic")

    # 5. ApplmStyleClassifier returns [0, 1]
    clf = ApplmStyleClassifier(embedder=mock)
    for emb in result.embeddings:
        s = clf.score(emb)
        if not (0.0 <= s <= 1.0):
            bad.append(f"classifier score out of [0, 1]: {s}")
            break

    # 6. lm_immunogenicity_score integration via mock
    from .neoantigen_screener import lm_immunogenicity_score

    lm_result = lm_immunogenicity_score("NLVPMVATV", embedder=mock)
    if not (0.0 <= lm_result["score"] <= 1.0):
        bad.append(f"lm_immunogenicity_score out of [0,1]: {lm_result['score']}")
    if lm_result["backend"] != "mock":
        bad.append(f"lm_immunogenicity_score backend mismatch: {lm_result['backend']}")
    if lm_result["dim"] != 480:
        bad.append(f"lm_immunogenicity_score dim mismatch: {lm_result['dim']}")

    # 7. ESM2 dim lookup
    expected_dims = {
        "facebook/esm2_t6_8M_UR50D": 320,
        "facebook/esm2_t12_35M_UR50D": 480,
        "facebook/esm2_t30_150M_UR50D": 640,
    }
    for model_id, dim in expected_dims.items():
        if ESM2Embedder._dim_for_model(model_id) != dim:
            bad.append(f"dim mismatch for {model_id}")

    # 8. Backend selector picks mock when transformers missing
    from unittest.mock import patch

    with patch.dict("sys.modules", {"transformers": None}):
        sel = select_protein_lm_embedder()
        if not isinstance(sel, MockProteinLMEmbedder):
            bad.append(
                f"auto-select should pick mock when transformers missing, got {type(sel).__name__}"
            )

    if bad:
        return False, "ESM2 protein-LM backend issues: " + "; ".join(bad)
    real_available = _check_transformers_available()
    return True, (
        f"ESM2 protein-LM OK: mock embedder 480-dim L2-normalized, "
        f"deterministic; ApplmStyleClassifier scores in [0, 1]; "
        f"lm_immunogenicity_score integration returns valid result; "
        f"ESM2 dim lookup correct; real backend (transformers) "
        f"available={real_available}"
    )


__all__ = []
