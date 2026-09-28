"""Foundation-model embedder for scRNA-seq cells.

Provides a deterministic, stdlib-only TF-IDF + truncated SVD cell embedder
that approximates what scGPT / Geneformer / UNI-RNA do for cell-type
identification — without requiring GPU or external model downloads.

If ``scgpt`` is installed (optional extra), the real model is used instead
via the plug point in ``embed_cells``.

References
----------
Cui et al., scGPT: toward building a foundation model for single-cell
multi-omics. *Nat Methods* 21, 1480–1491 (2024).
"""

from __future__ import annotations

import math
from typing import Callable

# ---------- stdlib TF-IDF + SVD embedder ------------------------------------


def _tfidf_svd_embed(
    matrix: list[list[float]],
    *,
    n_components: int = 32,
    seed: int = 42,
) -> list[list[float]]:
    """Deterministic TF-IDF + truncated-SVD cell embedder. No numpy needed.

    Steps:
      1. Treat each gene as a "term", each cell as a "document".
      2. Compute log-normalized TF (per cell, per gene).
      3. Compute IDF (per gene across cells).
      4. Project to ``n_components`` dimensions via a deterministic SVD-style
         dimensionality reduction using truncated power iteration on the
         gene × cell × gene^T cross-product.

    Output: ``(n_cells, n_components)`` dense float matrix.
    """
    n_cells = len(matrix)
    if n_cells == 0:
        return []
    n_genes = len(matrix[0]) if matrix[0] else 0
    if n_genes == 0:
        return [[0.0] * n_components for _ in range(n_cells)]

    # 1. log-normalize each cell (cell-library normalization). Clamp negatives
    #    to 0 first (gauss-distributed synthetic data can dip below zero).
    cell_totals = [max(1e-9, sum(max(0.0, v) for v in row)) for row in matrix]
    scale = 1e4
    norm = [[max(0.0, v) / cell_totals[i] * scale for v in row] for i, row in enumerate(matrix)]
    norm = [[math.log1p(v) for v in row] for row in norm]

    # 2. IDF: fraction of cells expressing each gene
    df = [0] * n_genes
    for row in norm:
        for j, v in enumerate(row):
            if v > 0:
                df[j] += 1
    idf = [math.log((1 + n_cells) / (1 + d)) + 1 for d in df]

    # 3. TF-IDF matrix
    tfidf = [[norm[i][j] * idf[j] for j in range(n_genes)] for i in range(n_cells)]

    # 4. Truncated SVD via deterministic power iteration on tfidf^T @ tfidf.
    #    Eigenvectors of tfidf^T @ tfidf (size n_genes x n_genes) are the
    #    right singular vectors of tfidf. Project tfidf onto the top-k
    #    eigenvectors to get cell embeddings.
    rng = __import__("random").Random(seed)
    k = min(n_components, n_genes, n_cells)

    # Form the cross-product C = tfidf^T @ tfidf
    # For memory efficiency, compute (C @ v) on the fly.
    def matvec(v: list[float]) -> list[float]:
        # tfidf @ (tfidf^T @ v) — equivalent to C @ v when C = tfidf^T @ tfidf
        t = [sum(tfidf[i][j] * v[j] for j in range(n_genes)) for i in range(n_cells)]
        out = [sum(tfidf[i][j] * t[i] for i in range(n_cells)) for j in range(n_genes)]
        return out

    components: list[list[float]] = []
    for _ in range(k):
        # Random init
        v = [rng.gauss(0, 1) for _ in range(n_genes)]
        # Power iteration
        for _it in range(20):
            w = matvec(v)
            norm_w = math.sqrt(sum(x * x for x in w)) or 1.0
            v = [x / norm_w for x in w]
        # Orthogonalize against previously found components
        for prev in components:
            dot = sum(v[j] * prev[j] for j in range(n_genes))
            v = [v[j] - dot * prev[j] for j in range(n_genes)]
        norm_v = math.sqrt(sum(x * x for x in v)) or 1.0
        v = [x / norm_v for x in v]
        components.append(v)

    # Project cells onto components: cell_embed[i][c] = sum_j tfidf[i][j] * components[c][j]
    embeddings = [
        [sum(tfidf[i][j] * components[c][j] for j in range(n_genes)) for c in range(k)]
        for i in range(n_cells)
    ]
    # Pad to n_components if we had fewer usable components
    while len(embeddings[0]) < n_components:
        for row in embeddings:
            row.append(0.0)
    return embeddings


# ---------- public API -----------------------------------------------------

_EMBED_FN: Callable[[list[list[float]]], list[list[float]]] | None = None


def embed_cells(
    matrix: list[list[float]],
    *,
    model: str = "tfidf-svd",
    n_components: int = 32,
    seed: int = 42,
) -> list[list[float]]:
    """Embed cells into a fixed-dimensional vector space.

    Parameters
    ----------
    model : str
        - ``"tfidf-svd"`` (default, stdlib-only): deterministic TF-IDF + truncated SVD.
        - ``"scgpt"``: use the real scGPT foundation model. Requires the
          ``scgpt`` package and a downloaded checkpoint.
        - ``"identity"``: per-cell mean expression (one dim).
    """
    if model == "scgpt":
        # Use the real scGPT foundation model via HuggingFace weights.
        # Falls back with a clear error if weights aren't downloaded.
        try:
            from .scgpt_integration import embed_with_scgpt, scgpt_available
        except ImportError as e:
            raise RuntimeError(
                f"scGPT backend requested but scgpt_integration not importable: {e}. "
                f"Use --model tfidf-svd (default) or pip install torch transformers."
            )
        if not scgpt_available():
            raise FileNotFoundError(
                "scGPT weights not found at ~/.cache/mrnavax/. "
                "Download from https://huggingface.co/perturblab/scgpt-human "
                "(best_model.pt, vocab.json, args.json)."
            )
        return embed_with_scgpt(matrix, max_cells=64)
    if model == "identity":
        return [[sum(row) / len(row) if row else 0.0] for row in matrix]
    if model == "tfidf-svd":
        return _tfidf_svd_embed(matrix, n_components=n_components, seed=seed)
    raise ValueError(f"unknown model: {model!r}")


__all__ = [
    # Public API.
    "embed_cells",
]

