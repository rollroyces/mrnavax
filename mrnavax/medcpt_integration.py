"""MedCPT integration plug-in.

MedCPT (Jin et al., *Nat Commun* 2024) is a contrastively-trained biomedical
sentence encoder for zero-shot clinical-trial retrieval. Two components
ship as HuggingFace models:

- ``ncbi/MedCPT-Query-Encoder`` — encodes queries (patient summaries)
- ``ncbi/MedCPT-Article-Encoder`` — encodes articles (trial text)

When both ``transformers`` and ``torch`` are installed AND the model
weights are cached locally, this module activates automatically and
overrides the stdlib TF-IDF retriever. Otherwise it stays out of the way
and the toolkit continues to run with the deterministic TF-IDF baseline.

Setup
-----
    pip install torch transformers
    python -c "from transformers import AutoModel; AutoModel.from_pretrained('ncbi/MedCPT-Query-Encoder')"
    python -c "from transformers import AutoModel; AutoModel.from_pretrained('ncbi/MedCPT-Article-Encoder')"

The first call downloads ~440 MB of weights and caches them in
``~/.cache/huggingface/``.

Reference
---------
Jin, Q., Leeman, R., Liu, Y., et al. MedCPT: Contrastive Pre-trained
Transformers for large-scale zero-shot biomedical information retrieval.
*Nat Commun* 15, 1785 (2024).
"""

from __future__ import annotations

import os
import threading

# Lazy-loaded module references. Imports happen on first encode call to
# avoid paying the ~5 s transformers-import cost for users who never call it.
_QUERY_MODEL = None
_ARTICLE_MODEL = None
_QUERY_TOKENIZER = None
_ARTICLE_TOKENIZER = None
_LOCK = threading.Lock()
_AVAILABLE: bool | None = None
_AVAIL_ERROR: str | None = None


def medcpt_available() -> tuple[bool, str]:
    """Return ``(True, '')`` if MedCPT can be loaded, else ``(False, reason)``.

    Cached after the first call. Cheap (no model download).
    """
    global _AVAILABLE, _AVAIL_ERROR
    if _AVAILABLE is not None:
        return _AVAILABLE, _AVAIL_ERROR or ""

    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
    except ImportError as e:
        _AVAILABLE = False
        _AVAIL_ERROR = f"missing dependency: {e}. Install torch + transformers."
        return _AVAILABLE, _AVAIL_ERROR

    # Allow the operator to skip the model check via env var (the CI runner
    # would otherwise download 440 MB every run).
    if os.environ.get("MRNA_AI_SKIP_MEDCPT_CHECK") == "1":
        _AVAILABLE = False
        _AVAIL_ERROR = "skipped via MRNA_AI_SKIP_MEDCPT_CHECK=1"
        return _AVAILABLE, _AVAIL_ERROR

    try:
        from transformers import AutoModel, AutoTokenizer  # noqa: F401
    except Exception as e:
        _AVAILABLE = False
        _AVAIL_ERROR = f"transformers import failed: {e}"
        return _AVAILABLE, _AVAIL_ERROR

    # Optional: confirm the model ID exists locally. We don't actually load
    # the weights here — that happens on the first encode call. Just verify
    # the identifier resolves.
    try:
        from huggingface_hub import HfApi  # type: ignore

        api = HfApi()
        # HEAD request on the model files — will raise if unreachable
        list(api.list_repo_files("ncbi/MedCPT-Query-Encoder"))
    except Exception as e:
        _AVAILABLE = False
        _AVAIL_ERROR = (
            f"MedCPT model not reachable on HuggingFace: {e}. "
            "Set MRNA_AI_SKIP_MEDCPT_CHECK=1 to force-disable this check."
        )
        return _AVAILABLE, _AVAIL_ERROR

    _AVAILABLE = True
    _AVAIL_ERROR = ""
    return _AVAILABLE, _AVAIL_ERROR


def _load_models() -> tuple[object, object]:
    """Lazy-load MedCPT encoders. Idempotent and thread-safe."""
    global _QUERY_MODEL, _ARTICLE_MODEL, _QUERY_TOKENIZER, _ARTICLE_TOKENIZER
    if _QUERY_MODEL is not None:
        return
    with _LOCK:
        if _QUERY_MODEL is not None:
            return
        from transformers import AutoModel, AutoTokenizer

        _QUERY_MODEL = AutoModel.from_pretrained("ncbi/MedCPT-Query-Encoder")
        _QUERY_TOKENIZER = AutoTokenizer.from_pretrained("ncbi/MedCPT-Query-Encoder")
        _ARTICLE_MODEL = AutoModel.from_pretrained("ncbi/MedCPT-Article-Encoder")
        _ARTICLE_TOKENIZER = AutoTokenizer.from_pretrained("ncbi/MedCPT-Article-Encoder")


def encode_query(text: str) -> list[float]:
    """Encode a patient summary into a dense MedCPT query vector."""
    _load_models()
    import torch

    inputs = _QUERY_TOKENIZER(text, return_tensors="pt", truncation=True, max_length=512)
    with torch.no_grad():
        vec = _QUERY_MODEL(**inputs).last_hidden_state[:, 0, :]  # CLS token
    return vec[0].cpu().tolist()


def encode_article(text: str) -> list[float]:
    """Encode a trial description into a dense MedCPT article vector."""
    _load_models()
    import torch

    inputs = _ARTICLE_TOKENIZER(text, return_tensors="pt", truncation=True, max_length=512)
    with torch.no_grad():
        vec = _ARTICLE_MODEL(**inputs).last_hidden_state[:, 0, :]
    return vec[0].cpu().tolist()


def retrieve_medcpt(patient_text: str, trial_texts: list[str]) -> list[float]:
    """Encode patient + trials with MedCPT and return cosine similarities.

    Raises ``RuntimeError`` if MedCPT is not available.
    """
    available, reason = medcpt_available()
    if not available:
        raise RuntimeError(f"MedCPT unavailable: {reason}")

    import math

    def cos(a: list[float], b: list[float]) -> float:
        na = math.sqrt(sum(x * x for x in a)) or 1e-9
        nb = math.sqrt(sum(x * x for x in b)) or 1e-9
        return sum(x * y for x, y in zip(a, b)) / (na * nb)

    p_vec = encode_query(patient_text)
    scores = []
    for t in trial_texts:
        t_vec = encode_article(t)
        scores.append(cos(p_vec, t_vec))
    return scores
