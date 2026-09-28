"""MedCPT-style dense retriever for trial matching.

Provides a deterministic TF-IDF + cosine-similarity retrieval baseline that
captures the same operational idea behind MedCPT (Jin et al., *Nat Commun*
2024): semantic matching between a patient summary and clinical-trial text.

This stdlib implementation:

1. Tokenizes with **biomedical synonym expansion** (lay → medical terms),
   so "skin cancer" matches "melanoma", "high blood pressure" matches
   "hypertension", etc.
2. Computes **TF-IDF vectors** for both the patient summary and each trial.
3. Returns trials ranked by **cosine similarity**.

If a real biomedical encoder (MedCPT, BioBERT, PubMedBERT) is installed,
swap in the dense embedder via the plug point.

Reference
---------
Jin et al., MedCPT: Contrastive Pre-trained Transformers for large-scale
zero-shot biomedical information retrieval. *Nat Commun* 15, 1785 (2024).
"""

from __future__ import annotations

import math
import re
from collections import Counter
from typing import Callable

# ---------- biomedical synonym table (subset, stdlib-only) ----------------
# Maps lay / non-technical terms to the canonical medical term(s).
# Used at tokenization time to expand the patient summary vocabulary so it
# matches the medical-terminology-heavy trial text.
SYNONYMS: dict[str, list[str]] = {
    "skin cancer": ["melanoma", "carcinoma", "neoplasm"],
    "cancer": ["neoplasm", "tumor", "carcinoma", "malignancy"],
    "tumor": ["neoplasm", "mass", "lesion"],
    "high blood pressure": ["hypertension"],
    "heart attack": ["myocardial", "infarction", "MI"],
    "stroke": ["cerebrovascular", "accident", "CVA"],
    "kidney": ["renal"],
    "liver": ["hepatic"],
    "lung": ["pulmonary", "respiratory"],
    "breast": ["mammary"],
    "bowel": ["colorectal", "colonic"],
    "blood clot": ["thrombosis", "thromboembolism", "embolism"],
    "diabetes": ["diabetic", "hyperglycemia"],
    "obesity": ["BMI", "overweight"],
    "weight loss": ["cachexia", "anorexia"],
    "chemo": ["chemotherapy"],
    "chemotherapy": ["cytotoxic"],
    "radiation": ["radiotherapy"],
    "surgery": ["resection", "surgical", "excision"],
    "stage 4": ["stage IV", "metastatic", "advanced"],
    "stage 3": ["stage III"],
    "stage 2": ["stage II"],
    "stage 1": ["stage I"],
    "spread": ["metastatic", "metastasis", "disseminated"],
    "remission": ["response", "regression"],
    "relapse": ["recurrence", "recurrent"],
    "biopsy": ["histology", "histopathologic"],
    "scan": ["imaging", "CT", "MRI", "PET"],
    "gene": ["mutation", "biomarker", "molecular"],
    "shot": ["vaccine", "vaccination", "immunization"],
    "joint pain": ["arthralgia"],
    "tired": ["fatigue", "lethargy"],
    "rash": ["dermatitis", "erythema"],
}


def _tokenize(s: str) -> list[str]:
    """Lowercase + alphanumeric tokenization."""
    return [w.lower() for w in re.findall(r"[A-Za-z][A-Za-z0-9*-]+", s)]


# Build a flat expansion: lay_term -> set of canonical expansions
_EXPANSIONS: dict[str, set[str]] = {}
for lay, canonicals in SYNONYMS.items():
    _EXPANSIONS.setdefault(lay, set()).update(canonicals)
    for c in canonicals:
        _EXPANSIONS.setdefault(c, set()).add(c)


def _tokenize(s: str) -> list[str]:
    """Lowercase + alphanumeric tokenization."""
    return [w.lower() for w in re.findall(r"[A-Za-z][A-Za-z0-9*-]+", s)]


def expand_tokens(tokens: list[str]) -> list[str]:
    """Expand a token list with biomedical synonyms.

    E.g. ``["skin", "cancer"]`` → ``["skin", "cancer", "melanoma", ...]``.
    """
    out = list(tokens)
    text = " ".join(tokens).lower()
    for phrase, expansions in _EXPANSIONS.items():
        if phrase in text:
            for e in expansions:
                if e not in out:
                    out.append(e)
    return out


# ---------- TF-IDF + cosine similarity ------------------------------------


def _tfidf(documents: list[list[str]]) -> tuple[list[list[float]], dict[str, int]]:
    """Compute TF-IDF vectors for a list of tokenized documents.

    Returns ``(vectors, vocab)`` where ``vectors[i]`` is a dense list of
    TF-IDF weights aligned to ``vocab``.
    """
    n = len(documents)
    vocab: dict[str, int] = {}
    for doc in documents:
        for t in doc:
            if t not in vocab:
                vocab[t] = len(vocab)
    # document frequency
    df: dict[str, int] = Counter()
    for doc in documents:
        for t in set(doc):
            df[t] += 1
    # TF-IDF vectors
    vectors: list[list[float]] = []
    for doc in documents:
        tf = Counter(doc)
        total = max(1, len(doc))
        v = [0.0] * len(vocab)
        for t, c in tf.items():
            idx = vocab[t]
            idf = math.log((1 + n) / (1 + df[t])) + 1
            v[idx] = (c / total) * idf
        vectors.append(v)
    return vectors, vocab


def _cosine(a: list[float], b: list[float]) -> float:
    na = math.sqrt(sum(x * x for x in a)) or 1e-9
    nb = math.sqrt(sum(x * x for x in b)) or 1e-9
    return sum(x * y for x, y in zip(a, b)) / (na * nb)


def retrieve_dense(
    patient_text: str,
    trials: list[dict],
    *,
    top_k: int = 20,
    field_weights: dict[str, float] | None = None,
) -> list[tuple[int, float]]:
    """Retrieve top-k trials by semantic similarity to the patient summary.

    Parameters
    ----------
    patient_text : str
        Free-text patient summary.
    trials : list[dict]
        Each trial must have ``title``, ``condition``, ``inclusion``,
        ``exclusion``, ``biomarkers`` keys.
    top_k : int
        Return at most this many trials.
    field_weights : dict, optional
        Weight per field when building the trial document. Defaults to:
        ``{"title": 2.0, "condition": 1.5, "inclusion": 1.0, "exclusion": 0.5,
        "biomarkers": 2.0}``.

    Returns
    -------
    list of ``(trial_index, score)`` sorted by score descending.
    """
    if not trials:
        return []
    fw = field_weights or {
        "title": 2.0,
        "condition": 1.5,
        "inclusion": 1.0,
        "exclusion": 0.5,
        "biomarkers": 2.0,
    }

    def _safe_text(v: object) -> str:
        if isinstance(v, str):
            return v
        if isinstance(v, (list, tuple)):
            return " ".join(str(x) for x in v)
        return str(v or "")

    # Build a single document per trial by concatenating weighted fields.
    # Repeat tokens according to weight (simple integer multiplier).
    trial_docs: list[list[str]] = []
    for t in trials:
        tokens: list[str] = []
        for field in ("title", "condition", "biomarkers"):
            text = _safe_text(t.get(field, ""))
            for tok in expand_tokens(_tokenize(text)):
                tokens.extend([tok] * int(round(fw.get(field, 1.0))))
        for inclusion in t.get("inclusion", []) or []:
            for tok in expand_tokens(_tokenize(_safe_text(inclusion))):
                tokens.extend([tok] * int(round(fw.get("inclusion", 1.0))))
        for exclusion in t.get("exclusion", []) or []:
            for tok in expand_tokens(_tokenize(_safe_text(exclusion))):
                tokens.extend([tok] * int(round(fw.get("exclusion", 1.0))))
        trial_docs.append(tokens)

    # Patient doc (no field weights)
    patient_doc = expand_tokens(_tokenize(patient_text))

    # Compute TF-IDF over (patient + all trials)
    all_docs = [patient_doc] + trial_docs
    vectors, _ = _tfidf(all_docs)
    p_vec = vectors[0]
    t_vecs = vectors[1:]

    scored = [(i, _cosine(p_vec, t_vecs[i])) for i in range(len(trials))]
    scored.sort(key=lambda x: -x[1])
    return scored[:top_k]


# ---------- plug point for real biomedical encoders ------------------------

_REAL_ENCODER: Callable[[list[str]], list[list[float]]] | None = None


def set_real_encoder(fn: Callable[[list[str]], list[list[float]]]) -> None:
    """Plug a real encoder (MedCPT, BioBERT, etc.) into the retriever.

    The function should take a list of strings and return one dense vector
    per string. Once set, ``retrieve_dense`` will use it instead of TF-IDF.
    """
    global _REAL_ENCODER
    _REAL_ENCODER = fn


def retrieve_real(patient_text: str, trial_texts: list[str]) -> list[float]:
    """Retrieve using a real encoder if one is plugged in. Otherwise raises."""
    if _REAL_ENCODER is None:
        raise RuntimeError(
            "no real encoder plugged in — install MedCPT or call "
            "mrnavax.medcpt_retriever.set_real_encoder(...)"
        )
    vectors = _REAL_ENCODER([patient_text] + trial_texts)
    p = vectors[0]
    scores = [_cosine(p, v) for v in vectors[1:]]
    return scores
