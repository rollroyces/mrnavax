"""Sim-ICL: sequence-similarity demonstration selection for biological LLM matching.

Implements the **Sim-ICL** strategy from Fung et al. 2026 (Genome
Biology, in press): when a downstream task is solved by in-context
learning, **selecting demonstrations by similarity to the query** —
not random sampling — yields competitive performance with protein-LM
classifiers, especially in low-shot regimes.

Reference
---------
Fung S.H., Zhang Z., Wang R., Miao C., Wong B.S.H., Li K.Y., Hong C.,
Zhou J., Yip K.Y.#, Tsui S.K.W.#, and Cao Q.#. (2026) A Systematic
Evaluation of In-Context Learning in Large Language Models for Antibody
Characterization. Genome Biology (in press).

Applied to mrnavax
--------------------------
The toolkit's :mod:`trial_llm` module already uses an LLM to judge
each clinical-trial criterion against a patient summary (the TrialGPT
pattern, Jin et al. 2024). Sim-ICL enhances this by **injecting
top-K prior patient-trial-eligibility triples** as few-shot examples
into the prompt, ranked by similarity to the current query.

For our task the natural similarity is **TF-IDF cosine similarity
over the concatenated (patient + trial) text** — stdlib-only,
deterministic, and matches the approach already used in
:mod:`medcpt_retriever`. The paper used sequence similarity (BLOSUM
or k-mer) for antibodies; here we extend the same idea to **mixed
text similarity** for clinical-trial eligibility.

Demo store
----------
We ship a small JSON of synthetic ``DemoCase`` triples
(``examples/simicl_demos.json``) so Sim-ICL has something to rank
against out of the box. Real users with access to labelled
patient-trial pairs can swap in their own store by either:

  - pointing ``MRNA_AI_SIMICL_DEMOS`` at a JSON file of the same
    shape, or
  - programmatically constructing a :class:`DemoStore` and passing it
    into :func:`build_simicl_prompt`.
"""

from __future__ import annotations

import json
import math
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# Default config (env-var overrideable)
# ---------------------------------------------------------------------------


DEFAULT_TOPK = 32  # paper default (Fung et al. 2026)
ENV_TOPK = "MRNA_AI_SIMICL_TOPK"
ENV_ENABLED = "MRNA_AI_SIMICL_ENABLED"
ENV_DEMOS = "MRNA_AI_SIMICL_DEMOS"


def _get_topk(default: int = DEFAULT_TOPK) -> int:
    """Read MRNA_AI_SIMICL_TOPK with validation. Falls back to default."""
    raw = os.environ.get(ENV_TOPK, str(default))
    try:
        k = int(raw)
        return k if k > 0 else default
    except (TypeError, ValueError):
        return default


def _get_enabled() -> bool:
    """Read MRNA_AI_SIMICL_ENABLED. Default True."""
    raw = os.environ.get(ENV_ENABLED, "1").strip().lower()
    return raw not in {"0", "false", "no", "off"}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DemoCase:
    """A single few-shot demonstration for Sim-ICL.

    Attributes
    ----------
    demo_id
        Stable identifier (e.g., ``"demo_001"``). Used for traceability
        in the prompt and audit logs.
    patient_text
        Free-text patient summary (same shape as the query).
    trial_nct
        Trial NCT identifier.
    trial_title
        Trial title.
    inclusion
        Bullet-pointed inclusion criteria.
    exclusion
        Bullet-pointed exclusion criteria.
    ground_truth_verdicts
        The known-eligibility label for this patient-trial pair, as a
        dict ``{"eligible": bool, "n_met": int, "n_total": int}``.
        Optional — when None, the case is used for similarity ranking
        only and the LLM produces its own verdicts.
    """

    demo_id: str
    patient_text: str
    trial_nct: str
    trial_title: str
    inclusion: tuple[str, ...]
    exclusion: tuple[str, ...]
    ground_truth_verdicts: dict | None = None

    def to_dict(self) -> dict:
        return {
            "demo_id": self.demo_id,
            "patient_text": self.patient_text,
            "trial_nct": self.trial_nct,
            "trial_title": self.trial_title,
            "inclusion": list(self.inclusion),
            "exclusion": list(self.exclusion),
            "ground_truth_verdicts": self.ground_truth_verdicts,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "DemoCase":
        return cls(
            demo_id=str(d["demo_id"]),
            patient_text=str(d["patient_text"]),
            trial_nct=str(d["trial_nct"]),
            trial_title=str(d["trial_title"]),
            inclusion=tuple(d.get("inclusion") or ()),
            exclusion=tuple(d.get("exclusion") or ()),
            ground_truth_verdicts=d.get("ground_truth_verdicts"),
        )

    @property
    def query_text(self) -> str:
        """Concatenated text used for similarity ranking.

        Joining patient + trial fields with newlines gives the TF-IDF
        ranker enough signal to distinguish e.g. (BRAF V600E patient,
        BRAF-inhibitor trial) from (KRAS G12D patient, NSCLC trial).
        """
        return "\n".join(
            [
                self.patient_text,
                self.trial_nct,
                self.trial_title,
                " ".join(self.inclusion),
                " ".join(self.exclusion),
            ]
        )


@dataclass
class DemoStore:
    """In-memory store of DemoCase triples.

    The store is intentionally lightweight: a list + a TF-IDF index
    rebuilt on demand. For >10k demos, swap in a real vector store
    (FAISS, Milvus) — the public API stays the same.

    Attributes
    ----------
    demos
        The list of :class:`DemoCase` triples.
    """

    demos: list[DemoCase] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.demos)

    def __iter__(self):
        return iter(self.demos)

    @classmethod
    def from_json(cls, path: str | Path) -> "DemoStore":
        """Load a store from a JSON file of the canonical shape.

        The JSON may be either a list of demo dicts, or a dict with
        a ``"demos"`` key wrapping that list.
        """
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(
                f"DemoStore JSON not found: {p}. Pass an absolute path or "
                f"set $MRNA_AI_SIMICL_DEMOS to a file with the canonical "
                f"shape (list of {{demo_id, patient_text, trial_nct, "
                f"trial_title, inclusion, exclusion, ground_truth_verdicts}})."
            )
        raw = json.loads(p.read_text())
        items = raw["demos"] if isinstance(raw, dict) and "demos" in raw else raw
        if not isinstance(items, list):
            raise ValueError(
                f"DemoStore JSON must be a list or a dict with a 'demos' "
                f"list; got {type(items).__name__} in {p}"
            )
        return cls(demos=[DemoCase.from_dict(x) for x in items])

    def add(self, demo: DemoCase) -> None:
        """Append a single demo (in-place)."""
        self.demos.append(demo)

    def rank(self, query: str, k: int | None = None) -> list[DemoCase]:
        """Return the top-``k`` demos by TF-IDF cosine similarity to ``query``.

        Parameters
        ----------
        query
            Free-text query — typically ``patient_text + " " + trial_title``.
        k
            Number of demos to return. Defaults to ``_get_topk()``.
            If ``k >= len(self.demos)``, returns all demos sorted.

        Returns
        -------
        list[DemoCase]
            Sorted by descending similarity score. Ties broken by
            ``demo_id`` for determinism. Returns an empty list if the
            store is empty.
        """
        if not self.demos:
            return []
        if k is None:
            k = _get_topk()
        if k >= len(self.demos):
            k = len(self.demos)
        scored = _rank_by_tfidf_cosine(query, self.demos)
        # Sort by score desc, then by demo_id asc (stable, deterministic).
        scored.sort(key=lambda x: (-x[1], x[0].demo_id))
        return [d for d, _ in scored[:k]]


# ---------------------------------------------------------------------------
# TF-IDF cosine ranker (stdlib-only)
# ---------------------------------------------------------------------------


_TOKEN_RE = re.compile(r"[a-z0-9]+")


def _tokenize(text: str) -> list[str]:
    """Lowercase alphanum tokenization. Same vocabulary as medcpt_retriever."""
    return _TOKEN_RE.findall(text.lower())


def _rank_by_tfidf_cosine(query: str, demos: Iterable[DemoCase]) -> list[tuple[DemoCase, float]]:
    """Rank ``demos`` against ``query`` by TF-IDF cosine similarity.

    Algorithm: standard textbook TF-IDF with cosine similarity.
        tf(t, d)  = count(t, d) / sum(count(*, d))  (per-doc normalization)
        idf(t)    = log((N + 1) / (df(t) + 1)) + 1   (smoothed, ≥1 always)
        cosine(q, d) = dot(q, d) / (||q|| * ||d||)

    Pure-stdlib, deterministic, fast enough for thousands of demos
    in CI. For >10k demos, swap in scikit-learn's TfidfVectorizer
    via an optional extra; the public API stays the same.
    """
    docs: list[list[str]] = [_tokenize(d.query_text) for d in demos]
    q_tokens = _tokenize(query)
    if not q_tokens:
        return [(d, 0.0) for d in demos]

    # Document frequencies
    df: Counter[str] = Counter()
    for toks in docs:
        for t in set(toks):
            df[t] += 1
    n = len(docs)

    def idf(t: str) -> float:
        return math.log((n + 1) / (df[t] + 1)) + 1.0

    # Query vector
    q_tf: Counter[str] = Counter(q_tokens)
    q_len = len(q_tokens)
    q_vec: dict[str, float] = {t: (q_tf[t] / q_len) * idf(t) for t in q_tf}
    q_norm = math.sqrt(sum(v * v for v in q_vec.values())) or 1.0

    # Doc vectors + dot products
    scored: list[tuple[DemoCase, float]] = []
    for demo, toks in zip(demos, docs):
        if not toks:
            scored.append((demo, 0.0))
            continue
        tf_d: Counter[str] = Counter(toks)
        d_len = len(toks)
        d_vec: dict[str, float] = {t: (tf_d[t] / d_len) * idf(t) for t in tf_d}
        # Dot product over intersection
        dot = sum(q_vec.get(t, 0.0) * v for t, v in d_vec.items())
        d_norm = math.sqrt(sum(v * v for v in d_vec.values())) or 1.0
        score = dot / (q_norm * d_norm)
        scored.append((demo, score))
    return scored


# ---------------------------------------------------------------------------
# Sim-ICL prompt construction
# ---------------------------------------------------------------------------


# Header prepended to the few-shot examples block in the prompt.
FEWSHOT_HEADER = """You are an oncology clinical trial eligibility screener.
For EACH example below, the prior verdict shows what the correct
eligibility outcome was for that (patient, trial) pair. Use these as
calibration when judging the new pair at the end of this prompt."""


def _format_demo_block(demo: DemoCase) -> str:
    """Render one demo as a few-shot example block."""
    gt = demo.ground_truth_verdicts or {}
    eligible = gt.get("eligible")
    if eligible is None:
        verdict_str = "(verdict not provided)"
    elif eligible:
        verdict_str = f"ELIGIBLE (n_met={gt.get('n_met', '?')}/{gt.get('n_total', '?')})"
    else:
        verdict_str = f"NOT ELIGIBLE (n_met={gt.get('n_met', '?')}/{gt.get('n_total', '?')})"
    inc = "\n".join(f"- {c}" for c in demo.inclusion) or "(none)"
    exc = "\n".join(f"- {c}" for c in demo.exclusion) or "(none)"
    return (
        f"=== Example [{demo.demo_id}] ===\n"
        f"Patient: {demo.patient_text.strip()}\n"
        f"Trial: {demo.trial_nct} — {demo.trial_title}\n"
        f"Inclusion criteria:\n{inc}\n"
        f"Exclusion criteria:\n{exc}\n"
        f"Prior verdict: {verdict_str}\n"
    )


def build_simicl_prompt(
    patient_text: str,
    nct_id: str,
    title: str,
    inclusion: list[str],
    exclusion: list[str],
    demos: list[DemoCase],
    *,
    base_prompt: str,
) -> str:
    """Prepend a few-shot examples block to ``base_prompt``.

    Parameters
    ----------
    patient_text, nct_id, title, inclusion, exclusion
        The query patient-trial pair (same as ``build_match_prompt``).
    demos
        Top-K demos ranked by similarity (output of ``DemoStore.rank``).
        When empty or Sim-ICL is disabled, returns ``base_prompt``
        unchanged.
    base_prompt
        The original TrialGPT prompt (output of ``build_match_prompt``).

    Returns
    -------
    str
        ``base_prompt`` with a few-shot block prepended, or
        ``base_prompt`` unchanged when ``demos`` is empty.
    """
    if not demos or not _get_enabled():
        return base_prompt
    fewshot_block = (
        FEWSHOT_HEADER
        + "\n\n"
        + "\n\n".join(_format_demo_block(d) for d in demos)
        + "\n\n=== Now your turn ===\n"
    )
    # Splice the few-shot block in immediately before the
    # "Return ONLY a JSON object" instruction, so the model sees
    # the calibration examples first, then the schema reminder.
    marker = "Return ONLY a JSON object"
    if marker in base_prompt:
        return base_prompt.replace(marker, fewshot_block + marker, 1)
    # Fallback: prepend at the top.
    return fewshot_block + base_prompt


# ---------------------------------------------------------------------------
# Convenience: load the bundled demo store
# ---------------------------------------------------------------------------


def load_default_demo_store() -> DemoStore:
    """Load the shipped demo store from ``examples/simicl_demos.json``.

    Returns an empty store (with a single sentinel ``DemoCase`` removed
    by the loader) when the file doesn't exist — never raises. CI
    tests rely on this behavior to avoid a hard dep on the example
    file.
    """
    candidates = [
        Path(__file__).parent.parent / "examples" / "simicl_demos.json",
        Path("/Users/hermes/mrnavax/examples/simicl_demos.json"),
    ]
    env = os.environ.get(ENV_DEMOS)
    if env:
        candidates.insert(0, Path(env))
    for c in candidates:
        if c.exists():
            try:
                return DemoStore.from_json(c)
            except (json.JSONDecodeError, KeyError, ValueError):
                continue
    return DemoStore(demos=[])


__all__ = [
    # Default config constants.
    "DEFAULT_TOPK",
    "ENV_TOPK",
    "ENV_ENABLED",
    "ENV_DEMOS",
    # Few-shot prompt header.
    "FEWSHOT_HEADER",
    # Dataclasses.
    "DemoCase",
    "DemoStore",
    # Public API.
    "build_simicl_prompt",
    "load_default_demo_store",
]

