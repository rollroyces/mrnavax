"""TrialGPT-style per-criterion LLM matching.

Implements the **TrialGPT-Matching** approach from Jin et al.
(Nature Communications 2024): per-criterion LLM reasoning where each
inclusion and exclusion criterion is judged independently against
the patient summary, producing a granular eligibility assessment.

References
----------
Jin, Qiao, et al. "Matching patients to clinical trials with large
language models." *Nature Communications* 15 (2024): 9074.
DOI: 10.1038/s41467-024-53081-z

Reported: 87.3% accuracy on 1,015 patient-criterion pairs, close to
expert performance.

Usage
-----
With an OPENAI_API_KEY set, the per-criterion matcher is used. Without
one, the system falls back to keyword overlap (handled in
``trial_matcher.match()``).

The :func:`score_trial_with_llm` function calls ``llm_json`` once per
trial with all criteria in a single prompt — cheaper than per-criterion
calls while preserving the granular per-criterion verdict structure.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from typing import Any

# TrialGPT-aligned prompt template. The model is asked to label each
# criterion with one of {"met", "unmet", "uncertain"} and provide a
# short evidence snippet from the patient summary.
MATCH_PROMPT_TEMPLATE = """You are an oncology clinical trial eligibility screener.
For EACH criterion below, decide if the patient meets it. Use only the
information in the patient summary; do not infer or fabricate.

Patient summary:
{patient}

Trial: {nct} — {title}

Inclusion criteria (judge each independently):
{inclusion}

Exclusion criteria (judge each independently):
{exclusion}

Return ONLY a JSON object with this exact shape:
{{
  "inclusion": [
    {{"criterion": "<text>", "verdict": "met"|"unmet"|"uncertain", "evidence": "<one short sentence>"}}
  ],
  "exclusion": [
    {{"criterion": "<text>", "verdict": "met"|"unmet"|"uncertain", "evidence": "<one short sentence>"}}
  ]
}}

Definitions:
    - "met"      : the patient clearly meets this criterion
    - "unmet"    : the patient clearly does NOT meet this criterion
    - "uncertain": cannot determine from the patient summary alone

Do not include any prose, markdown, or text outside the JSON object."""


# Verdict values accepted as inputs to the scoring.
VERDICT_MET = "met"
VERDICT_UNMET = "unmet"
VERDICT_UNCERTAIN = "uncertain"
VALID_VERDICTS = {VERDICT_MET, VERDICT_UNMET, VERDICT_UNCERTAIN}


@dataclass
class CriterionVerdict:
    """One criterion's LLM verdict.

    Attributes
    ----------
    criterion
        The criterion text.
    verdict
        One of ``"met"``, ``"unmet"``, ``"uncertain"``.
    evidence
        Short sentence explaining the verdict (from the patient summary).
    """

    criterion: str
    verdict: str
    evidence: str = ""

    def to_dict(self) -> dict:
        return {"criterion": self.criterion, "verdict": self.verdict, "evidence": self.evidence}


@dataclass
class TrialMatchResult:
    """Aggregate result of TrialGPT-style matching for one trial.

    Attributes
    ----------
    nct_id
        Trial identifier.
    inclusion_verdicts
        Per-inclusion-criterion verdicts.
    exclusion_verdicts
        Per-exclusion-criterion verdicts.
    n_met_inclusion
        Count of inclusion criteria labelled ``"met"``.
    n_total_inclusion
        Total inclusion criteria.
    n_unmet_exclusion
        Count of exclusion criteria labelled ``"unmet"`` (good — patient
        does NOT have the exclusion, so still eligible).
    n_total_exclusion
        Total exclusion criteria.
    eligibility_score
        Aggregate score, 0..1. Computed as:
        ``(n_met_inclusion / n_total_inclusion) *
        (n_unmet_exclusion / n_total_exclusion)``
        when both are > 0; falls back to one of the two when the other
        is empty.
    notes
        Free-form debug notes (e.g., "schema-mismatch fallback").
    raw_llm_response
        The raw JSON the LLM returned (for audit / debugging).
    """

    nct_id: str
    inclusion_verdicts: list[CriterionVerdict] = field(default_factory=list)
    exclusion_verdicts: list[CriterionVerdict] = field(default_factory=list)
    n_met_inclusion: int = 0
    n_total_inclusion: int = 0
    n_unmet_exclusion: int = 0
    n_total_exclusion: int = 0
    eligibility_score: float = 0.0
    notes: list[str] = field(default_factory=list)
    raw_llm_response: str = ""

    def to_dict(self) -> dict:
        return {
            "nct_id": self.nct_id,
            "n_met_inclusion": self.n_met_inclusion,
            "n_total_inclusion": self.n_total_inclusion,
            "n_unmet_exclusion": self.n_unmet_exclusion,
            "n_total_exclusion": self.n_total_exclusion,
            "eligibility_score": round(self.eligibility_score, 3),
            "notes": list(self.notes),
            "inclusion": [v.to_dict() for v in self.inclusion_verdicts],
            "exclusion": [v.to_dict() for v in self.exclusion_verdicts],
        }


def build_match_prompt(
    patient_text: str,
    nct_id: str,
    title: str,
    inclusion: list[str],
    exclusion: list[str],
) -> str:
    """Construct the TrialGPT-style per-criterion prompt."""
    inc_lines = "\n".join(f"- {c}" for c in inclusion) if inclusion else "(none)"
    exc_lines = "\n".join(f"- {c}" for c in exclusion) if exclusion else "(none)"
    return MATCH_PROMPT_TEMPLATE.format(
        patient=patient_text.strip(),
        nct=nct_id,
        title=title,
        inclusion=inc_lines,
        exclusion=exc_lines,
    )


def _parse_criteria(raw: Any, expected: list[str]) -> list[CriterionVerdict]:
    """Coerce an LLM response field into a list of CriterionVerdict.

    If the LLM did not follow the schema, fall back to one
    ``"uncertain"`` verdict per expected criterion.
    """
    out: list[CriterionVerdict] = []
    if not isinstance(raw, list):
        return [CriterionVerdict(criterion=c, verdict=VERDICT_UNCERTAIN) for c in expected]
    seen: set[str] = set()
    for item in raw:
        if not isinstance(item, dict):
            continue
        text = str(item.get("criterion", "")).strip()
        verdict = str(item.get("verdict", VERDICT_UNCERTAIN)).lower().strip()
        if verdict not in VALID_VERDICTS:
            verdict = VERDICT_UNCERTAIN
        evidence = str(item.get("evidence", "")).strip()
        # Normalize for dedup: lowercase
        key = text.lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(CriterionVerdict(criterion=text, verdict=verdict, evidence=evidence))
    # Backfill any expected criteria that the LLM omitted.
    for c in expected:
        if c.lower() not in seen:
            out.append(CriterionVerdict(criterion=c, verdict=VERDICT_UNCERTAIN))
    return out


def _eligibility_score(met_inc: int, total_inc: int, unmet_exc: int, total_exc: int) -> float:
    """Compute aggregate eligibility score 0..1.

    Multiplicative: must clear both inclusion AND not have any
    exclusion triggered. When a trial has no exclusion criteria, only
    inclusion counts.
    """
    inc_pct = (met_inc / total_inc) if total_inc else 0.0
    if total_exc == 0:
        return inc_pct
    exc_pct = unmet_exc / total_exc
    return inc_pct * exc_pct


def score_trial_with_llm(
    patient_text: str,
    nct_id: str,
    title: str,
    inclusion: list[str],
    exclusion: list[str],
    *,
    backend: str | None = None,
    demo_store: Any | None = None,
    use_simicl: bool | None = None,
) -> TrialMatchResult:
    """Run TrialGPT-style per-criterion LLM matching.

    Parameters
    ----------
    patient_text
        Free-text patient summary (diagnosis, stage, biomarkers,
        prior therapies, performance status).
    nct_id
        Trial NCT identifier.
    title
        Trial title.
    inclusion
        List of inclusion criteria.
    exclusion
        List of exclusion criteria.
    backend
        Optional LLM backend override. ``None`` = auto-detect from
        ``OPENAI_API_KEY`` (uses real OpenAI) or fall back to ``mock``.
    demo_store
        Optional :class:`~mrnavax.trial_similar.DemoStore`. When
        ``None``, loads the bundled store from
        ``examples/simicl_demos.json`` via
        :func:`~mrnavax.trial_similar.load_default_demo_store`.
        Pass an explicit ``DemoStore(demos=[])`` to disable.
    use_simicl
        Override ``$MRNA_AI_SIMICL_ENABLED``. When ``True``, the top-K
        demos (per ``$MRNA_AI_SIMICL_TOPK``, default 32) are injected
        into the prompt as few-shot examples before asking the LLM
        for verdicts. When ``False``, behaves as plain TrialGPT.

    Returns
    -------
    TrialMatchResult with per-criterion verdicts and aggregate score.
    On any failure, all verdicts are set to ``"uncertain"`` and a note
    is added.

    Notes
    -----
    Sim-ICL integration: this function implements the **Sim-ICL**
    demonstration-selection strategy from Fung et al. 2026 (Genome
    Biology, in press). The few-shot examples are chosen by TF-IDF
    cosine similarity between (patient+trial) query and the demo
    store, rather than random sampling. This matches the paper's
    finding that sequence-similar demonstrations yield competitive
    performance with protein-LM classifiers in low-shot regimes.
    """
    from .llm import llm_json  # late import to avoid circular deps
    from .trial_similar import (
        _get_enabled,
        build_simicl_prompt,
        load_default_demo_store,
    )

    # Sim-ICL: resolve config + load demos
    if demo_store is None:
        demo_store = load_default_demo_store()
    simicl_on = _get_enabled() if use_simicl is None else bool(use_simicl)

    prompt = build_match_prompt(patient_text, nct_id, title, inclusion, exclusion)
    notes: list[str] = []
    raw = ""
    parsed: dict[str, Any] = {}

    if simicl_on and len(demo_store) > 0:
        query_text = "\n".join(
            [patient_text, nct_id, title, " ".join(inclusion), " ".join(exclusion)]
        )
        top_demos = demo_store.rank(query_text)
        prompt = build_simicl_prompt(
            patient_text,
            nct_id,
            title,
            inclusion,
            exclusion,
            top_demos,
            base_prompt=prompt,
        )
        notes.append(f"simicl-k{len(top_demos)}")
        # Record which demos were used, for auditability
        notes.append("simicl-demo-ids=" + ",".join(d.demo_id for d in top_demos))

    try:
        parsed = llm_json(prompt, backend=backend)
        raw = json.dumps(parsed)
    except Exception as e:
        notes.append(f"llm-fallback: {e}")
        parsed = {}

    inclusion_v = _parse_criteria(parsed.get("inclusion"), inclusion)
    exclusion_v = _parse_criteria(parsed.get("exclusion"), exclusion)

    n_met_inc = sum(1 for v in inclusion_v if v.verdict == VERDICT_MET)
    n_total_inc = len(inclusion_v)
    n_unmet_exc = sum(1 for v in exclusion_v if v.verdict == VERDICT_UNMET)
    n_total_exc = len(exclusion_v)
    score = _eligibility_score(n_met_inc, n_total_inc, n_unmet_exc, n_total_exc)
    if not notes and (n_met_inc + n_unmet_exc == 0):
        notes.append("all-uncertain (LLM returned no confident verdicts)")

    return TrialMatchResult(
        nct_id=nct_id,
        inclusion_verdicts=inclusion_v,
        exclusion_verdicts=exclusion_v,
        n_met_inclusion=n_met_inc,
        n_total_inclusion=n_total_inc,
        n_unmet_exclusion=n_unmet_exc,
        n_total_exclusion=n_total_exc,
        eligibility_score=score,
        notes=notes,
        raw_llm_response=raw,
    )


def llm_matching_available() -> bool:
    """True iff an LLM backend can run without an explicit OPENAI_API_KEY.

    The mock backend is always available; the openai backend requires
    the env var. Used by ``backends.py`` to skip the LLM check when
    no key is present and we want to avoid a real API call.
    """
    if os.environ.get("OPENAI_API_KEY"):
        return True
    # Mock backend always available for offline runs.
    if os.environ.get("MRNA_AI_LLM_BACKEND", "auto").lower() in {"mock", "auto"}:
        return True
    return False


__all__ = [
    # Verdict string constants.
    "VERDICT_MET",
    "VERDICT_UNMET",
    "VERDICT_UNCERTAIN",
    "VALID_VERDICTS",
    # Dataclasses.
    "CriterionVerdict",
    "TrialMatchResult",
    # Public API.
    "build_match_prompt",
    "score_trial_with_llm",
    "llm_matching_available",
]

