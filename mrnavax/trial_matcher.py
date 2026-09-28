"""TrialGPT-style patient-to-trial matcher.

End-to-end stub of the three-stage pipeline from
Jin et al. *Nature Communications* 15, 9074 (2024):

  1. TrialGPT-Retrieval — keyword generation + candidate filtering.
  2. TrialGPT-Matching — criterion-level eligibility with explanations.
  3. TrialGPT-Ranking  — aggregate criterion scores to a trial-level rank.

This is a *demonstrator* (no fine-tuned retriever, no production medical use).
The LLM does the matching + ranking; the retrieval stage is a simple keyword
overlap so the demo is fast and offline-tunable.

References
----------
TrialGPT:           https://www.nature.com/articles/s41467-024-53081-z
Biomarker LLM match: https://www.nature.com/articles/s41746-025-01673-4
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

# ---------- data shapes -----------------------------------------------------


@dataclass
class Trial:
    nct_id: str
    title: str
    condition: str
    phase: str
    inclusion: list[str]
    exclusion: list[str]
    biomarkers: list[str] = field(default_factory=list)  # e.g. ["BRAF V600E", "PD-L1>=50%"]

    @classmethod
    def from_dict(cls, d: dict) -> "Trial":
        return cls(
            nct_id=d["nct_id"],
            title=d["title"],
            condition=d.get("condition", ""),
            phase=d.get("phase", ""),
            inclusion=d.get("inclusion", []),
            exclusion=d.get("exclusion", []),
            biomarkers=d.get("biomarkers", []),
        )


@dataclass
class RankedTrial:
    nct_id: str
    title: str
    score: float
    eligibility_pct: float
    n_met: int
    n_total: int
    reasons: list[str]


# ---------- stage 1: retrieval ---------------------------------------------

STOPWORDS = set(
    "a an the of for with to and or in on at is are was were be been patient patients".split()
)


def _tokenize(s: str) -> list[str]:
    return [
        w.lower() for w in re.findall(r"[A-Za-z][A-Za-z0-9*-]+", s) if w.lower() not in STOPWORDS
    ]


def retrieve_candidates(
    patient_text: str,
    trials: list[Trial],
    *,
    top_k: int = 20,
) -> list[Trial]:
    """Keyword-overlap retrieval. Replaceable with a dense encoder in prod."""
    patient_tokens = set(_tokenize(patient_text))
    scored: list[tuple[float, Trial]] = []
    for t in trials:
        text = " ".join([t.title, t.condition] + t.inclusion + t.biomarkers)
        tokens = set(_tokenize(text))
        if not tokens:
            continue
        overlap = len(patient_tokens & tokens) / (len(tokens) ** 0.5)
        scored.append((overlap, t))
    scored.sort(key=lambda x: -x[0])
    return [t for _, t in scored[:top_k]]


# ---------- stage 2 + 3: matching + ranking --------------------------------

_MATCH_PROMPT = """You are an oncology trial eligibility screener. Decide
whether the patient meets EACH criterion. Return a JSON object with one
entry per inclusion criterion and one per exclusion criterion.

Patient summary:
{patient}

Trial NCT: {nct}
Inclusion criteria (one per item):
{inclusion}

Exclusion criteria (one per item):
{exclusion}

For each criterion output an object:
{{"criterion": "...", "verdict": "met"|"unmet"|"uncertain", "evidence": "..."}}

Return a JSON object: {{"inclusion":[...], "exclusion":[...]}}. No prose."""


def _match_one(
    trial: Trial,
    patient_text: str,
    *,
    backend: str | None = None,
) -> tuple[dict, list[str]]:
    from .llm import llm_json

    prompt = _MATCH_PROMPT.format(
        patient=patient_text,
        nct=trial.nct_id,
        inclusion="\n".join(f"- {c}" for c in trial.inclusion),
        exclusion="\n".join(f"- {c}" for c in trial.exclusion),
    )
    notes: list[str] = []
    out: dict = {"inclusion": [], "exclusion": []}
    try:
        candidate = llm_json(prompt, backend=backend)
        # Validate shape; fall back if LLM didn't follow the schema
        if (
            isinstance(candidate, dict)
            and isinstance(candidate.get("inclusion"), list)
            and isinstance(candidate.get("exclusion"), list)
        ):
            out = candidate
        else:
            notes.append("schema-mismatch fallback")
    except Exception as e:
        notes.append(f"heuristic-fallback: {e}")

    # If after LLM call inclusion/exclusion are empty, fall back to keyword overlap
    if not out["inclusion"] and not out["exclusion"]:
        inc: list[dict] = []
        patient_tokens = set(_tokenize(patient_text))
        for c in trial.inclusion:
            t = set(_tokenize(c))
            verdict = "met" if (t & patient_tokens) else "uncertain"
            inc.append({"criterion": c, "verdict": verdict, "evidence": "keyword overlap fallback"})
        exc = [
            {"criterion": c, "verdict": "uncertain", "evidence": "no LLM/fallback signal"}
            for c in trial.exclusion
        ]
        out = {"inclusion": inc, "exclusion": exc}
        if not notes:
            notes.append("keyword fallback")
    return out, notes


def rank(
    trial: Trial,
    match_result: dict,
) -> RankedTrial:
    inc = match_result.get("inclusion", [])
    exc = match_result.get("exclusion", [])
    n_met_inc = sum(1 for c in inc if c.get("verdict") == "met")
    n_total_inc = len(inc)
    n_met_exc = sum(1 for c in exc if c.get("verdict") == "met")
    elig_pct = (n_met_inc / n_total_inc * 100) if n_total_inc else 0.0
    score = elig_pct / 100.0 - 0.2 * n_met_exc
    reasons = [c["criterion"] for c in inc if c.get("verdict") == "met"][:3]
    return RankedTrial(
        nct_id=trial.nct_id,
        title=trial.title,
        score=round(score, 3),
        eligibility_pct=round(elig_pct, 1),
        n_met=n_met_inc,
        n_total=n_total_inc,
        reasons=reasons,
    )


def match(
    patient_text: str,
    trials: list[Trial],
    *,
    top_k: int = 20,
    backend: str | None = None,
    retriever: str = "keyword",
    matcher: str = "auto",
) -> tuple[list[RankedTrial], list[dict]]:
    """End-to-end pipeline: retrieve → match → rank.

    Parameters
    ----------
    retriever : str
        - ``"keyword"``: simple keyword overlap (default, fastest)
        - ``"dense"``: TF-IDF + biomedical synonym expansion (MedCPT-style)
        - ``"auto"``: dense if available, else keyword
    matcher : str
        - ``"auto"`` (default): use ``TrialGPT-style`` per-criterion LLM
          matcher when ``OPENAI_API_KEY`` is set, else keyword fallback
        - ``"trialgpt"``: force the per-criterion LLM matcher
        - ``"keyword"``: force keyword-only matching (no LLM)
    """
    if retriever == "dense" or retriever == "auto":
        # Try MedCPT first (real semantic encoder), then TF-IDF dense,
        # then keyword fallback.
        from .medcpt_integration import medcpt_available

        medcpt_ok, _ = medcpt_available()
        if medcpt_ok:
            try:
                from .medcpt_integration import retrieve_medcpt

                trial_dicts = [
                    {
                        "title": t.title,
                        "condition": t.condition,
                        "inclusion": list(t.inclusion),
                        "exclusion": list(t.exclusion),
                        "biomarkers": list(t.biomarkers),
                    }
                    for t in trials
                ]
                # Build a single concatenated article string per trial
                article_strings = [
                    " ".join([d["title"], d["condition"]] + d["inclusion"] + d["biomarkers"])
                    for d in trial_dicts
                ]
                scores = retrieve_medcpt(patient_text, article_strings)
                ranked = sorted(zip(trials, scores), key=lambda x: -x[1])[:top_k]
                candidates = [t for t, _ in ranked]
            except Exception:
                # MedCPT loaded but encode failed — fall through
                candidates = _dense_fallback(
                    patient_text,
                    trials,
                    top_k,
                    from_medcpt=True,
                )
        else:
            candidates = _dense_fallback(patient_text, trials, top_k)
    else:
        candidates = retrieve_candidates(patient_text, trials, top_k=top_k)
    out: list[RankedTrial] = []
    debug: list[dict] = []
    # Determine matcher backend
    use_trialgpt = False
    use_simicl = False
    if matcher == "trialgpt":
        use_trialgpt = True
    elif matcher == "trialgpt-simicl":
        use_trialgpt = True
        use_simicl = True
    elif matcher == "auto":
        # Use TrialGPT when an OpenAI key is present or when caller
        # explicitly forced an LLM backend. Otherwise fall back to
        # keyword matching.
        import os

        use_trialgpt = bool(os.environ.get("OPENAI_API_KEY") or (backend and backend != "mock"))

    for t in candidates:
        if use_trialgpt:
            from .trial_llm import score_trial_with_llm

            try:
                llm_result = score_trial_with_llm(
                    patient_text,
                    t.nct_id,
                    t.title,
                    list(t.inclusion),
                    list(t.exclusion),
                    backend=backend,
                    use_simicl=use_simicl,
                )
                match_result = llm_result.to_dict()
                # Trim to the shape that rank() expects
                match_result = {
                    "inclusion": llm_result.to_dict()["inclusion"],
                    "exclusion": llm_result.to_dict()["exclusion"],
                }
                notes = list(llm_result.notes)
            except Exception as e:
                # Fall back to keyword matcher on failure
                match_result, fallback_notes = _match_one(t, patient_text, backend=backend)
                notes = ["trialgpt-fallback"] + list(fallback_notes) + [str(e)]
        else:
            match_result, notes = _match_one(t, patient_text, backend=backend)
        ranked = rank(t, match_result)
        out.append(ranked)
        debug.append({"nct": t.nct_id, "match": match_result, "notes": notes})
    out.sort(key=lambda r: -r.score)
    return out, debug


def _dense_fallback(
    patient_text: str, trials: list[Trial], top_k: int, from_medcpt: bool = False
) -> list[Trial]:
    """TF-IDF dense retriever, with keyword fallback."""
    try:
        from .medcpt_retriever import retrieve_dense

        trial_dicts = [
            {
                "title": t.title,
                "condition": t.condition,
                "inclusion": list(t.inclusion),
                "exclusion": list(t.exclusion),
                "biomarkers": list(t.biomarkers),
            }
            for t in trials
        ]
        scored = retrieve_dense(patient_text, trial_dicts, top_k=top_k)
        return [trials[i] for i, _ in scored]
    except Exception:
        return retrieve_candidates(patient_text, trials, top_k=top_k)


# ---------- IO --------------------------------------------------------------


def load_trials_jsonl(path: str | Path) -> list[Trial]:
    trials: list[Trial] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            trials.append(Trial.from_dict(json.loads(line)))
    return trials


# ---------- CLI ------------------------------------------------------------


def _run_cli(argv: list[str]) -> int:
    import argparse
    import json as _json
    import os

    p = argparse.ArgumentParser(prog="mrnavax trial")
    p.add_argument("--patient", required=True)
    p.add_argument("--trials", required=True, help="JSONL of trials")
    p.add_argument("--top-k", type=int, default=10)
    p.add_argument("--backend", choices=["auto", "mock", "openai"], default="auto")
    p.add_argument(
        "--retriever",
        choices=["keyword", "dense", "auto"],
        default="keyword",
        help="retrieval method: keyword (fast) or dense (MedCPT-style semantic)",
    )
    p.add_argument(
        "--matcher",
        choices=["auto", "trialgpt", "trialgpt-simicl", "keyword"],
        default="auto",
        help="matching method: trialgpt (per-criterion LLM), "
        "trialgpt-simicl (per-criterion LLM + Sim-ICL few-shot demos "
        "from $MRNA_AI_SIMICL_DEMOS or examples/simicl_demos.json), "
        "or keyword (no LLM)",
    )
    p.add_argument("--out")
    args = p.parse_args(argv)

    patient_text = Path(args.patient).read_text() if Path(args.patient).exists() else args.patient
    trials = load_trials_jsonl(args.trials)
    backend = None if args.backend == "auto" else args.backend
    if backend is None and os.environ.get("MRNA_AI_FORCE_MOCK") == "1":
        backend = "mock"
    ranked, debug = match(
        patient_text,
        trials,
        top_k=args.top_k,
        backend=backend,
        retriever=args.retriever,
        matcher=args.matcher,
    )
    out_obj = {"ranked": [asdict(r) for r in ranked], "n_candidates_screened": len(debug)}
    out_text = _json.dumps(out_obj, indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_run_cli(sys.argv[1:]))


__all__ = [
    # Trial data dataclasses.
    "Trial",
    "RankedTrial",
    # Pipeline stages.
    "retrieve_candidates",
    "rank",
    "match",
    # IO.
    "load_trials_jsonl",
]

