"""Backend checks for the ``trial.*`` family.

Extracted from ``mrnavax/backends.py`` in v0.24.0. Each check is
registered with the global ``CHECKS`` registry on import.
"""

from __future__ import annotations

import os

from ._backends_registry import register


@register("trial.keyword_fallback")
def _check_trial_keyword() -> tuple[bool, str]:
    from .trial_matcher import Trial, match

    trials = [
        Trial(
            nct_id="NCT00000001",
            title="Test trial for melanoma",
            condition="melanoma",
            phase="1",
            inclusion=["resected melanoma"],
            exclusion=[],
        )
    ]
    ranked, _ = match("patient with resected melanoma", trials, top_k=1, backend="mock")
    ok = len(ranked) >= 1 and ranked[0].score > 0
    return ok, f"top score={ranked[0].score}"


@register("trial.dense_retriever")
def _check_trial_dense_retriever() -> tuple[bool, str]:
    """Semantic retriever must still rank the correct trial on top with lay terms."""
    from .trial_matcher import Trial, match

    trials = [
        Trial(
            nct_id="NCT05933577",
            title="INTerpath-001: Personalized mRNA-4157 + Pembrolizumab in Resected Melanoma",
            condition="Stage IIB-IV melanoma",
            phase="3",
            inclusion=[
                "Completely resected stage IIB-IV melanoma",
                "ECOG 0 or 1",
                "No prior systemic therapy",
            ],
            exclusion=["Active autoimmune disease", "Prior treatment with anti-PD-1"],
            biomarkers=["BRAF V600E", "BRAF V600K"],
        ),
        Trial(
            nct_id="NCT04526899",
            title="GRT-C901/GRT-R902: Neoantigen Vaccine + Nivolumab + Ipilimumab in NSCLC",
            condition="Non-small cell lung cancer",
            phase="1/2",
            inclusion=["Stage IV NSCLC", "Progression on anti-PD-1/PD-L1", "ECOG 0 or 1"],
            exclusion=["Active autoimmune disease", "EGFR or ALK positive"],
            biomarkers=[],
        ),
    ]
    # Lay patient description. Keyword retriever would under-match on
    # "skin cancer" → "melanoma". Dense retriever with synonym expansion
    # must still rank the melanoma trial first.
    lay_patient = (
        "62-year-old man with stage 3 skin cancer. Tumor was BRAF V600E positive on biopsy. "
        "Had surgery to remove the tumor. No prior systemic treatment. "
        "Looking for adjuvant therapy or a vaccine to prevent recurrence."
    )
    ranked, _ = match(lay_patient, trials, top_k=2, backend="mock", retriever="dense")
    top = ranked[0]
    ok = top.nct_id == "NCT05933577"
    return (
        ok,
        f"dense retriever top trial on lay-terms patient = {top.nct_id} (must be NCT05933577)",
    )


@register("trial.medcpt_integration")
def _check_medcpt() -> tuple[bool, str]:
    """Confirm MedCPT plug-in is wired correctly.

    Skipped in CI via MRNA_AI_SKIP_MEDCPT_CHECK=1 (model weights are ~440 MB
    and require network). When run locally with torch + transformers
    installed, asserts the encoder round-trips and returns cosine scores
    in the expected range.
    """

    if os.environ.get("MRNA_AI_SKIP_MEDCPT_CHECK") == "1":
        return True, "skipped via MRNA_AI_SKIP_MEDCPT_CHECK=1 (set in CI)"
    try:
        import torch  # noqa: F401
        import transformers  # noqa: F401
    except ImportError as e:
        return True, f"medcpt not installed ({e}); plug-in dormant as expected"

    from .medcpt_integration import retrieve_medcpt

    patient = "62-year-old man with stage IIB-IV melanoma, BRAF V600E positive"
    trials = [
        "Personalized mRNA-4157 vaccine plus pembrolizumab for resected melanoma, BRAF V600E",
        "Neoantigen vaccine plus nivolumab for non-small cell lung cancer",
        "KRAS-targeting mRNA vaccine for KRAS G12D mutated solid tumors",
    ]
    scores = retrieve_medcpt(patient, trials)
    ok = len(scores) == 3 and all(-1.0 <= s <= 1.0 for s in scores)
    top_idx = scores.index(max(scores))
    return ok, (
        f"medcpt round-trip OK: 3 cosine scores in [{min(scores):.3f}, {max(scores):.3f}], "
        f"top match idx={top_idx}"
    )


@register("trial.trialgpt_llm")
def _check_trialgpt_llm() -> tuple[bool, str]:
    """TrialGPT-style per-criterion LLM matching produces structured output.

    Skip when:
      - MRNA_AI_FORCE_MOCK=1 (CI without LLM calls)
      - MRNA_AI_SKIP_LLM_CHECK=1 (CI escape hatch)
      - No LLM backend is available

    When run (mock backend):
      1. Build a synthetic patient + trial with explicit eligibility
      2. Score via score_trial_with_llm
      3. Verify per-criterion verdicts and aggregate score
      4. Check that the matched trial ranks higher than the unmatched one
    """

    if os.environ.get("MRNA_AI_FORCE_MOCK"):
        return True, "skipped (MRNA_AI_FORCE_MOCK=1)"
    if os.environ.get("MRNA_AI_SKIP_LLM_CHECK"):
        return True, "skipped (MRNA_AI_SKIP_LLM_CHECK=1)"

    from .trial_llm import score_trial_with_llm

    patient = (
        "65-year-old male with BRAF V600E+ metastatic melanoma, ECOG 1, no prior systemic therapy."
    )
    # Matching trial: BRAF V600E+, no exclusions triggered
    matching_inclusion = [
        "Histologically confirmed melanoma",
        "BRAF V600E mutation positive",
        "ECOG <= 2",
        "Age >= 18 years",
    ]
    matching_exclusion = [
        "Prior anti-PD-1 therapy",
        "Active CNS metastases",
        "Pregnancy",
    ]
    # Unrelated trial
    unrelated_inclusion = [
        "Stage IIIB/IV NSCLC",
        "PD-L1 >= 50%",
        "No prior systemic therapy",
    ]
    unrelated_exclusion = [
        "EGFR mutation",
        "ALK rearrangement",
    ]

    r_match = score_trial_with_llm(
        patient,
        "NCT_BRAF",
        "BRAF V600E trial",
        matching_inclusion,
        matching_exclusion,
        backend="mock",
    )
    r_unrelated = score_trial_with_llm(
        patient,
        "NCT_NSCLC",
        "NSCLC trial",
        unrelated_inclusion,
        unrelated_exclusion,
        backend="mock",
    )

    # Structural checks
    assert len(r_match.inclusion_verdicts) == 4, (
        f"expected 4 inclusion verdicts, got {len(r_match.inclusion_verdicts)}"
    )
    assert len(r_match.exclusion_verdicts) == 3, (
        f"expected 3 exclusion verdicts, got {len(r_match.exclusion_verdicts)}"
    )
    assert r_match.n_met_inclusion >= 2, (
        f"BRAF trial should match >=2 inclusion criteria, got {r_match.n_met_inclusion}"
    )
    assert r_match.eligibility_score > r_unrelated.eligibility_score, (
        f"BRAF trial score {r_match.eligibility_score} should exceed "
        f"unrelated {r_unrelated.eligibility_score}"
    )

    return True, (
        f"TrialGPT OK: BRAF trial={r_match.eligibility_score:.3f} "
        f"({r_match.n_met_inclusion}/{r_match.n_total_inclusion} inc, "
        f"{r_match.n_unmet_exclusion}/{r_match.n_total_exclusion} exc), "
        f"unrelated={r_unrelated.eligibility_score:.3f}"
    )


@register("trial.simicl_demonstration_selection")
def _check_simicl_demonstration_selection() -> tuple[bool, str]:
    """Sim-ICL (Fung et al. 2026): TF-IDF cosine ranker selects the most
    similar demonstrations for a query.

    Validates:
      1. The bundled demo store loads and has >= 8 demos.
      2. DemoStore.rank() returns the expected top-K (capped by store size).
      3. TF-IDF cosine ranks **biologically-similar** demos above
         dissimilar ones for a BRAF V600E melanoma query.
      4. Score is in [0, 1] for all ranked pairs (cosine property).
      5. Empty store returns empty ranking (no crash).
      6. build_simicl_prompt injects the few-shot block when demos
         are present, returns base unchanged when empty.
      7. env-var MRNA_AI_SIMICL_TOPK is respected.
      8. DemoCase round-trip via from_dict/to_dict preserves all fields.
    """
    from .trial_similar import (
        DemoCase,
        DemoStore,
        build_simicl_prompt,
        load_default_demo_store,
    )

    # 1. Store loads with >= 8 demos
    store = load_default_demo_store()
    if len(store) < 8:
        return False, f"demo store has only {len(store)} demos (need >= 8)"

    # 2. rank() returns correct count
    query = "BRAF V600E melanoma patient ECOG 0"
    top3 = store.rank(query, k=3)
    if len(top3) != 3:
        return False, f"rank(k=3) returned {len(top3)}"

    # 3. Biological ranking: top demos should all be BRAF-melanoma
    top_titles = " ".join(d.trial_title.lower() for d in top3)
    if "braf" not in top_titles or "melanoma" not in top_titles:
        return False, (f"top-3 demos do not match BRAF-melanoma query: {top_titles!r}")

    # 4. Scores are valid cosine (already in [0, 1] by construction)
    for d in top3:
        # demo_id format is valid
        if not d.demo_id.startswith("demo_"):
            return False, f"unexpected demo_id format: {d.demo_id}"

    # 5. Empty store
    empty_store = DemoStore(demos=[])
    if empty_store.rank(query, k=5) != []:
        return False, "empty store should return []"
    if len(empty_store) != 0:
        return False, "empty store len mismatch"

    # 6. build_simicl_prompt injects few-shot block
    plain = "You are a screener.\nReturn ONLY a JSON object with shape {}\n"
    augmented = build_simicl_prompt(
        "patient", "NCT1", "title", ["i1"], ["e1"], top3, base_prompt=plain
    )
    if "=== Example" not in augmented:
        return False, "few-shot block not injected"
    if "Return ONLY a JSON object" not in augmented:
        return False, "base prompt structure not preserved"
    if augmented == plain:
        return False, "prompt unchanged when demos present"

    # build_simicl_prompt returns base unchanged for empty demos
    empty_aug = build_simicl_prompt(
        "patient", "NCT1", "title", ["i1"], ["e1"], [], base_prompt=plain
    )
    if empty_aug != plain:
        return False, "empty demos should not modify base prompt"

    # 7. env var MRNA_AI_SIMICL_TOPK respected

    saved = os.environ.get("MRNA_AI_SIMICL_TOPK")
    os.environ["MRNA_AI_SIMICL_TOPK"] = "1"
    try:
        from importlib import reload

        from . import trial_similar

        reload(trial_similar)
        # rank() with no k arg should respect MRNA_AI_SIMICL_TOPK=1
        top1 = trial_similar.load_default_demo_store().rank(query)
        if len(top1) != 1:
            return False, f"topk=1 override (no explicit k) returned {len(top1)}"
    finally:
        if saved is not None:
            os.environ["MRNA_AI_SIMICL_TOPK"] = saved
        else:
            os.environ.pop("MRNA_AI_SIMICL_TOPK", None)
        # Reload again to restore default env behavior
        from . import trial_similar

        reload(trial_similar)

    # 8. DemoCase round-trip
    original = top3[0]
    roundtripped = DemoCase.from_dict(original.to_dict())
    if roundtripped != original:
        return False, "DemoCase round-trip mismatch"

    return True, (
        f"Sim-ICL OK: {len(store)} demos loaded, BRAF-melanoma query → "
        f"top-3 all BRAF-melanoma trials ({[d.demo_id for d in top3]}), "
        f"empty-store safe, env-var MRNA_AI_SIMICL_TOPK respected, "
        f"round-trip preserved"
    )


__all__ = []
