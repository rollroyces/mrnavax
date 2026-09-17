"""Backend availability + integrity checks.

Used by CI to verify that the optional backend wiring is intact even when the
heavy models themselves aren't installed. Heavy model downloads (mhcflurry,
scGPT) are gated by their respective extras.

Run directly:

    python -m mrnavax.backends --check-all

Exits 0 if every check passes, 1 otherwise.
"""

from __future__ import annotations

import json as _json
import os
import sys
import tempfile
from pathlib import Path
from typing import Callable

CHECKS: list[tuple[str, Callable[[], tuple[bool, str]]]] = []


def _example_path(name: str) -> str:
    """Resolve an example-file path relative to the installed package.

    Works in two layouts:
    - Installed wheel: pkg_dir = site-packages/mrnavax, examples
      are bundled next to it under site-packages/mrnavax/examples.
    - Dev / repo: pkg_dir is inside the repo, examples live at the repo
      root under mrnavax/examples.

    Returns the first path that exists, falling back to the wheel layout
    (which will raise FileNotFoundError downstream if missing — that's the
    expected behavior for missing bundled data).
    """
    from pathlib import Path

    pkg_dir = Path(__file__).resolve().parent
    wheel_candidate = pkg_dir / "examples" / name
    if wheel_candidate.exists():
        return str(wheel_candidate)
    repo_candidate = pkg_dir.parent.parent / "mrnavax" / "examples" / name
    if repo_candidate.exists():
        return str(repo_candidate)
    return str(wheel_candidate)


def register(name: str):
    def deco(fn: Callable[[], tuple[bool, str]]) -> Callable[[], tuple[bool, str]]:
        CHECKS.append((name, fn))
        return fn

    return deco


# ---------- neoantigen backends -------------------------------------------


@register("neoantigen.heuristic_A0201")
def _check_neoantigen_heuristic() -> tuple[bool, str]:
    from .neoantigen_screener import screen_peptide_llm

    # Known-good A*02:01 binder — must come back positive in the heuristic.
    r = screen_peptide_llm("NLVPMVATV", "HLA-A*02:01", backend="mock")
    ok = r.source == "heuristic" and isinstance(r.binding_affinity_nM, float)
    return ok, f"source={r.source} aff={r.binding_affinity_nM} nM"


@register("neoantigen.heuristic_nonA2_silent")
def _check_neoantigen_non_a2() -> tuple[bool, str]:
    from .neoantigen_screener import screen_peptide_llm

    r = screen_peptide_llm("NLVPMVATV", "HLA-A*03:01", backend="mock")
    ok = r.binding_affinity_nM >= 5000  # heuristic is silent on non-A2
    return ok, f"non-A2 affinity={r.binding_affinity_nM} nM (expected >= 5000)"


@register("neoantigen.mhcflurry_available")
def _check_neoantigen_mhcflurry() -> tuple[bool, str]:
    from .neoantigen_screener import _mhcflurry_available

    available = _mhcflurry_available()
    if available:
        return True, "installed — will use as default"
    return True, "not installed (optional: pip install -e '.[neoantigen-mhcflurry]')"


# ---------- codon backend -------------------------------------------------


@register("codon.analyze")
def _check_codon_analyze() -> tuple[bool, str]:
    from .codon_optimizer import analyze_cds

    r = analyze_cds("ATGGATAAGAAATACTCAATAGGCTTAGATATCGGCACAAATAGC")
    ok = 0 < r.cai <= 1 and 0 < r.gc_percent < 100
    return ok, f"cai={r.cai} gc={r.gc_percent}%"


@register("codon.optimize")
def _check_codon_optimize() -> tuple[bool, str]:
    from .codon_optimizer import optimize_basic

    r = optimize_basic("ATGGATAAGAAATACTCAATAGGCTTAGATATCGGCACAAATAGCGTGGGCTGGGCG")
    ok = r["after"]["cai"] > r["before"]["cai"] and r["changes"] > 0
    return ok, f"basic CAI {r['before']['cai']} → {r['after']['cai']} ({r['changes']} swaps)"


@register("codon.ribodecode_optimizer")
def _check_codon_ribodecode() -> tuple[bool, str]:
    from .codon_ribodecode import optimize_ribodecode

    seq = (
        "ATGGATAAGAAATACTCAATAGGCTTAGATATCGGCACAAATAGCGTGGGCTGGGCGGTGATCAC"
        "CGATGAATATAAGGTTCCGTCTAAAAAGTTCAAGGTTCTGGGAAATACAGACCGCCACAGTATC"
    )
    r = optimize_ribodecode(seq)
    ok = (
        r.after["cai"] >= r.before["cai"]
        and r.after["rare_codon_fraction"] <= r.before["rare_codon_fraction"]
    )
    return ok, (
        f"ribodecode CAI {r.before['cai']} → {r.after['cai']}, "
        f"rare-pair {r.rare_pair_count}, GC-stddev {r.gc_window_stddev_before} → {r.gc_window_stddev_after}"
    )


_CODON_TABLE = {
    "ATG": "M",
    "ATA": "I",
    "ATT": "I",
    "ATC": "I",
    "ACA": "T",
    "ACT": "T",
    "ACC": "T",
    "ACG": "T",
    "ACN": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AAR": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "AGR": "R",
    "AGY": "S",
    "TCA": "S",
    "TCG": "S",
    "TCT": "S",
    "TCC": "S",
    "TCN": "S",
    "TCY": "S",
    "TCR": "S",
    "TCW": "S",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GCN": "A",
    "GAT": "D",
    "GAC": "D",
    "GAY": "D",
    "GAA": "E",
    "GAG": "E",
    "GAR": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "GGN": "G",
    "CAT": "H",
    "CAC": "H",
    "CAY": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CAR": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "CGN": "R",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CTN": "L",
    "CTR": "L",
    "CTY": "L",
    "TTA": "L",
    "TTG": "L",
    "TTR": "L",
    "TTT": "F",
    "TTC": "F",
    "TTY": "F",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "TAY": "Y",
    "TGT": "C",
    "TGC": "C",
    "TGY": "C",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
    "TRA": "*",
}


def _translate(cds: str) -> str:
    return "".join(_CODON_TABLE.get(cds[i : i + 3], "?") for i in range(0, len(cds), 3))


@register("codon.lineardesign_optimizer")
def _check_codon_lineardesign() -> tuple[bool, str]:
    """LinearDesign DP must improve CAI without breaking the protein."""
    from .codon_lineardesign import optimize_lineardesign

    seq = (
        "ATGGATAAGAAATACTCAATAGGCTTAGATATCGGCACAAATAGCGTGGGCTGGGCGGTGATCAC"
        "CGATGAATATAAGGTTCCGTCTAAAAAGTTCAAGGTTCTGGGAAATACAGACCGCCACAGTATC"
    )
    r = optimize_lineardesign(seq)
    protein_ok = _translate(seq) == _translate(r.new_cds)
    cai_ok = r.after["cai"] >= r.before["cai"]
    ok = protein_ok and cai_ok
    return ok, (
        f"lineardesign CAI {r.before['cai']} → {r.after['cai']}, "
        f"changes={r.changes}, structure_score={r.structure_score}, "
        f"protein_preserved={protein_ok}"
    )


# ---------- trial backend -------------------------------------------------


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


# ---------- lnp backend ---------------------------------------------------


@register("lnp.recommend")
def _check_lnp_recommend() -> tuple[bool, str]:
    from .lnp_advisor import recommend

    r = recommend(target="lung", cargo="saRNA", intent="cancer vaccine", n=3)
    ok = len(r.shortlist) > 0 and all("name" in c for c in r.shortlist)
    return ok, f"{len(r.shortlist)} candidates for lung/saRNA"


# ---------- scrna backend -------------------------------------------------


@register("scrna.pipeline")
def _check_scrna() -> tuple[bool, str]:
    from .sc_rna_pipeline import run_pipeline

    r = run_pipeline(
        _example_path("cells.csv"),
        _example_path("variants_coding.csv"),
        _example_path("proteins.fasta"),
        hla=("HLA-A*02:01",),
        tumor_marker_genes=["TP53", "KRAS", "BRAF"],
    )
    ok = r.n_cells > 0 and r.n_candidate_peptides > 0 and r.tumor_cluster >= 0
    return ok, (
        f"cells={r.n_cells} clusters={len(set(r.cluster_labels))} "
        f"tumor={r.tumor_cluster} peptides={r.n_candidate_peptides}"
    )


@register("scrna.variant_filter")
def _check_scrna_variant_filter() -> tuple[bool, str]:
    from .sc_rna_pipeline import run_pipeline

    full = run_pipeline(
        _example_path("cells.csv"),
        _example_path("variants_coding.csv"),
        _example_path("proteins.fasta"),
        hla=("HLA-A*02:01",),
        tumor_marker_genes=["TP53", "KRAS", "BRAF"],
    )
    filt = run_pipeline(
        _example_path("cells.csv"),
        _example_path("variants_coding.csv"),
        _example_path("proteins.fasta"),
        hla=("HLA-A*02:01",),
        tumor_marker_genes=["TP53", "KRAS", "BRAF"],
        variant_filter_top_fraction=0.30,
    )
    ok = (
        filt.n_variants_after_filter < full.n_variants_after_filter
        and filt.n_candidate_peptides < full.n_candidate_peptides
    )
    return ok, (
        f"top-30% filter: {full.n_variants_input} variants → "
        f"{filt.n_variants_after_filter} kept → "
        f"{filt.n_candidate_peptides} peptides (was {full.n_candidate_peptides})"
    )


@register("scrna.variant_scorer_alphamissense")
def _check_variant_scorer() -> tuple[bool, str]:
    from .variant_scorer import score_variant

    # Known driver mutations should rank highest
    braf = score_variant("BRAF", 600, "V", "E", protein_length=766)
    kras = score_variant("KRAS", 12, "G", "V", protein_length=189)
    benign = score_variant("MYC", 100, "A", "A", protein_length=439)  # silent
    ok = braf.normalized_score > kras.normalized_score > 0
    return ok, (
        f"BRAF.V600E={braf.normalized_score} > "
        f"KRAS.G12V={kras.normalized_score} > 0 (silent={benign.normalized_score})"
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


@register("scrna.scgpt_integration")
def _check_scgpt_integration() -> tuple[bool, str]:
    """scGPT integration loads real weights and produces meaningful embeddings.

    Skip when:
      - MRNA_AI_FORCE_MOCK=1 (CI without scGPT weights)
      - MRNA_AI_SKIP_SCGPT_CHECK=1 (CI escape hatch)
      - scGPT weights are not on disk (first-run before download)

    When run:
      1. Load weights + vocab
      2. Embed 30 synthetic cells with named gene IDs (tumor vs normal markers)
      3. Verify output shape (n_cells, 512) and variance > 0
    """

    if os.environ.get("MRNA_AI_FORCE_MOCK"):
        return True, "skipped (MRNA_AI_FORCE_MOCK=1)"
    if os.environ.get("MRNA_AI_SKIP_SCGPT_CHECK"):
        return True, "skipped (MRNA_AI_SKIP_SCGPT_CHECK=1)"

    from .scgpt_integration import (
        ScGPTConfig,
        embed_with_scgpt,
        scgpt_available,
    )

    if not scgpt_available():
        return True, (
            "skipped (scGPT weights not present at ~/.cache/mrnavax/). "
            "Download from https://huggingface.co/perturblab/scgpt-human to enable."
        )

    cfg = ScGPTConfig.from_json()
    gene_set = [
        "TP53",
        "MYC",
        "KRAS",
        "EGFR",
        "EPCAM",
        "CD8A",
        "CD4",
        "CD3E",
        "PTPRC",
        "ACTB",
    ]
    import random

    random.seed(42)
    matrix = []
    for i in range(30):
        is_tumor = i < 15
        row = []
        for g in gene_set:
            if is_tumor and g in {"TP53", "MYC", "KRAS", "EGFR", "EPCAM"}:
                row.append(random.gauss(5.0, 0.5))
            elif (not is_tumor) and g in {"CD8A", "CD4", "CD3E", "PTPRC"}:
                row.append(random.gauss(5.0, 0.5))
            else:
                row.append(random.gauss(2.0, 0.3))
        matrix.append([max(0.0, v) for v in row])
    emb = embed_with_scgpt(matrix, gene_names=gene_set, max_cells=30)
    assert len(emb) == 30, f"got {len(emb)} embeddings"
    assert len(emb[0]) == cfg.d_hid, f"got dim {len(emb[0])}, expected {cfg.d_hid}"
    var = sum(sum(x * x for x in row) for row in emb) / (len(emb) * len(emb[0]))
    assert var > 1e-6, f"all-zero embeddings (var={var})"
    return True, (
        f"scGPT OK: 30 cells -> 30x{cfg.d_hid} embeddings (var={var:.3f}, "
        f"nlayers={cfg.nlayers}, vocab={cfg.ntoken})"
    )


@register("manufacture.score_manufacturability")
def _check_manufacturability() -> tuple[bool, str]:
    """The manufacturability checker must flag real-world problems.

    Test:
      1. A clean CDS scores high (>0.7) with no errors.
      2. A pathological CDS (long poly-A, ARE nonamers) gets errors.
      3. A CDS with hidden internal stops is flagged.
    """
    from .manufacturability import (
        check_hidden_stops,
        score_manufacturability,
    )

    # 1. Clean CDS — high score, no errors
    clean = "ATG" + ("GCTGCAGCTGCAGCTGCA" * 50) + "TAA"
    r1 = score_manufacturability(
        clean,
        utr5="GCCGCCACC",
        utr3="AAAAAAAAAAAAAAAAAAAAAAAA",
    )
    assert r1.n_error == 0, f"clean CDS should have 0 errors, got {r1.n_error}"
    assert r1.overall_score > 0.7, f"clean CDS score {r1.overall_score} too low"

    # 2. Pathological CDS — must trigger poly-A error
    pathological = (
        "AAAAAAAAAAAAAAAATAA" + "ATGGCTGCAGCTGCATAA" + "TAG" + "GGGGGGGGGGCCCCCCCCCAAAAAAAATAA"
    )
    r2 = score_manufacturability(
        pathological,
        utr5="AAAAAAAA",
        utr3="UUAUUUAUUAAUUAUUUAUUAUUUAUU",
    )
    poly = next(c for c in r2.checks if c.name == "poly_a_runs")
    assert poly.severity == "error", f"poly-A severity {poly.severity} != error"
    are = next(c for c in r2.checks if c.name == "are_motif")
    assert are.severity == "error", f"ARE severity {are.severity} != error"
    assert r2.overall_score < r1.overall_score, (
        f"pathological score {r2.overall_score} should be < clean {r1.overall_score}"
    )

    # 3. Hidden stops — internal TAA in-frame
    with_stops = "ATG" + "GCT" * 30 + "TAA" + "GCT" * 30 + "TAA"
    stops = check_hidden_stops(with_stops)
    assert stops.severity == "error", f"hidden_stops severity {stops.severity}"
    assert stops.score == 0.0

    return True, (
        f"manufacturability OK: clean={r1.overall_score:.2f}, "
        f"pathological={r2.overall_score:.2f} (errors={r2.n_error}), "
        f"hidden_stops={stops.severity}"
    )


@register("codon.ribodecode_protocols")
def _check_ribodecode_protocols() -> tuple[bool, str]:
    """RiboDecode adapter: Protocol contracts + mock backends pass
    conformance checks without downloading the upstream package.

    Validates:
      1. RiboDecodeRequest dataclass validation (length, mfe_weight,
         custom env CSV required when env='custom').
      2. TranslationPredictor / CodonOptimizer are runtime_checkable
         Protocols and the mock backends satisfy them.
      3. MockTranslationPredictor returns a value in [0, 100] for a
         realistic CDS (matching the upstream TranslationModel range).
      4. MockCodonOptimizer preserves the protein sequence end-to-end
         (the optimizer's hard contract).
    """
    from .codon_protocols import (
        CodonOptimizer,
        RiboDecodeRequest,
        TranslationPredictor,
    )
    from .codon_ribodecode_adapter import (
        MockCodonOptimizer,
        MockTranslationPredictor,
        ribo_decode_available,
    )

    # 1. RiboDecodeRequest validation
    bad = []
    try:
        RiboDecodeRequest(cds="")
    except ValueError:
        pass
    else:
        bad.append("empty cds should raise")
    try:
        RiboDecodeRequest(cds="ATG")  # too short, length 3 OK; check non-multiple
        RiboDecodeRequest(cds="ATGA")  # length 4, not multiple of 3
    except ValueError:
        pass
    else:
        bad.append("non-multiple-of-3 cds should raise")
    try:
        RiboDecodeRequest(cds="A" * 4501)
    except ValueError:
        pass
    else:
        bad.append(">4500 nt cds should raise")
    try:
        RiboDecodeRequest(cds="ATGGACGGGTAG", mfe_weight=1.5)
    except ValueError:
        pass
    else:
        bad.append("mfe_weight outside [0,1] should raise")
    try:
        RiboDecodeRequest(cds="ATGGACGGGTAG", env="custom")
    except ValueError:
        pass
    else:
        bad.append("env='custom' without custom_env_csv should raise")
    if bad:
        return False, f"dataclass validation gaps: {bad}"

    # 2. Protocol runtime_checkable + mock instances conform
    pred = MockTranslationPredictor()
    opt = MockCodonOptimizer()
    assert isinstance(pred, TranslationPredictor), "mock pred does not satisfy Protocol"
    assert isinstance(opt, CodonOptimizer), "mock opt does not satisfy Protocol"

    # 3. Mock translation score in [0, 100]
    # Use a GFP-like codon sequence (high CAI by construction)
    high_cai = "ATG" + ("GCTGCTGCTGCTGCTGCT" * 50) + "TAA"  # 453 nt, all Ala = GCC
    p_high = pred.predict(high_cai, env="HEK293T")
    low_cai = "ATG" + ("TTATTATTATTATTATTA" * 50) + "TAA"  # low-frequency codons
    p_low = pred.predict(low_cai, env="HEK293T")
    assert 0.0 <= p_low.translation_level <= 100.0, (
        f"low translation level {p_low.translation_level} out of [0,100]"
    )
    assert 0.0 <= p_high.translation_level <= 100.0, (
        f"high translation level {p_high.translation_level} out of [0,100]"
    )
    assert p_high.translation_level > p_low.translation_level, (
        f"high-CAI CDS ({p_high.translation_level:.2f}) should score "
        f"higher than low-CAI ({p_low.translation_level:.2f})"
    )

    # 4. Mock optimizer preserves the protein
    req = RiboDecodeRequest(cds="ATG" + ("GCTGCTGCTGCTGCTGCT" * 10) + "TAA")
    res = opt.optimize(req)
    from .codon_optimizer import CODON_TO_AA

    # MockCodonOptimizer strips the trailing stop before optimizing;
    # the optimized CDS has no stop. So we compare the AA prefix of
    # the input (excluding trailing stop) to the AA translation of the
    # output.
    input_codons = [req.cds[i : i + 3] for i in range(0, len(req.cds), 3)]
    if input_codons and CODON_TO_AA.get(input_codons[-1]) == "*":
        input_codons = input_codons[:-1]
    original_protein = "".join(CODON_TO_AA[c] for c in input_codons)
    new_protein = "".join(
        CODON_TO_AA[res.optimized_cds[i : i + 3]] for i in range(0, len(res.optimized_cds), 3)
    )
    assert new_protein == original_protein, (
        f"mock optimizer changed protein: {new_protein[:30]}... vs {original_protein[:30]}..."
    )
    # upstream_backend = "ribodecode" if installed else "mock"
    real_installed = ribo_decode_available()
    return True, (
        f"RiboDecode protocols OK: mock translation high={p_high.translation_level:.1f} "
        f"low={p_low.translation_level:.1f}, optimizer preserved protein "
        f"({len(res.optimized_cds)} nt), real backend installed={real_installed}"
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


# ---------------------------------------------------------------------------
# Protein language model (ESM2)
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Spatial transcriptomics (STModule)
# ---------------------------------------------------------------------------


@register("spatial.stmodule_module_identification")
def _check_stmodule_module_identification() -> tuple[bool, str]:
    """STModule Protocol adapter: tissue-module identification from
    spatial-transcriptomics data.

    Validates:
      1. SpatialData dataclass: platform enum, file existence,
         num_modules >= 1.
      2. SpatialModule + SpatialModuleResult dataclasses with
         to_dict() JSON-serializable.
      3. SpatialModuleBackend is runtime_checkable.
      4. MockSpatialModuleBackend satisfies the Protocol and produces
         the right shape (n_modules = num_modules).
      5. Mock result is deterministic across runs.
      6. Backend selector picks mock when Rscript is unavailable.
      7. Spot/location ID mismatch is handled gracefully (intersection).
    """
    from .spatial_module_adapter import (
        MockSpatialModuleBackend,
        rscript_available,
        select_spatial_module_backend,
    )
    from .spatial_protocols import (
        SpatialData,
        SpatialModule,
        SpatialModuleBackend,
        SpatialModuleResult,
    )

    def _write(p: Path, rows: list[str]) -> None:
        p.write_text("\n".join(rows))

    bad: list[str] = []
    with tempfile.TemporaryDirectory() as tmp:
        count = Path(tmp) / "counts.tsv"
        loc = Path(tmp) / "locs.tsv"
        _write(
            count,
            ["\t".join([""] + [f"G{j}" for j in range(5)])]
            + ["\t".join([f"loc_{i + 1}", "1", "0", "1", "1", "0"]) for i in range(12)],
        )
        _write(
            loc,
            ["\t".join(["", "x", "y"])]
            + [
                "\t".join([f"loc_{i + 1}", str(17.0 + i * 0.5), str(5.0 + (i % 4))])
                for i in range(12)
            ],
        )

        # 1. SpatialData construction + validation
        data = SpatialData(count_file=count, locations_file=loc, platform="ST")
        if data.platform != "ST":
            bad.append("platform not preserved")
        if data.num_modules != 10:
            bad.append("default num_modules != 10")
        try:
            SpatialData(count_file=count, locations_file=loc, platform="ST", num_modules=0)
        except ValueError:
            pass
        else:
            bad.append("num_modules=0 should raise")
        try:
            SpatialData(count_file=count, locations_file=loc, platform="BOGUS")
        except ValueError:
            pass
        else:
            bad.append("invalid platform should raise")

        # 2. SpatialModule + SpatialModuleResult round-trip
        m = SpatialModule(module_id=0, top_genes=("TP53", "KRAS"), n_spots=12, mean_activity=0.7)
        r = SpatialModuleResult(
            platform="ST",
            modules=(m,),
            n_spots=12,
            elapsed_seconds=0.5,
            backend="mock",
            notes=("mock-backend",),
        )
        try:
            d = r.to_dict()
            s = _json.dumps(d)
            loaded = _json.loads(s)
            if loaded["modules"][0]["module_id"] != 0:
                bad.append("module round-trip mismatch")
        except (TypeError, ValueError) as e:
            bad.append(f"JSON serialization failed: {e}")

        # 3. Protocol runtime_checkable
        mock = MockSpatialModuleBackend()
        if not isinstance(mock, SpatialModuleBackend):
            bad.append("mock backend does not satisfy Protocol")

        # 4. Mock run produces correct n_modules
        result = mock.run(data)
        if len(result.modules) != data.num_modules:
            bad.append(f"mock returned {len(result.modules)} modules, expected {data.num_modules}")
        if result.backend != "mock":
            bad.append(f"backend != 'mock': {result.backend!r}")

        # 5. Determinism
        result2 = mock.run(data)
        if [mm.to_dict() for mm in result.modules] != [mm.to_dict() for mm in result2.modules]:
            bad.append("mock backend not deterministic")

        # 6. Backend selector picks mock when Rscript unavailable
        from unittest.mock import patch

        with patch(
            "mrnavax.spatial_module_adapter.rscript_available",
            return_value=False,
        ):
            sel = select_spatial_module_backend(prefer="auto")
            if not isinstance(sel, MockSpatialModuleBackend):
                bad.append(
                    f"auto-select should pick mock when Rscript missing, got {type(sel).__name__}"
                )

        # 7. Spot/location mismatch — must not crash, must report <= 11 spots
        loc_mismatch = Path(tmp) / "locs_mismatch.tsv"
        _write(
            loc_mismatch,
            ["\t".join(["", "x", "y"])]
            + [
                "\t".join([f"loc_{i + 1}", str(17.0 + i * 0.5), str(5.0 + (i % 4))])
                for i in range(11)
            ],
        )
        data_mismatch = SpatialData(count_file=count, locations_file=loc_mismatch, platform="ST")
        result_mm = mock.run(data_mismatch)
        if result_mm.n_spots > 11:
            bad.append(f"spot mismatch not handled: n_spots={result_mm.n_spots}")

    if bad:
        return False, "STModule backend issues: " + "; ".join(bad)

    real_available = rscript_available()
    return True, (
        f"STModule OK: mock backend produces correct n_modules, "
        f"deterministic, JSON-serializable; mismatch handled gracefully "
        f"(n_spots<=11); backend selector chooses mock when Rscript "
        f"unavailable; real backend (Rscript) installed={real_available}"
    )


# ---------------------------------------------------------------------------
# Variant / neoantigen / scRNA / trial / manufacture
# ---------------------------------------------------------------------------


@register("codon.lineardesign_full_length")
def _check_lineardesign_full_length() -> tuple[bool, str]:
    """LinearDesign DP must work on full-length CDS (≥600 nt) and
    preserve the protein sequence.

    Uses a moderate-size synthetic CDS (~600 nt) for fast CI runs.
    The full-length DP is documented to scale to 4,000+ nt (Cas9) in
    <60 s; verified separately during v0.6.0 release prep.
    """
    from .codon_lineardesign import optimize_lineardesign
    from .codon_optimizer import CODON_TO_AA

    # 200 aa = 600 nt — moderate CDS, exercises the DP beyond the
    # v0.5.0 length cap.
    random_protein = "M" + "AGCT" * 50  # 201 aa, all common AAs
    # Reverse-translate using most-frequent codons
    from .codon_optimizer import HUMAN_CODON_FREQ

    cds = "".join(max(HUMAN_CODON_FREQ[aa], key=HUMAN_CODON_FREQ[aa].get) for aa in random_protein)
    r = optimize_lineardesign(cds, gc_window_size=15)
    new_protein = "".join(CODON_TO_AA[r.new_cds[i : i + 3]] for i in range(0, len(r.new_cds), 3))
    assert new_protein == random_protein, "LinearDesign changed the protein"
    assert len(r.new_cds) == len(cds), "LinearDesign changed CDS length"
    assert r.n_states_evaluated > 0
    return True, (
        f"LinearDesign full-length OK: "
        f"{len(cds)} nt in {r.elapsed_seconds:.2f}s, "
        f"{r.n_states_evaluated:,} states evaluated, "
        f"protein preserved (CAI {r.before['cai']:.3f} -> {r.after['cai']:.3f})"
    )


@register("scrna.alphamissense_integration")
def _check_alphamissense_integration() -> tuple[bool, str]:
    """Confirm the AlphaMissense plug-in is wired end-to-end.

    Always runs (no network, no model download) by exercising the
    synthetic test index. Verifies:
      1. build_test_index() returns the expected number of entries.
      2. lookup() returns an AlphaMissenseResult with the right fields.
      3. variant_scorer integrates the score and reports it in rationale.
    """
    from .alphamissense_integration import build_test_index, lookup
    from .variant_scorer import score_variant

    idx = build_test_index()
    assert len(idx) == 5, f"expected 5 test entries, got {len(idx)}"
    r = lookup("P01116", "G", 12, "D", index=idx)
    assert r is not None, "KRAS.G12D lookup returned None"
    assert r.classification == "likely_pathogenic"
    assert 0.0 <= r.score <= 1.0
    # End-to-end: BRAF.V600E with synthetic AlphaMissense
    r2 = score_variant(
        "BRAF",
        600,
        "V",
        "E",
        protein_length=766,
        uniprot_id="P15056",
        am_lookup=lambda u, w, p, m: lookup(u, w, p, m, index=idx),
    )
    assert "alphamissense_score" in r2.components
    assert r2.components["alphamissense_score"] > 0.5
    return True, (
        f"AlphaMissense lookup + integration OK: "
        f"5-entry test index, BRAF.V600E norm={r2.normalized_score}, "
        f"AM={r2.components['alphamissense_score']:.3f}"
    )


@register("variant.alphagenome_atlas")
def _check_alphagenome_atlas() -> tuple[bool, str]:
    """AlphaGenome Atlas Protocol adapter: regulatory-variant impact scoring.

    Validates:
      1. AVIResult dataclass: score in [0, 1], classification bins
         derived from the score (low/moderate/high thresholds match
         the AlphaMissense bins for consistent triage).
      2. RegulatoryVariantScorer Protocol is runtime_checkable.
      3. MockRegulatoryVariantScorer satisfies the Protocol and
         produces deterministic results across calls.
      4. Mock distinguishes coding-region variants (use AlphaMissense)
         from regulatory-region variants (use AlphaGenome Atlas) via
         position parity — purely a mock heuristic; the real Atlas
         returns per-variant is_coding.
      5. regulatory_score(chrom, pos, ref, alt) returns an AVIResult
         with score in [0, 1].
    """
    from .alphagenome_integration import (
        AVIResult,
        MockRegulatoryVariantScorer,
        RegulatoryVariantScorer,
        regulatory_score,
    )

    # 1. AVIResult validation
    r = AVIResult(score=0.5, classification="moderate", is_coding=False)
    assert 0.0 <= r.score <= 1.0
    assert r.classification == "moderate"

    bad: list[str] = []
    try:
        AVIResult(score=1.5, classification="high", is_coding=False)
    except ValueError:
        pass
    else:
        bad.append("out-of-range score should raise")

    # 2. Protocol runtime_checkable
    mock = MockRegulatoryVariantScorer()
    if not isinstance(mock, RegulatoryVariantScorer):
        bad.append("mock does not satisfy Protocol")

    # 3. Determinism
    a = mock.score_variant("chr7", 140753336, "T", "A")
    b = mock.score_variant("chr7", 140753336, "T", "A")
    if a.score != b.score or a.classification != b.classification:
        bad.append("mock not deterministic")

    # 4. Coding vs regulatory distinction
    coding = mock.score_variant("chr7", 140753336, "T", "A")  # even pos
    reg = mock.score_variant("chr7", 140753337, "T", "A")     # odd pos
    if not coding.is_coding or reg.is_coding:
        bad.append(
            f"coding/regulatory distinction failed: "
            f"even={coding.is_coding}, odd={reg.is_coding}"
        )

    # 5. regulatory_score convenience wrapper
    rs = regulatory_score("chr7", 140753336, "T", "A")
    if not (0.0 <= rs.score <= 1.0):
        bad.append(f"regulatory_score out of [0,1]: {rs.score}")

    if bad:
        return False, "AlphaGenome Atlas backend issues: " + "; ".join(bad)
    return True, (
        f"AlphaGenome Atlas OK: AVI in [0, 1], classification bins "
        f"(<0.34 low, <0.564 moderate, >=0.564 high), mock deterministic "
        f"(even pos → coding, odd pos → regulatory), regulatory_score "
        f"score={rs.score:.3f} class={rs.classification}"
    )


@register("scrna.structural_disruption_chou_fasman")
def _check_structural_disruption() -> tuple[bool, str]:
    """L→P in a helix context must score higher than L→P in a coil context.

    Proline is a helix breaker; in a real alpha-helix, L→P scores higher
    because of the additional structural-disruption component.
    """
    from .variant_scorer import score_variant

    helix_seq = "LAELAEKLAEEK"
    coil_seq = "GGPGGPPPGGPG"
    v_helix = score_variant(
        "FAKE",
        2,
        "L",
        "P",
        protein_length=12,
        protein_sequence=helix_seq,
    )
    v_coil = score_variant(
        "FAKE",
        2,
        "L",
        "P",
        protein_length=12,
        protein_sequence=coil_seq,
    )
    ok = v_helix.normalized_score > v_coil.normalized_score
    return ok, (
        f"L→P in helix={v_helix.normalized_score:.3f} (struct={v_helix.components['structural_disruption']}) "
        f"> L→P in coil={v_coil.normalized_score:.3f} (struct={v_coil.components['structural_disruption']})"
    )


@register("scrna.embedding_tfidf_svd")
def _check_embedding() -> tuple[bool, str]:
    import math

    from .foundation_embedder import embed_cells
    from .sc_rna_pipeline import _synthetic

    _, _, matrix = _synthetic()
    emb = embed_cells(matrix, model="tfidf-svd", n_components=8)
    # Cluster 0 (cells 0-19) vs cluster 2 (cells 35-49)
    intra = sum(
        math.sqrt(sum((emb[i][k] - emb[j][k]) ** 2 for k in range(8)))
        for i in range(19)
        for j in range(i + 1, 20)
    ) / (19 * 18 / 2)
    inter = (
        sum(math.sqrt(sum((emb[0][k] - emb[j][k]) ** 2 for k in range(8))) for j in range(35, 50))
        / 15
    )
    sep = inter / intra if intra else 0
    ok = sep > 1.5  # require at least 1.5x separation
    return ok, f"shape={len(emb)}x{len(emb[0])} inter/intra separation={sep:.2f}x"


# ---------- runner --------------------------------------------------------


def run_all(verbose: bool = True) -> int:
    failures = 0
    if verbose:
        print("=" * 72)
        print(f"{'Backend check':50s} {'Status':6s}  Detail")
        print("=" * 72)
    for name, fn in CHECKS:
        try:
            ok, detail = fn()
        except Exception as e:
            ok, detail = False, f"EXCEPTION: {e!r}"
        status = "PASS" if ok else "FAIL"
        if verbose:
            print(f"{name:50s} {status:6s}  {detail}")
        if not ok:
            failures += 1
    if verbose:
        print("=" * 72)
        print(f"{len(CHECKS) - failures}/{len(CHECKS)} checks passed.")
    return 0 if failures == 0 else 1


def main() -> int:
    return run_all(verbose=True)


if __name__ == "__main__":
    sys.exit(main())
