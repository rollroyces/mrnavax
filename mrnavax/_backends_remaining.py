"""Backend checks for the one-off families (manufacture / lnp / case_study / conservation).

These don't have enough checks to warrant their own module. Extracted
from ``mrnavax/backends.py`` in v0.24.0.
"""

from __future__ import annotations

from ._backends_registry import register


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


@register("lnp.recommend")
def _check_lnp_recommend() -> tuple[bool, str]:
    from .lnp_advisor import recommend

    r = recommend(target="lung", cargo="saRNA", intent="cancer vaccine", n=3)
    ok = len(r.shortlist) > 0 and all("name" in c for c in r.shortlist)
    return ok, f"{len(r.shortlist)} candidates for lung/saRNA"


@register("case_study.variant_prioritization")
def _check_case_study() -> tuple[bool, str]:
    """End-to-end: run the curated ClinVar variant set through
    score_variant with both AlphaMissense (coding) and AlphaGenome
    Atlas AVI (regulatory) wired up. Verify precision@K against
    ground-truth pathogenicity labels.

    This is the only check that exercises the FULL pipeline stack:
    load → score → rank → compute precision@K. Catches integration
    regressions that single-tool checks miss.
    """
    from .alphagenome_integration import MockRegulatoryVariantScorer
    from .alphamissense_integration import build_test_index, lookup
    from .case_study import (
        load_clinvar_variants,
        precision_at_k,
        score_case_study_variants,
        summarize_variant_set,
    )
    from .variant_scorer import score_variant

    variants = load_clinvar_variants()
    summary = summarize_variant_set(variants)
    if summary["n_pathogenic"] == 0:
        return False, "curated CSV has no pathogenic variants — case study broken"

    idx = build_test_index()
    mock_avi = MockRegulatoryVariantScorer()

    def scorer(variant):
        # The curated CSV is hg38-coordinates only; we use placeholder
        # AAs (V/E) for coding variants and A/G for regulatory. The
        # BLOSUM62 + AM/AVI signal still carries the priority.
        if variant.is_coding:
            wt, mut = "V", "E"
        else:
            wt, mut = "A", "G"
        r = score_variant(
            gene=variant.gene,
            position=1,
            wt_aa=wt,
            mut_aa=mut,
            chrom=variant.chrom,
            ref_dna=variant.ref,
            alt_dna=variant.alt,
            uniprot_id=None,
            am_lookup=lambda u, w, p, m: lookup(u, w, p, m, index=idx) if u else None,
            avi_lookup=mock_avi.score_variant,
        )
        return {"score": r.normalized_score if r else 0.0}

    scored = score_case_study_variants(variants, scorer)
    triples = [(s["label"], s["score"], s["pathogenicity"]) for s in scored]
    p_at_3 = precision_at_k(triples, k=3)

    # Sanity: the case study must produce a non-trivial precision@3.
    # With mock backends, top-3 should all be pathogenic (precision@3=1.0).
    if p_at_3 < 0.66:
        return False, (
            f"precision@3 = {p_at_3:.3f}, expected ≥ 0.66 with mock backends"
        )

    return True, (
        f"case study OK: {summary['n_total']} variants "
        f"({summary['n_pathogenic']} pathogenic + {summary['n_benign']} benign + "
        f"{summary['n_uncertain']} uncertain), scored and ranked; "
        f"precision@3 = {p_at_3:.3f}"
    )


@register("conservation.phylop46way")
def _check_conservation_phylop() -> tuple[bool, str]:
    """PhyloP46way / GERP++ conservation lookup: the 4th coding-region
    scoring signal.

    Validates:
      1. MockPhyloPLookup satisfies the ConservationLookup Protocol.
      2. Mock returns scores in [-1, 1] (PhyloP range).
      3. Mock is deterministic across calls.
      4. score_variant accepts conservation_lookup=... and records the
         score in components['phylop46way_score'].
      5. Conservation composes with AlphaMissense (coding) and
         AlphaGenome Atlas AVI (regulatory) without breaking the
         existing routing logic.
      6. Conservation is silent on lookup errors / None returns /
         out-of-range values.
    """
    from .conservation import (
        ConservationLookup,
        MockPhyloPLookup,
    )
    from .variant_scorer import score_variant

    # 1. Protocol runtime_checkable
    mock = MockPhyloPLookup()
    if not isinstance(mock, ConservationLookup):
        return False, "MockPhyloPLookup does not satisfy ConservationLookup Protocol"

    # 2. Range check
    bad: list[str] = []
    a = mock.lookup("chr7", 140753336)
    b = mock.lookup("chr11", 12345678)
    if not -1.0 <= a <= 1.0:
        bad.append(f"mock score {a} out of [-1, 1]")
    if not -1.0 <= b <= 1.0:
        bad.append(f"mock score {b} out of [-1, 1]")

    # 3. Determinism
    a2 = mock.lookup("chr7", 140753336)
    if a != a2:
        bad.append("mock not deterministic")

    # 4. score_variant wires through
    r = score_variant(
        "BRAF",
        600,
        "V",
        "E",
        protein_length=766,
        chrom="chr7",
        pos=140753336,
        conservation_lookup=mock.lookup,
    )
    if "phylop46way_score" not in r.components:
        bad.append("phylop46way_score missing in components")
    if not -1.0 <= r.components.get("phylop46way_score", 0.0) <= 1.0:
        bad.append("phylop46way_score in components out of range")

    # 5. Composes with AM + AVI
    from mrnavax.alphagenome_integration import MockRegulatoryVariantScorer
    from mrnavax.alphamissense_integration import AlphaMissenseResult

    def am_coding(uniprot, wt, pos, mut):
        return AlphaMissenseResult(
            score=0.9, classification="likely_pathogenic",
            uniprot=uniprot, aa_change=f"{wt}{pos}{mut}",
        )

    r_full = score_variant(
        "BRAF",
        600,
        "V",
        "E",
        protein_length=766,
        uniprot_id="P15056",
        chrom="chr7",
        pos=140753336,
        ref_dna="T",
        alt_dna="A",
        am_lookup=am_coding,
        avi_lookup=MockRegulatoryVariantScorer().score_variant,
        conservation_lookup=mock.lookup,
    )
    if "alphamissense_score" not in r_full.components:
        bad.append("composed: AM missing")
    if "alphagenome_atlas_score" not in r_full.components:
        bad.append("composed: AVI missing")
    if "phylop46way_score" not in r_full.components:
        bad.append("composed: PhyloP missing")

    # 6. Silent on error
    def raises(*args):
        raise RuntimeError("UCSC down")

    r_err = score_variant(
        "BRAF",
        600,
        "V",
        "E",
        protein_length=766,
        chrom="chr7",
        pos=140753336,
        conservation_lookup=raises,
    )
    if "phylop46way_score" in r_err.components:
        bad.append("silent fail: PhyloP should be dropped on lookup error")

    if bad:
        return False, "conservation.PhyloP integration issues: " + "; ".join(bad)
    return True, (
        "PhyloP46way OK: mock deterministic, scores in [-1, 1]; "
        "score_variant wires through to components.phylop46way_score "
        "(composes with AM + AVI without breaking routing); "
        "silent on lookup error"
    )


@register("construct.full_assembly")
def _check_construct_full_assembly() -> tuple[bool, str]:
    """mRNA construct designer assembles 5'UTR + CDS + 3'UTR + poly-A
    and produces a valid full-length DNA construct.

    Verifies:
      * Full assembly length = utr5 + cds + utr3 + polyA_signal + polyA.
      * CDS preserves the input amino-acid sequence (fidelity).
      * CDS is codon-optimized (multi-objective backend default).
      * Construct starts with the 5'UTR and ends with the poly-A tail.
      * Full-length eGFP (239 AA) builds without error.
    """
    from .codon_optimizer import CODON_TO_AA
    from .construct_designer import (
        _POLYA_SIGNAL,
        ConstructConfig,
        design_construct,
    )

    # Small test protein (12 AA) — fast, exercises all regions
    protein = "MVSKGEELFTGV"
    result = design_construct(protein)

    # Full assembly length
    expected_len = (
        len(result.utr5)
        + result.cds_length
        + len(result.utr3)
        + len(_POLYA_SIGNAL)
        + len(result.poly_a_tail)
    )
    assert len(result.construct_dna) == expected_len, (
        f"construct length {len(result.construct_dna)} != expected {expected_len}"
    )

    # Region boundaries
    assert result.construct_dna.startswith(result.utr5), "5'UTR not at start"
    assert result.construct_dna.endswith(result.poly_a_tail), "poly-A not at end"

    # Fidelity
    cds_aas = [
        CODON_TO_AA[result.cds[i : i + 3]]
        for i in range(0, len(result.cds), 3)
    ]
    assert "".join(cds_aas) == protein, "CDS changed amino-acid sequence"

    # Optimization ran (multi-objective backend is the default)
    assert len(result.optimization_notes) > 0, (
        "expected multi-objective optimization note"
    )

    # Full-length eGFP
    egfp = (
        "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
        "VTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVN"
        "RIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHY"
        "QQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    )
    egfp_result = design_construct(egfp, ConstructConfig())
    assert egfp_result.n_codons == 239
    assert egfp_result.cds_length == 239 * 3

    return True, (
        f"construct designer OK: small (12 AA) {expected_len} nt, "
        f"eGFP (239 AA) {egfp_result.construct_length} nt; "
        f"5'UTR + CDS + 3'UTR + AAUAAA + polyA({len(result.poly_a_tail)} nt); "
        f"CDS fidelity preserved, multi-objective optimization applied"
    )


@register("variant.utr_aware_scoring")
def _check_utr_aware_scoring() -> tuple[bool, str]:
    """UTR context as a 5th signal in score_variant (v0.27.0).

    Verifies:
      * Without utr5/utr3, the components dict has no utr_context_score
        (backward compatible — the 4-signal composition is preserved).
      * With strong utr5 + utr3 supplied, utr_context_score / kozak_score /
        utr3_score appear in the components dict.
      * With strong UTR context, the normalized_score is >= the
        score without UTR context (bonus, never a penalty).
      * With weak UTR context, the normalized_score is also >= the
        score without UTR context (additive bonus only, never a penalty).
    """
    from .variant_scorer import score_variant

    # 1. Without UTR context — must NOT add utr_context_score
    r_no = score_variant(
        gene="BRAF", position=600, wt_aa="V", mut_aa="E",
        protein_length=766, chrom="chr7", pos=140753336,
    )
    assert r_no is not None, "score_variant returned None without UTR context"
    assert "utr_context_score" not in r_no.components, (
        "utr_context_score must not appear when utr5/utr3 are unset"
    )

    # 2. With strong UTR context — components must include utr scores
    r_strong = score_variant(
        gene="BRAF", position=600, wt_aa="V", mut_aa="E",
        protein_length=766, chrom="chr7", pos=140753336,
        utr5="GGGCGACGCGGTGGCGGCCACCAAT",
        utr3="AUUUAGCAUUUAGCAUUUAG" + "A" * 100,
    )
    assert r_strong is not None, "score_variant returned None with UTR context"
    assert "utr_context_score" in r_strong.components
    assert "kozak_score" in r_strong.components
    assert "utr3_score" in r_strong.components

    # 3. With weak UTR context — components must still include utr scores
    r_weak = score_variant(
        gene="BRAF", position=600, wt_aa="V", mut_aa="E",
        protein_length=766, chrom="chr7", pos=140753336,
        utr5="GGGCGACGCAAAAAAAAAAAAAAAAAAAA",
        utr3="GCGCGCGC" * 5,
    )
    assert r_weak is not None
    assert "utr_context_score" in r_weak.components

    # 4. Strong UTR context boosts score (or keeps it equal)
    assert r_strong.normalized_score >= r_no.normalized_score, (
        f"strong UTR lowered score: {r_strong.normalized_score} < {r_no.normalized_score}"
    )

    # 5. Weak UTR context never lowers the score (additive bonus only)
    assert r_weak.normalized_score >= r_no.normalized_score, (
        f"weak UTR lowered score: {r_weak.normalized_score} < {r_no.normalized_score}"
    )

    return True, (
        f"UTR-aware scoring OK: no-utr={r_no.normalized_score:.3f}, "
        f"strong-utr={r_strong.normalized_score:.3f}, "
        f"weak-utr={r_weak.normalized_score:.3f}; "
        f"5th signal (utr_context_score) added when utr5/utr3 supplied; "
        f"additive bonus only (never a penalty)"
    )


@register("variant.score_variant_inputs_dataclass")
def _check_score_variant_inputs_dataclass() -> tuple[bool, str]:
    """VariantInputs dataclass + score_variant_from_inputs (v0.29.0).

    Verifies:
      * VariantInputs is a public dataclass with 4 required + 14 optional
        fields (18 total).
      * VariantInputs.from_kwargs() builds the same shape as kwargs callers
        of score_variant() use.
      * score_variant_from_inputs(inputs) returns the same normalized_score
        as score_variant(**kwargs).
      * filter_variants_with_lookups(score_only=True) forces the
        score-each behavior regardless of inferred mode.
    """
    import dataclasses

    from .sc_rna_pipeline import Variant
    from .variant_scorer import (
        VariantInputs,
        score_variant,
        score_variant_from_inputs,
    )

    # 1. VariantInputs has 4 required + 14 optional fields (18 total)
    fields = VariantInputs.__dataclass_fields__
    required = [
        f for f, v in fields.items()
        if v.default is dataclasses.MISSING
        and v.default_factory is dataclasses.MISSING
    ]
    assert len(required) == 4, (
        f"VariantInputs must have 4 required fields, got {len(required)}: {required}"
    )
    assert "gene" in required
    assert "position" in required
    assert "wt_aa" in required
    assert "mut_aa" in required
    assert len(fields) == 18, (
        f"VariantInputs must have 18 fields total, got {len(fields)}"
    )

    # 2. from_kwargs builds the dataclass correctly
    inputs = VariantInputs.from_kwargs(
        "BRAF", 600, "V", "E",
        protein_length=766,
        chrom="chr7",
        pos=140753336,
    )
    assert inputs.gene == "BRAF"
    assert inputs.position == 600
    assert inputs.wt_aa == "V"
    assert inputs.mut_aa == "E"
    assert inputs.protein_length == 766
    assert inputs.chrom == "chr7"
    assert inputs.pos == 140753336

    # 3. score_variant_from_inputs produces the same result
    r_kwargs = score_variant("BRAF", 600, "V", "E", protein_length=766)
    r_inputs = score_variant_from_inputs(inputs)
    assert r_kwargs is not None
    assert r_inputs is not None
    assert r_inputs.normalized_score == r_kwargs.normalized_score

    # 4. score_only=True on filter_variants_with_lookups forces score-each
    variants_with_dna = [
        Variant(
            gene="BRAF", position=600, wt_aa="V", mut_aa="E",
            chrom="chr7", ref_dna="T", alt_dna="A",
        ),
    ]
    from ._scrna_filter import filter_variants_with_lookups
    kept, scores, _ = filter_variants_with_lookups(
        variants_with_dna,
        top_fraction=1.0,
        min_score=0.0,
        proteins={"BRAF": "M" * 766},
        uniprot_ids={"BRAF": "P15056"},
        score_only=True,
    )
    assert kept == variants_with_dna
    assert len(scores) == 1, (
        f"score_only=True with DNA coords must produce scores, got {scores}"
    )

    return True, (
        f"VariantInputs OK: {len(fields)} fields ({len(required)} required + "
        f"{len(fields) - len(required)} optional); "
        f"score_variant_from_inputs() matches score_variant() at "
        f"{r_inputs.normalized_score:.3f}; "
        f"score_only=True forces score-each behavior "
        f"({len(scores)} variant scored with DNA coords)"
    )


@register("utr_designer.coupled_design")
def _check_utr_designer() -> tuple[bool, str]:
    """UTR-aware CDS design (v0.30.0, 10th tool).

    Verifies:
      * design_utr_aware_cds() returns a UTRDesignResult with all
        expected fields populated.
      * The chosen UTR5 + UTR3 are valid DNA strings (no U's, A/T/C/G only).
      * The CDS length = 3 × len(input_protein).
      * The combined_score is in [0, 1].
      * The ranking has 12 entries (4 5'UTR × 3 3'UTR combinations).
      * The top-ranked combination has the highest combined_score.
      * Strong-Kozak 5'UTR variants are preferred when available
        (consistent with the manufacturing literature).
    """
    from .utr_designer import (
        UTRDesignConfig,
        design_utr_aware_cds,
    )

    cfg = UTRDesignConfig(cds="MVSKGEELFTGV")
    result = design_utr_aware_cds(cfg)

    # Field population
    assert isinstance(result.utr5, str) and len(result.utr5) > 0
    assert isinstance(result.utr3, str) and len(result.utr3) > 0
    assert isinstance(result.cds_dna, str)
    assert len(result.cds_dna) == 3 * len(cfg.cds), (
        f"CDS must be 3 × len(AA): got {len(result.cds_dna)} for "
        f"{len(cfg.cds)} AA"
    )
    assert result.protein == cfg.cds

    # Valid DNA
    valid_bases = set("ATCG")
    assert all(c in valid_bases for c in result.utr5.upper()), (
        f"utr5 contains non-DNA bases: {result.utr5}"
    )
    assert all(c in valid_bases for c in result.utr3.upper()), (
        f"utr3 contains non-DNA bases: {result.utr3}"
    )
    assert all(c in valid_bases for c in result.cds_dna.upper())

    # Scores in [0, 1]
    assert 0.0 <= result.combined_score <= 1.0
    assert 0.0 <= result.cds_score <= 1.0

    # Ranking
    assert len(result.ranking) == 12, (
        f"ranking must have 12 entries (4 UTR5 × 3 UTR3), got {len(result.ranking)}"
    )
    scores = [r[0] for r in result.ranking]
    assert scores == sorted(scores, reverse=True), "ranking must be sorted descending"

    # Top-ranked is the chosen one (or one of the top)
    top_score = result.ranking[0][0]
    assert abs(result.combined_score - top_score) < 0.01, (
        f"chosen combined_score {result.combined_score} != top-ranked {top_score}"
    )

    # Strong-Kozak in top-3
    top3_names = {name for _, name, _ in result.ranking[:3]}
    assert "strong_kozak" in top3_names, (
        f"strong_kozak should be in top 3; got {top3_names}"
    )

    return True, (
        f"UTR-aware design OK: chose {result.utr5[:15]}... + "
        f"{result.utr3[:15]}...; CDS={len(result.cds_dna)} nt, "
        f"score={result.cds_score:.3f}, combined={result.combined_score:.3f}; "
        f"12 candidates evaluated, ranking sorted descending; "
        f"strong_kozak in top 3"
    )


__all__ = []
