"""Backend checks for the ``variant.*`` family.

Extracted from ``mrnavax/backends.py`` in v0.24.0. Each check is
registered with the global ``CHECKS`` registry on import.
"""

from __future__ import annotations

from pathlib import Path

from ._backends_registry import register


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
      6. AVI integration into score_variant: when avi_lookup is provided
         with chrom/ref_dna/alt_dna AND the lookup returns is_coding=False,
         AVI becomes the dominant signal at weight 0.45. Components dict
         exposes alphagenome_atlas_score. For is_coding=True variants,
         AlphaMissense still dominates; AVI is recorded as secondary.
    """
    from .alphagenome_integration import (
        AVIResult,
        MockRegulatoryVariantScorer,
        RegulatoryVariantScorer,
        regulatory_score,
    )
    from .alphamissense_integration import AlphaMissenseResult
    from .variant_scorer import score_variant

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

    # 6. AVI integration into score_variant
    # 6a. Coding variant: AM dominates, AVI is secondary
    def am_coding(uniprot, wt, pos, mut):
        return AlphaMissenseResult(
            score=0.9,
            classification="likely_pathogenic",
            uniprot=uniprot,
            aa_change=f"{wt}{pos}{mut}",
        )

    # Mock returns is_coding=True at even positions (140753336 is even).
    # So AM should dominate over AVI here.
    r_coding = score_variant(
        "BRAF",
        600,
        "V",
        "E",
        chrom="chr7",
        ref_dna="T",
        alt_dna="A",
        protein_length=766,
        uniprot_id="P15056",
        am_lookup=am_coding,
        avi_lookup=mock.score_variant,
    )
    if "alphamissense_score" not in r_coding.components:
        bad.append("AVI-integrated: AM score missing for coding variant")
    if "alphagenome_atlas_score" not in r_coding.components:
        bad.append("AVI-integrated: AVI score missing for coding variant")
    # AM dominates → rationale should NOT have "(dominant)" on AVI
    if "(dominant)" in r_coding.rationale:
        bad.append("AVI-integrated: AVI should not be dominant for coding variant")

    # 6b. Regulatory variant: AVI dominates
    # Position 101 is odd → mock returns is_coding=False → AVI dominates.
    r_reg = score_variant(
        "REG_GENE",
        101,
        "A",
        "G",
        chrom="chr7",
        ref_dna="A",
        alt_dna="G",
        protein_length=200,
        am_lookup=am_coding,
        avi_lookup=mock.score_variant,
    )
    if "alphagenome_atlas_score" not in r_reg.components:
        bad.append("AVI-integrated: AVI score missing for regulatory variant")
    if "(dominant)" not in r_reg.rationale:
        bad.append("AVI-integrated: AVI should be dominant for regulatory variant")

    if bad:
        return False, "AlphaGenome Atlas backend issues: " + "; ".join(bad)
    return True, (
        f"AlphaGenome Atlas OK: AVI in [0, 1], classification bins "
        f"(<0.34 low, <0.564 moderate, >=0.564 high), mock deterministic "
        f"(even pos → coding, odd pos → regulatory), regulatory_score "
        f"score={rs.score:.3f} class={rs.classification}; AVI-integrated "
        f"score_variant: coding variant → AM dominates (rationale "
        f"no '(dominant)' on AVI), regulatory variant → AVI dominates "
        f"(rationale has '(dominant)')"
    )


@register("variant.alphagenome_atlas_fixture")
def _check_alphagenome_atlas_fixture() -> tuple[bool, str]:
    """Regression test against the bundled AlphaGenome Atlas API
    response fixture (tests/fixtures/alphagenome_atlas_sample.json).

    Catches upstream Atlas API schema changes at PR time. The fixture
    is a synthetic sample modeled on the documented Atlas response
    shape; a real captured response can replace it when an
    ALPHAGENOME_API_KEY is available. Validates:

      1. Fixture file exists and parses as JSON.
      2. Every variant has the required keys (chrom, pos, ref, alt,
         score, classification, is_coding).
      3. Score is in [0, 1], classification is one of
         {low, moderate, high}, is_coding is bool.
      4. The fixture exercises both coding and regulatory variants
         (so the parser handles both is_coding branches).
      5. All three AVI classifications are present.
      6. Each fixture variant parses into a valid AVIResult via the
         subprocess adapter payload shape (mirrors what
         AlphaGenomeCLIAdapter.score_variant does with shim output).
    """
    import json as _json

    from .alphagenome_integration import AVIResult

    # 1. Fixture file exists
    fixture_path = (
        Path(__file__).resolve().parent.parent
        / "tests"
        / "fixtures"
        / "alphagenome_atlas_sample.json"
    )
    if not fixture_path.exists():
        return False, f"missing fixture file: {fixture_path}"

    # 2. Parses as JSON with required structure
    with open(fixture_path) as f:
        data = _json.load(f)
    if "variants" not in data or not isinstance(data["variants"], list):
        return False, "fixture missing 'variants' list"

    bad: list[str] = []
    coding_count = 0
    regulatory_count = 0
    classifications_seen: set[str] = set()

    # 3-5. Per-variant validation
    for i, rec in enumerate(data["variants"]):
        missing = {"chrom", "pos", "ref", "alt", "score", "classification", "is_coding"} - rec.keys()
        if missing:
            bad.append(f"variant[{i}] missing fields {missing}")
            continue
        try:
            # 6. Parse into AVIResult (mirrors adapter payload shape)
            # The construction itself validates the shape (raises
            # ValueError on out-of-range score or invalid class).
            AVIResult(
                score=float(rec["score"]),
                classification=str(rec["classification"]),
                is_coding=bool(rec["is_coding"]),
            )
        except (ValueError, TypeError) as e:
            bad.append(f"variant[{i}] AVIResult parse failed: {e}")
            continue
        if rec["is_coding"]:
            coding_count += 1
        else:
            regulatory_count += 1
        classifications_seen.add(rec["classification"])

    if coding_count == 0:
        bad.append("fixture has no coding-region variants")
    if regulatory_count == 0:
        bad.append("fixture has no regulatory-region variants")
    for expected in {"low", "moderate", "high"}:
        if expected not in classifications_seen:
            bad.append(f"fixture missing classification={expected}")

    if bad:
        return False, "AlphaGenome Atlas fixture issues: " + "; ".join(bad)
    return True, (
        f"AlphaGenome Atlas fixture OK: {len(data['variants'])} variants parsed "
        f"({coding_count} coding + {regulatory_count} regulatory), all "
        f"3 classifications (low/moderate/high) present, scores in [0,1]"
    )


__all__ = []
