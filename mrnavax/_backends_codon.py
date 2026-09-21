"""Backend checks for the ``codon.*`` family.

Extracted from ``mrnavax/backends.py`` in v0.24.0. Each check is
registered with the global ``CHECKS`` registry on import. To run:

    python -m mrnavax.backends --check-all
"""

from __future__ import annotations

from ._backends_registry import register

_CODON_TABLE: dict[str, str] = {
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
    """Translate a CDS to its amino-acid sequence (with '*' for stops)."""
    return "".join(_CODON_TABLE.get(cds[i : i + 3], "?") for i in range(0, len(cds), 3))


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


__all__ = []
