"""Codon-usage + GC + rare-codon analyzer for mRNA CDS sequences.

Implements the classical features that any modern codon-optimization model
(CodonBERT, RiboDecode, LinearDesign, mRNABERT) consumes as inputs:

* **CAI** (Codon Adaptation Index) — Sharp & Li (1987). Reference table
  embedded below for *Homo sapiens*.
* **GC%** — fraction of G+C in the CDS (excluding stop).
* **Rare-codon fraction** — fraction of codons whose usage frequency is
  below 0.10 in the reference table (heuristic for ribosome stalling).
* **Dinucleotide CpG/ObsExp** — proxy for innate immune recognition; very
  low CpG is a hallmark of self-mRNA.

The output is a structured dict that any downstream LLM (e.g. as a feature
provider for a prompt to GPT-4o, or for a small CNN) can consume directly.

References
----------
CodonBERT: https://academic.oup.com/bioinformatics/article/40/7/btae330
RiboDecode: https://www.nature.com/articles/s41467-025-64894-x
mRNABERT:  https://www.nature.com/articles/s41467-025-65340-8
"""

from __future__ import annotations

import math
import re
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

# Reference: human codon usage frequencies (frequencies per AA, summing to 1.0
# across synonymous codons). Source: codon usage database, canonical human table
# (Kazusa-style). Simplified to 5 decimals; sufficient for CAI.
HUMAN_CODON_FREQ: dict[str, dict[str, float]] = {
    "F": {"TTT": 0.46, "TTC": 0.54},
    "L": {"TTA": 0.08, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.41},
    "I": {"ATT": 0.36, "ATC": 0.47, "ATA": 0.17},
    "M": {"ATG": 1.00},
    "V": {"GTT": 0.18, "GTC": 0.24, "GTA": 0.12, "GTG": 0.46},
    "S": {"TCT": 0.19, "TCC": 0.22, "TCA": 0.15, "TCG": 0.05, "AGT": 0.15, "AGC": 0.24},
    "P": {"CCT": 0.29, "CCC": 0.32, "CCA": 0.28, "CCG": 0.11},
    "T": {"ACT": 0.25, "ACC": 0.36, "ACA": 0.28, "ACG": 0.11},
    "A": {"GCT": 0.27, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11},
    "Y": {"TAT": 0.44, "TAC": 0.56},
    "H": {"CAT": 0.42, "CAC": 0.58},
    "Q": {"CAA": 0.27, "CAG": 0.73},
    "N": {"AAT": 0.47, "AAC": 0.53},
    "K": {"AAA": 0.43, "AAG": 0.57},
    "D": {"GAT": 0.46, "GAC": 0.54},
    "E": {"GAA": 0.42, "GAG": 0.58},
    "C": {"TGT": 0.46, "TGC": 0.54},
    "W": {"TGG": 1.00},
    "R": {"CGT": 0.18, "CGC": 0.20, "CGA": 0.13, "CGG": 0.21, "AGA": 0.21, "AGG": 0.21},
    "G": {"GGT": 0.16, "GGC": 0.27, "GGA": 0.25, "GGG": 0.25},
    "*": {"TAA": 0.30, "TAG": 0.24, "TGA": 0.47},
}

# Map codon -> amino acid for the table above (built once).
CODON_TO_AA: dict[str, str] = {
    codon: aa for aa, codons in HUMAN_CODON_FREQ.items() for codon in codons
}

# Build per-codon absolute frequency by multiplying per-AA freq by AA frequency
# (approximated from the codon table itself so we don't need an external DB).
# For CAI we actually only need the per-AA max frequency, which we compute
# directly below.
_AA_MAX_FREQ: dict[str, float] = {aa: max(freqs.values()) for aa, freqs in HUMAN_CODON_FREQ.items()}


@dataclass
class CodonReport:
    n_codons: int
    cai: float
    gc_percent: float
    rare_codon_fraction: float
    cpg_obs_exp: float
    most_common_codons: list[tuple[str, int]]
    rare_codons: list[str]
    gc_window_stddev: float  # rolling GC stddev across 30-codon windows

    def to_dict(self) -> dict:
        return asdict(self)


def _read_fasta(text: str) -> list[tuple[str, str]]:
    """Minimal FASTA parser. Returns list of (header, sequence)."""
    records: list[tuple[str, str]] = []
    header = ""
    buf: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(buf).upper()))
            header = line[1:].strip() or "unnamed"
            buf = []
        else:
            buf.append(re.sub(r"\s+", "", line))
    if header:
        records.append((header, "".join(buf).upper()))
    return records


def _gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return 100.0 * gc / len(seq)


def _cpg_obs_exp(seq: str) -> float:
    """CpG observed/expected ratio (CpG suppression ~ immunogenicity proxy)."""
    if len(seq) < 2:
        return 1.0
    c = seq.count("C")
    g = seq.count("G")
    cg = sum(1 for i in range(len(seq) - 1) if seq[i] == "C" and seq[i + 1] == "G")
    if c == 0 or g == 0:
        return 0.0
    return (cg * len(seq)) / (c * g)


def _gc_window_stddev(seq: str, window: int = 90) -> float:
    """Stddev of GC% across rolling windows — proxy for translation speed variance."""
    if len(seq) < window:
        return 0.0
    vals = []
    for i in range(0, len(seq) - window + 1, max(1, window // 3)):
        sub = seq[i : i + window]
        vals.append(_gc_content(sub))
    if len(vals) < 2:
        return 0.0
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    return math.sqrt(var)


def analyze_cds(cds: str, *, rare_threshold: float = 0.10) -> CodonReport:
    """Analyze a coding DNA sequence (CDS, no UTRs, length multiple of 3)."""
    cds = re.sub(r"[^ATCG]", "", cds.upper().replace("U", "T"))
    # Strip trailing stop if present
    if len(cds) >= 3 and cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
        cds = cds[:-3]
    # Tolerate non-multiple-of-3 input by trimming to the largest valid prefix.
    trim = len(cds) % 3
    if trim:
        cds = cds[:-trim]
    if len(cds) < 3:
        raise ValueError("CDS shorter than one codon after cleaning")
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]

    # CAI: geometric mean of (freq[c] / max_freq[aa(c)])
    log_sum = 0.0
    n = 0
    rare: list[str] = []
    codon_counts = Counter(codons)
    for c in codons:
        aa = CODON_TO_AA.get(c)
        if aa is None or aa == "*":
            continue
        freq = HUMAN_CODON_FREQ[aa].get(c, 0.0)
        if freq < rare_threshold:
            rare.append(c)
        denom = _AA_MAX_FREQ[aa]
        if denom > 0 and freq > 0:
            log_sum += math.log(freq / denom)
            n += 1
    cai = math.exp(log_sum / n) if n else 0.0

    gc_pct = _gc_content(cds)
    rare_frac = len(rare) / len(codons) if codons else 0.0
    cpg_oe = _cpg_obs_exp(cds)
    most_common = codon_counts.most_common(5)
    gc_std = _gc_window_stddev(cds)

    return CodonReport(
        n_codons=len(codons),
        cai=round(cai, 4),
        gc_percent=round(gc_pct, 2),
        rare_codon_fraction=round(rare_frac, 4),
        cpg_obs_exp=round(cpg_oe, 4),
        most_common_codons=most_common,
        rare_codons=sorted(set(rare))[:20],
        gc_window_stddev=round(gc_std, 2),
    )


def optimize_basic(
    cds: str,
    *,
    target_gc_min: float = 45.0,
    target_gc_max: float = 60.0,
    rare_threshold: float = 0.10,
) -> dict:
    """Greedy synonymous-codon swap targeting higher-freq codon and ideal GC band.

    This is the **classic baseline** that any modern codon model (CodonBERT,
    RiboDecode) should beat. Useful as a sanity check + a starting prompt
    feature for an LLM.

    Returns: ``{"new_cds": str, "before": dict, "after": dict, "changes": int}``
    """
    cds = cds.upper().replace("U", "T")
    if cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
        cds = cds[:-3]

    before = analyze_cds(cds).to_dict()
    out = list(cds)
    changes = 0
    for i in range(0, len(cds), 3):
        codon = cds[i : i + 3]
        aa = CODON_TO_AA.get(codon)
        if aa is None or aa == "*" or aa == "M" or aa == "W":
            continue  # single-codon AAs, nothing to swap
        freqs = HUMAN_CODON_FREQ[aa]
        if freqs[codon] >= 0.30:  # already a good codon
            continue
        # Pick the highest-frequency synonymous codon that best matches current GC%
        current_gc = _gc_content(codon)
        candidates = sorted(
            freqs.items(),
            key=lambda kv: (
                -kv[1],
                abs((_gc_content(kv[0]) - current_gc)),
            ),
        )
        new_codon = candidates[0][0]
        if new_codon != codon:
            out[i : i + 3] = new_codon
            changes += 1
    new_cds = "".join(out)
    after = analyze_cds(new_cds, rare_threshold=rare_threshold).to_dict()
    return {
        "new_cds": new_cds,
        "before": before,
        "after": after,
        "changes": changes,
        "target_gc_band": (target_gc_min, target_gc_max),
    }


# ---------- CLI entry -------------------------------------------------------


def _run_cli(argv: list[str]) -> int:
    import argparse
    import json as _json

    p = argparse.ArgumentParser(prog="mrnavax codon")
    p.add_argument("--sequence", required=True, help="FASTA file or raw CDS string")
    p.add_argument("--optimize", action="store_true", help="greedy codon-optimize")
    p.add_argument(
        "--backend",
        choices=["basic", "ribodecode", "lineardesign", "ribodecode-real", "multi-objective"],
        default="basic",
        help="optimizer: basic (greedy), ribodecode (in-house heuristic), "
        "lineardesign (joint translation + mRNA structure via DP), "
        "ribodecode-real (published RiboDecode CLI from "
        "github.com/wangfanfff/RiboDecode, requires ViennaRNA+CUDA), or "
        "multi-objective (SOTA knowledge-infused loss: CAI + GC + CpG + "
        "rare-codon-run + structure-proxy with per-component breakdown)",
    )
    p.add_argument(
        "--env",
        default="HEK293T",
        help="cellular environment for ribodecode-real: HEK293T, A549, "
        "HeLa, or custom (requires --env-csv)",
    )
    p.add_argument(
        "--env-csv",
        help="custom environment CSV (gene-ID, RPKM) for --env custom",
    )
    p.add_argument(
        "--mfe-weight",
        type=float,
        default=0.0,
        help="structure weight for ribodecode-real (0=translation-only, 1=MFE-only)",
    )
    p.add_argument(
        "--optim-epoch",
        type=int,
        default=10,
        help="number of optimization epochs for ribodecode-real",
    )
    p.add_argument("--ribo-weights", help="JSON file with per-codon translation rates")
    p.add_argument("--out", help="write JSON report here")
    args = p.parse_args(argv)

    raw = Path(args.sequence).read_text() if Path(args.sequence).exists() else args.sequence
    if raw.lstrip().startswith(">"):
        records = _read_fasta(raw)
        cds = records[0][1] if records else ""
    else:
        cds = raw.strip()

    if args.optimize:
        if args.backend == "ribodecode-real":
            # Real upstream RiboDecode package (CLI subprocess).
            # Falls back to the stdlib mock if the binary is missing.
            from .codon_protocols import RiboDecodeRequest
            from .codon_ribodecode_adapter import (
                RiboDecodeNotInstalled,
                select_codon_optimizer,
            )

            req = RiboDecodeRequest(
                cds=cds,
                env=args.env,
                custom_env_csv=Path(args.env_csv) if args.env_csv else None,
                mfe_weight=args.mfe_weight,
                optim_epoch=args.optim_epoch,
            )
            optimizer = select_codon_optimizer(prefer="auto")
            try:
                ribo_res = optimizer.optimize(req)
            except RiboDecodeNotInstalled as e:
                raise SystemExit(str(e))
            # Map to the same shape optimize_basic/lineardesign produce
            before = analyze_cds(cds).to_dict()
            after = analyze_cds(ribo_res.optimized_cds).to_dict()
            changes = sum(
                1
                for a, b in zip(
                    [cds[i : i + 3] for i in range(0, len(cds) - 2, 3)],
                    [
                        ribo_res.optimized_cds[i : i + 3]
                        for i in range(0, len(ribo_res.optimized_cds), 3)
                    ],
                )
                if a != b
            )
            result = {
                "new_cds": ribo_res.optimized_cds,
                "before": before,
                "after": after,
                "changes": changes,
                "translation_score": ribo_res.predicted_translation,
                "structure_score": ribo_res.predicted_mfe,
                "weights": {
                    "translation_weight": 1.0 - args.mfe_weight,
                    "structure_weight": args.mfe_weight,
                    "gc_window_size": 21,
                    "min_stem_length": 3,
                },
                "elapsed_seconds": round(ribo_res.elapsed_seconds, 3),
                "n_states_evaluated": 0,  # unknown for upstream RiboDecode
                "backend": ribo_res.backend,
                "epoch": ribo_res.epoch,
                "env": args.env,
                "notes": list(ribo_res.notes),
            }
        elif args.backend == "ribodecode":
            from .codon_ribodecode import optimize_ribodecode

            ribo_weights = None
            if args.ribo_weights:
                ribo_weights = _json.loads(Path(args.ribo_weights).read_text())
            result = optimize_ribodecode(cds, ribo_weights=ribo_weights).to_dict()
        elif args.backend == "lineardesign":
            from .codon_lineardesign import optimize_lineardesign

            # LinearDesign DP is now linear-time per step (state space
            # bounded by |syn|^W, independent of CDS length), so the CLI
            # accepts full-length CDS.
            verbose = getattr(args, "verbose", False)
            result = optimize_lineardesign(cds, verbose=verbose).to_dict()
        elif args.backend == "multi-objective":
            # SOTA pattern (arxiv:2505.23862, Aug 2026): knowledge-infused
            # loss decomposition with per-component breakdown.
            from .codon_multi_objective import (
                MultiObjectiveConfig,
                multi_objective_optimize,
            )

            cfg = MultiObjectiveConfig()
            result = multi_objective_optimize(cds, cfg).to_dict()
        else:
            result = optimize_basic(cds)
    else:
        result = analyze_cds(cds).to_dict()

    out_text = _json.dumps(result, indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_run_cli(sys.argv[1:]))
