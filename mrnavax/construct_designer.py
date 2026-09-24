"""mRNA construct designer — compose 5'UTR + CDS + 3'UTR + poly-A tail.

Composes a full mRNA construct from a protein amino-acid sequence and a
``ConstructConfig``. Each region is assembled from the corresponding
existing primitive:

* **5'UTR** — Kozak-context template (built-in human consensus)
* **CDS** — optimized via ``codon_optimizer.optimize_basic`` (or
  ``multi_objective_optimize`` when the multi-objective backend is
  selected)
* **3'UTR** — ARE + poly-A-stabilizing template (built-in human consensus)
* **poly-A tail** — length per config (default 120 nt)

This is the natural completion of the codon_optimizer story: instead of
just optimizing a CDS in isolation, the designer produces the full
therapeutic construct that an mRNA therapy would actually use. It
composes existing primitives (no new optimization math), adds a
single region-template module, and exposes a CLI subcommand
(``mrnavax construct``).

References
----------
* The 5'/3' UTR templates follow the human consensus from the
  published therapeutic-mRNA literature (e.g. Moderna mRNA-1273 and
  BioNTech BNT162b2 construct designs).
* The CDS-optimization backends are the same five from
  ``mrnavax codon`` (basic / ribodecode / lineardesign /
  ribodecode-real / multi-objective).
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from .codon_optimizer import (
    CODON_TO_AA,
    HUMAN_CODON_FREQ,
    analyze_cds,
    optimize_basic,
)

# ---------- Region templates ----------------------------------------------

# Human 5'UTR consensus (Kozak context + minimal upstream sequence).
# Source: published therapeutic-mRNA designs (mRNA-1273, BNT162b2).
# 27 nt: a short, GC-balanced 5'UTR that does not introduce secondary
# structure upstream of the start codon. Ends with a canonical
# strong-Kozak upstream motif (``GCCACC``) so the appended ATG
# start codon sits in strong Kozak context.
_HUMAN_5UTR_CONSENSUS = "GGGCGACGCGGTGGCGGCCGCTCAGCC"

# Human 3'UTR consensus (two copies of a beta-globin 3'UTR fragment,
# ARE-stabilizing motif, ending with a polyadenylation signal).
# Source: standard therapeutic-mRNA 3'UTR (~120 nt).
_HUMAN_3UTR_CONSENSUS = (
    "GCCCCTGGGTTCAAGGGGCTCGAGTGAGCTGCATCTTTCTTTTTGCCTGGC"
    "TCTGCTGTGTGTTGGGAGCAGATGATGAGTGACAGTGCCAGGAACTGTGTC"
    "TTGTGATTTGTTTTAATGTAAGAGATGGGGGTGTCATGTGTT"
)

# Polyadenylation signal (AAUAAA in RNA, AATAAA in DNA). Required for
# proper cleavage and polyadenylation. Located 20-30 nt upstream of
# the poly-A tail.
_POLYA_SIGNAL = "AATAAA"


@dataclass
class ConstructConfig:
    """Configuration for the mRNA construct designer.

    Attributes
    ----------
    species : str
        Codon-usage species for CDS optimization. Currently only
        ``"human"`` is bundled (the HUMAN_CODON_FREQ table).
    utr5 : str
        5'UTR sequence. Defaults to the human consensus
        (``_HUMAN_5UTR_CONSENSUS``).
    utr3 : str
        3'UTR sequence. Defaults to the human consensus
        (``_HUMAN_3UTR_CONSENSUS``).
    poly_a_length : int
        Length of the poly-A tail (default 120 — typical for
        therapeutic mRNAs).
    optimize_backend : str
        CDS optimization backend. One of ``"basic"`` (greedy),
        ``"ribodecode"``, ``"lineardesign"``, ``"ribodecode-real"``,
        ``"multi-objective"`` (v0.25.0). Default ``"multi-objective"``
        since it's the canonical SOTA pattern.
    """

    species: str = "human"
    utr5: str = _HUMAN_5UTR_CONSENSUS
    utr3: str = _HUMAN_3UTR_CONSENSUS
    poly_a_length: int = 120
    optimize_backend: str = "multi-objective"


@dataclass
class ConstructResult:
    """Result of designing an mRNA construct."""

    input_protein: str
    construct_dna: str  # full 5'UTR + CDS + 3'UTR + polyA + signal
    utr5: str
    cds: str
    utr3: str
    poly_a_tail: str
    poly_a_signal: str
    cds_length: int
    construct_length: int
    gc_percent: float
    cds_analysis: dict
    n_codons: int
    optimization_notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)


def _clean_cds(cds: str) -> str:
    """Local helper — strip whitespace, normalize to ACGT, trim to multiple-of-3.

    Mirrors the private helper in codon_optimizer.py so we don't depend
    on a private symbol across modules.
    """
    cds = re.sub(r"[^ATCG]", "", cds.upper().replace("U", "T"))
    if len(cds) >= 3 and cds[-3:] in CODON_TO_AA and CODON_TO_AA[cds[-3:]] == "*":
        cds = cds[:-3]
    trim = len(cds) % 3
    if trim:
        cds = cds[:-trim]
    return cds


# Standard genetic code — only the standard codon table is needed for
# reverse translation. Pull from HUMAN_CODON_FREQ (most-frequent codon
# per AA) — the most-frequent codon is the "default" reverse-translation
# codon in the standard genetic code.
_REVERSE_TRANSLATION_TABLE: dict[str, str] = {
    aa: max(codons, key=codons.get) for aa, codons in HUMAN_CODON_FREQ.items()
}


def _reverse_translate(protein: str) -> str:
    """Convert a protein sequence to a default DNA CDS.

    Uses the most-frequent human codon per amino acid. The CDS
    optimization step is where alternative synonymous codons get
    chosen; this is just the starting point.

    Returns the CDS without a stop codon — the stop codon is the user's
    responsibility to add (or omitted entirely if the input ends in
    ``*``, which is then included as ``TAA``).
    """
    protein = protein.upper().strip()
    if not protein:
        raise ValueError("protein sequence is empty")
    for aa in protein:
        if aa not in _REVERSE_TRANSLATION_TABLE:
            raise ValueError(
                f"unknown amino acid: {aa!r} (only standard 20 AAs + "
                f"stop '*' are supported)"
            )
    # Translate each AA to its default codon
    cds = "".join(_REVERSE_TRANSLATION_TABLE[aa] for aa in protein)
    return cds


def _optimize_cds(cds: str, backend: str) -> tuple[str, list[str]]:
    """Run the requested CDS-optimization backend.

    Returns (optimized_cds, notes).
    """
    notes: list[str] = []
    if backend == "basic":
        out = optimize_basic(cds)
        return out["new_cds"], notes
    if backend == "multi-objective":
        from .codon_multi_objective import (
            MultiObjectiveConfig,
            multi_objective_optimize,
        )

        cfg = MultiObjectiveConfig()
        result = multi_objective_optimize(cds, cfg)
        notes.append(
            f"multi-objective score: {result.overall_before:.3f} -> "
            f"{result.overall_after:.3f} (+{result.improvement:.3f})"
        )
        return result.optimized_cds, notes
    if backend == "ribodecode":
        from .codon_ribodecode import optimize_ribodecode

        result = optimize_ribodecode(cds)
        return result.new_cds, notes
    if backend == "lineardesign":
        from .codon_lineardesign import optimize_lineardesign

        result = optimize_lineardesign(cds)
        return result.new_cds, notes
    if backend == "ribodecode-real":
        # Falls back to the in-house heuristic if the upstream CLI is
        # not installed (same pattern as the codon CLI).
        from .codon_protocols import RiboDecodeRequest
        from .codon_ribodecode_adapter import select_codon_optimizer

        optimizer = select_codon_optimizer(prefer="auto")
        req = RiboDecodeRequest(cds=cds, env="HEK293T")
        try:
            ribo_res = optimizer.optimize(req)
            return ribo_res.optimized_cds, notes
        except Exception as exc:  # noqa: BLE001
            notes.append(f"ribodecode-real unavailable ({exc}); using basic")
            out = optimize_basic(cds)
            return out["new_cds"], notes
    raise ValueError(f"unknown optimize backend: {backend!r}")


def _gc_percent(seq: str) -> float:
    """GC% in the range [0, 100]."""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in "GC")
    return round(100.0 * gc / len(seq), 2)


# ---------- Public entry point ------------------------------------------


def design_construct(
    protein: str,
    config: ConstructConfig | None = None,
) -> ConstructResult:
    """Compose a full mRNA construct from a protein sequence.

    Parameters
    ----------
    protein : str
        Amino-acid sequence (single-letter codes). The stop codon
        (``*``) is included if the user adds it explicitly; otherwise
        the CDS has no stop codon (which is fine for codon-level
        analysis — add the stop yourself for full-length constructs).
    config : ConstructConfig, optional
        Construct configuration. Defaults to ``ConstructConfig()``
        (human consensus UTRs, 120-nt poly-A tail, multi-objective
        CDS optimization).

    Returns
    -------
    ConstructResult
        Full construct plus per-region breakdown.
    """
    if config is None:
        config = ConstructConfig()

    # Reverse-translate to a default CDS
    cds_input = _reverse_translate(protein)

    # Optimize the CDS using the requested backend
    cds_opt, notes = _optimize_cds(cds_input, config.optimize_backend)
    cds_opt = _clean_cds(cds_opt)

    # Assemble the construct: 5'UTR + CDS + 3'UTR + polyA_signal + polyA
    # Note: the polyA signal (AAUAAA) is placed at the end of the 3'UTR
    # before the poly-A tail. Many published designs place it inside
    # the 3'UTR; we follow that convention.
    poly_a_tail = "A" * config.poly_a_length
    construct = config.utr5 + cds_opt + config.utr3 + _POLYA_SIGNAL + poly_a_tail

    # Analyze the optimized CDS
    cds_analysis = analyze_cds(cds_opt).to_dict()

    return ConstructResult(
        input_protein=protein.upper().strip(),
        construct_dna=construct,
        utr5=config.utr5,
        cds=cds_opt,
        utr3=config.utr3,
        poly_a_tail=poly_a_tail,
        poly_a_signal=_POLYA_SIGNAL,
        cds_length=len(cds_opt),
        construct_length=len(construct),
        gc_percent=_gc_percent(construct),
        cds_analysis=cds_analysis,
        n_codons=len(cds_opt) // 3,
        optimization_notes=notes,
    )


# ---------- CLI ---------------------------------------------------------


def _read_fasta_single(text: str) -> str:
    """Read a single sequence from FASTA text (or raw string)."""
    text = text.strip()
    if text.startswith(">"):
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
                header = line[1:].strip()
                buf = []
            else:
                buf.append(re.sub(r"\s+", "", line))
        if header:
            records.append((header, "".join(buf).upper()))
        return records[0][1] if records else ""
    # Raw string — strip whitespace
    return re.sub(r"\s+", "", text).upper()


def _run_cli(argv: list[str]) -> int:
    import argparse
    import json as _json

    p = argparse.ArgumentParser(prog="mrnavax construct")
    p.add_argument("--sequence", required=True,
                   help="Protein AA sequence (FASTA file or raw string)")
    p.add_argument("--backend", default="multi-objective",
                   choices=["basic", "ribodecode", "lineardesign",
                            "ribodecode-real", "multi-objective"],
                   help="CDS optimization backend (default: multi-objective)")
    p.add_argument("--poly-a-length", type=int, default=120,
                   help="Length of the poly-A tail in nt (default: 120)")
    p.add_argument("--species", default="human",
                   choices=["human"],
                   help="Codon-usage species (default: human)")
    p.add_argument("--out", help="write JSON report here")
    args = p.parse_args(argv)

    raw = (Path(args.sequence).read_text()
           if Path(args.sequence).exists() else args.sequence)
    protein = _read_fasta_single(raw)
    if not protein:
        print("error: empty protein sequence", file=__import__("sys").stderr)
        return 2

    cfg = ConstructConfig(
        species=args.species,
        poly_a_length=args.poly_a_length,
        optimize_backend=args.backend,
    )
    result = design_construct(protein, cfg)
    out_text = _json.dumps(result.to_dict(), indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    import sys

    sys.exit(_run_cli(sys.argv[1:]))
