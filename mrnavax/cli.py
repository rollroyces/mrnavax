"""Unified CLI: `python -m mrnavax.cli <tool> ...`"""

from __future__ import annotations

import argparse
import json as _json
import sys
from pathlib import Path

from .alphagenome_integration import score_variants_from_csv
from .codon_optimizer import _run_cli as codon_run
from .lnp_advisor import _run_cli as lnp_run
from .manufacturability import score_manufacturability
from .neoantigen_screener import _run_cli as neo_run
from .sc_rna_pipeline import _run_cli as scrna_run
from .spatial_module_adapter import (
    STModuleNotInstalled,
    select_spatial_module_backend,
)
from .spatial_protocols import SpatialData
from .trial_matcher import _run_cli as trial_run


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="mrnavax",
        description="mRNA × AI toolkit "
        "(codon / neoantigen / trial / lnp / scrna / manufacture / spatial / construct)",
    )
    sub = p.add_subparsers(dest="tool", required=True)

    sub.add_parser("codon", help="codon-usage analysis & greedy optimization")
    sub.add_parser("neoantigen", help="peptide×HLA neoantigen screen")
    sub.add_parser("trial", help="patient-to-trial eligibility matching")
    sub.add_parser("lnp", help="LNP composition recommender")
    sub.add_parser("scrna", help="scRNA-seq → neoantigen handoff pipeline")
    sub.add_parser("manufacture", help="wet-lab manufacturability checks on a CDS")
    sub.add_parser(
        "spatial",
        help="spatial-transcriptomics tissue-module identification "
        "(STModule CLI or stdlib mock)",
    )
    sub.add_parser(
        "variant-regulatory",
        help="AlphaGenome Atlas regulatory-variant impact scoring "
        "(AVI score for non-coding regions; AlphaMissense handles coding)",
    )
    sub.add_parser(
        "construct",
        help="compose a full mRNA construct (5'UTR + CDS + 3'UTR + poly-A) "
        "from a protein amino-acid sequence",
    )

    args, rest = p.parse_known_args(argv)
    runners = {
        "codon": codon_run,
        "neoantigen": neo_run,
        "trial": trial_run,
        "lnp": lnp_run,
        "scrna": scrna_run,
    }
    if args.tool == "manufacture":
        return _manufacture_run(rest)
    if args.tool == "spatial":
        return _spatial_run(rest)
    if args.tool == "variant-regulatory":
        return _variant_regulatory_run(rest)
    if args.tool == "construct":
        return _construct_run(rest)
    return runners[args.tool](rest)


def _manufacture_run(argv: list[str]) -> int:
    """CLI for the manufacturability checker."""
    p = argparse.ArgumentParser(prog="mrnavax manufacture")
    p.add_argument(
        "--cds",
        required=True,
        help="FASTA file or raw CDS string (DNA, multiples of 3)",
    )
    p.add_argument(
        "--utr5",
        default="",
        help="optional 5' UTR (DNA) for Kozak scoring",
    )
    p.add_argument(
        "--utr3",
        default="",
        help="optional 3' UTR (DNA) for AU-rich element detection",
    )
    p.add_argument(
        "--out",
        default=None,
        help="output JSON file (default: stdout)",
    )
    args = p.parse_args(argv)

    cds = Path(args.cds).read_text() if Path(args.cds).exists() else args.cds

    if not cds.strip():
        print("error: empty CDS", file=sys.stderr)
        return 2
    report = score_manufacturability(
        cds,
        utr5=args.utr5 if args.utr5 else None,
        utr3=args.utr3 if args.utr3 else None,
    )
    out_text = _json.dumps(report.to_dict(), indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


def _spatial_run(argv: list[str]) -> int:
    """CLI for spatial-transcriptomics tissue-module identification."""
    import json as _json_spatial

    p = argparse.ArgumentParser(prog="mrnavax spatial")
    p.add_argument(
        "--count-file",
        required=True,
        help="TSV count matrix (spots × genes). Header row + spot IDs in "
        "first unnamed column.",
    )
    p.add_argument(
        "--locations-file",
        required=True,
        help="TSV spatial coordinates. Header row + spot IDs in first "
        "unnamed column, columns 'x' and 'y'.",
    )
    p.add_argument(
        "--platform",
        choices=["ST", "Visium", "SlideSeqV2", "StereoSeq", "Other"],
        default="ST",
        help="SRT platform — drives preprocessing and max_iter defaults "
        "in upstream STModule",
    )
    p.add_argument(
        "--num-modules",
        type=int,
        default=10,
        help="number of tissue modules to identify (default: 10)",
    )
    p.add_argument(
        "--backend",
        choices=["auto", "mock", "stmodule"],
        default="auto",
        help="backend: auto picks stmodule when Rscript on $PATH, "
        "else mock",
    )
    p.add_argument("--out", default=None)
    args = p.parse_args(argv)

    try:
        data = SpatialData(
            count_file=Path(args.count_file),
            locations_file=Path(args.locations_file),
            platform=args.platform,
            num_modules=args.num_modules,
        )
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    if args.backend == "auto":
        backend = select_spatial_module_backend()
    elif args.backend == "mock":
        from .spatial_module_adapter import MockSpatialModuleBackend

        backend = MockSpatialModuleBackend()
    else:
        from .spatial_module_adapter import STModuleCLIAdapter

        backend = STModuleCLIAdapter()

    try:
        result = backend.run(data)
    except STModuleNotInstalled as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    out_text = _json_spatial.dumps(result.to_dict(), indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


def _variant_regulatory_run(argv: list[str]) -> int:
    """CLI for AlphaGenome Atlas regulatory-variant impact scoring.

    Reads a CSV with columns: chrom,pos,ref,alt,label (label optional)
    and emits a JSON document with one entry per variant: AVI score,
    classification, and is_coding flag. Variants are sorted by AVI
    score descending (most-impactful first).

    Without an API key (or with --backend=mock), uses a deterministic
    stdlib mock. With --backend=alphagenome + ALPHAGENOME_API_KEY env
    var, calls the official AlphaGenome Atlas API via subprocess.
    """
    import json as _json_vr

    p = argparse.ArgumentParser(prog="mrnavax variant-regulatory")
    p.add_argument(
        "--csv",
        required=True,
        help="CSV with columns chrom,pos,ref,alt,label (label optional)",
    )
    p.add_argument(
        "--backend",
        choices=["auto", "mock", "alphagenome"],
        default="auto",
        help="auto picks alphagenome when ALPHAGENOME_API_KEY is set, "
        "else mock",
    )
    p.add_argument("--out", default=None, help="output JSON file (default: stdout)")
    args = p.parse_args(argv)

    if args.backend == "mock":
        from .alphagenome_integration import MockRegulatoryVariantScorer

        backend = MockRegulatoryVariantScorer()
    elif args.backend == "alphagenome":
        from .alphagenome_integration import AlphaGenomeCLIAdapter

        try:
            backend = AlphaGenomeCLIAdapter()
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
    else:
        from .alphagenome_integration import select_regulatory_scorer

        backend = select_regulatory_scorer()

    try:
        results = score_variants_from_csv(args.csv, backend=backend)
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    out_text = _json_vr.dumps(
        {
            "tool": "variant-regulatory",
            "backend": type(backend).__name__,
            "n_variants": len(results),
            "results": [
                {
                    "score": r.score,
                    "classification": r.classification,
                    "is_coding": r.is_coding,
                }
                for r in results
            ],
        },
        indent=2,
    )
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


def _construct_run(argv: list[str]) -> int:
    """CLI for the mRNA construct designer (5'UTR + CDS + 3'UTR + poly-A)."""
    from .construct_designer import (
        ConstructConfig,
        _read_fasta_single,
        design_construct,
    )

    p = argparse.ArgumentParser(prog="mrnavax construct")
    p.add_argument(
        "--sequence",
        required=True,
        help="Protein AA sequence (FASTA file or raw string)",
    )
    p.add_argument(
        "--backend",
        choices=["basic", "ribodecode", "lineardesign",
                 "ribodecode-real", "multi-objective"],
        default="multi-objective",
        help="CDS optimization backend (default: multi-objective, the SOTA pattern)",
    )
    p.add_argument(
        "--poly-a-length", type=int, default=120,
        help="Length of the poly-A tail in nt (default: 120)",
    )
    p.add_argument(
        "--species", default="human", choices=["human"],
        help="Codon-usage species (default: human)",
    )
    p.add_argument(
        "--out", default=None,
        help="Write JSON report here (default: stdout)",
    )
    args = p.parse_args(argv)

    raw = Path(args.sequence).read_text() if Path(args.sequence).exists() else args.sequence
    protein = _read_fasta_single(raw)
    if not protein:
        print("error: empty protein sequence", file=sys.stderr)
        return 2

    cfg = ConstructConfig(
        species=args.species,
        poly_a_length=args.poly_a_length,
        optimize_backend=args.backend,
    )
    try:
        result = design_construct(protein, cfg)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2

    out_text = _json.dumps(result.to_dict(), indent=2)
    if args.out:
        Path(args.out).write_text(out_text + "\n")
    print(out_text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
