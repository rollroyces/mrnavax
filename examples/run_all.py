#!/usr/bin/env python3
"""End-to-end mrnavax demo.

Runs all 11 tools on the example files in `examples/` and prints a
one-line summary per tool. Lets you verify your install is healthy
in under 30 seconds (mock mode, no model downloads).

For real Atlas calls (`mrnavax predict`), you need:

    pip install -e '.[variant-alphagenome]'
    echo "AIza..." > ~/projects/alphagenome-work/.alphagenome_key

Usage:
    python examples/run_all.py            # mock mode
    MRNA_AI_FORCE_MOCK=1 python examples/run_all.py   # explicit
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = REPO_ROOT / "examples"


def _run(tool: str, *args: str) -> dict | None:
    """Invoke a tool via the CLI and parse its JSON output."""
    cmd = [sys.executable, "-m", "mrnavax.cli", tool, *args]
    env = {**os.environ, "MRNA_AI_FORCE_MOCK": "1"}
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=60, cwd=REPO_ROOT, env=env
        )
    except subprocess.TimeoutExpired:
        return {"error": "timeout"}
    if result.returncode != 0:
        # Some tools log warnings to stdout but still return valid JSON.
        try:
            return json.loads(result.stdout)
        except json.JSONDecodeError:
            return {"error": result.stderr.strip().splitlines()[-1] if result.stderr else "no output"}
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        return {"raw": result.stdout.strip()[:200]}


def _summary(label: str, value: str) -> None:
    print(f"  {label:<28} {value}")


def main() -> int:
    """Run all 11 tools and print a one-line summary per tool."""
    print("mrnavax end-to-end demo (mock mode)\n")

    # 1. codon
    out = _run("codon", "--sequence", str(EXAMPLES / "cas9.fasta"))
    if out and "cai" in out:
        _summary("codon", f"CAI={out['cai']:.3f}, n_codons={out['n_codons']}")
    else:
        _summary("codon", f"FAIL: {out}")

    # 2. neoantigen
    out = _run(
        "neoantigen",
        "--variants",
        str(EXAMPLES / "tp53_variants.csv"),
        "--hla",
        "HLA-A*02:01",
    )
    if out and "candidates" in out:
        n_top = sum(1 for c in out["candidates"] if c.get("binder"))
        _summary("neoantigen", f"{len(out['candidates'])} candidates, {n_top} binders")
    else:
        _summary("neoantigen", f"FAIL: {out}")

    # 3. trial
    out = _run(
        "trial",
        "--patient",
        str(EXAMPLES / "patient_summary.txt"),
        "--trials",
        str(EXAMPLES / "trials.jsonl"),
        "--top-k",
        "5",
    )
    if out and "ranked" in out:
        _summary("trial", f"top-1 = {out['ranked'][0].get('nct_id', '?')}")
    else:
        _summary("trial", f"FAIL: {out}")

    # 4. lnp
    out = _run("lnp", "--target", "lung", "--cargo", "saRNA", "--intent", "cancer vaccine")
    if out and "cargo" in out:
        _summary("lnp", f"target={out['target']}, cargo={out['cargo']}")
    else:
        _summary("lnp", f"FAIL: {out}")

    # 5. scrna
    out = _run(
        "scrna",
        "--expression",
        str(EXAMPLES / "cells.csv"),
        "--variants",
        str(EXAMPLES / "variants_coding.csv"),
        "--proteins",
        str(EXAMPLES / "proteins.fasta"),
        "--tumor-markers",
        "TP53,KRAS,BRAF",
    )
    if out and "n_cells" in out:
        n_clusters = len(set(out.get("cluster_labels", [])))
        _summary(
            "scrna",
            f"{out['n_cells']} cells × {out['n_genes']} genes, "
            f"{n_clusters} clusters (tumor={out.get('tumor_cluster', '?')})",
        )
    else:
        _summary("scrna", f"FAIL: {out}")

    # 6. manufacture
    out = _run("manufacture", "--cds", str(EXAMPLES / "cds_gfp.json"))
    if out and "overall_score" in out:
        n_total = (out.get("n_pass", 0) + out.get("n_warn", 0) + out.get("n_error", 0))
        _summary(
            "manufacture",
            f"score={out['overall_score']:.3f}, "
            f"{out.get('n_pass', '?')} pass / {n_total} total",
        )
    else:
        _summary("manufacture", f"FAIL: {out}")

    # 7. spatial
    out = _run(
        "spatial",
        "--count-file",
        str(EXAMPLES / "st_bc2_count_matrix.tsv"),
        "--locations-file",
        str(EXAMPLES / "st_bc2_locations.tsv"),
        "--platform",
        "ST",
        "--num-modules",
        "10",
    )
    if out and "modules" in out:
        _summary("spatial", f"{len(out['modules'])} tissue modules")
    else:
        _summary("spatial", f"FAIL: {out}")

    # 8. construct (uses --sequence, not --protein)
    out = _run(
        "construct",
        "--sequence",
        "ATGGATAAGCTACCGTATCGTATGAATAA",
    )
    if out and "construct_dna" in out:
        cds_analysis = out.get("cds_analysis", {})
        cai = cds_analysis.get("cai", 0.0) if isinstance(cds_analysis, dict) else 0.0
        _summary(
            "construct",
            f"{out['construct_length']} nt, CAI={cai:.3f}",
        )
    else:
        _summary("construct", f"FAIL: {out}")

    # 9. utr-design
    out = _run("utr-design", "--protein", "MVSKGEELFTGVVPILVELDGDVNGHKFS")
    if out and "combined_score" in out:
        _summary("utr-design", f"combined={out['combined_score']:.3f}")
    else:
        _summary("utr-design", f"FAIL: {out}")

    # 10. variant-regulatory
    out = _run("variant-regulatory", "--csv", str(EXAMPLES / "regulatory_variants.csv"))
    if out and "results" in out:
        n_high = sum(1 for r in out["results"] if r.get("classification") == "high")
        _summary(
            "variant-regulatory",
            f"{len(out['results'])} variants, {n_high} high-impact",
        )
    else:
        _summary("variant-regulatory", f"FAIL: {out}")

    # 11. predict (only runs in live mode if ALPHAGENOME_API_KEY is set)
    api_key = os.environ.get("ALPHAGENOME_API_KEY") or (
        (Path.home() / "projects" / "alphagenome-work" / ".alphagenome_key").exists()
    )
    if api_key:
        out = _run(
            "predict",
            "--chrom",
            "chr7",
            "--pos",
            "140753336",
            "--ref",
            "T",
            "--alt",
            "A",
            "--name",
            "BRAF_V600E",
        )
        results = out.get("results", []) if isinstance(out, dict) else []
        if results:
            r0 = results[0]
            _summary(
                "predict",
                f"score={r0['score']:.3f}, "
                f"classification={r0.get('classification', '?')}",
            )
        elif out and "error" not in out:
            _summary("predict", f"unexpected shape: {list(out.keys()) if isinstance(out, dict) else type(out).__name__}")
        else:
            _summary("predict", f"FAIL: {out}")
    else:
        _summary(
            "predict",
            "skipped (no API key — see docs/operations/atlas-activation.md)",
        )

    print("\nDone. Sample outputs committed under examples/sample_outputs/.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
