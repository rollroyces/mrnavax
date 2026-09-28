"""Live Atlas regulatory-variant prediction demo.

Tests 6 well-characterized non-coding regulatory variants from the
GWAS catalog. Atlas is designed for regulatory impact prediction —
these are the canonical use case.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SHIM = REPO_ROOT / "mrnavax" / "_shims" / "alphagenome_cli.py"

# Non-coding regulatory variants with hg38 coordinates. The "context"
# describes why each variant is interesting from a regulatory-impact
# perspective.
REG_VARIANTS = [
    {
        "name": "T2D_rsjp_TCF7L2",
        "chrom": "chr10",
        "pos": 114757988,
        "ref": "G",
        "alt": "A",
        "context": "T2D GWAS risk SNP (rs7903146) — TCF7L2 intron; cis-eQTL",
    },
    {
        "name": "LPA_CAD_rskb",
        "chrom": "chr6",
        "pos": 161137102,
        "ref": "G",
        "alt": "A",
        "context": "CAD risk SNP (rs10455872) — LPA intronic",
    },
    {
        "name": "BC_FGFR2_rs298",
        "chrom": "chr10",
        "pos": 121592059,
        "ref": "C",
        "alt": "T",
        "context": "Breast cancer risk (rs2981578) — FGFR2 intron",
    },
    {
        "name": "RA_TYK2_protect",
        "chrom": "chr19",
        "pos": 10365682,
        "ref": "G",
        "alt": "C",
        "context": "RA protective coding (rs34536443) — TYK2 P1104P",
    },
    {
        "name": "T1D_IL2RA_enh",
        "chrom": "chr10",
        "pos": 6099047,
        "ref": "C",
        "alt": "T",
        "context": "Type-1 diabetes — IL2RA enhancer region",
    },
    {
        "name": "SCZ_cand_rs107",
        "chrom": "chr11",
        "pos": 113906690,
        "ref": "C",
        "alt": "T",
        "context": "Schizophrenia candidate region",
    },
]


def predict_one(variant: dict) -> dict | None:
    payload = json.dumps({k: v for k, v in variant.items() if k != "context"})
    proc = subprocess.run(
        [sys.executable, str(SHIM)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0:
        return None
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None


def main():
    print("=" * 78)
    print("Live Atlas REGULATORY-variant prediction: 6 GWAS SNPs")
    print("=" * 78)
    print()
    print(f"{'Variant':<22} {'Score':>6}  {'Class':<8}  Context")
    print("-" * 78)

    for v in REG_VARIANTS:
        result = predict_one(v)
        if result is None:
            print(f"{v['name']:<22} FAILED")
            continue
        print(
            f"{v['name']:<22} {result['score']:>6.3f}  "
            f"{result['classification']:<8}  {v['context']}"
        )
    print("-" * 78)
    print()
    print("Interpretation:")
    print("  Atlas reports regulatory impact (AVI score) on the surrounding")
    print("  16kb window. High scores mean the variant disrupts regulatory")
    print("  tracks strongly; low scores mean the variant is benign on")
    print("  chromatin / accessibility / expression tracks.")
    print()
    print("  Atlas does NOT claim disease causation — that requires")
    print("  genetic + clinical interpretation beyond the model's scope.")
    print("  Use these scores to PRIORITIZE variants for downstream")
    print("  functional validation, not to make diagnostic claims.")


if __name__ == "__main__":
    main()
