"""Live Atlas prediction demo — score canonical cancer variants.

Run with:
    .venv-am/bin/python scripts/predict_demo.py

This exercises the full prediction stack on 6 canonical cancer variants
spanning 5 genes (KRAS, BRAF, EGFR, TP53, MYC). Each variant has a
known clinical interpretation; comparing the live Atlas prediction
against the literature is a meaningful end-to-end test.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SHIM = REPO_ROOT / "mrnavax" / "_shims" / "alphagenome_cli.py"


# Canonical cancer variants with their hg38 coordinates and known
# clinical interpretations. Coordinates cross-checked against ClinVar
# / COSMIC v100.
VARIANTS = [
    {
        "name": "BRAF V600E",
        "chrom": "chr7",
        "pos": 140753336,
        "ref": "T",
        "alt": "A",
        "expected_class": "high",
        "is_coding_expected": True,
        "context": "BRAF V600E — canonical melanoma activating mutation",
    },
    {
        "name": "KRAS G12D",
        "chrom": "chr12",
        "pos": 25245350,
        "ref": "C",
        "alt": "T",
        "expected_class": "high",
        "is_coding_expected": True,
        "context": "KRAS G12D — pancreatic / colorectal oncogenic",
    },
    {
        "name": "EGFR L858R",
        "chrom": "chr7",
        "pos": 55191822,
        "ref": "T",
        "alt": "G",
        "expected_class": "high",
        "is_coding_expected": True,
        "context": "EGFR L858R — NSCLC activating mutation",
    },
    {
        "name": "TP53 R175H",
        "chrom": "chr17",
        "pos": 7675088,
        "ref": "C",
        "alt": "T",
        "expected_class": "high",
        "is_coding_expected": True,
        "context": "TP53 R175H — Li-Fraumeni hotspot loss-of-function",
    },
    {
        "name": "MYC T58A",
        "chrom": "chr8",
        "pos": 127736664,
        "ref": "A",
        "alt": "G",
        "expected_class": "moderate",
        "is_coding_expected": True,
        "context": "MYC T58A — Burkitt lymphoma mutation",
    },
    {
        "name": "EGFR T790M",
        "chrom": "chr7",
        "pos": 55181378,
        "ref": "C",
        "alt": "T",
        "expected_class": "high",
        "is_coding_expected": True,
        "context": "EGFR T790M — TKI resistance mutation",
    },
]


def predict_one(variant: dict) -> dict | None:
    """Run the live Atlas call on one variant. Returns the result dict or None."""
    payload = json.dumps({k: v for k, v in variant.items() if k not in {"name", "expected_class", "is_coding_expected", "context"}})
    proc = subprocess.run(
        [sys.executable, str(SHIM)],
        input=payload,
        capture_output=True,
        text=True,
        timeout=120,
    )
    if proc.returncode != 0:
        print(f"  [ERROR] {variant['name']}: {proc.stderr.strip()}")
        return None
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError as exc:
        print(f"  [ERROR] {variant['name']}: bad JSON: {exc}")
        return None


def main():
    print("=" * 70)
    print("Live Atlas prediction: 6 canonical cancer variants")
    print("=" * 70)
    print()
    print(f"{'Variant':<14} {'Score':>6}  {'Class':<8}  {'Cod':<3}  {'nSc':>4}  {'Match':<5}  Context")
    print("-" * 70)

    matches = 0
    for v in VARIANTS:
        result = predict_one(v)
        if result is None:
            print(f"{v['name']:<14} FAILED")
            continue
        score = result["score"]
        classification = result["classification"]
        is_coding = result["is_coding"]
        n_scorers = result["n_scorers"]
        # Naive match check: AVI class agrees with expected + is_coding agrees
        match = (
            classification == v["expected_class"]
            and is_coding == v["is_coding_expected"]
        )
        if match:
            matches += 1
        match_marker = "✅" if match else "❌"
        print(
            f"{v['name']:<14} {score:>6.3f}  {classification:<8}  "
            f"{'Y' if is_coding else 'N':<3}  {n_scorers:>4}  "
            f"{match_marker:<5}  {v['context']}"
        )

    print("-" * 70)
    print(f"Match: {matches} / {len(VARIANTS)}")
    if matches == len(VARIANTS):
        print("ALL VARIANTS MATCH EXPECTED → live Atlas is calibrated correctly.")
    elif matches >= len(VARIANTS) * 0.8:
        print("Most variants match — Atlas is operational; small miscalibrations")
        print("are expected since Atlas reports regulatory impact, not pathogenicity.")
    else:
        print("Several mismatches — Atlas calibration may have drifted;")
        print("inspect the live results above and consider re-running the weekly")
        print("regression check to refresh the baseline.")
    return 0 if matches >= len(VARIANTS) * 0.8 else 1


if __name__ == "__main__":
    sys.exit(main())
