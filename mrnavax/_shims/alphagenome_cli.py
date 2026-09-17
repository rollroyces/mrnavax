#!/usr/bin/env python3
"""alphagenome-cli shim — small wrapper around the official alphagenome
Python package. Reads a JSON payload from stdin, calls
``alphagenome.models.dna_client.DnaClient.predict_variant``, emits a
single JSON document to stdout, and exits 0.

Heavy dep: ``alphagenome`` Python package (``pip install alphagenome``).
This shim keeps the adapter module stdlib-only.

This shim is bundled inside the ``mrnavax`` package at
``mrnavax/_shims/alphagenome_cli.py`` and shipped via the
``[variant-alphagenome]`` extra. The real adapter
(``mrnavax.alphagenome_integration.AlphaGenomeCLIAdapter``) calls this
shim via subprocess so that the heavy upstream dependency stays
opt-in.

Usage:
    echo '{"api_key": "...", "chrom": "chr7", "pos": 140753336, "ref": "T", "alt": "A"}' \\
        | python alphagenome_cli.py
    # -> {"score": 0.72, "classification": "high", "is_coding": true}

Exit codes:
    0 — success
    1 — upstream call failed (see stderr)
    2 — ``alphagenome`` package not installed
"""

from __future__ import annotations

import json
import sys


def _classify(score: float) -> str:
    """Map an AVI score to its bin label."""
    if score < 0.34:
        return "low"
    if score < 0.564:
        return "moderate"
    return "high"


def main() -> int:
    try:
        payload = json.loads(sys.stdin.read())
    except json.JSONDecodeError as exc:
        sys.stderr.write(f"invalid JSON payload: {exc}\n")
        return 2

    api_key = payload.get("api_key", "")
    chrom = payload.get("chrom", "")
    pos = int(payload.get("pos", 0))
    ref = payload.get("ref", "")
    alt = payload.get("alt", "")

    try:
        from alphagenome.data import genome
        from alphagenome.models import dna_client
    except ImportError:
        sys.stderr.write(
            "alphagenome package not installed; "
            "pip install -e '.[variant-alphagenome]'\n"
        )
        return 2

    try:
        model = dna_client.create(api_key)
        interval = genome.Interval(chrom, pos - 1, pos)
        result = model.predict_variant(
            interval=interval,
            variant=genome.Variant(chrom, pos, ref, alt),
            organism=genome.Organism.HOMO_SAPIENS,
        )
        # The Atlas exposes a single combined ``AVI score`` per variant.
        avi = float(getattr(result, "avi_score", 0.0))
        is_coding = bool(getattr(result, "is_coding", False))
        sys.stdout.write(json.dumps({
            "score": avi,
            "classification": _classify(avi),
            "is_coding": is_coding,
        }))
        return 0
    except Exception as exc:  # noqa: BLE001 — adapter boundary catches all
        sys.stderr.write(f"alphagenome error: {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
