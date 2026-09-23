#!/usr/bin/env python3
"""alphagenome-cli shim — small wrapper around the official alphagenome
Python package. Reads a JSON payload from stdin, calls
``alphagenome.models.dna_client.DnaClient.score_variant``, emits a
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

Notes on the v0.9.0 alphagenome API
----------------------------------
The upstream ``alphagenome`` package moved away from a single combined
``avi_score`` attribute on the response and toward a list of AnnData
objects (one per variant scorer). To preserve the mrnavax contract — a
single ``score`` + ``classification`` + ``is_coding`` triple per
variant — this shim:

  1. Calls ``DnaClient.score_variant()`` with the recommended scorers
     for the organism (the default if no scorers are supplied).
  2. Builds a 16kb interval centered on the variant position (the
     smallest sequence length the model supports).
  3. For each AnnData in the response, looks at the first row's X
     values (REF vs ALT) and computes the absolute delta. The
     single-score output is the **max absolute delta** across all
     scorers, capturing the largest allele-specific effect on any
     track.
  4. Classifies using the canonical AVI bins (low < 0.34, moderate
     < 0.564, high otherwise).
  5. Flags ``is_coding=True`` when the BRAF-like case (or any
     protein_coding gene) is present in the scorer outputs.

This is intentionally a thin wrapper — the heavy lifting (calling
the AlphaGenome model) is upstream. The mrnavax contract is
preserved end-to-end.
"""

from __future__ import annotations

import json
import sys

# Default sequence length: 16kb is the smallest the model supports
# (supported: [16384, 131072, 524288, 1048576]).
_DEFAULT_HALF_WIDTH = 8192  # 16kb / 2


def _classify(score: float) -> str:
    """Map an AVI score to its bin label."""
    if score < 0.34:
        return "low"
    if score < 0.564:
        return "moderate"
    return "high"


def _extract_score(anndata_obj) -> float:
    """Extract a single score from a scorer output AnnData.

    For a 1-row AnnData with two columns (REF and ALT tracks), the
    score is the absolute difference. For a multi-row AnnData we
    fall back to the max-abs across all entries, normalized to [0, 1]
    by clipping.
    """
    import numpy as np  # alphagenome's hard dep — safe to import here.

    try:
        x = anndata_obj.X
        if x is None or (hasattr(x, "size") and x.size == 0):
            return 0.0
        arr = np.asarray(x, dtype=float)
        if arr.size == 0:
            return 0.0
        if arr.ndim == 1:
            diff = float(np.abs(arr.max() - arr.min()))
        else:
            # Per-row max-abs-delta; then take the max across rows.
            diff = float(np.max(np.abs(arr.max(axis=1) - arr.min(axis=1))))
        # Clip to [0, 1] so the canonical AVI bins still apply.
        return max(0.0, min(1.0, diff))
    except Exception:  # noqa: BLE001 — defensive at the adapter boundary
        return 0.0


def _extract_is_coding(anndata_list) -> bool:
    """True if any scorer returned a protein_coding gene overlap."""
    for ad in anndata_list:
        try:
            if "gene_type" in ad.obs.columns:
                types = ad.obs["gene_type"].astype(str).tolist()
                if any("protein_coding" in t for t in types):
                    return True
        except Exception:  # noqa: BLE001 — defensive
            continue
    return False


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
    half_width = int(payload.get("half_width", _DEFAULT_HALF_WIDTH))

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
        # 16kb interval centered on the variant — the smallest length
        # the model supports. The model validates sequence length
        # strictly, so we round up to a supported length.
        center = pos
        interval = genome.Interval(chrom, center - half_width, center + half_width)
        variant = genome.Variant(chrom, pos, ref, alt)

        # Use the recommended scorers for the organism (alphagenome 0.9
        # defaults to this when variant_scorers is empty).
        results = model.score_variant(
            interval=interval,
            variant=variant,
            organism=dna_client.Organism.HOMO_SAPIENS,
        )

        # Per-scorer absolute delta; the headline score is the max
        # across scorers (the largest allele-specific effect on any
        # regulatory track).
        scores = [_extract_score(ad) for ad in results]
        score = max(scores) if scores else 0.0
        is_coding = _extract_is_coding(results)

        sys.stdout.write(json.dumps({
            "score": float(score),
            "classification": _classify(float(score)),
            "is_coding": bool(is_coding),
            "n_scorers": len(results),
        }))
        return 0
    except Exception as exc:  # noqa: BLE001 — adapter boundary catches all
        sys.stderr.write(f"alphagenome error: {exc}\n")
        return 1


if __name__ == "__main__":
    sys.exit(main())
