# `predict` — live AlphaGenome Atlas regulatory-variant prediction

The **11th tool** in the mrnavax toolkit. Calls the real AlphaGenome
Atlas API to score the regulatory impact (AVI score) of one or more
variants. This is the only tool in the toolkit that makes **live
network calls** by default; everything else is local.

## Why a separate CLI tool?

The v0.25.1 release added the helper-file auto-detection for
`ALPHAGENOME_API_KEY`. With that in place, scoring a variant
end-to-end via the live Atlas API becomes a one-liner. The `predict`
tool wraps that one-liner with:

- Single-variant and CSV batch modes
- Automatic API key resolution (env var / helper file)
- Robust subprocess error handling
- JSON output suitable for downstream pipelines

## Quick start

### Single variant

```bash
mrnavax predict --chrom chr7 --pos 140753336 --ref T --alt A --name BRAF_V600E
```

### Batch from CSV

```csv
name,chrom,pos,ref,alt
BRAF_V600E,chr7,140753336,T,A
KRAS_G12D,chr12,25245350,C,T
```

```bash
mrnavax predict --input variants.csv
```

### Output

```json
{
  "n_variants": 1,
  "results": [
    {
      "name": "BRAF_V600E",
      "chrom": "chr7",
      "pos": 140753336,
      "ref": "T",
      "alt": "A",
      "score": 1.0,
      "classification": "high",
      "is_coding": true,
      "n_scorers": 19
    }
  ]
}
```

## Output fields

| Field | Type | Meaning |
|---|---|---|
| `score` | float in [0, 1] | Max-abs delta across all 19 Atlas scorers, clipped to [0, 1]. Higher = larger regulatory impact. |
| `classification` | str in {"low", "moderate", "high"} | Canonical AVI bins (low < 0.34, moderate < 0.564, high otherwise). |
| `is_coding` | bool | True if the variant position overlaps a protein_coding gene in the 16kb window. Best-effort; use as a hint, not a precise annotation. |
| `n_scorers` | int | Number of Atlas scorers that returned a response (typically 19). |

## API key resolution

The shim auto-resolves the API key in this order:

1. `api_key=...` explicit (in the JSON payload — for direct programmatic use)
2. `ALPHAGENOME_API_KEY` env var (CI path)
3. `~/projects/alphagenome-work/.alphagenome_key` helper file (local dev)
4. `~/.alphagenome_key` alternate location

If none of these resolves, the shim returns exit code 1 with
`"API key not configured"` in stderr.

## CLI vs. programmatic

For programmatic use, prefer the Python API:

```python
import json, subprocess
from pathlib import Path

shim = Path("/path/to/mrnavax/_shims/alphagenome_cli.py")
proc = subprocess.run(
    [sys.executable, str(shim)],
    input=json.dumps({"chrom": "chr7", "pos": 140753336, "ref": "T", "alt": "A"}),
    capture_output=True, text=True, timeout=180,
)
result = json.loads(proc.stdout)
print(result["score"], result["classification"])
```

The `predict` CLI is just a thin wrapper around this subprocess call.

## Performance

- Single variant: ~3–5 seconds (after API warmup). Cold-start can be
  ~30 seconds (first call).
- Batch: sequential; ~3–5 seconds per variant. No rate-limit handling
  is built in; if you need thousands of variants, batch through the
  Python API directly with concurrency control.

## Cost / quota

The Atlas API is free for academic use (preview); check the current
Google AI Studio Atlas terms before running large batches. The
toolkit makes no claims about quota; the shim passes through whatever
the upstream API returns.

## What `predict` is NOT

- **Not a pathogenicity predictor**: AVI score measures regulatory
  impact, not disease causation. Use AlphaMissense (for coding
  variants) + clinical genetics interpretation for diagnostic claims.
- **Not a structural-impact predictor**: Use AlphaMissense for protein
  structure / function effects.
- **Not a polygenic risk score**: Atlas reports single-variant impact;
  PRS requires aggregating many variants.

## Reference

- Avsec et al., "AlphaGenome: AI for genomics at scale", *Nature*
  2026 (DOI placeholder until publication).
- https://deepmind.google.com/science/alphagenome — official API docs.