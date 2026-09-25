# `utr-design` — coupled 5'UTR + CDS + 3'UTR design for max expression

The **10th tool** in the mrnavax toolkit. Couples the v0.27.0 UTR
context scorer (Kozak + 3'UTR quality) with the v0.25.0 multi-objective
CDS optimizer (RNop knowledge-infused loss pattern) to pick the best
combination for a given protein amino-acid sequence.

## Why "coupled"?

The v0.26.0 `construct` tool composes a 5'UTR + CDS + 3'UTR + poly-A
construct from fixed UTR templates — useful for assembly, but the UTR
choice is hardcoded. The v0.27.0 UTR context scorer can rank UTRs by
their impact on expression, but doesn't help choose one.

`utr-design` is the **bridge**: a bounded grid search over 12
candidate UTR combinations (4 5'UTR variants × 3 3'UTR variants),
scoring each via the v0.27.0 UTR context scorer, optimizing the CDS
once via the v0.25.0 multi-objective optimizer, and picking the
combination with the highest **joint score** (0.4 × UTR context + 0.6
× CDS).

## Quick start

```python
from mrnavax.utr_designer import UTRDesignConfig, design_utr_aware_cds

cfg = UTRDesignConfig(cds="MVSKGEELFTGV")  # 12-AA eGFP fragment
result = design_utr_aware_cds(cfg)

print(f"Chosen UTR5: {result.utr5}")
print(f"Chosen UTR3: {result.utr3}")
print(f"Optimized CDS: {result.cds_dna}")
print(f"Combined score: {result.combined_score:.3f}")
```

## CLI

```bash
mrnavax utr-design --protein "MVSKGEELFTGV"
mrnavax utr-design --protein my_protein.fasta --prefer-kozak 0.85
mrnavax utr-design --protein my_protein.fasta --library minimal
```

## Candidate UTR library

The default `library="all"` evaluates 12 combinations:

| 5'UTR variant | Source | Kozak score |
|---|---|---|
| `strong_kozak` | 9-nt slice `GCCACCAAT` (R=A purine) | ~0.9 |
| `moderate_kozak` | 9-nt slice `GCATGG` (deviant) | ~0.5 |
| `weak_kozak` | 9-nt slice `AAAAAA..` (no GCC) | ~0.1 |
| `mrna1273_like` | mRNA-1273 published 5'UTR | ~0.5 |

| 3'UTR variant | Source | ARE burden |
|---|---|---|
| `short_constitutive` | 32 nt GC-balanced, no ARE | ~1.0 |
| `mrna1273_like` | mRNA-1273 published 3'UTR | ~0.7 |
| `long_conservative` | 60 nt GC-balanced | ~1.0 |

The `library="minimal"` setting uses only the strong-Kozak +
short-constitutive pair (1 combination).

## Custom thresholds

By default the designer requires Kozak ≥ 0.7 and 3'UTR ≥ 0.5. If no
combination meets both thresholds, it falls back to the highest-ranked
overall combination.

```python
# Be strict about Kozak but lenient about 3'UTR
cfg = UTRDesignConfig(
    cds="MVSKGEELFTGV",
    prefer_kozak=0.85,
    prefer_utr3=0.3,
)
```

## Choosing a CDS backend

By default the CDS is optimized via the multi-objective backend
(v0.25.0 RNop knowledge-infused loss). The `basic` backend is also
available:

```python
cfg = UTRDesignConfig(cds="MVSKGEELFTGV", backend="basic")
```

The multi-objective backend is recommended because it accounts for
CAI + GC + CpG + rare-run + structure-proxy in a single weighted
loss, giving a more nuanced CDS choice than CAI+GC alone.

## Output

The result is a `UTRDesignResult` dataclass with all the information
you need to assemble the final construct:

```python
@dataclass
class UTRDesignResult:
    utr5: str
    utr3: str
    cds_dna: str
    protein: str
    utr_context: UTRContextResult
    cds_score: float
    combined_score: float
    candidates_evaluated: int
    ranking: list[tuple[float, str, str]]  # (combined_score, utr5_name, utr3_name)
```

Use `result.utr5 + result.cds_dna + result.utr3` as the
CDS-flanked-by-UTRs portion of a full mRNA construct.

## What's NOT optimized

- **5'UTR / 3'UTR sequence content**: we sample from a hand-curated
  library of 4 + 3 = 7 variants. For sequence-level UTR optimization
  (which is a different problem — gradient-based design over a
  differentiable expression model), use a dedicated tool like
  UTR-Function or 5UTR-Tailor. The library approach is more honest
  about uncertainty and avoids overfitting.
- **CDS per-position choice across UTR pairs**: the CDS is optimized
  once, then reused across all 12 UTR pairs. This is correct because
  the CDS is independent of the UTR choice, but it means we don't
  search over (CDS, UTR) jointly. In practice the dominant signal is
  the CDS itself, so this is fine.
- **Poly-A tail length**: see the `construct` tool for full poly-A
  assembly. `utr-design` focuses on the UTR × CDS coupling.