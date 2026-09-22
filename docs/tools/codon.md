# `codon` — sequence analysis & optimization

Computes the canonical codon-usage features that every modern mRNA design model
(CodonBERT, RiboDecode, LinearDesign, mRNABERT) consumes. Ships four
optimization backends:

| Backend | Algorithm | Setup |
|---|---|---|
| `basic` (default) | Greedy per-codon frequency swap, GC% band 45–60% | — |
| `ribodecode` | Context-aware hill-climb matching RiboDecode's published algorithm (Li et al., *Nat Commun* 16, 9957, 2025) | — (stdlib) |
| `lineardesign` | Joint translation × MFE DP, no length cap | — (stdlib) |
| `ribodecode-real` | Subprocess to upstream `ribo-decode` CLI | `pip install ribodecode-1.3.0-py3-none-any.whl` + Rscript on `$PATH` |
| `multi-objective` | SOTA knowledge-infused loss decomposition (CAI + GC + CpG + rare-codon-run + structure-proxy) with per-component breakdown. Stdlib-only. See [Multi-objective optimization](#multi-objective-knowledge-infused-loss). | — (stdlib) |

## Usage

```bash
# Codon analysis (no optimization)
mrnavax codon --sequence mrnavax/examples/cas9.fasta

# Greedy optimization
mrnavax codon --sequence mrnavax/examples/cas9.fasta --optimize

# LinearDesign (joint translation + mRNA structure, O(L) DP)
mrnavax codon --sequence mrnavax/examples/cas9.fasta \
    --optimize --backend lineardesign

# RiboDecode heuristic (context-aware hill-climb)
mrnavax codon --sequence mrnavax/examples/cas9.fasta \
    --optimize --backend ribodecode

# Real RiboDecode CLI (requires upstream R package)
mrnavax codon --sequence mrnavax/examples/cas9.fasta \
    --optimize --backend ribodecode-real \
    --env HEK293T --mfe-weight 0.3 --optim-epoch 10

# LinearDesign works on full-length CDS (Cas9 4.1 kb in ~9 s)
mrnavax codon --sequence mrnavax/examples/cas9.fasta \
    --optimize --backend lineardesign
```

## Output schema

```json
{
  "n_codons": 210,
  "cai": 0.6954,                        // Codon Adaptation Index (0–1)
  "gc_percent": 37.3,                   // overall GC%
  "rare_codon_fraction": 0.0429,        // fraction of codons below 0.10 freq
  "cpg_obs_exp": 1.1331,                // CpG O/E ratio — proxy for innate immunity
  "most_common_codons": [["GAT", 12], ...],
  "rare_codons": ["CTA", "TTA"],
  "gc_window_stddev": 4.75               // rolling GC stddev (translation speed proxy)
}
```

## `--optimize` backends

### `basic` (default)

Greedy synonymous-codon swap that maximizes per-codon usage frequency
while keeping the GC% in the 45–60% band.

Expected CAI improvement on a bacterial gene expressed in human cells: ~0.2
absolute (e.g. 0.7 → 0.93 for the Cas9 example).

This is the **classical baseline** that any modern codon model should beat.

### `ribodecode` (heuristic)

Context-aware hill-climb: at each position, choose the synonym with the
best ribosome-profiling-weighted frequency. Stdlib-only — no upstream
dependency. Suitable for CI and offline runs.

### `lineardesign` (DP)

Joint translation × MFE dynamic programming. O(L) per step (state space
bounded by |syn|^W for window size W). **No length cap** — works on
full-length mRNA constructs up to and beyond 4 kb (Cas9 in ~9 s on
Apple Silicon). Protein sequence preserved exactly.

### `ribodecode-real`

Subprocess to the upstream `ribo-decode` R package. Requires:

```bash
# Install per the upstream vignette (R + Seurat not actually needed
# for the optimizer, just the .whl)
pip install TranslationModel-1.1.0-py3-none-any.whl
pip install ribodecode-1.3.0-py3-none-any.whl

# Optional: custom cellular environment via RPKM CSV
mrnavax codon --sequence gfp.fasta --optimize --backend ribodecode-real \
    --env HEK293T \
    --csv env_hek293t.csv \
    --mfe-weight 0.3 --optim-epoch 10
```

The adapter shells out to the upstream CLI with `--env {HEK293T,A549,HeLa,custom}`,
`--mfe-weight` (0=translation-only, 1=MFE-only), `--optim-epoch N`. When the
upstream binary isn't on `$PATH`, the adapter raises
`RiboDecodeNotInstalled` with a clear remediation message including the
Google Drive `.whl` download links.

## Multi-objective optimization (knowledge-infused loss)

The `multi-objective` backend implements the SOTA pattern from arXiv:2505.23862
(RNop, *mRNA Design and Optimization with Deep Knowledge-Infused Approach*,
Aug 2026), where each biologically-motivated objective is decomposed
into its own score so the user can see **which** component drove each
codon swap. This addresses the "impossible triangle" of mRNA
optimization (fidelity, multi-objective, efficiency) the way the
top-tier projects (RNop / LinearDesign / RiboDecode) all converge on:
fidelity is preserved (only synonymous swaps), multi-objective is
explicit (a weighted loss vector), and efficiency is stdlib-Python.

Five normalized components:

| Component | Range | What it measures |
|---|---|---|
| `cai` | (0, 1] | Codon Adaptation Index — translation efficiency |
| `gc_score` | [0, 1] | GC% inside target band — transcript stability proxy |
| `cpg_score` | [0, 1] | CpG Obs/Exp avoidance — innate-immune proxy |
| `rare_run_penalty` | [0, 1] | Longest rare-codon run — ribosome-stalling proxy |
| `structure_proxy` | [0, 1] | GC-window stddev uniformity — local-structure proxy |

The overall score is a weighted sum of the components (with
`rare_run_penalty` subtracted as a penalty) clamped to `[0, 1]`. The
default weights follow the canonical SOTA pattern: CAI is the
dominant signal (~40%), structure-proxy is next (~25%), then GC and
CpG, with rare-codon-run counted as a soft penalty.

```bash
mrnavax codon --sequence gfp.fasta \
    --optimize --backend multi-objective
```

Output includes `before` and `after` per-component breakdowns plus
`overall_before`, `overall_after`, `improvement`, and `n_changes`,
making the contribution of each biological axis fully auditable.

**Customize the loss vector** — `MultiObjectiveConfig` exposes all 5
weights (and the GC band, rare threshold, structure window). Wire your
own config in Python:

```python
from mrnavax.codon_multi_objective import (
    MultiObjectiveConfig, multi_objective_optimize,
)

cfg = MultiObjectiveConfig(
    cai_weight=0.30,        # less CAI, more on structure
    gc_weight=0.20,
    cpg_weight=0.15,
    rare_run_weight=0.25,   # emphasize avoiding rare-codon runs
    structure_weight=0.10,
    target_gc_min=50.0,     # tighter GC band
    target_gc_max=58.0,
)
result = multi_objective_optimize(cds, cfg)
print(result.improvement, result.n_changes)
```

## RNop reference

Gong, Z., Jiang, Z., Gao, W., Wang, Y., Cai, Z., Deng, Z., Ma, L.
(2026). *mRNA Design and Optimization with Deep Knowledge-Infused
Approach.* arXiv:2505.23862v2. The "impossible triangle" framing
(fidelity, multi-objective, efficiency) and the knowledge-infused loss
decomposition are the conceptual basis for the `multi-objective`
backend's component breakdown.

## RiboDecode references

- Li, Y., Wang, F., Yang, J., et al. (2025). *Deep generative optimization
  of mRNA codon sequences for enhanced mRNA translation and therapeutic
  efficacy.* Nat Commun 16, 9957. DOI: 10.1038/s41467-025-64894-x
- GitHub: [wangfanfff/RiboDecode](https://github.com/wangfanfff/RiboDecode)

## LinearDesign reference

The DP algorithm follows the joint translation × structure DP from
LinearDesign (Zhang et al., 2023) — O(L) per step via suffix-state
pruning, parent-pointer backtrack, L2-normalized translation scores.
Protein preservation invariant verified across all test CDSs.
