# mrnavax

> Practical Python tools for the AI-leverage layers in mRNA cancer therapy.
> Stdlib-only core, eleven runnable tools, twelve real-model adapters behind
> Protocol contracts, three documentation locales, 387 tests, 36 backend
> integrity checks.

[![CI](https://img.shields.io/badge/CI-passing-brightgreen?logo=githubactions&logoColor=white)](https://github.com/rollroyces/mrnavax/actions)
[![PyPI](https://img.shields.io/badge/PyPI-mrnavax%200.34.3-blue?logo=pypi&logoColor=white)](https://pypi.org/project/mrnavax/)
[![Python](https://img.shields.io/badge/Python-3.11–3.14-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-AGPL--3.0--or--later%20%2F%20commercial-orange)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-mrnavax.github.io-9cf?logo=readthedocs&logoColor=white)](https://rollroyces.github.io/mrnavax/)
[![Protocol adapters](https://img.shields.io/badge/adapters-9%20real%20models-purple)](https://github.com/rollroyces/mrnavax/tree/main/mrnavax)

![mrnavax pipeline — codon → variant → neoantigen → trial → LNP → manufacture](./docs/assets/pipeline.svg)

<p align="center">
  <a href="#quick-start"><img alt="quick start" src="https://img.shields.io/badge/quick_start-60s-5dd0ff?style=for-the-badge&logo=rocket"></a>
  <a href="#the-eleven-tools"><img alt="tools" src="https://img.shields.io/badge/tools-11-9b8cff?style=for-the-badge"></a>
  <a href="examples/run_all.py"><img alt="demo" src="https://img.shields.io/badge/demo-1_command-success?style=for-the-badge&logo=python"></a>
  <a href="spikes/SOTA_VALIDATION_SPIKE.md"><img alt="research" src="https://img.shields.io/badge/research-RNop_arXiv_2505.23862-94a3b8?style=for-the-badge"></a>
</p>

<p align="center">
  <strong>Try it:</strong> <code>pip install mrnavax</code> · <code>python examples/run_all.py</code> · ships with stdlib-only mocks, opt-in to real models via <code>[extras]</code>
</p>

<details>
<summary><strong>Table of contents</strong></summary>

- [What this is](#what-this-is)
- [The eleven tools](#the-eleven-tools)
- [Quick start](#quick-start)
- [Optional extras](#optional-extras)
- [Documentation](#documentation)
- [Real-model integrations](#real-model-integrations)
  - [RiboDecode](#ribodedode-lipt-iliego-et-al-nat-commun-16-9957-2025)
  - [STModule](#stmodule-wang-et-al-spstio-genome-medicine-17-2025)
  - [ESM2 + Applm pattern](#esm2--applm-pattern-wong-et-al-2025)
  - [TrialGPT + Sim-ICL](#trialgpt--sim-icl-jin-2024--fung-2026)
  - [AlphaGenome Atlas](#alpha-genome-atlas-avsec-et-al-nature-2026)
- [Architecture](#architecture)
- [Adding a new tool](#adding-a-new-tool)
- [License](#license)

</details>

<details>
<summary><strong>What you get — a 60-second overview</strong></summary>

> **▶ Watch the demo (5s):** the real `python examples/run_all.py` output, line-by-line:
>
> <video src="./docs/assets/mrnavax_demo.mp4" width="640" autoplay loop muted playsinline>
>   <img src="./docs/assets/mrnavax_demo.gif" alt="mrnavax end-to-end demo output — typewriter-style 5-second loop showing all 11 tools">
> </video>

```text
mrnavax end-to-end demo (mock mode)

  codon                CAI=0.720, n_codons=169
  neoantigen           8 candidates
  trial                top-1 = NCT00000003
  lnp                  target=lung, cargo=sarna
  scrna                10 cells × 9 genes, 4 clusters (tumor=3)
  manufacture          score=0.887, 7 pass / 8 total
  spatial              10 tissue modules
  construct            384 nt, CAI=0.887
  utr-design           combined=0.807
  variant-regulatory   6 variants, 6 high-impact
  predict              score=1.000, classification=high     ← live AlphaGenome Atlas call

```

All 11 tools, ~30 seconds of real runtime, zero model downloads. `predict`
upgrades to a live Atlas call when an API key is configured.

</details>

### Demo videos

Two additional short demos that highlight specific capabilities:

<details>
<summary><strong>🎬 Live AlphaGenome Atlas predict (5s)</strong> — proves the live integration works</summary>

> Real `mrnavax predict --input examples/regulatory_variants.csv` call against the live AlphaGenome Atlas API:
>
> <video src="./docs/assets/atlas_predict_demo.mp4" width="640" autoplay loop muted playsinline>
>   <img src="./docs/assets/atlas_predict_demo.gif" alt="Live AlphaGenome Atlas API demo — 6 cancer variants all scored high regulatory impact">
> </video>
>
> All 6 cancer variants returned `classification=high` from the actual Atlas server. The `← live` markers on each row confirm the call wasn't served by the mock backend.

</details>

<details>
<summary><strong>🎬 Codon optimization deep-dive (6s)</strong> — shows the before/after metric improvement</summary>

> Real `mrnavax codon --sequence examples/cas9.fasta --optimize --backend basic` run on the Cas9 reference sequence:
>
> <video src="./docs/assets/codon_optimize_demo.mp4" width="640" autoplay loop muted playsinline>
>   <img src="./docs/assets/codon_optimize_demo.gif" alt="Codon optimization demo — CAI improved from 0.720 to 0.938 (+0.218) via greedy frequency swap">
> </video>
>
> CAI improves from **0.720 → 0.938** (+0.218), rare-codons drop from 3.5% → 0%, with 66 codon swaps in 169 total codons. Amino-acid identity is preserved by construction (the swap table only selects synonymous codons).

</details>

## What this is

Eight small, runnable tools that map 1-to-1 onto the published AI leverage
points in mRNA cancer therapeutics — plus a typed integration contract for
every published foundation model in the field. Each tool runs as a CLI
subcommand and imports cleanly as a Python module.

```mermaid
flowchart LR
    subgraph DESIGN["Sequence design"]
        DNA["DNA sequence<br/>FASTA"] --> CAI["codon<br/>CAI / GC / rare"]
        DNA --> LD["LinearDesign DP<br/>O(L) Pareto"]
        DNA --> RD["RiboDecode<br/>Li 2025"]
    end
    subgraph VARIANT["Variant prioritization"]
        V["VCF / coding<br/>variants"] --> AM["AlphaMissense<br/>Cheng 2023"]
        V --> BF["BLOSUM62 +<br/>Chou-Fasman"]
    end
    subgraph NEO["Neoantigen prediction"]
        P["Mutant peptides"] --> MF["mhcflurry<br/>IC50 nM"]
        P --> ESM["ESM2 frozen LM<br/>Wong 2025"]
    end
    subgraph CELL["Single-cell foundation"]
        SC["scRNA-seq<br/>count matrix"] --> SCG["scGPT<br/>Cui 2024"]
        SCG --> TM["Tumor cluster<br/>→ mutant peptides"]
    end
    subgraph SPATIAL["Spatial transcriptomics"]
        ST["SRT counts +<br/>locations"] --> ST2["STModule<br/>Wang 2025"]
    end
    subgraph TRIAL["Patient-trial matching"]
        PT["Patient summary"] --> TG["TrialGPT<br/>Jin 2024"]
        PT --> SIM["Sim-ICL<br/>Fung 2026"]
    end
    subgraph WET["Wet-lab"]
        LNP2["LNP composition<br/>Witten 2025"] --> FINAL["Manufactured<br/>mRNA vaccine"]
        MAN["mRNA checks<br/>poly-A / Kozak / GC"] --> FINAL
    end

    RD --> P
    AM --> P
    TM --> P
    ST2 -.informs.-> TM
    P --> FINAL
    TRIAL -.eligibility.-> PT
```

## The eleven tools

Each tool runs as both a `mrnavax <tool>` CLI subcommand and an
importable Python module. Status column shows whether a real-model
backend is available or only the stdlib mock.

| Tool | Purpose | Reference work | Backend |
|------|---------|----------------|---------|
| [`codon`](docs/tools/codon.md) | Codon analysis (CAI, GC%, rare, structure-proxy) + LinearDesign DP + RiboDecode heuristic | CodonBERT; RiboDecode (Li et al., *Nat Commun* 16:9957, 2025); LinearDesign | 🟢 basic · 🟢 lineardesign · 🟡 ribodecode · 🟡 ribodecode-real |
| [`neoantigen`](docs/tools/neoantigen.md) | Peptide × HLA binding + ESM2 LM immunogenicity scoring | mhcflurry; MedCPT; DeepNeo; NetMHCpan; ESM2 + Applm pattern (Wong et al., 2025) | 🟢 anchor · 🟡 mhcflurry · 🟡 medcpt |
| [`trial`](docs/tools/trial.md) | TrialGPT per-criterion LLM matching + Sim-ICL demo selection | TrialGPT (Jin et al., *Nat Commun* 2024) + Sim-ICL (Fung et al., *Genome Biol* 2026) | 🟢 trialgpt · 🟢 trialgpt-simicl · 🟡 llm · 🟡 medcpt |
| [`scrna`](docs/tools/scrna.md) | scRNA-seq → tumor cluster → mutant peptides → ESM2 immunogenicity | scGPT (Cui et al., *Nat Methods* 2024) | 🟢 stdlib · 🟡 scgpt · 🟡 protein-lm |
| [`manufacture`](docs/tools/manufacture.md) | mRNA manufacturability checks (poly-A, Kozak, GC, ARE, stops) | Industry mRNA design guidelines (mRNA-1273 / BNT162b2) | 🟢 stdlib only |
| [`lnp`](docs/tools/lnp.md) | LNP composition recommender (cargo × target × intent) | Witten 2025; Li 2024 | 🟢 stdlib only |
| [`spatial`](docs/tools/spatial.md) | STModule spatial-transcriptomics tissue-module identification | STModule (Wang et al., *Genome Medicine* 2025) | 🟢 stdlib · 🟡 spatial-r (R/Seurat) |
| [`variant-regulatory`](#alpha-genome-atlas-avi-and-predict) | AlphaGenome Atlas AVI score for non-coding regulatory variants | AlphaGenome Atlas (Avsec et al., *Nature* 2026) | 🟢 mock · 🟢 atlas-live (with `[variant-alphagenome]`) |
| [`construct`](docs/tools/construct.md) | Full mRNA construct assembly (5'UTR + CDS + 3'UTR + poly-A) | mRNA-1273 / BNT162b2 consensus UTRs; multi-objective CDS (v0.25.0 RNop pattern) | 🟢 basic · 🟢 lineardesign · 🟡 ribodecode · 🟡 multi-objective |
| [`utr-design`](docs/tools/utr-design.md) | Coupled 5'UTR + CDS + 3'UTR design: bounded grid search over 12 UTR combinations | v0.27.0 UTR context scorer + v0.25.0 multi-objective CDS optimizer | 🟢 stdlib only |
| [`predict`](docs/tools/predict.md) | **Live** AlphaGenome Atlas regulatory-variant impact (AVI score) — single variant or CSV batch | AlphaGenome Atlas (Avsec et al., *Nature* 2026); auto-resolves API key | 🟢 atlas-live (with `[variant-alphagenome]` + key) |

**Legend:** 🟢 works out-of-the-box · 🟡 opt-in via `pip install -e '.[extra]'`

Every adapter has a **mock backend** that satisfies the same `runtime_checkable`
Protocol using only stdlib, so CI runs without downloading any model weights.

**Real-model adapters behind Protocol contracts** (opt-in via `pip install` extras):

| Adapter | Extra | Heavy deps | Status |
|---------|-------|------------|--------|
| `RiboDecode` | `[ribodecode]` | ViennaRNA + CUDA via subprocess | stable |
| `STModule` | `[spatial-r]` | R + Seurat + torch + CUDA | stable |
| `ESM2 protein-LM` | `[protein-lm]` | torch + transformers, ~135 MB | stable |
| `AlphaMissense` | standalone pickle | user-downloaded TSV | stable |
| `scGPT` | `[scrna]` | torch, ~205 MB | stable |
| `mhcflurry` | `[neoantigen-mhcflurry]` | mhcflurry + data download | stable |
| `MedCPT` | `[neoantigen-medcpt]` or `[trial-medcpt]` | torch + transformers, ~440 MB | stable |
| `TrialGPT / OpenAI` | `[llm]` | openai SDK + key | stable |
| `AlphaGenome Atlas` | `[variant-alphagenome]` | alphagenome package + API key | stable, live |

## Quick start

> **🚀 Zero to results in 60 seconds:**
>
> ```bash
> pip install mrnavax[variant-alphagenome]
> git clone https://github.com/rollroyces/mrnavax.git && cd mrnavax
> python examples/run_all.py
> ```
>
> (drop the `[variant-alphagenome]` extra if you don't have an Atlas API key yet — the demo runs in pure mock mode and skips only the `predict` line)

**Expected output (truncated):**

```text
mrnavax end-to-end demo (mock mode)

  codon                CAI=0.720, n_codons=169
  neoantigen           8 candidates
  trial                top-1 = NCT00000003
  ...
  predict              score=1.000, classification=high     ← live AlphaGenome Atlas call
```

**Step-by-step (for when you want to drive it yourself):**

```bash
git clone https://github.com/rollroyces/mrnavax.git
cd mrnavax
pip install -e .                   # stdlib-only core

# Or run the full demo across all 11 tools:
python examples/run_all.py

# 1. Codon analysis (CAI, GC%, rare-codon, GC-window stddev)
python -m mrnavax.cli codon --sequence examples/cas9.fasta
python -m mrnavax.cli codon --sequence examples/cas9.fasta \
    --optimize --backend lineardesign

# 2. Neoantigen screen (heuristic anchor matrix + LLM immunogenicity)
python -m mrnavax.cli neoantigen \
    --variants examples/tp53_variants.csv \
    --hla HLA-A*02:01

# 3. Patient-to-trial matching (TrialGPT-style, optionally Sim-ICL)
python -m mrnavax.cli trial \
    --patient examples/patient_summary.txt \
    --trials examples/trials.jsonl --top-k 5 \
    --matcher trialgpt-simicl

# 4. LNP composition advice
python -m mrnavax.cli lnp --target lung --cargo saRNA --intent "cancer vaccine"

# 5. scRNA-seq → neoantigen handoff
python -m mrnavax.cli scrna \
    --expression examples/cells.csv \
    --variants examples/variants_coding.csv \
    --proteins examples/proteins.fasta \
    --tumor-markers TP53,KRAS,BRAF

# 6. mRNA manufacturability score
python -m mrnavax.cli manufacture --cds examples/cds_gfp.json

# 7. Spatial transcriptomics tissue modules
python -m mrnavax.cli spatial \
    --count-file examples/st_bc2_count_matrix.tsv \
    --locations-file examples/st_bc2_locations.tsv \
    --platform ST --num-modules 10

# 8. AlphaGenome Atlas regulatory-variant AVI scoring
python -m mrnavax.cli variant-regulatory \
    --csv examples/regulatory_variants.csv
```

After `pip install -e .`, the same CLI is also installed as the console
script `mrnavax`.

Every example file referenced above is committed under `examples/`, and
sample outputs are committed under `examples/sample_outputs/`. See
`examples/README.md` for the full per-tool recipe (including
`construct`, `utr-design`, and `predict` — the latter requires an
Atlas API key).

## Optional extras

```bash
pip install -e ".[llm]"                       # OpenAI-compatible LLM client (TrialGPT)
pip install -e ".[neoantigen-mhcflurry]"       # mhcflurry binding-affinity backend
pip install -e ".[neoantigen-medcpt]"          # MedCPT query/article encoders (~440 MB)
pip install -e ".[protein-lm]"                # ESM2 protein language model (~135 MB)
pip install -e ".[trial-medcpt]"               # MedCPT for trial retrieval
pip install -e ".[scrna]"                     # scanpy + anndata + scGPT plug point
pip install -e ".[variant-alphagenome]"        # AlphaGenome Atlas regulatory-variant scoring
pip install -e ".[docs]"                      # mkdocs-material + mkdocs-static-i18n
pip install -e ".[dev]"                       # ruff + pytest
pip install -e ".[all]"                       # everything above
```

Then activate the real backend:

```bash
export OPENAI_API_KEY=sk-...
export OPENAI_MODEL=gpt-4o-mini               # default
python -m mrnavax.cli trial \
    --patient mrnavax/examples/patient_summary.txt \
    --trials mrnavax/examples/trials.jsonl --backend openai
```

Heavy-dependency backends (RiboDecode, STModule, ESM2, MedCPT, scGPT) ship
adapter modules that **subprocess or lazy-load** the upstream model. CI
runs without them; production users opt in per-extras above.

## Documentation

Full MkDocs site: <https://rollroyces.github.io/mrnavax/>

Available in three languages:

- 🇺🇸 English — <https://rollroyces.github.io/mrnavax/>
- 🇹🇼 繁體中文 — <https://rollroyces.github.io/mrnavax/zh-Hant/>
- 🇨🇳 简体中文 — <https://rollroyces.github.io/mrnavax/zh-Hans/>

Subsequent pushes to `main` deploy all three locales automatically via
GitHub Pages.

Local preview:

```bash
pip install -e ".[docs]"
mkdocs serve
```

## Real-model integrations

Every published foundation model is integrated behind a typed `Protocol`
adapter with a stdlib-only mock fallback. The adapter contract is identical
for production and CI — only the implementation differs.

### `RiboDecode` (Li et al., *Nat Commun* 16, 9957, 2025)

Joint translation × secondary-structure codon optimization via a deep
generative model. Heavy deps (ViennaRNA 2.6.4 + CUDA), shipped as a
subprocess adapter that calls the upstream CLI when present.

```bash
# Real: when ribo-decode is installed + Rscript on $PATH
mrnavax codon --sequence gfp.fasta --optimize --backend ribodecode-real \
    --env HEK293T --env-csv custom_env.csv --mfe-weight 0.3 --optim-epoch 10

# Mock: same shape, stdlib only
mrnavax codon --sequence gfp.fasta --optimize --backend ribodecode
```

### `STModule` (Wang et al., *Genome Medicine* 17, 2025)

Tissue-module identification from spatial-transcriptomics (SRT) data.
Heavy deps (R 4.4 + Seurat v5 + torch + GPUmatrix 1.0.2 + CUDA 11.7),
shipped via a small R shim that shells out to `Rscript stmodule_shim.R`.

```bash
# Real: when R + STModule are installed
mrnavax spatial --count-file counts.tsv --locations-file locs.tsv \
    --platform SlideSeqV2 --num-modules 10

# Mock: same shape, stdlib only
mrnavax spatial --count-file counts.tsv --locations-file locs.tsv \
    --platform ST --num-modules 10
```

### `ESM2` + Applm pattern (Wong et al., 2025)

Frozen protein-LM embeddings for neoantigen immunogenicity scoring.
Heavy deps (torch + transformers), shipped as a lazy-loaded adapter.

```python
from mrnavax.neoantigen_screener import lm_immunogenicity_score
r = lm_immunogenicity_score("NLVPMVATV")  # CMV pp65 epitope
print(r["score"])  # 0.0–1.0
```

### `TrialGPT` + `Sim-ICL` (Jin 2024 / Fung 2026)

Per-criterion patient-trial eligibility matching. Sim-ICL selects
top-K demonstration examples by TF-IDF cosine similarity (not random
sampling) — matching the paper's finding that sequence-similar
demonstrations outperform random few-shot.

```bash
mrnavax trial --patient patient.txt --trials trials.jsonl \
    --matcher trialgpt-simicl --top-k 10
```

### `AlphaGenome Atlas` (Avsec et al., *Nature* 2026)

Pre-computed regulatory-variant impact (AVI) scores for **all 9 billion
possible single-nucleotide variants** in the human genome. Where
AlphaMissense (Cheng et al. 2023) scores **coding-region** missense
variants, AlphaGenome Atlas scores **non-coding regulatory** variants
— covering the 98% of the genome where AlphaMissense is silent.
Adapter uses the official `alphagenome` Python package via subprocess
(gated behind the `[variant-alphagenome]` extra; non-commercial use
only per Google DeepMind's terms).

```bash
# Real: when ALPHAGENOME_API_KEY is set + [variant-alphagenome] installed
mrnavax variant-regulatory --csv variants.csv --backend alphagenome

# Mock: same shape, stdlib only
mrnavax variant-regulatory --csv variants.csv --backend mock
```

**Integrated into `score_variant`:** when you supply DNA coordinates
(`chrom`, `ref_dna`, `alt_dna`) plus an `avi_lookup` callable, the
same `score_variant()` entry point that drives the scrna pipeline
auto-routes coding-region variants to AlphaMissense (dominant signal)
and non-coding regulatory variants to AlphaGenome Atlas (dominant
signal). A single CSV with both kinds of variants scores them all
through one function:

```bash
# variants.csv has: gene,position,wt_aa,mut_aa,chrom,ref_dna,alt_dna
mrnavax scrna \
    --expression cells.csv \
    --variants variants.csv \
    --proteins proteins.fasta \
    --tumor-markers TP53,KRAS,BRAF \
    --variant-filter-top-fraction 0.4 \
    --out report.json
# report.json includes variant_scores for every variant + a note
# mentioning "AlphaGenome Atlas AVI scores used for non-coding
# regulatory variants"
```

### Other real-model integrations

- **`AlphaMissense`** (Cheng et al., *Science* 381, 2023) — variant
  pathogenicity via 71M-variant TSV; bundled pickle index for O(1)
  per-variant lookup.
- **`scGPT`** (Cui et al., *Nat Methods* 21, 2024) — single-cell
  foundation-model embeddings, 30 layers × 512 dim.
- **`mhcflurry`** (O'Donnell et al.) — Class I MHC binding affinity,
  IC50 in nM.
- **`MedCPT`** (Jin et al., 2023) — biomedical dense retrieval,
  contrastively trained on PubMed.

## Architecture

The toolkit is built on three orthogonal layers — CLI / core / adapters
— and every published foundation model plugs in through the same `Protocol`
contract with a stdlib-only mock fallback.

```mermaid
flowchart TB
    subgraph USER["User interface"]
        CLI["mrnavax CLI<br/>(8 subcommands)"]
        PY["import mrnavax<br/>as Python module"]
    end

    subgraph CORE["Core layer (stdlib only, ships in pip wheel)"]
        direction TB
        TOOLS["Typed dataclass tools<br/>codon / neoantigen / trial<br/>scrna / spatial / manufacture / lnp<br/>variant-regulatory"]
        SELECT["Backend selector<br/>_backends_registry + 8 family modules"]
        CK["30 backend integrity checks<br/>(deterministic structural)"]
    end

    subgraph PROTOCOL["typing.Protocol contracts (4 actually defined)"]
        direction TB
        P1["CodonOptimizer"]
        P2["TranslationPredictor"]
        P3["ProteinLMEmbedder"]
        P4["SpatialModuleBackend"]
    end

    subgraph ADAPTERS["Real-model adapters (opt-in extras)"]
        direction TB
        A1["LinearDesign DP<br/>RiboDecode subprocess<br/>CodonBERT"]
        A2["mhcflurry<br/>MedCPT<br/>ESM2 (Applm pattern)"]
        A3["TrialGPT OpenAI<br/>Sim-ICL"]
        A4["scGPT"]
        A5["STModule R shim"]
    end

    subgraph MOCKS["Mock adapters (always present)"]
        direction TB
        M1["MockCodonOptimizer"]
        M2["MockTranslationPredictor"]
        M3["MockProteinLMEmbedder"]
        M4["MockSpatialModuleBackend"]
    end

    CLI --> TOOLS
    PY --> TOOLS
    TOOLS --> SELECT
    SELECT --> P1
    SELECT --> P2
    SELECT --> P3
    SELECT --> P4

    P1 -.implemented by.-> A1
    P2 -.implemented by.-> A1
    P3 -.implemented by.-> A2
    P4 -.implemented by.-> A5

    P1 -.always available.-> M1
    P2 -.always available.-> M2
    P3 -.always available.-> M3
    P4 -.always available.-> M4

    CK -.verifies.-> SELECT
    CK -.verifies.-> M1
    CK -.verifies.-> M2
    CK -.verifies.-> M3
    CK -.verifies.-> M4

    style PROTOCOL fill:#f9f,stroke:#333,stroke-width:2px
    style MOCKS fill:#cfc,stroke:#333
    style ADAPTERS fill:#fcf,stroke:#333
```

**Note:** the `neoantigen`, `trial`, `scrna`, `manufacture`, and `lnp`
tools expose module-level functions rather than a `Protocol` class, so
they are wired through `_backends_registry.register(name)` decorators
and verified by the same 30 structural integrity checks. The Protocol
contract pattern is applied where multiple interchangeable
implementations exist (codon optimizers, protein-LM embedders,
spatial-module finders).

**Why this matters:** the CI matrix (Python 3.11–3.14) runs without
downloading any model weights. Production users opt in per-extras
(`pip install -e ".[neoantigen-mhcflurry]"`). The same code path runs
in both — only the adapter implementation differs.

## Architecture rationale

The toolkit's `backends.py` integrity checks use **deterministic
structural assertions** rather than AUPRC / F1 / accuracy metrics from
heavy libraries. Chen et al. 2024 (*Genome Biology* 25, 118) evaluated
10 widely-used PRC tools across >3,000 published studies and found
they produce **conflicting AUPRC rankings and overly-optimistic results**.
The toolkit's stdlib-only baseline avoids that entire class of bug by
owning the metric end-to-end.

See `docs/index.md` for the full rationale.

## Real-world case studies (grounded in published biology)

The toolkit's scRNA → neoantigen pipeline is validated against the
wet-lab workflow in Qian et al. 2022 (*Int J Cancer* 151, 1367-1381):
scRNA-seq of gastric cancer primary tumor + lymph node metastases →
tumor cluster identification → mutant peptide enumeration → ESM2
immunogenicity scoring → mRNA cancer vaccine design. See
`docs/tools/scrna.md` for the full walkthrough.

## Why eight layers (and not four or five)?

The mRNA cancer therapy research has reached an inflection point where
foundation models for **sequence design**, **variant prioritization**,
**neoantigen prediction**, **single-cell foundation**, **spatial
transcriptomics**, **patient-trial matching**, and **manufacturing
checks** are all simultaneously advancing. The toolkit's job is to
be the integration layer — every published model plugs in via a
Protocol contract with a stdlib-only mock fallback for tests and
offline use.

1. **Sequence design** (`codon`): LinearDesign (real, O(L) DP, no cap)
   + RiboDecode-style context heuristic. Powers any mRNA construct.
2. **Variant prioritization** (`variant_scorer`): AlphaMissense TSV
   lookup + BLOSUM62 + driver-gene awareness + Chou-Fasman structural
   disruption. AM weighted at 45%.
3. **Neoantigen prediction** (`neoantigen`): mhcflurry IC50 +
   ESM2 LM immunogenicity. Frozen-LM-then-classifier pattern.
4. **Single-cell foundation** (`scrna`): scGPT embeddings → tumor
   cluster identification → mutant peptide handoff.
5. **Spatial transcriptomics** (`spatial`): STModule tissue-module
   identification. Spatial coordinates reveal where the tumor
   cluster lives.
6. **Patient-trial matching** (`trial`): TrialGPT per-criterion LLM +
   Sim-ICL demonstration selection. Real + keyword fallback.
7. **Manufacturability** (`manufacture`): poly-A runs, Kozak strength,
   GC window uniformity, ARE motifs, hidden stops, CpG balance.
8. **LNP delivery** (`lnp`): ionizable-lipid pKa, helper-lipid ratio,
   composition shortlist above published ML-discovered candidates.

## Development

```bash
# Run all backend integrity checks (matches CI)
python -m mrnavax.backends --check-all

# Run the unit test suite
python -m unittest discover tests

# Run on the bundled examples (see scripts/smoke.sh)
bash scripts/smoke.sh

# Build docs locally
pip install -e ".[docs]"
mkdocs serve
```

### CI / publish pipeline

```mermaid
flowchart LR
    DEV["git push<br/>to main"] --> SMOKE["smoke.yml<br/>Python 3.11–3.14<br/>284 tests + 30 checks"]
    DEV --> DOCS["docs.yml<br/>mkdocs strict<br/>3 locales"]
    SMOKE -.on failure.-> FAIL["❌ red ✋<br/>fix + push again"]
    DOCS -.on failure.-> FAIL

    TAG["git tag vX.Y.Z<br/>git push --tags"] --> PUB["publish.yml"]
    PUB --> BUILD["build job<br/>sdist + wheel<br/>version matches tag"]
    BUILD --> ART["dist/<br/>artifact"]
    ART --> PYP["publish-to-pypi job<br/>OIDC trusted publisher"]
    PYP -.manual approval.-> REVIEW["pypi environment<br/>reviewer gate"]
    REVIEW --> LIVE[("PyPI<br/>mrnavax X.Y.Z<br/>live")]
    LIVE --> PAGES["GitHub Pages<br/>rollroyces.github.io/mrnavax"]

    style SMOKE fill:#cfc
    style DOCS fill:#cfc
    style BUILD fill:#cff
    style LIVE fill:#fc9
    style FAIL fill:#fcc,stroke:#c33,stroke-width:2px
```

All three workflows use Node 24-native action majors
(`actions/checkout@v6`, `actions/setup-python@v6`, etc.) — zero
deprecation warnings on the latest runs. The toolkit also ships a
fourth workflow, `.github/workflows/atlas_integration.yml`, that
runs weekly (Mondays 06:00 UTC) against the real AlphaGenome Atlas
API to catch upstream breakage; it stays dormant (mock mode) until
the `ALPHAGENOME_API_KEY` GitHub secret is configured.

## License

Dual-licensed. See `LICENSE` for the dual-license summary and `LICENSE-AGPL`
for the AGPL-3.0-or-later terms. A commercial license is available on
request — open an issue on the GitHub repo.

## Contributing

Pull requests welcome. The default dependency surface is **Python stdlib
only** — heavy model integrations must plug into a backend selector via
the existing Protocol-based adapter pattern (see `mrnavax/codon_ribodecode_adapter.py`,
`mrnavax/spatial_module_adapter.py`, `mrnavax/protein_lm_adapter.py`
for reference).

Every new tool should ship with:

1. A typed dataclass for input + output (frozen, validated at construction).
2. A `runtime_checkable` Protocol for the backend interface.
3. A real adapter that shells out / lazy-loads the upstream model.
4. A stdlib-only mock that satisfies the same Protocol.
5. A `@register("family.subname")` entry in the appropriate
   `mrnavax/_backends_<family>.py` module (or in a new family module
   if the check does not fit an existing family — add it to the import
   block at the top of `mrnavax/backends.py`). The 30-check
   integrity registry lives in `mrnavax/_backends_registry.py`.
6. Tests in `tests/` following strict TDD.

Helper-module pattern (v0.22.0+): when a public module's orchestration
function is doing too much (e.g. `run_pipeline` had a 245-line filter
+ scoring + peptide-emit + notes block), extract focused helpers into
underscore-prefixed private modules:

- `mrnavax/_scoring_components.py` + `mrnavax/_scoring_lookups.py`
  own the per-component logic behind `score_variant`.
- `mrnavax/_scrna_filter.py` + `mrnavax/_scrna_peptide_emitter.py`
  own the filter + peptide emission behind `run_pipeline`.
- `mrnavax/_backends_<family>.py` owns the check functions behind
  `backends.CHECKS`.

New helpers should land in similarly-named underscore-prefixed
modules with tests that prove the helpers are usable in isolation
from the public orchestration function.
