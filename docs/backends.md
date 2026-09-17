# Backends

Seven tools, each with a typed `Protocol` adapter and a stdlib mock
fallback. Heavy deps are opt-in per `pip install` extras — CI runs
without them.

## Backend matrix

| Backend | Tools | What it does | Setup |
|---|---|---|---|
| `mock` | all | Stdlib stub. Always available. | — |
| `openai` | `neoantigen`, `trial` | Calls OpenAI Chat Completions (or any OpenAI-compatible endpoint) | `OPENAI_API_KEY=...` |
| `mhcflurry` | `neoantigen` | Pan-allele MHC-I binding affinity | `pip install -e ".[neoantigen-mhcflurry]"` |
| `MedCPT` | `neoantigen`, `trial` | Dense biomedical retrieval (PubMed contrastive) | `pip install -e ".[neoantigen-medcpt]"` or `[trial-medcpt]` |
| `scGPT` | `scrna` | Single-cell foundation-model embeddings | `pip install -e ".[scrna]"` |
| `AlphaMissense` | `variant_scorer` | Pathogenicity via 71M-variant TSV | Standalone TSV (CC BY-NC-SA) |
| `AlphaGenome Atlas` | `variant-regulatory`, `variant_scorer`, `scrna` | Regulatory-variant impact (AVI) for all 9B possible SNVs | `pip install -e \".[variant-alphagenome]"` + `ALPHAGENOME_API_KEY` |
| `ESM2` | `neoantigen` | Frozen protein-LM embeddings for immunogenicity | `pip install -e ".[protein-lm]"` |
| `RiboDecode` (real) | `codon` | Joint translation × MFE codon optimization | `pip install ribodecode-1.3.0-py3-none-any.whl` |
| `RiboDecode` (heuristic) | `codon` | Stdlib RiboDecode-style hill-climb | — |
| `LinearDesign` | `codon` | Joint translation × MFE DP | — (stdlib) |
| `STModule` (real) | `spatial` | Tissue-module identification from SRT data | `pip install STModule` + Rscript on `$PATH` |
| `STModule` (mock) | `spatial` | Stdlib stub with per-platform gene universes | — |
| `TrialGPT` | `trial` | Per-criterion LLM eligibility matching | (uses `openai` backend) |
| `Sim-ICL` | `trial` | Top-K demo selection by TF-IDF cosine | — (stdlib) |

## How Protocol adapters work

Each real-model backend is wrapped in a `runtime_checkable` Protocol:

```python
@runtime_checkable
class TranslationPredictor(Protocol):
    def predict(self, cds: str, env: str = "HEK293T",
                custom_env_csv: Path | None = None) -> TranslationPrediction: ...
```

The toolkit ships both a real adapter (e.g. `TranslationModelCLIAdapter`
subprocess to `pred-translation`) and a mock (`MockTranslationPredictor`
using CAI-derived score). Consumers see the same Protocol; the backend
selector picks real-or-mock based on `$PATH` + installed deps.

## Backend selectors

| Tool | Selector | Real-or-mock dispatch |
|---|---|---|
| `codon` | (CLI flag `--backend`) | `ribodecode-real` if `ribo-decode` on `$PATH` else `ribodecode` (heuristic) |
| `neoantigen` | `select_translation_predictor` / `select_codon_optimizer` | mhcflurry if installed else OpenAI if key set else mock |
| `trial` | `select_spatial_module_backend` / `select_protein_lm_embedder` | OpenAI if key set else mock |
| `scrna` | `embed_with_foundation_model` | scGPT if installed else identity (no-op fallback) |
| `manufacture` | — | All stdlib (no LLM) |
| `lnp` | — | All stdlib (no LLM) |
| `spatial` | `select_spatial_module_backend` | `Rscript` on `$PATH` else mock |

## Auto-detection

When `--backend auto` (the default), the tool picks the strongest
available backend in this order:

1. **Heavy upstream binary** (e.g. `Rscript` for STModule, `pred-translation`
   for RiboDecode) — best accuracy, requires user setup.
2. **Installed Python package** (mhcflurry, transformers, OpenAI) — good
   accuracy, opt-in.
3. **Mock** — always works, deterministic, ~0 ms, stdlib-only.

This makes local dev painless and lets production deployments override
via env var or CLI flag.

## Recorded-response fixtures

For backends that talk to external APIs whose schemas could break
between releases, we bundle a JSON fixture + a backend check that
parses it. The fixture locks the parser shape so CI catches schema
changes at PR time, not at user runtime.

| Fixture | Backend check | When to refresh |
|---|---|---|
| `tests/fixtures/alphagenome_atlas_sample.json` | `variant.alphagenome_atlas_fixture` | When Google's AlphaGenome Atlas response shape changes |

To refresh a fixture: capture a real API response (with appropriate
auth), save it under the existing filename, run
`python -m mrnavax.backends --check-all`, and commit. The synthetic
fixtures ship today are sufficient to lock the parser shape; a real
captured response can replace them whenever convenient.
