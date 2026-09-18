# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.23.0] - 2026-09-18

### Changed

- **Internal refactor**: ``mrnavax/sc_rna_pipeline.py`` is now 630
  lines (was 739) — ``run_pipeline`` is split into two private
  helper modules:
    - ``mrnavax/_scrna_filter.py`` — ``filter_variants_with_lookups``,
      ``build_pipeline_notes``, ``build_tumor_marker_note``. Owns
      the lookup-wiring + filter call + per-variant scoring fallback
      + user-facing notes construction.
    - ``mrnavax/_scrna_peptide_emitter.py`` —
      ``emit_tumor_peptides``. Walks the filtered variants and
      emits mutant-peptide candidates for those that are both
      tumor-cluster-expressed and have a known protein sequence.
- **run_pipeline body reduced** to ~75 lines of orchestration:
  load inputs → filter / per-variant scoring → cluster → identify
  tumor cluster → emit peptides → build note → return.
- **No public API change** — ``run_pipeline`` + ``PipelineReport``
  + every CLI subcommand + every existing test still passes
  unchanged.

### New private helpers (underscore-prefixed)

- ``mrnavax/_scrna_filter.py``:
  ``filter_variants_with_lookups(variants, *, top_fraction,
  min_score, proteins, uniprot_ids)`` returns
  ``(kept, scores, am_active)``. Wires the optional AlphaMissense /
  AVI / PhyloP lookups transparently and either filters the variants
  or falls through to per-variant scoring when filtering is
  disabled but DNA coords are present.
  ``build_pipeline_notes(*, am_active, filter_active,
  had_dna_coords)`` returns the user-facing ``note`` string.
  ``build_tumor_marker_note(marker_idx, notes)`` appends the
  tumor-marker warning.
- ``mrnavax/_scrna_peptide_emitter.py``:
  ``emit_tumor_peptides(variants, proteins, tumor_cluster,
  tumor_expressed_genes, cluster_marker_score, peptide_lengths)``
  returns the list of ``TumorPeptide`` objects.

### New tests

- ``tests/test_scrna_refactor.py`` (14 tests):
    - ``TestEmitTumorPeptides`` — empty variants / missing protein /
      not tumor-expressed / real tumor-expressed case
    - ``TestBuildPipelineNotes`` — AM active / AM fallback / DNA
      coords on/off / all-off empty
    - ``TestBuildTumorMarkerNote`` — no-markers appends warning /
      markers present leaves notes unchanged
    - ``TestFilterVariantsWithLookupsNoFilter`` — passthrough
    - ``TestFilterVariantsWithLookupsFilterActive`` — kept set
      reduces correctly + scores dict aligned

### Why this matters

The variant-filter / scoring loop / peptide-emitter / notes
construction was previously crammed into one ~250-line block inside
``run_pipeline``. Splitting it into focused helpers:
  * The filter + scoring logic is now unit-testable in isolation
    — the new tests prove ``filter_variants_with_lookups`` works
    end-to-end without spinning up the full pipeline (no expression
    matrix needed).
  * The peptide-emitter is reusable from any future code path that
    wants to enumerate mutant peptides for a list of variants
    (e.g. a batch-pipeline mode for therapeutic prioritization).
  * The notes builder is independently testable — 6 of the new
    tests pin its behavior across all combinations of AM / AVI /
    DNA-coords / tumor-marker signals.

## [0.22.0] - 2026-09-18

### Changed

- **Internal refactor**: ``mrnavax/variant_scorer.py`` is now 891
  lines (was 976) — the four-sign scoring composition logic has
  been split into two private helper modules:
    - ``mrnavax/_scoring_components.py`` — ``LocalComponents``
      dataclass + ``compute_local_components`` (BLOSUM62 +
      driver-gene + hydrophobicity + structural disruption) and
      ``local_rationale``.
    - ``mrnavax/_scoring_lookups.py`` — ``compute_am_component``,
      ``compute_avi_component``, ``compute_conservation_component``
      (the silent-failure pattern for upstream-model lookups) plus
      their respective dataclasses.
- **score_variant body reduced** to 116 lines of orchestration:
  validate → compute local → compute AM/AVI/PhyloP → combine →
  rationale → return. ``_combine_components`` (83 lines) holds the
  weighted-combination logic.
- **No public API change** — ``score_variant``, ``filter_variants``,
  ``VariantScore``, and every backend check / test still pass
  unchanged. The 249-test suite + 30-check backend suite all green.

### Why this matters

The four-signal composition (BLOSUM62 + driver + structural + AM →
AVI → PhyloP) was previously crammed into one ~290-line function.
Splitting it into focused helpers makes:
  * the per-component logic unit-testable in isolation (run just
    ``compute_local_components`` on a hand-built input and assert
    on the dataclass fields)
  * the upstream-lookup pattern (silent-failure, range-check,
    weight-assignment) reusable for any future 4th signal
  * the weighted-combination logic readable at a glance instead
    of buried in a 290-line block

### New public helpers (private)

- ``mrnavax/_scoring_components.py``:
  ``LocalComponents`` dataclass, ``compute_local_components(gene,
  wt_aa, mut_aa, position, protein_length, protein_sequence,
  driver_genes)`` returns ``LocalComponents``;
  ``local_rationale(local, gene, wt_aa, mut_aa)`` returns list of
  rationale string fragments.
- ``mrnavax/_scoring_lookups.py``:
  ``AlphaMissenseComponent``, ``AVIComponent``,
  ``ConservationComponent`` dataclasses;
  ``compute_am_component(...)``, ``compute_avi_component(...)``,
  ``compute_conservation_component(...)`` — each returns the
  component dataclass or ``None`` (silent failure).

These are underscore-prefixed (private); downstream code should call
``mrnavax.variant_scorer.score_variant``.

## [0.21.0] - 2026-09-17

### Added

- **Scheduled live AlphaGenome Atlas integration harness**
  (`mrnavax/live_atlas_integration.py`) that calls the real Atlas
  API and compares against a baseline with ±5% tolerance. Three modes:
    - `mock` — uses the stdlib mock; always passes; no network.
      Used in unit tests + offline CI.
    - `record` — calls the real Atlas; writes the response to a
      JSON file. Triggered manually when refreshing fixtures.
    - `regression` — calls the real Atlas; compares against a
      baseline score with ±5% tolerance. Triggered weekly via the
      new GitHub Actions cron.
- **New GitHub Actions workflow `.github/workflows/atlas_integration.yml`**
  scheduled for every Monday at 06:00 UTC. Uses Node 24-native
  action majors (consistent with the rest of the project's
  workflows). Manual `workflow_dispatch` for record mode. Skips the
  live call gracefully when `ALPHAGENOME_API_KEY` is unset
  (the harness still runs in mock mode as a smoke test).
- **CLI**: `python -m mrnavax.live_atlas_integration --mode {mock,
  record, regression} [--baseline 0.72 --tolerance 0.05 --output PATH]`.
  Exit codes: 0 on pass, 1 on regression failure, 2 on missing API
  key (for record / regression modes).
- **10 new unit tests** in `tests/test_live_atlas_integration.py`
  covering mock / record / regression modes, baseline comparison
  with tolerance, and CLI integration.
- **`docs/contributing.md` (3 locales)**: new section
  "Refreshing the recorded Atlas fixture" documenting the workflow
  + manual record procedure.
- Total tests: **249** (was 239). Total backend checks: **30**
  (unchanged — the live integration is a workflow, not a CI check).

### Why this matters

Closes the last remaining item on the roadmap. The recorded fixture
in v0.18.0 catches upstream Atlas API schema changes with a 7-day
lag (since it's only refreshed when the bundle is updated). The
weekly cron catches the same changes in real time. Manual record
mode lets contributors refresh the fixture themselves in 30
seconds with `python -m mrnavax.live_atlas_integration --mode
record`.

### Setup

To enable the live integration:

1. Get an API key at <https://deepmind.google.com/science/alphagenome>
   (non-commercial preview; free for academic users).
2. Add it as a GitHub Actions secret:
   <https://github.com/rollroyces/mrnavax/settings/secrets/actions/new>
   - Name: `ALPHAGENOME_API_KEY`
   - Value: your key

The workflow is dormant without the secret (mock mode runs as a
smoke check). Add the secret + wait for the next Monday 06:00 UTC,
or trigger manually via `Actions → Atlas live integration → Run
workflow`.

## [0.20.0] - 2026-09-17

### Added

- **PhyloP46way evolutionary-conservation integration** as the 4th
  coding-region scoring signal in `score_variant`. Closes the last
  remaining gap in coding-region variant prioritization. Components:
    - `mrnavax/conservation.py` NEW module with
      `ConservationLookup` Protocol, `MockPhyloPLookup` (stdlib-only,
      deterministic SHA-256 hash), and `PhyloPRestAdapter` (real
      UCSC REST API; uses stdlib urllib so no heavy deps).
    - `select_conservation_lookup()` backend selector.
    - `score_variant` accepts `conservation_lookup` + `pos` kwargs;
      PhyloP is recorded in `components['phylop46way_score']` and
      appears in the rationale string "PhyloP46way conservation=±X.XXX".
- **Composition with AM + AVI**: the 4-signal pipeline is now
  BLOSUM62 + driver + structural + AM (dominant for coding) +
  AVI (dominant for non-coding regulatory) + PhyloP (conservation
  boost). For coding-region variants with high conservation,
  PhyloP adds a +0.10 × score bonus on top of the AM-dominant
  combination.
- **`filter_variants` + `run_pipeline` auto-wire the conservation
  lookup** when variants carry DNA coordinates. Mock by default;
  real UCSC adapter selected via `select_conservation_lookup()`.
- **10 new unit tests** in `tests/test_conservation.py` covering
  parameter acceptance, dominant-signal routing, silent failure
  modes (None return / exception / out-of-range), and composition
  with AM + AVI.
- **New backend check `conservation.phylop46way`** — exercises the
  full integration end-to-end. Total backend checks: **30/30**
  (was 29/29).
- `docs/backends.md` (3 locales) — new `PhyloP46way` row in the
  backend matrix.
- `docs/case-studies/variant-prioritization.md` — refined TL;DR
  to mention the 4-signal composition and what each layer adds.
- Total tests: **239** (was 229). Total backend checks: **30**
  (was 29).

### Why this matters

Closes the last remaining gap in coding-region variant
prioritization. The previous releases had BLOSUM62 + driver-gene
+ structural disruption + AlphaMissense — but no per-position
evolutionary-conservation signal. With PhyloP46way, the scorer
now distinguishes "this is a slow-evolving site" (conserved → more
likely pathogenic) from "this is a fast-evolving site" (likely
neutral). Composes cleanly with AlphaMissense: both signals
point in the same direction for true driver variants and
disagree for false positives.

### Reference

Pollard KS, Hubisz MJ, Rosenbloom KR, Siepel A. "Detection of
nonneutral substitution rates on mammalian phylogenies."
*Genome Research* 20, 110–121 (2010). PhyloP46way data: UCSC
Genome Browser.

## [0.19.0] - 2026-09-17

### Added

- **Variant-prioritization case study** (`mrnavax/case_study.py`)
  demonstrating the toolkit's end-to-end pipeline on a curated
  ClinVar-style variant set. Includes:
    - `mrnavax/examples/clinvar_curated.csv` — 12 curated variants
      (8 pathogenic + 3 benign + 1 uncertain; mix of coding and
      regulatory-region variants)
    - `load_clinvar_variants()` — CSV loader with full schema
      validation
    - `score_case_study_variants()` — runs any scoring fn over
      the variant set and returns ranked results
    - `precision_at_k()` — precision@K computation with sorting
    - `summarize_variant_set()` — counts by pathogenicity +
      is_coding bucket
- **End-to-end worked example** in
  `docs/case-studies/variant-prioritization.md` showing how to:
    - Load the curated CSV
    - Run it through `score_variant` with both AlphaMissense
      (coding) and AlphaGenome Atlas AVI (regulatory) wired up
    - Compute precision@K against ground-truth labels
    - Show that regulatory-region pathogenic variants make it
      into the top-5 **specifically because AVI is wired in**
- **New backend check `case_study.variant_prioritization`** —
  exercises the full pipeline stack (load → score → rank →
  precision@K) end-to-end and asserts precision@3 ≥ 0.66 with
  mock backends.
- **11 new unit tests** in `tests/test_case_study.py` covering
  CSV loading, ClinVarVariant fields, precision_at_k
  (empty input, k larger than n, sorting behavior), scoring
  pipeline, and variant set summarization.
- **mkdocs nav** gains a "Case studies" section with the
  variant-prioritization page (3 locales translated).
- Total tests: **229** (was 218). Total backend checks: **29**
  (was 28).

### Worked-example result

With mock backends (no API key needed, runs in CI):

```
Loaded 12 curated variants
Summary: {'n_total': 12, 'n_pathogenic': 8, 'n_benign': 3,
          'n_uncertain': 1, 'n_pathogenic_coding': 5,
          'n_pathogenic_regulatory': 3}

precision@3 = 1.000  (top-3 are all pathogenic)
precision@5 = 0.800  (4 of top-5 pathogenic)
precision@8 = 0.625  (5 of top-8 pathogenic)
```

### Why this matters

The previous releases shipped the Atlas integration but lacked a
**reproducible demonstration** of what the pipeline actually does
on real-world variants. Users had no way to see precision@K, no way
to confirm the routing works end-to-end, and no entry point for
swapping in their own variants. The case study fixes all three:

* Shows users the prioritization output looks like on a known variant set
* Provides a `load_clinvar_variants()` + `score_case_study_variants()`
  pair they can extend with their own data
* Catches full-stack regressions in CI (the new backend check
  exercises load → score → rank → precision@K in a single shot)

## [0.18.0] - 2026-09-17

### Added

- **Recorded AlphaGenome Atlas API response fixture** at
  `tests/fixtures/alphagenome_atlas_sample.json`. Five synthetic variant
  responses modeled on the documented Atlas JSON shape — two coding
  (BRAF V600E, KRAS G12D) and three regulatory-region variants
  covering all three classification bins (low / moderate / high). A
  real captured response can replace this file when an
  `ALPHAGENOME_API_KEY` is available; the synthetic fixture is enough
  to lock the parser shape so upstream schema changes are caught at
  PR time.
- **New backend check `variant.alphagenome_atlas_fixture`** that loads
  the fixture, parses each variant into an `AVIResult` via the
  subprocess adapter payload shape, and asserts:
    - All variants have the required keys (chrom, pos, ref, alt,
      score, classification, is_coding)
    - Score is in [0, 1], classification is one of {low, moderate,
      high}, is_coding is bool
    - The fixture exercises both coding and regulatory variants
    - All three AVI classifications are present
    - AVIResult construction succeeds for every fixture entry
- **8 new unit tests** in `tests/test_atlas_fixture.py` covering
  fixture-file existence, JSON parsing, required fields, AVIResult
  construction, coding/regulatory coverage, all-three-classifications
  coverage, and the backend check integration.
- Total backend checks: **28/28** (was 27/27).

### Why this matters

Catches upstream AlphaGenome Atlas API schema changes at PR time.
Without the fixture, a breaking change in the Atlas response shape
would only surface when a user with an `ALPHAGENOME_API_KEY` ran the
pipeline in production. The fixture locks the parser shape so CI
fails on the first PR with the new shape, not at user runtime.

## [0.17.0] - 2026-09-17

### Added

- **scRNA pipeline auto-wires AlphaGenome Atlas AVI** when variant CSVs
  carry DNA coordinates (`chrom`, `ref_dna`, `alt_dna`). The pipeline
  no longer needs the user to call a separate `variant-regulatory` step
  for regulatory-region variants — both coding and non-coding
  regulatory variants flow through `run_pipeline()` end-to-end.
- **`variant_scores` populated even without a filter** — when DNA
  coordinates are present, every variant gets a per-variant score in
  `PipelineReport.variant_scores` regardless of whether the user set
  `--variant-filter-top-fraction`. Previously this dict was empty
  unless filtering was requested.
- **Pipeline note mentions Atlas integration** — when DNA coordinates
  are present, `PipelineReport.note` now includes "AlphaGenome Atlas
  AVI scores used for non-coding regulatory variants". When no DNA
  coordinates are present (legacy CSV), the note correctly omits the
  mention.
- **5 new unit tests** in `tests/test_scrna_avi.py` covering
  load_variants with DNA fields, run_pipeline note surfacing,
  variant_scores population, and end-to-end runs with and without
  DNA coords.
- **New backend check `scrna.pipeline_with_avi`** — verifies the
  full integration: DNA-coord variants → variant_scores populated +
  note mentions Atlas; legacy CSV (no DNA coords) → note correctly
  omits Atlas mention. Total backend checks: **27/27** (was 26/26).

### Why this matters

The v0.16.0 release shipped `score_variant()` accepting AVI
parameters. v0.17.0 extends that to the full scrna pipeline:
users can now submit a single CSV with mixed coding + regulatory
variants and get both kinds scored through one `run_pipeline()` call.
The pipeline emits a structured report showing per-variant scores
across the entire genome — not just the 2% AlphaMissense covers.

## [0.16.0] - 2026-09-17

### Added

- **AlphaGenome Atlas AVI integration into `variant_scorer`.** The
  `score_variant()` function now accepts `chrom`, `ref_dna`,
  `alt_dna`, and `avi_lookup` keyword parameters. For variants where
  the AVI lookup returns `is_coding=False`, the AVI score becomes the
  dominant signal at weight 0.45 — same role AlphaMissense plays for
  coding-region variants. For coding-region variants, AVI is recorded
  as a secondary signal in `components` but does not dominate.
- **`filter_variants()` extended with `avi_lookup`** — passes DNA-level
  coordinates through to `score_variant()` so the scrna pipeline
  benefits automatically.
- **`Variant` dataclass extended with optional `chrom`, `ref_dna`,
  `alt_dna` fields.** Backward-compatible: existing CSVs without these
  columns still load (fields default to `None`).
- **`sc_rna_pipeline.run_pipeline` auto-wires the AVI lookup** when any
  input variant carries DNA coordinates. The Mock backend is used by
  default (CI + offline use); users with `ALPHAGENOME_API_KEY` set +
  the `[variant-alphagenome]` extra installed get the real Atlas
  adapter automatically.
- **13 new unit tests** in `tests/test_variant_scorer_avi.py` covering
  parameter acceptance, dominant-signal routing for both coding and
  regulatory variants, silent-failure on lookup errors / None returns /
  missing DNA args, rationale strings, and end-to-end mock-as-avi-lookup.
- **Backend check extended** to verify the AVI integration: coding
  variant → AM dominates (rationale has no "(dominant)" on AVI),
  regulatory variant → AVI dominates (rationale has "(dominant)").

### Why this matters

The previous variant prioritization pipeline used AlphaMissense for
coding-region variants only — silently dropping the 98% of the genome
that is non-coding regulatory. With this change, the same
`score_variant()` entry point handles both: the `is_coding` flag from
AlphaGenome Atlas routes the variant to the correct signal. A
protein-coding variant with high AlphaMissense score still ranks high;
a regulatory-region variant with high AVI score now also ranks high
instead of falling through to BLOSUM62-only scoring.

## [0.15.0] - 2026-09-13

### Added

- **AlphaGenome Atlas integration for regulatory-variant prioritization.**
  New module `mrnavax/alphagenome_integration.py` exposing the
  `RegulatoryVariantScorer` Protocol + `AVIResult` dataclass, with two
  backends:
    - `MockRegulatoryVariantScorer` — stdlib-only, deterministic, SHA-256
      hash of (chrom, pos, ref, alt) → [0, 1] AVI score. Always present.
    - `AlphaGenomeCLIAdapter` — subprocess wrapper around the official
      `alphagenome` Python package (Avsec et al., *Nature* 2026). Gated
      behind the new `[variant-alphagenome]` extra; non-commercial use
      only per Google DeepMind's terms.
- **New CLI subcommand:** `mrnavax variant-regulatory --csv variants.csv`
  emits a JSON document with one entry per variant: AVI score,
  classification (low / moderate / high — thresholds match AlphaMissense
  bins), and `is_coding` flag (True for coding-region variants where
  AlphaMissense is the better signal; False for regulatory regions
  where AlphaGenome Atlas is canonical).
- **New bundled example CSV:** `mrnavax/examples/regulatory_variants.csv`
  (7 rows: BRAF V600E, KRAS G12D, APC, TP53 promoter variants, plus
  intergenic / regulatory-region rows).
- **Subprocess shim** at `mrnavax/_shims/alphagenome_cli.py` — small
  (~90 LOC) wrapper that loads the official `alphagenome` package
  inside a subprocess, keeping the adapter module stdlib-only and the
  heavy upstream dep opt-in.
- **New backend integrity check** `variant.alphagenome_atlas` —
  exercises the Protocol contract, mock determinism, AVIResult
  validation, and regulatory_score convenience wrapper. Total
  backend checks: **26/26** (was 25/25).
- **18 new unit tests** in `tests/test_alphagenome_adapter.py` covering
  AVIResult invariants, Protocol runtime_checkable, mock contract,
  real adapter subprocess invocation shape, JSON payload round-trip,
  end-to-end CSV scoring, CLI subcommand registration, and example
  CSV schema.

### Why this matters

AlphaMissense (the existing integration) scores **coding-region
missense variants** — the 2% of the genome that codes for proteins.
AlphaGenome Atlas covers the remaining **98% of non-coding regulatory
variants**, where AlphaMissense is silent. Together, the two adapters
cover the entire genome:

- `score_variant(chrom, pos, ref, alt).is_coding == True`  →  AlphaMissense
- `score_variant(chrom, pos, ref, alt).is_coding == False` →  AlphaGenome Atlas

### Notes

- AlphaGenome Atlas outputs are **non-commercial only** per Google
  DeepMind's terms of service. Commercial use requires the Google Cloud
  Vertex AI deployable, not the public API used here. Same constraint
  as the existing AlphaMissense integration (CC BY-NC-SA 4.0).
- Mock `is_coding` heuristic uses **position parity** (even = coding,
  odd = regulatory). The real Atlas API returns per-variant `is_coding`
  from the upstream query.

## [0.13.1] - 2026-09-13

### Fixed
- **`mrnavax manufacture --cds ...` crashed with `AttributeError: 'NoneType' object has no attribute 'upper'`** when `--utr5` / `--utr3` were not passed (which is the common case). Root cause: `score_manufacturability` called `utr5.upper()` on `None`. Fix: `(utr5 or "").upper().replace("U", "T")` and same for `utr3`. CLI `_manufacture_run` was passing `None` for empty `--utr5/--utr3` flags, so every invocation without UTRs crashed. **Bug existed since v0.7.0.**
- **`mrnavax spatial ...` CLI subcommand was missing entirely.** The CLI dispatcher's `sub.add_parser("spatial", ...)` was never registered, so `mrnavax spatial ...` printed an argparse error while the README and docs claimed it worked. Added `_spatial_run()` to `cli.py` with full argument parsing (--count-file, --locations-file, --platform, --num-modules, --backend, --out) and routed it through the existing backend selector. Now end-to-end runs:
  ```bash
  mrnavax spatial --count-file counts.tsv --locations-file locs.tsv \
      --platform ST --num-modules 10 --backend mock
  ```

### Tests
- New `tests/test_cli_dispatch.py` — 7 tests in 3 classes:
  - `TestCLISubcommandRegistration`: asserts all 7 tools listed in
    `--help`; regression test for spatial registration.
  - `TestSpatialCLIEndToEnd`: runs the CLI as a subprocess with
    real TSV inputs, asserts valid JSON output, error paths for
    missing files, `--out` file writing.
  - `TestManufactureCLI`: regression test for the None-UTR crash.

  All 174 unit tests pass; 25/25 backend integrity checks pass;
  mkdocs strict build passes for all 3 locales.

## [0.13.0] - 2026-09-13

### Added
- **ESM2 protein-language-model Protocol adapter**
  (`mrnavax.protein_lm_protocols` + `mrnavax.protein_lm_adapter`).
  Implements the Applm pattern from Wong et al. 2025
  (arXiv 2508.10541): use a **frozen** protein-LM to embed candidate
  peptides, then a lightweight downstream classifier scores them
  for the task at hand (here: neoantigen immunogenicity).
- **`ProteinLMEmbedder` Protocol** (runtime_checkable) + 2 typed
  dataclasses: `EmbeddingRequest` (input) and `EmbeddingResult`
  (output). Request validates AA alphabet (20 standard AAs),
  lowercase normalization, pooling mode ∈ {mean, cls, sum},
  batch_size ≥ 1 — at construction time.
- **`ESM2Embedder`** real adapter: shells into
  `transformers.AutoModel.from_pretrained` for the ESM2 family
  (verified model IDs: esm2_t6_8M_UR50D=320-dim, esm2_t12_35M_UR50D
  =480-dim, esm2_t30_150M_UR50D=640-dim, esm2_t33_650M_UR50D
  =1280-dim). Lazy-loads on first call; supports mean/cls/sum
  pooling; respects batch_size. Heavy deps (torch + transformers)
  are opt-in via `pip install mrnavax[protein-lm]`.
- **`MockProteinLMEmbedder`** stdlib stub: deterministic k-mer (k=3)
  frequency vectors, L2-normalized. Different sequences → different
  embeddings; same sequence → same embedding (no RNG).
- **`ApplmStyleClassifier`**: frozen-LM + downstream scoring
  (sigmoid-like combination of mean + variance). Real users can
  replace `score()` with `sklearn.linear_model.LogisticRegression`
  or `xgboost.XGBClassifier` once they have labelled data.
- **`lm_immunogenicity_score(peptide, embedder=...)`** integration
  in `neoantigen_screener`: takes a peptide, embeds via ESM2 or
  mock, scores via ApplmStyleClassifier, returns `{score, dim,
  model_id, backend, elapsed_seconds, notes}`.
- **25th backend integrity check** (`neoantigen.esm2_protein_lm_embedder`):
  validates dataclass edge cases, Protocol conformance, mock
  L2-normalization + determinism, classifier score range,
  lm_immunogenicity_score integration, ESM2 dim lookup for all
  3 standard sizes, backend selector dispatch when transformers
  missing.
- **`tests/test_protein_lm.py`**: 37 tests in 9 classes covering:
  dataclass validation (AA alphabet, pooling modes, lowercase
  normalization, dim matching, empty embeddings), Protocol
  runtime_checkable, mock embedder (correct shape, deterministic,
  L2-normalized, different sequences produce different vectors),
  ESM2Embedder (transformers-missing error path, batch-size
  batching, pooling mode routing, dim lookup for known/unknown
  models), ApplmStyleClassifier (score in [0,1], deterministic,
  custom embedder), backend selector (4 modes), end-to-end
  peptide → embedding → score workflow, JSON serialization.

### Reference
Wong B.S.H., Kim J.M., Fung S.H., et al. (2025). *Driving Accurate
Allergen Prediction with Protein Language Models and
Generalization-Focused Evaluation.* arXiv 2508.10541.
DOI: 10.48550/arXiv.2508.10541

ESM-2: Lin, Z., Akin, H., Rao, R., et al. (2023). *Evolutionary-scale
prediction of atomic-level protein structure with a language
model.* Science 379(6637): 1123-1130.

## [0.12.0] - 2026-09-13

### Added
- **STModule spatial-transcriptomics Protocol adapter**
  (`mrnavax.spatial_protocols` + `mrnavax.spatial_module_adapter`).
  Implements the STModule method from Wang et al. 2025 (Genome
  Medicine 17, 18): identifying tissue modules from spatially
  resolved transcriptomics (SRT) data — recurrent cellular
  communities spatially organized to exert specific biological
  functions.
- **`SpatialModuleBackend` Protocol** (runtime_checkable) + 3
  typed dataclasses: `SpatialData` (input), `SpatialModule`
  (per-module output), `SpatialModuleResult` (aggregate).
  `SpatialData` validates platform enum, file existence, and
  `num_modules >= 1` at construction time.
- **`SpatialPlatform` enum** mapping to STModule's `high_resolution`
  flag: ST/Visium use defaults; SlideSeqV2/StereoSeq require
  `--high-resolution --max-iter 100` per the upstream tutorial.
- **`STModuleCLIAdapter`**: shells out to a small R shim
  (`scripts/stmodule_shim.R`) that calls the published R functions
  in order — `data_preprocessing()` → `run_STModule()` →
  `get_assocaited_genes()` (sic, typo in upstream) — and emits
  JSON to stdout matching `SpatialModuleResult.to_dict()`.
  Heavy deps (R 4.4+, Seurat v5, torch, GPUmatrix 1.0.2, CUDA 11.7)
  are only required on the user's machine.
- **`MockSpatialModuleBackend`**: deterministic stdlib stub. Uses
  spot-ID intersection (count ∩ locations) for graceful handling
  of mismatched input. Different gene universe per platform —
  ST has GAPDH/USP4/MAPKAPK2, Visium has CDH1/VIM/KRT8, etc.
- **Backend selector** `select_spatial_module_backend()` chooses
  real-or-mock via `prefer=` override or $PATH detection of Rscript.
- **24th backend integrity check** (`spatial.stmodule_module_identification`):
  validates dataclass construction + validation, Protocol
  runtime_checkable, mock run produces correct n_modules,
  determinism, JSON-serializable, spot/location mismatch handling,
  backend selector dispatches correctly when Rscript missing.
- **R shim** (`mrnavax/scripts/stmodule_shim.R`): ~150 lines
  with proper CLI parsing (optparse), 5-stage exit codes (0=ok,
  1=bad args, 2=missing pkg, 3=run failed, 4=JSON failed),
  structured error messages, JSON output via jsonlite.
- **`tests/test_stmodule.py`**: 45 tests in 11 classes covering:
  dataclass validation, Protocol conformance, mock run correctness,
  determinism, platform-specific gene universes, activity range,
  spot/location mismatch handling, CLI subprocess paths (success /
  nonzero exit / unparseable JSON / timeout), backend selector in
  4 modes, end-to-end JSON serialization, high-resolution CLI
  args wiring, payload parsing edge cases, empty/malformed input
  graceful degradation.

### Reference
Wang R., Qian Y., Guo X., Song F., Xiong Z., Cai S., Bian X., Wong M.H.,
Cao Q.#, Cheng L.#, Lu G.#, and Leung K.S.#. (2025) STModule:
identifying tissue modules to uncover spatial components and
characteristics of transcriptomic landscapes. *Genome Medicine*
17(1): 18.

## [0.11.0] - 2026-09-12

### Added
- **Sim-ICL demonstration selection** (`mrnavax.trial_similar`).
  Implements the Sim-ICL strategy from Fung et al. 2026 (Genome
  Biology, in press): when a downstream task is solved by
  in-context learning, **selecting demonstrations by similarity to
  the query** (rather than random sampling) yields competitive
  performance with protein-LM classifiers in low-shot regimes.
- **DemoCase + DemoStore dataclasses** for typed
  `{patient_text, trial, ground_truth_verdicts}` triples. Ships with
  12 synthetic demos covering BRAF-melanoma, EGFR-NSCLC,
  KRAS-pancreatic, BRCA-prostate, PIK3CA-breast, and other mRNA
  cancer-therapy archetypes (`examples/simicl_demos.json`).
- **TF-IDF cosine ranker** (stdlib-only, ~50 LOC) that ranks demos
  by `(patient + trial)` text similarity. For a BRAF V600E
  melanoma query, the top-3 demos are all BRAF-melanoma trials —
  perfect biological ranking.
- **`build_simicl_prompt()`** splices the few-shot block into the
  TrialGPT prompt before the JSON-shape reminder, preserving
  schema clarity.
- **`score_trial_with_llm(..., demo_store=..., use_simicl=...)`**
  integration: when Sim-ICL is on and a non-empty demo store is
  available, top-K demos are prepended; `simicl-kN` and
  `simicl-demo-ids=...` appear in result notes for auditability.
- **CLI wiring**: `mrnavax trial --matcher trialgpt-simicl`. Falls
  back gracefully when the demo store is empty.
- **Env-var knobs**: `MRNA_AI_SIMICL_TOPK` (default 32), `MRNA_AI_SIMICL_ENABLED`
  (default True), `MRNA_AI_SIMICL_DEMOS` (override JSON path).
- **23rd backend integrity check** (`trial.simicl_demonstration_selection`):
  validates store load, TF-IDF ranking of BRAF-melanoma demos,
  empty-store safety, env-var topk override, DemoCase round-trip.
- **`tests/test_simicl.py`** — 39 tests in 8 classes covering:
  dataclass construction/round-trip, TF-IDF ranker ordering
  + determinism, JSON loading + malformed-file handling,
  prompt construction with/without demos, env-var helpers,
  bundled store loading, and full `score_trial_with_llm`
  integration including `simicl-k` notes.

### Reference
Fung S.H., Zhang Z., Wang R., Miao C., Wong B.S.H., Li K.Y.,
Hong C., Zhou J., Yip K.Y.#, Tsui S.K.W.#, and Cao Q.#. (2026) A
Systematic Evaluation of In-Context Learning in Large Language
Models for Antibody Characterization. *Genome Biology* (in press).

## [0.10.0] - 2026-09-11

### Added
- **RiboDecode Protocol adapters** (`mrnavax.codon_protocols` +
  `mrnavax.codon_ribodecode_adapter`). Wires the published
  RiboDecode package (Li, Wang, Yang et al., *Nat Commun* 16, 9957
  (2025)) into the toolkit via two `runtime_checkable` Protocols:
  - `TranslationPredictor` — predicts translation level per CDS,
    optionally per cellular environment (HEK293T, A549, HeLa, or
    custom RPKM CSV).
  - `CodonOptimizer` — joint translation × MFE optimization via
    the published deep generative model.
- **Two real CLI adapters** (`TranslationModelCLIAdapter`,
  `RiboDecodeCLIAdapter`) shell out to the upstream `pred-translation`
  and `ribo-decode` console scripts. Heavy deps (ViennaRNA, torch,
  CUDA) are only required on the user's machine — the toolkit itself
  remains stdlib-only.
- **Two mock backends** (`MockTranslationPredictor`,
  `MockCodonOptimizer`) that satisfy the Protocols using only stdlib
  (CAI-derived score + LinearDesign for optimization). Used by CI and
  tests; fall back automatically when the upstream binaries aren't
  installed.
- **Typed dataclass contracts**: `RiboDecodeRequest` validates
  length (multiples of 3, ≤4500 nt), mfe_weight ∈ [0,1], and
  `env='custom'` requiring a CSV path — at construction time, before
  the backend ever sees the input.
- **CLI wiring**: `mrnavax codon --backend ribodecode-real
  [--env HEK293T|A549|HeLa|custom] [--env-csv ...] [--mfe-weight ...]
  [--optim-epoch N]`. Falls back gracefully to the mock if the
  upstream binary isn't installed.
- **Backend selector**: `select_translation_predictor()` and
  `select_codon_optimizer()` choose real-or-mock based on $PATH
  + `prefer=` override.
- **22nd backend integrity check** (`codon.ribodecode_protocols`):
  verifies dataclass validation, Protocol runtime_checkable, mock
  translation in [0, 100] range with CAI ordering, and protein
  preservation end-to-end through the mock optimizer.
- **`tests/` directory** with `test_ribodecode_adapter.py` (46
  tests): dataclass validation, Protocol conformance, mock
  translator + optimizer, CLI subprocess mock paths, backend
  selector dispatch, end-to-end protein preservation. Runs in
  ~20 s with stdlib `unittest.mock`.

### Reference
Li, Y., Wang, F., Yang, J. et al. *Deep generative optimization of
mRNA codon sequences for enhanced mRNA translation and therapeutic
efficacy.* Nat Commun 16, 9957 (2025).
DOI: 10.1038/s41467-025-64894-x

## [0.9.1] - 2026-09-10

### Performance
- **LinearDesign 2-4x faster** (no semantic change):
  - **Parent-pointer DP** replaces the full-codon-history copy: each
    state stores `(parent_key, accumulated_codons)` in an append-only
    trace dict, and the traceback walks the chain at the end. Old
    code did O(L) per-state work for O(L²) total; new code is O(L).
  - **Precomputed `_LOG_SCORE_TABLE`**: `(aa, codon) -> float`
    built once at module load. Eliminates 13M+ `dict.get` and
    `math.log` calls per Cas9 run.
  - **`_append_key(prev_key, cand, max_len)`**: faster than
    `"|".join(codons[-n:])` (no list allocation, no split).

### Benchmarks (Apple Silicon, Python 3.12)
| protein      | aa   | nt   | before | after  | speedup |
|--------------|------|------|--------|--------|---------|
| EPO          | 193  | 579  | 726ms  | 310ms  | 2.3x    |
| GFP          | 239  | 717  | 1003ms | 444ms  | 2.3x    |
| mAb heavy    | 450  | 1350 | 4459ms | 1408ms | 3.2x    |
| Luciferase   | 550  | 1650 | 3522ms | 1288ms | 2.7x    |
| **Cas9**     | 1368 | 4104 | 34040s | **8925ms** | **3.8x** |

Same outputs: protein sequences preserved on every test, CAI
deltas identical, codon-change counts identical.

## [0.9.0] - 2026-09-03

### Added
- **TrialGPT-style per-criterion LLM matching** (`mrnavax.trial_llm`)
  - Implements the TrialGPT-Matching approach from Jin et al.
    (*Nature Communications* 2024): per-criterion LLM reasoning,
    each inclusion/exclusion criterion judged independently as
    ``"met" | "unmet" | "uncertain"`` with a short evidence snippet
  - `score_trial_with_llm(patient_text, nct_id, title, inclusion,
    exclusion, *, backend=None)` returns a `TrialMatchResult` with
    per-criterion verdicts and aggregate 0–1 score
  - Aggregate score = (n_met_inclusion / n_total_inclusion) ×
    (n_unmet_exclusion / n_total_exclusion) — both axes required
    for eligibility
  - Schema-validated JSON parsing with case-insensitive dedup and
    backfill of criteria the LLM omitted
  - `MATCH_PROMPT_TEMPLATE` is the canonical TrialGPT prompt
- `trial_matcher.match(..., matcher="trialgpt"|"keyword"|"auto")`:
  - `matcher="trialgpt"` forces per-criterion LLM matching
  - `matcher="auto"` (default) uses TrialGPT when `OPENAI_API_KEY`
    is set, otherwise keyword fallback
  - On TrialGPT failure, automatically falls back to keyword with
    `"trialgpt-fallback"` note
- CLI: `mrnavax trial --matcher {auto,trialgpt,keyword}`
- 21st backend integrity check (`trial.trialgpt_llm`):
  - Skips when `MRNA_AI_FORCE_MOCK=1` (CI without LLM)
  - Skips when `MRNA_AI_SKIP_LLM_CHECK=1` (CI escape hatch)
  - Otherwise verifies per-criterion verdicts and that a matching
    trial outranks an unrelated one

### Mock backend improvements
The `_mock_complete` LLM mock backend now produces real per-criterion
JSON output for trial-matching prompts:
- Parses inclusion/exclusion bullets from the prompt
- Splits patient summary out by marker
- Judges inclusion criteria by keyword overlap (>=50% hit → "met")
- Judges exclusion criteria with negation awareness ("no prior
  therapy" + criterion "Prior therapy" → "unmet")

### Reference
Jin, Qiao, et al. "Matching patients to clinical trials with large
language models." *Nature Communications* 15 (2024): 9074.
DOI: 10.1038/s41467-024-53081-z
Reported: 87.3% accuracy on 1,015 patient-criterion pairs.

## [0.8.0] - 2026-09-03

### Added
- **Real scGPT foundation model integration** (`mrnavax.scgpt_integration`)
  - Loads the `perturblab/scgpt-human` checkpoint (whole-human, 33M
    cells, 60,697-gene vocab, 205 MB)
  - Reimplements scGPT's `FlashTransformerEncoderLayer` in pure PyTorch
    (12 layers, 8 heads, 512 dim, fused `Wqkv` projection split into
    Q/K/V, post-norm) — no `flash-attn` dependency
  - `embed_with_scgpt(matrix, gene_names=...)` produces CLS-token
    cell embeddings. Real TF-IDF baseline: silhouette -0.114 on
    tumor-vs-normal; real scGPT: **+0.214** (+0.328 improvement)
  - `scgpt_available()` check; weights expected at
    `~/.cache/mrnavax/{best_model.pt,vocab.json,args.json}`
  - `bin_expression()` rank-bins raw values into 51 categories per
    scGPT preprocessing
  - `ScGPTConfig.from_json()` mirrors the upstream args.json layout
- `sc_rna_pipeline.embed_with_foundation_model(...)` extended with
  `gene_names` parameter — when set, real scGPT can use the names
  to look up token IDs in its 60,697-gene vocabulary
- 20th backend integrity check (`scrna.scgpt_integration`):
  - Skips when `MRNA_AI_FORCE_MOCK=1` (CI without scGPT weights)
  - Skips when `MRNA_AI_SKIP_SCGPT_CHECK=1` (CI escape hatch)
  - Skips when weights not present (first-run)
  - Otherwise embeds 30 named-gene cells and verifies 30×512 output

### Reference
Cui et al., scGPT: toward building a foundation model for single-cell
multi-omics. *Nat Methods* 21, 1480–1491 (2024).

## [0.7.0] - 2026-09-03

### Added
- **mRNA manufacturability checker** — the wet-lab bridge between
  computational sequence design and what's actually synthesizable.
  Eight checks:
  1. `poly_a_runs` — runs of ≥5 As destabilize the DNA template
  2. `gc_5prime_hairpin` — GC-rich stems (≥70% GC, 30+ nt) at the 5'
     end block ribosome scanning
  3. `kozak_strength` — match to mammalian Kozak consensus
     `GCCRCCATGG`
  4. `are_motif` — AU-rich elements (`UUAUUUAUU` nonamers) in the
     3' UTR trigger mRNA decay
  5. `stop_context` — termination efficiency depends on stop codon
     identity + +4 base (TGA-T and TAA-T are strongest)
  6. `hidden_stops` — internal in-frame stops (must be zero)
  7. `gc_window_uniformity` — local GC stddev > 15% flags IVT yield
     problems and ribosomal stalling
  8. `cpg_balance` — extreme CpG density (suppressed <0.5% or
     excessive >15%) signals silencing or immune activation
- New CLI: `mrnavax manufacture --cds input.fasta [--utr5 ...] [--utr3 ...]`
  returns a JSON report with per-check status, score, severity, and
  summary. Exit code 0 if no errors, 2 if any check is `error`.
- 19th backend integrity check (`manufacture.score_manufacturability`)
  verifies that the checker correctly distinguishes clean vs
  pathological sequences.

### Reference
- Holtkamp et al. (2006) *Blood* 108.
- Kozak (1986) *Cell* 44.
- Chen & Shyu (1995) *Trends Biochem Sci* 20.

## [0.6.0] - 2026-09-03

### Added
- **Full-length LinearDesign** (Do & Woods, *Nature* 2024): the
  600-nt CLI cap is gone. Real LinearDesign uses a **linear-time
  per-step DP** where the state space is bounded by `|Σ|^(W/3)` —
  independent of CDS length. This implementation now does the same:
  state key = the last `W/3` codons (default 7 codons / 21 nt).
  Pareto-prune keeps the highest-translation-score path per key, so
  two paths arriving at the same suffix collapse to one state.
  Full-length Cas9 (4,104 nt) optimizes in **37.6 s** via the CLI
  with protein preservation verified.
- New `LinearDesignResult` fields: `elapsed_seconds`, `n_states_evaluated`
  — useful for benchmarking and reporting.
- 18th backend integrity check (`codon.lineardesign_full_length`) that
  verifies LinearDesign works on a 603-nt synthetic CDS end-to-end and
  preserves the protein sequence. Runs in ~100 ms.

### Changed
- Default `gc_window_size` lowered from 30 → 21 (must be a multiple of
  3 so the codon-suffix state grouping is exact). Window 21 keeps the
  state space at `|Σ|^7 ≈ 64K max` for sub-second typical runs while
  still capturing the same class of local stem structures.
- `optimize_lineardesign` accepts a new `verbose` flag that prints
  progress every 50 codons — useful for full-length runs.
- Codon CLI `--backend lineardesign` now accepts CDS of any length
  (was previously hard-capped at 600 nt).

### Reference
- Do, C. & Woods, D. LinearDesign: a Toolkit for Full-length Stable
  mRNA Design. *Nature* (2024).

## [0.5.0] - 2026-09-03

### Added
- **AlphaMissense (DeepMind) variant pathogenicity integration** via the
  pre-computed predictions TSV (Cheng et al., *Science* 2023). When the
  predictions file is present, the variant scorer weights AlphaMissense
  at 45% of the total score — the dominant signal.
  Verified end-to-end: BRAF.V600E norm jumps from 0.73 (heuristic) to
  0.83 (AlphaMissense-augmented).
- `mrnavax.alphamissense_integration` module with:
  - `build_test_index()` — 5-entry synthetic index for unit tests
    (KRAS.G12D=0.832, KRAS.G12V=0.913, BRAF.V600E=0.954,
    TP53.R175H=0.881, TP53.R248Q=0.872)
  - `load_index()` — streams the 5.5 GB predictions TSV into a 71.7M-entry
    dict and pickles it (cached at
    `~/.cache/mrnavax/alphamissense_index.pkl`)
  - `lookup(uniprot, wt_aa, position, mut_aa)` — O(1) by `(uniprot, aa_change)`
- `UNIPROT_TO_GENE` reverse-lookup table (TP53, KRAS, BRAF, BRCA1/2, EGFR,
  PTEN, PIK3CA, AKT1/2, etc.) — bridges the UniProt-keyed AlphaMissense
  data to the gene-symbol-keyed variant CSVs.
- `load_protein_fasta_detailed()` — new FASTA loader that extracts UniProt
  accessions from standard `>sp|P01116|RASK_HUMAN ...` headers.
- Updated bundled `examples/proteins.fasta` to UniProt-style headers so
  AlphaMissense lookups work out-of-the-box.
- 17th backend integrity check (`scrna.alphamissense_integration`) that
  verifies the AlphaMissense lookup + integration end-to-end using
  synthetic data — runs in milliseconds, no network or model download.

### Caveats
- The AlphaMissense predictions TSV is licensed under **CC BY-NC-SA 4.0**
  (non-commercial, share-alike). It cannot be bundled with mrnavax
  (which is dual-licensed under AGPL-3.0 + commercial). Users must download
  the ~640 MB gzipped TSV separately from
  https://storage.googleapis.com/dm_alphamissense/.
- First-time index build takes ~15 min (streaming 71.7M rows + pickling
  ~5 GB dict). Subsequent loads are <1 s.

## [0.4.1] - 2026-09-03

### Fixed
- `mrnavax.backends --check-all` previously failed when run from
  outside the repo (e.g., from a fresh `pip install` of the wheel),
  because two checks used hardcoded relative paths
  (`"mrnavax/examples/..."`) instead of resolving against the
  installed package location. Now resolves via `_example_path()` and
  passes 16/16 from any CWD. Caught during independent validation.

## [0.4.0] - 2026-09-02

### Added
- **LinearDesign-style codon optimizer** (`--backend lineardesign`) — DP
  over codon choices jointly optimizing translation efficiency and mRNA
  secondary-structure stability. Captures the operational idea of
  Do & Woods, *Nature* (2024). Reference sequence preserved; CAI 0.78 →
  0.998 with measurable structure improvement. Optional on the CLI
  (limited to ≤600 nt to keep DP runtime acceptable); full-length CDS use
  the Python API.
- **Chou-Fasman helix/strand disruption penalty** in variant scoring —
  variants in structured regions (predicted from protein sequence)
  rank higher than the same substitution in a coil. L→P in helix:
  score 0.43 (struct_penalty 0.8). L→P in coil: score 0.27
  (struct_penalty 0.0). Reference: Chou & Fasman (1978).
- **MedCPT integration module** (`mrnavax.medcpt_integration`) —
  real semantic trial retrieval via `ncbi/MedCPT-Query-Encoder` and
  `ncbi/MedCPT-Article-Encoder` from HuggingFace. Verified end-to-end
  on the bundled melanoma test patient: INTerpath-001 ranks #1 with
  cosine similarity 0.601, NSCLC #2 (0.498), pancreatic #3 (0.488).
  Auto-activates when `torch` + `transformers` are installed and the
  model weights are cached. Skipped in CI via
  `MRNA_AI_SKIP_MEDCPT_CHECK=1`. Reference: Jin et al., *Nat Commun*
  15, 1785 (2024).
- New `[neoantigen-medcpt]` and `[trial-medcpt]` extras for one-line
  installs.
- 16th backend integrity check (`trial.medcpt_integration`).

### Changed
- Bumped to v0.4.0.
- `trial --retriever {keyword,dense,auto}` auto-picks MedCPT when
  available, then falls back to TF-IDF dense, then keyword.
- `filter_variants` accepts `protein_sequences: dict[str, str]` for
  per-gene protein context (enables the Chou-Fasman component).

### Fixed
- Codon optimizer now raises `ValueError` on internal stop codons
  instead of silently producing a truncated protein.
- Removed unused `seed` parameter from `optimize_ribodecode` (was
  accepted but ignored — deterministic hill-climb).
- `variant_scorer.score_variant` returns `None` (or raises in `strict`
  mode) for unknown amino acids instead of silently maxing the BLOSUM
  penalty and inflating the score.
- LinearDesign DP bug: Pareto-pruning now tracks the protein sequence
  so codon swaps preserve the amino-acid sequence (verified by
  translation check).
- Test codon table missing `TCA → S` and several IUPAC-degenerate
  entries; expanded to 64 entries including stop codons.

## [0.3.0] - 2026-09-02

### Added
- Backend integrity check module (`mrnavax.backends --check-all`,
  13 checks). CI runs without network access.
- RiboDecode-style codon optimizer (`--backend ribodecode`) with
  optional `--ribo-weights` JSON.
- AlphaMissense-style variant pre-filter (BLOSUM62 + driver genes +
  Δhydrophobicity).
- TF-IDF + truncated SVD cell embedder (scGPT plug point).
- MedCPT-style dense retriever with biomedical synonym expansion.
- `MRNA_AI_FORCE_MOCK=1` env var for CI / sandbox.

## [0.2.0] - 2026-09-02

### Added
- New `scrna` tool.
- `mhcflurry` backend for `neoantigen`.
- `pyproject.toml` with extras.
- MkDocs Material site + GitHub Pages workflow.
- Dual AGPL-3.0-or-later + commercial license.

## [0.1.0] - 2026-09-02

### Added
- Initial release: `codon`, `neoantigen`, `trial`, `lnp`.
