# Contributing

Pull requests welcome.

## Constraints

- **Stdlib-only core.** Tools that ship in the default install
  (`pip install -e .`) must depend only on the Python standard
  library. New deps go behind an `optional-dependencies` extra in
  `pyproject.toml`.
- **Determinism.** Tools must produce the same output for the same
  input (no timestamps in default outputs). Mock backends are always
  deterministic.
- **Schema validation.** LLM-backed tools must validate the response
  shape and fall back to the heuristic if the LLM returns malformed
  JSON.
- **Small CLI surface.** Each tool's `_run_cli` should accept `--out`
  and print to stdout if unset.
- **Protocol-based pluggability.** Heavy model integrations
  (RiboDecode, STModule, ESM2, scGPT, mhcflurry) must live behind a
  `runtime_checkable` Protocol or abstract adapter class installed
  via optional extras (e.g. `pip install mrnavax[sota]`).
  See `mrnavax/codon_ribodecode_adapter.py`,
  `mrnavax/spatial_module_adapter.py`,
  `mrnavax/protein_lm_adapter.py` for reference.
- **TDD for new features.** Every new public function must have a test
  written before implementation (see *Development* below).

## Development

```bash
git clone https://github.com/rollroyces/mrnavax.git
cd mrnavax
pip install -e ".[dev,llm]"

# Run the 25 backend integrity checks
python -m mrnavax.backends --check-all

# Run the unit test suite (167 tests)
python -m unittest discover tests

# Smoke test
bash scripts/smoke.sh

# Docs locally
pip install -e ".[docs]"
mkdocs serve
```

## Pull request process

1. Open an issue first for non-trivial changes.
2. **TDD cycle** for new features (per the `test-driven-development`
   skill):
   - **RED**: write a failing test that exercises the wished-for API.
     Run it and confirm it fails for the right reason.
   - **GREEN**: write the minimal implementation to pass.
   - **REFACTOR**: clean up duplication, names, helpers.
3. Follow the **api-integration-verify** discipline when wrapping a
   third-party model: probe the upstream docs / GitHub first, verify
   parameter names + return shapes from real URLs, then design the
   Protocol adapter.
4. Add a `register()` entry in `backends.py` for any new backend.
5. CI must be green on Python 3.11–3.14 before merge.
6. Squash-merge with a `[verified]` commit message once
   `bash scripts/smoke.sh` + `python -m mrnavax.backends --check-all`
   + `python -m unittest discover tests` all pass.

## Adding a new tool

Every new tool should ship with:

1. A typed `dataclass` for input + output (frozen, validated at
   construction).
2. A `runtime_checkable` Protocol for the backend interface.
3. A real adapter that shells out / lazy-loads the upstream model.
4. A stdlib-only mock that satisfies the same Protocol.
5. A `register()` entry in `backends.py` for CI integrity.
6. Tests in `tests/` following strict TDD.
7. A `docs/tools/<name>.{md,zh-Hant.md,zh-Hans.md}` page describing
   the tool, CLI usage, backend matrix, and reference paper.

## Refreshing the recorded Atlas fixture

The bundled fixture at
`tests/fixtures/alphagenome_atlas_sample.json` locks the parser
shape against upstream changes. To refresh it with a real captured
response (once you have an `ALPHAGENOME_API_KEY`):

```bash
# 1. Export your API key (https://deepmind.google.com/science/alphagenome)
export ALPHAGENOME_API_KEY=***

# 2. Trigger the GitHub Actions workflow manually:
#    Actions → Atlas live integration → Run workflow
#    Inputs: mode=record, leave baseline/tolerance as default
#    This runs the real Atlas call and commits the response as
#    tests/fixtures/alphagenome_atlas_live_<timestamp>.json.
#
# 3. Review the PR, copy the captured response over the bundled
#    fixture (or commit alongside), and update the expected-score
#    assertions in any integration tests.
#
# 4. Locally you can also do:
python -m mrnavax.live_atlas_integration --mode record \
    --output tests/fixtures/alphagenome_atlas_live.json
```

The weekly cron (`.github/workflows/atlas_integration.yml`) runs
the same harness in `--mode regression` and compares the live
score against a baseline with ±5% tolerance — failing the workflow
if the score drifts outside tolerance. This catches upstream Atlas
API breakage in real time (the recorded fixture catches it with
a 7-day lag).

## Considered but not built

The following features have been **spiked and intentionally
deferred**. Each entry documents the question, the verdict, and the
conditions that would justify revisiting it. Future contributors
should read this list before proposing equivalent features — the
"obvious" version of each has already been considered and rejected
with reasons.

### `docs/live-atlas.md` — auto-updating weekly results page (v0.24.0 spike)

**Question:** Should the weekly Atlas cron also publish its result
to the public docs site as a "last live run + history" page?

**Verdict:** PARTIAL → **deferred** (smallest viable change is
~150 LOC of workflow + script + page; cost outweighs payoff until
the API key is configured).

**Negative case (why deferred):**

1. The API key isn't set on this repo yet, so the page would
   display its own disabled state to the world for the foreseeable
   future. We're not delivering visible value.
2. The 7-day cadence + the v0.18.0 recorded fixture already catch
   Atlas breakage. The page adds *visibility* only when something
   breaks — but during normal weeks the table is ~52 identical
   "✅ within ±0.05" rows.
3. Cost: one new script, one new workflow job, ongoing
   maintenance, and a `contents: write` permission on the cron.

**Revisit when:** (a) the `ALPHAGENOME_API_KEY` GitHub secret is
configured, AND (b) we want a user-visible freshness signal that
distinguishes "alive" from "abandoned" tooling. The spike's full
design is preserved at `.considerations/live-atlas-docs-page.md`
in the repo root (or search the git history for the v0.24.0 spike
commit).

### More entries to come

When a feature is rejected via this pattern, document the verdict
here in the same commit that closes the spike. Future readers
should be able to grep CONTRIBUTING for "Considered but not built"
and see the full history of feature requests that were spiked but
not shipped.

## License

By contributing, you agree that your contributions will be licensed
under the project's dual-license terms (AGPL-3.0-or-later for
open-source use; commercial license available on request).
