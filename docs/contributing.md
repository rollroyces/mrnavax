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
