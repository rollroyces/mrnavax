# Atlas API key activation

The toolkit ships with real-data AlphaGenome Atlas integration, but the
Google AI Studio API key needs to be configured before the live weekly
cron can pull real data. This page covers **three paths** to activate
the key:

## Why three paths?

The key can live in (at most) three places:

  1. `ALPHAGENOME_API_KEY` **environment variable** — preferred for CI
  2. **`~/projects/alphagenome-work/.alphagenome_key`** — preferred for
     local development (auto-detected; added in v0.25.1)
  3. **`~/.alphagenome_key`** — alternate home location (added in v0.25.1)

The loader resolves them in this order: explicit `api_key=...` arg >
`ALPHAGENOME_API_KEY` env var > helper file at
`~/projects/alphagenome-work/.alphagenome_key` > `~/.alphagenome_key`.

## Path 1: Local development (recommended)

If you cloned this repo and want Atlas calls to "just work" from your
laptop, write the 39-char Google AI Studio key (starts with `AIza`) to:

```bash
mkdir -p ~/projects/alphagenome-work
echo -n 'YOUR_39_CHAR_KEY_HERE' > ~/projects/alphagenome-work/.alphagenome_key
chmod 600 ~/projects/alphagenome-work/.alphagenome_key
```

The loader auto-detects this file the next time `score_variant()` or
`live_atlas_integration.py` runs.

## Path 2: GitHub Actions weekly cron (recommended for CI)

To enable the weekly Atlas integration test:

  1. Get a Google AI Studio API key from
     https://aistudio.google.com/apikey (free tier available).
  2. Open https://github.com/rollroyces/mrnavax/settings/secrets/actions/new
  3. Name: `ALPHAGENOME_API_KEY`
  4. Value: paste the 39-char key
  5. Click "Add secret"

The next scheduled cron (every Monday 06:00 UTC) or a manual `gh
workflow run atlas_integration.yml` will hit the live Atlas API and
record the result against the fixture at
`tests/fixtures/alphagenome_atlas_live.json`. If the upstream API
contract changes, the integration test catches it within a week.

## Path 3: CI / scripting (env var)

For ephemeral environments (Docker, CI runners without the helper
file):

```bash
export ALPHAGENOME_API_KEY='YOUR_39_CHAR_KEY_HERE'
```

## Verify the activation

After configuring any path, run the verification script to confirm
end-to-end Atlas connectivity:

```bash
python scripts/check_atlas_activation.py
```

Expected output (BRAF V600E from hg38):

```
[helper-file] OK: /Users/hermes/projects/alphagenome-work/.alphagenome_key (40 bytes, prefix matches Google AI key: True)
[package] alphagenome installed
[key] OK: resolved 39-char key (prefix matches Google AI: True)
[live-call] score=1.000, classification='high', is_coding=True, n_scorers=19
[fixture] recorded=1.000/'high', live=1.000/'high', matches=True
[result] PASS: live Atlas matches recorded fixture
[result] Activation verified — weekly Atlas cron will work
```

If the live call fails, the script prints the upstream error message
and exits non-zero.

## What "active" means

The Atlas API key is "active" when **all three** are true:

  1. ✅ The key resolves successfully (via any of the three paths above)
  2. ✅ The alphagenome package is installed (`pip install 'mrnavax[variant-alphagenome]'`)
  3. ✅ The live BRAF V600E call returns `score≈1.0`, `classification='high'`, `is_coding=True` (matches the recorded fixture)

When all three are true, the GitHub Actions weekly cron will exercise
real Atlas calls instead of skipping with the "ALPHAGENOME_API_KEY not
configured" notice.

## Privacy

The API key is never logged, never echoed, and never written to disk
by the toolkit. The verification script prints only its length and
prefix-match status, not the key itself.