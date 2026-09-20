# Spike: live Atlas weekly results docs page — design

**Status:** spike (no code yet)
**Date:** 2026-09-21
**Goal:** Decide whether it's worth shipping an auto-updating "Atlas weekly results" page on https://rollroyces.github.io/mrnavax/ that pulls from the existing weekly GitHub Actions cron.

---

## 1. The shape of the existing pieces

Verified by reading the files directly:

- **`mrnavax/live_atlas_integration.py`** (462 LOC, stdlib only) — `LiveAtlasReport` dataclass with `to_dict()`; `run_live_atlas_check(mode, baseline, tolerance, …)` returns a report; CLI prints JSON to stdout and exits 0/1. Already does **the only thing we need** when called in `regression` mode: captures `score`, `classification`, `is_coding`, `captured_at`, `delta`, `baseline`, `within_tolerance`, `error`. No code change required.
- **`.github/workflows/atlas_integration.yml`** (87 LOC) — weekly Monday 06:00 UTC cron + `workflow_dispatch` for manual `record|regression` runs. Required secret `ALPHAGENOME_API_KEY`. Has a `if: env.ALPHAGENOME_API_KEY != ''` guard that **already gracefully degrades** to a notice when the secret is unset. Currently no `upload-artifact` step — stdout JSON goes nowhere.
- **`docs.yml`** — runs on push to `main`, builds with `mkdocs build --strict`, uploads `site/` to GitHub Pages. **Never runs on a schedule.** That's the key constraint: the weekly cron runs at 06:00 UTC Monday; the Pages site only rebuilds when something is pushed to `main`.
- **`mkdocs.yml`** — `mkdocs-static-i18n` plugin with `docs_structure: suffix`. The `case-studies/variant-prioritization.md` page exists only in English and the build still passes — confirmed missing translations **fall back to the default-locale page** at the locale-prefixed URL. So we do **not** need to write three copies upfront.
- **Existing locale set:** en (default), zh-Hant, zh-Hans. New pages get en copy first; zh-Hant/zh-Hans can be added as separate `.zh-Hant.md` files later without blocking the ship.

**Key consequence of "docs only rebuilds on push":** the cron can't directly write to the docs site. We need a hand-off — the cron writes an artifact, and a separate push to `main` triggers a rebuild. Two viable shapes:

| Approach | Pros | Cons |
|---|---|---|
| **A. Cron → artifact → docs workflow `download-artifact` on push** | No new workflow; reuses `docs.yml`. | Requires `main` push to fire; docs stale until next push. |
| **B. Cron → commit back to `main` → `docs.yml` rebuilds** | Always-fresh site on Monday morning. | Adds a `contents: write` permission to the cron; risks commit thrash / merge conflicts if docs change concurrently. |
| **C. Cron → artifact → docs workflow runs on `workflow_run` trigger after cron** | Decoupled; always up-to-date; no `contents: write`. | `workflow_run` is a known but supported trigger; minor added complexity. |

**Recommendation: C.** It keeps the docs freshness contract (page reflects the most recent weekly run, not the most recent `main` push), avoids `contents: write`, and uses only well-supported GitHub features.

---

## 2. Concrete change set

### 2.1 Atlas cron: save + upload the report (~12 lines added)

**File:** `.github/workflows/atlas_integration.yml`

After the `Run live Atlas integration` step, in both the key-set branch and the skip-notice branch, write the JSON to a file and upload it. Two changes:

```yaml
      # Existing step unchanged except: capture stdout to a file in addition to the default log.
      - name: Run live Atlas integration
        if: env.ALPHAGENOME_API_KEY != ''
        env:
          ALPHAGENOME_API_KEY: ${{ secrets.ALPHAGENOME_API_KEY }}
        run: |
          MODE="${{ github.event.inputs.mode || 'regression' }}"
          BASELINE="${{ github.event.inputs.baseline || '0.72' }}"
          TOLERANCE="${{ github.event.inputs.tolerance || '0.05' }}"
          ARGS="--mode ${MODE} --baseline ${BASELINE} --tolerance ${TOLERANCE}"
          if [ "${MODE}" = "record" ]; then
            ARGS="${ARGS} --output tests/fixtures/alphagenome_atlas_live.json"
          fi
          mkdir -p live-atlas
          python -m mrnavax.live_atlas_integration ${ARGS} | tee live-atlas/report.json
          echo "captured_at=$(jq -r '.captured_at' live-atlas/report.json)" >> $GITHUB_OUTPUT

      - name: Write skipped-run stub when ALPHAGENOME_API_KEY is not set
        if: env.ALPHAGENOME_API_KEY == ''
        run: |
          mkdir -p live-atlas
          python -c "import json,datetime; json.dump({'skipped': True, 'reason': 'ALPHAGENOME_API_KEY not configured', 'captured_at': datetime.datetime.now(datetime.timezone.utc).isoformat()}, open('live-atlas/report.json','w'))"

      - name: Upload live Atlas report
        if: always()
        uses: actions/upload-artifact@v5
        with:
          name: live-atlas-report
          path: live-atlas/report.json
          retention-days: 90
          if-no-files-found: warn
```

The `if: always()` ensures the artifact is uploaded even when the live step fails (so the docs page can show "last successful run + last attempted run" honestly).

**Estimated diff:** ~12 added lines, 0 removed. Permissions unchanged (`contents: read` is sufficient — we only upload an artifact, no `contents: write`).

### 2.2 Docs workflow: consume the artifact on `workflow_run` (~25 lines added)

**File:** `.github/workflows/docs.yml`

Add a second job `build` keyed off `workflow_run` of `atlas_integration.yml`. It downloads the artifact, regenerates the live-results page from a tiny Python template, and commits the generated markdown to the repo, which triggers the existing push-driven `build` job.

```yaml
  refresh-live-atlas-page:
    if: ${{ github.event.workflow_run.conclusion == 'success' || github.event.workflow_run.conclusion == 'failure' }}
    runs-on: ubuntu-latest
    permissions:
      contents: write
    steps:
      - uses: actions/checkout@v6
        with:
          token: ${{ secrets.GITHUB_TOKEN }}

      - uses: actions/setup-python@v7
        with: { python-version: "3.12" }

      - name: Install
        run: pip install -e .

      - name: Download latest Atlas report artifact
        uses: actions/download-artifact@v5
        with:
          name: live-atlas-report
          path: live-atlas/
          github-token: ${{ secrets.GITHUB_TOKEN }}
          run-id: ${{ github.event.workflow_run.id }}

      - name: Render live Atlas page
        run: python scripts/render_live_atlas_page.py \
               --report live-atlas/report.json \
               --workflow-run "${{ github.event.workflow_run.html_url }}" \
               --out docs/live-atlas.md

      - name: Commit regenerated page
        run: |
          git config user.name "github-actions[bot]"
          git config user.email "41898282+github-actions[bot]@users.noreply.github.com"
          git add docs/live-atlas.md
          git diff --cached --quiet || git commit -m "docs: refresh live Atlas results page (run ${{ github.event.workflow_run.id }})"
          git push
```

Concurrency matters — the existing `build` job already has `concurrency: group: pages, cancel-in-progress: true`. The `refresh-live-atlas-page` job should serialize against itself so concurrent cron firings don't race:

```yaml
concurrency:
  group: live-atlas-refresh
  cancel-in-progress: false   # never cancel a partially-written page commit
```

**Estimated diff:** ~30 added lines.

### 2.3 Tiny stdlib template renderer (~80 lines, new file)

**File:** `scripts/render_live_atlas_page.py`

Reads `live-atlas/report.json` (already in the `LiveAtlasReport.to_dict()` shape), pulls the last 8 runs from `docs/live-atlas-history.json` (auto-appended), and emits a single markdown file. Stdlib only (`json`, `argparse`, `pathlib`, `datetime`). Renders a table like:

```
| Run | Captured (UTC) | Variant | AVI | Δ vs baseline | Bin | Status |
|---|---|---|---|---|---|---|
| #47 | 2026-09-15 06:00 | chr7:140753336 T>A | 0.718 | −0.002 | high | ✅ within ±0.05 |
| #46 | 2026-09-08 06:00 | chr7:140753336 T>A | 0.724 | +0.004 | high | ✅ within ±0.05 |
…
```

The renderer appends the new run to `docs/live-atlas-history.json` (committed by the same job) so the page shows history, not just the latest snapshot. Cap history at ~52 rows (one year of weekly runs).

**Estimated size:** ~80 LOC, one new file. Pure stdlib.

### 2.4 New page (~40 lines per locale, en only at ship time)

**File:** `docs/live-atlas.md` (en, canonical)

Header: "Last run: {{captured_at}} UTC — [View run #{{n}}]({{workflow_run_url}})". Then a status badge line ("✅ all green" / "⚠️ API call failed — see run #{{n}}"). Then the historical table. Then a brief footer explaining what this is and how to interpret it.

**File:** `docs/live-atlas.zh-Hant.md` and `docs/live-atlas.zh-Hans.md`: defer. The mkdocs-static-i18n fallback shows the English page at `/zh-Hant/live-atlas/` and `/zh-Hans/live-atlas/` automatically. Add localized versions as a follow-up PR — this is the same pattern `case-studies/variant-prioritization.md` already uses.

**File:** `mkdocs.yml` — add one nav line under a sensible section:

```yaml
  - Backends: backends.md
  - Live Atlas results: live-atlas.md          # ← new
  - Case studies:
      - Variant prioritization: case-studies/variant-prioritization.md
```

**Estimated diff:** ~45 added lines (40 page + 1 nav), 0 removed.

### 2.5 Summary of touched files

| File | Action | Lines added | Lines removed | Net |
|---|---|---:|---:|---:|
| `.github/workflows/atlas_integration.yml` | modify | ~12 | 0 | +12 |
| `.github/workflows/docs.yml` | modify | ~30 | 0 | +30 |
| `scripts/render_live_atlas_page.py` | create | ~80 | 0 | +80 |
| `docs/live-atlas.md` | create | ~45 | 0 | +45 |
| `mkdocs.yml` | modify | ~2 | 0 | +2 |
| **Total** | | **~170** | **0** | **~170** |

Plus auto-managed artifacts:
- `docs/live-atlas-history.json` (52-row rolling window, regenerated each run)

No Python source files in `mrnavax/` need to change. **Zero impact on the library's public surface.** All deps remain stdlib (the workflows use `actions/checkout@v6`, `actions/setup-python@v7`, `actions/upload-artifact@v5`, `actions/download-artifact@v5`, `jq` — all already standard in the repo's other workflows or provided by `ubuntu-latest`).

---

## 3. Graceful degradation under no-API-key

This is the most important "must work" requirement and it's already handled for free:

- The atlas cron has the existing `if: env.ALPHAGENOME_API_KEY == ''` branch that writes a stub JSON (`{skipped: True, reason: "..."}`) and now uploads it as an artifact.
- The docs workflow's `render_live_atlas_page.py` reads the report, sees `skipped: True`, and renders a page like:

  ```
  > ⚠️ **Live Atlas runs are disabled.** This repository does not have an
  > `ALPHAGENOME_API_KEY` configured. The page is updated weekly when
  > the key is set; until then it shows historical mock-mode snapshots.
  > [Set the secret →](https://github.com/rollroyces/mrnavax/settings/secrets/actions)
  ```

- The page builds, deploys, and shows honest content whether or not the key is set. **No code path crashes on missing data.**

The same graceful-degradation logic covers: API call fails, JSON malformed, classification empty, baseline missing — all those states map to a "⚠️ status" row in the table without breaking the page.

---

## 4. Risks & mitigations

| Risk | Likelihood | Mitigation |
|---|---|---|
| **Artifact download from `workflow_run` doesn't find the file** (artifact retention, zip layout) | Medium | Use `actions/download-artifact@v5` with the same `name` we uploaded; the workflow file path is stable. Add a `--report` existence check in `render_live_atlas_page.py` that errors clearly if the artifact is empty (so we see it in the run log). |
| **`workflow_run` doesn't fire when the source workflow is skipped** | Low | `workflow_run` fires for all conclusions including `skipped`; we render regardless and show the stub. |
| **Page commit clobbers concurrent human edits to `live-atlas.md`** | Low (only the cron edits it) | The render script owns the file; add a banner comment at the top: `<!-- AUTOGENERATED by scripts/render_live_atlas_page.py — do not edit -->` and add `live-atlas.md` to the merge-blocklist in `contributing.md`. |
| **History file grows unbounded** | Low | Cap at 52 entries (one year). Script trims oldest on append. |
| **The `captured_at` is the time of the *cron attempt*, not the API response time** | Low (cosmetic) | We use the report's `captured_at` field (set in `live_atlas_integration.py:319` *after* the API returns), which is what we want. |
| **`mkdocs build --strict` rejects a page missing from a nav section** | Low | Add the page to `mkdocs.yml` nav. The `--strict` mode also rejects broken internal links — render script must use `mkdocs build --strict`-friendly relative links only. |
| **`contents: write` permission expands the docs workflow's blast radius** | Low | Limited to the `refresh-live-atlas-page` job only; the `build`/`deploy` jobs keep `contents: read`. |
| **Three-locale fall-back confusion** ("why is this page in English on the Chinese site?") | Low | Footer note: "This page is maintained in English. Translations are tracked in issue #XYZ." (Or just don't explain — the rest of the site does the same fallback.) |
| **`gh api` path: tempting but not used** | — | We use GitHub's first-class artifact download via `actions/download-artifact@v5` rather than `gh api repos/.../actions/artifacts/.../zip` because the latter requires a token and is more brittle. |

---

## 5. Verdict

**PARTIAL** — ship a subset, not the full design.

**Ship:**
- §2.1 (cron writes + uploads artifact)
- §2.3 (renderer)
- §2.4 en page + §2.4 nav entry
- Defer §2.2 `refresh-live-atlas-page` job for the first iteration.

**Why the partial:** the spike target is "smallest viable change." The first iteration can be the cron artifact upload + the renderer + the page generated **locally** by the cron workflow itself, committed in-place. That is: drop the `contents: write` job and just have the atlas cron commit to `main` directly. This is uglier (cron has `contents: write`) but it's ~20 fewer lines of YAML and avoids the `workflow_run` complexity. Decide based on taste.

Specifically:
- **Option C-revised (recommended for first ship):** atlas workflow gets `contents: write`, runs the renderer, commits `docs/live-atlas.md` + `docs/live-atlas-history.json` to `main`. Push triggers `docs.yml` automatically. Total new code: ~150 LOC, one new script, no `workflow_run` indirection.
- **Option C (full design above):** ship if the cron committing to `main` feels wrong (e.g., the team dislikes non-PR-driven commits to `main`). Add ~30 LOC.

**Do not ship:**
- The three-locale versions of the page at first. Fallback handles it; add them in a follow-up PR when someone asks.

---

## 6. Honest negative case — why this might be wrong

1. **The whole point is undermined by the 7-day cadence.** The Atlas is stable. If we measure the same BRAF V600E variant weekly, we will see ~52 nearly-identical rows of "✅ within ±0.05" with the same number. **Real value of the page only emerges when something breaks.** Until then, the page is a vanity dashboard. *Mitigation: the page is cheap (150 LOC) and trivially removable if unused; the asset's value compounds only when breakage hits, which is exactly when you want it.*

2. **The API key isn't set on this repo (evidenced by the existing skip-notice branch).** That means for the foreseeable future the page will show the "⚠️ live runs disabled" stub. We will be deploying a page that displays its own disabled state to the world. *Mitigation: the page is still useful documentation (it explains what the harness does) and degrades honestly; the team can enable it whenever the key is configured. The 150 LOC investment is recoverable cost.*

3. **Workflow_run + artifact download has bitten repos before.** GitHub Actions artifact download APIs occasionally race with retention policy (90-day default; we set 90 explicitly). If a user lands on the site 91 days after the last live run, the artifact is gone and the page renders the stub. *Mitigation: the stub already handles this case; the in-repo `docs/live-atlas-history.json` is the persistent cache.*

4. **The atlas cron running weekly means a broken Atlas API generates a broken page on Monday morning.** A failed run still uploads the artifact (we set `if: always()`), the renderer still emits a page with a "⚠️" badge, and the site deploys — but a regression in CI is now visible to users before the team notices. *Mitigation: add a `if: failure()` step in the cron that posts a GitHub issue / opens a Slack ping — that's a separate, easy follow-up.*

5. **We are committing generated content to a public repo.** Every API response (even just the score) lands in git history. Today's payload is non-sensitive (it's a publicly-known BRAF V600E score), but if someone later expands the harness to capture more variants, the privacy story changes. *Mitigation: cap the captured payload at the existing `LiveAtlasReport` shape; add a CI lint that rejects PRs adding new fields to `to_dict()`. Cheap to add now, expensive later.*

6. **The expected payoff is "we notice breakage faster."** If no one reads the public docs site regularly, this is a dead feature. *Mitigation: add a "Last run" link from `backends.md` (where Atlas already lives) so visitors who care about the live integration see the freshness signal at the natural point. This costs ~1 line and is independent of the rest.*

---

## 7. Recommended next steps (if you proceed)

1. Decide between **Option C-revised** (cron commits to `main`, ~150 LOC) and **Option C** (full design, ~170 LOC + the `refresh-live-atlas-page` job). My recommendation: **C-revised** for first ship.
2. Add a `--mode mock` invocation to the docs workflow so the page can be generated in PR builds without an API key (and so the page renders meaningfully in `mkdocs build --strict` on PRs).
3. Add the GitHub issue templates / Slack ping for cron failures (deferred per risk #4).
4. Add zh-Hant + zh-Hans translations in a follow-up PR after the English page has lived for a week or two.
