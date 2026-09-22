# CLAUDE.md

Project-level instructions for Claude Code.

## Project Overview

NiceGUI-based web application for genetics cohort data visualization. Uses TanStack tables for dataframes, Polars for data processing, and ECharts for plots. Managed with `uv`.

Two deliberate exceptions to "ECharts for plots": `stats_panel.py` and the cohort PGS violins in `cohort_pgs_tab.py` use Plotly. ECharts has no violin or box series, and the two live as sibling tabs on the cohort page.

## Commands

- **Run app**: `uv run genetics-viz /path/to/config.yaml`
- **Run tests**: `uv run pytest`
- **Lint**: `uv run ruff check .`
- **Format**: `uv run ruff format .`
- **Compile check**: `uv run python -m py_compile src/genetics_viz/path/to/file.py`
- **Always use `uv run`** to execute Python commands, never bare `python`.
- **Always use `uv add`** to add dependencies, never bare `pip install`.

## Project Structure

```
src/genetics_viz/
  cli.py                          # Typer CLI entry point (accepts YAML config file)
  app.py                          # NiceGUI app init, config loading, static files
  config_model.py                 # YAML config dataclasses, load/save, password helpers
  models.py                       # Cohort, Family, Sample dataclasses + DataStore
  components/
    tanstack_table.py             # DataTable: TanStack table bridge (JS/CSS injection)
    column_selector.py            # Column visibility dialog with presets
    filters.py                    # Reusable filter menu components
    header.py                     # Shared page header with data dir dropdown + user menu
    validation_loader.py          # Validation TSV loading and badge logic
    variant_dialog.py             # IGV variant inspection dialog
    sv_dialog.py                  # IGV SV inspection dialog
    diagnostic_dialog.py          # Diagnostic review dialog
  config/                         # YAML config files (loaded at module level)
    ancestry.yaml                 # Ancestry region colors + PGS z-score color scale
    clinvar_colors.yaml           # ClinVar significance term -> color
    column_names.yaml             # Column display names, groups, sorting, drop flags
    cytobands_hg38.tsv            # Cytoband data for ideogram rendering
    continuous_scores.yaml        # Continuous score color ranges
    pipeline_steps.yaml           # ghfc-ngs step labels/order + run status chips
    vep_consequences.yaml         # VEP consequence terms, impacts, colors
    view_presets.yaml             # Column visibility presets
  pages/
    login.py                      # Login page (no auth required)
    profile.py                    # User profile (view role, change password)
    search.py                     # Cohort-wide variant search with individual filters
    admin/
      directories.py              # Admin: manage data directories
      users.py                    # Admin: manage users
    cohort/
      home.py                     # Home page with cohort cards
      cohort.py                   # Cohort overview page
      family.py                   # Family detail page with tabs
      components/
        ancestry_tab.py           # Ancestry PCA tab (per-family)
        pgs_tab.py                # Polygenic scores tab (per-family)
        cohort_ancestry_tab.py    # Ancestry PCA tab (cohort-wide)
        cohort_pgs_tab.py         # Polygenic score violins (cohort-wide)
        docs_panel.py             # ghfc-ngs params reference, fetched at runtime
        params_panel.py           # Workflow params file, shown verbatim
        status_panel.py           # Pipeline run status (Workflow tab + home card)
        wombat_tab.py             # WOMBAT analysis tab (per-family)
        dnm_tab.py                # DNM analysis tab (per-family)
        svs_tab.py                # SV analysis tab (per-family)
        stats_panel.py            # Carrier stats box/bar plots
    validation/
      file.py                     # Per-file validation page
      all.py                      # Aggregated validation page
      statistics.py               # Validation statistics page
      waves.py                    # WAVES validation page
      wave.py                     # Individual wave validation page
    diagnostic/
      all.py                      # Aggregated diagnostic page
      statistics.py               # Diagnostic statistics page
  utils/
    ancestry.py                   # Ancestry/PGS file discovery, reference panel, colors
    ancestry_plots.py             # Shared ECharts builders for the PCA scatters
    auth.py                       # Auth helpers: check_auth, can_write, get_current_user
    clinvar.py                    # ClinVar color/display utilities (from YAML)
    cohort_state.py               # .ghfc-ngs.state.json parsing, drift, step rows
    cytobands.py                  # Cytoband, chromosome, ideogram constants
    data.py                       # Multi-store registry, per-user data dir selection
    gene_scoring.py               # Gene scoring and color coding
    ghfc_ngs_docs.py              # Upstream params.md fetch + link rewriting
    pipeline_params.py            # <NAME>.params.yml discovery and parsing
    score_colors.py               # Continuous score color ranges
    vep.py                        # VEP consequence utilities (from YAML)
    pedigree_labels.py            # Sex/phenotype code normalisation and labels
    sharding.py                   # Two-level sharded directory path resolution
    view_presets.py               # View preset loading with reload support
    column_names.py               # Column name/group/sorting/schema utilities
  static/
    css/data_table.css            # TanStack table styles
    js/data_table.js              # TanStack table JS bridge
```

## Key Patterns

### State Management
- Filter state uses mutable dicts for pass-by-reference in closures: `{"value": [...]}` for lists, `{"value": False}` for booleans.
- `@ui.refreshable` functions support forward references (callbacks defined after the UI that calls `.refresh()`).

### Closure Factory
Use `make_handler(param)` returning inner `handler(e)` for loop variable capture:
```python
for name in names:
    def make_handler(n):
        def handler(e):
            # n is captured correctly
        return handler
    ui.checkbox(name, on_change=make_handler(name))
```

### Button Visual State
Toggle button appearance with Quasar props:
```python
button.props(remove="outline", add="unelevated color=green")  # active
button.props(remove="unelevated color=green", add="outline")  # inactive
```

### Config Utilities
Config files (YAML) are loaded once at module level.

### Authentication & Authorization
- `check_auth()` returns `RedirectResponse("/login")` if not authenticated, `None` otherwise. Use at top of every page handler: `if redirect := check_auth(): return redirect`.
- `get_current_user()` / `get_current_role()` read from `app.storage.user`.
- `can_write()` returns `True` if role is `"curator"` or `"administrator"`.
- `is_admin()` returns `True` if role is `"administrator"`.
- Write operations (save validation/diagnostic) are gated with `can_write()` both at UI level (hide save button) and backend level (guard in save function).

### Multi-Data-Directory
- `get_data_store()` reads `app.storage.user['data_dir']` for per-user data directory selection. Same signature as before — all calling files work unchanged.
- `get_static_prefix()` returns the per-user URL prefix (e.g., `/data-0`) for IGV static file URLs.
- Static files are registered centrally in `app.py` via `nicegui_app.add_static_files()` — NOT in individual pages.

### Sharded Directories
- `samples/` and `families/` directories may use two-level sharding: `samples/<shard1>/<shard2>/<id>/`.
- Shard keys: strip `-`, `.`, `_` from entity ID, then shard1 = last char (uppercased), shard2 = second-to-last char (uppercased).
- Auto-detected per data directory (if all children of `samples/` are single-char dirs → sharded). Cached.
- Use `get_sample_path(data_dir, id)` / `get_family_path(data_dir, id)` for filesystem paths.
- Use `get_sample_url(data_dir, id)` / `get_family_url(data_dir, id)` for URL segments.
- All functions in `utils/sharding.py`. Never construct `samples/{id}` or `families/{id}` paths manually.

### Ancestry / Polygenic Scores
- Per-family results live in `<family>/ancestry/`, named `<family_id>.apgs_b<bundle>_<filters>.<kind>.tsv` (kinds: `pcs`, `ancestry`, `pgs_zscore`). The bundle version contains dots, so discovery regexes need a lazy `(.+?)` capture, not `([^.]+)`.
- The PCs are projections onto a **versioned** reference panel, so the panel a family is plotted against is the one named in its own filename tag. Config supplies only the root (`ancestry_reference_dir`); the `v<bundle>/bundle/labels/reference_pcs.tsv.gz` path below it is derived per family. Never read the reference root from the `.qc.json` files — they record the compute-cluster path.
- A missing reference panel is a degraded state, not an error: the Ancestry tab still plots the family's own samples.
- `utils/ancestry.py` holds discovery, probes, the cached reference panel, region colours, and the shared `NULL_VALUES` — the PGS z-score file writes the literal string `nan`, which must be nulled or polars types the column as a string and sorting breaks.
- **Cohort-level results** live in `cohorts/<name>/ancestry/<name>.apgs_b<bundle>_<tag>.<kind>.tsv` — same naming, cohort name as the prefix, plus an extra `family_id` column. `find_ancestry_files_in(directory, prefix)` serves both layouts; `find_cohort_ancestry_files(cohort)` and `find_ancestry_files(data_dir, family_id)` are the wrappers. Not every cohort has them.
- **Several filter-tag variants of one bundle coexist** (`dp6gq15`, `dp10gq20` — minimum depth and GQ). `_tag_strictness` ranks them so the strictest wins by default; a plain string sort picks `dp6gq15` because "6" > "1" at the third character. Pages with more than one variant show a selector.
- Plot builders shared by the family and cohort tabs live in `utils/ancestry_plots.py`. They are pure functions returning ECharts option/series dicts, so neither page component has to import the other's privates.
- Normalise pedigree sex/phenotype through `utils/pedigree_labels.py` before grouping — some cohorts write the codes float-stringified (`-9.0`, `1.0`), which a bare comparison against `"-9"` misses.

### Cohort Discovery
- **A cohort is a directory under `cohorts/`.** The directory name is the cohort name.
- `DataStore.cohorts` holds only cohorts with a readable `<NAME>.pedigree.tsv`, so the header dropdown, validation pages and search keep their existing meaning. Directories with no pedigree — or one that will not parse — land in `DataStore.incomplete_cohorts` as `CohortStub`s, carded on the home page but not navigable. `CohortStub.reason` is `"missing"` or `"unreadable"`; the two need different fixes and do not share a message.
- Discovery lives in **one** place, `_iter_cohort_dirs()`. It used to be hand-copied into `load()`, `take_snapshot()` and `reload()`, which drifted the moment the rules changed.
- `take_snapshot()` watches the pedigree **and** `<NAME>.params.yml` (the newer of the two mtimes, `0.0` when absent) and keys every cohort directory, so a directory appearing or gaining a pedigree registers as a change.
- `create_header(cohort_name)` must never be handed a name outside `store.cohorts` as the select's value — Quasar raises on it, which used to turn every unknown cohort in the URL into a 500 and left the pages' own "not found" branches unreachable.

### Cohort Run State
Two files may sit in a cohort directory, both written by the ghfc-ngs workflow: `<NAME>.params.yml` (parameters, including the `steps:` array) and `.ghfc-ngs.state.json` (the run record — schema documented in `COHORT_STATE.md` of bourgeron-lab/ghfc-ngs). The Parameters and Status tabs on the cohort page exist only when the params file does.

- **genetics-viz never measures completion.** The record already carries per-step `{done, total, pct}` with denominators from the pedigree. Re-deriving it would duplicate the pipeline's planning logic and drift from it, and it is unaffordable anyway: a cold existence check costs ~5.6 ms over the CIFS mount, so a 2700-sample cohort is about a minute of stat calls.
- **A null step is unmeasured, not 0%.** `ancestry` is null unless requested, `extractor` always is. Render "not measured", never an empty bar.
- **Check `schema_version` before trusting any field.** Anything but `1` is refused whole rather than half-interpreted.
- **A missing state file means "not run since run recording shipped"**, not "never run". The Status view then lists the params file's steps, all unmeasured.
- **`status: running` is not a liveness signal.** There is no heartbeat; a killed run leaves the record behind until the next run closes it out. `RunRecord.is_stale_running` flags one older than 24 h.
- **The `path` fields in the record are compute-cluster paths** (`/pasteur/helix/...`) and do not resolve locally, so drift detection hashes `cohort.pedigree_file` / `cohort.params_file`. Same trap as the ancestry `.qc.json` files.
- The record is deliberately **not** in the change-monitor snapshot: the pipeline rewrites it on every run and `reload()` re-parses every pedigree in the data directory. `read_cohort_state` caches on the file's own `(path, mtime, size)` and picks up rewrites by itself.
- `config/pipeline_steps.yaml` holds step labels, display order and status chip colours only. A step the pipeline reports but the file does not know is appended under its raw name.

### Cohort Page Tabs
The page has two levels of tabs, and they are separate mechanisms.

- **Page-level, full width**: `General` (the individuals table plus the right-hand panels that follow its filters) and `Workflow`. `Workflow` is *created only when the cohort has a params file* — absent, not disabled — and its own lazy loader is the only `async` one, because it warms the run-record and document caches off-thread before rendering.
- **Right-hand panel, inside `General`**: `Statistics` / `Ancestry` / `Polygenic Scores`, nested inside the General tab panel. Anything describing the pipeline run belongs on `Workflow` instead: it is a property of the cohort, not of the current filter, and does not fit a 400px column.
- Both levels use the same lazy pattern — a `{"loaded": False}` dict, a `@ui.refreshable` body with `_loading(...)`, and `ui.timer(0.1, load, once=True)` on tab change. The inner level keys `tab_state` and `lazy_tabs` by the exact tab label and `tab_state[e.args]` is an unguarded lookup, so a label must be added to both dicts or removed from both.
- `render_status_panel(cohort, title=None)` suppresses the panel's own heading, because the Workflow tab's expansion header already carries the title, status chip and step summary.

### ghfc-ngs Documentation
The Workflow tab shows the pipeline's own parameters reference beside the cohort's `params.yml`. It is **fetched at runtime** from the public raw URL (`utils/ghfc_ngs_docs.py`), not vendored: the document changes often and a committed copy would go stale silently. stdlib `urllib` with a 10s timeout, so no new dependency; cached one hour on success and one minute on failure; a failure is a degraded state the panel reports, never an exception.

- **The document is read outside its own repository**, so its repo-relative links (`../README.md`, `../nextflow.config`, `../params_example/`) are rewritten to absolute GitHub URLs. A trailing slash means a directory, which GitHub serves under `/tree/main/`, not `/blob/main/`. Rewritten links become raw HTML anchors with `target="_blank"`, and backticked link text has to become `<code>` by hand — markdown2 does not process markdown inside an inline HTML span.
- **Its 73 in-page anchors work only because of the `header-ids` markdown2 extra**, whose slugs happen to match GitHub's for every anchor this document uses. Do not drop it from `MARKDOWN_EXTRAS`.
- **An anchor click is intercepted.** Left to the browser it scrolls the panel *and* the whole page (measured: `window.scrollY` 0 → 995) and appends a fragment to the URL. The handler moves only the container's `scrollTop`; `scrollIntoView` is wrong here because it scrolls every scrollable ancestor.
- **Styles and scripts for a lazy tab must go through `ui.run_javascript`, not `ui.add_css` / `ui.add_head_html`.** The panel is built by a timer callback, long after the page head was sent, so a head injection is silently dropped — verified: the rule never reached the document. `components/tanstack_table.py` records the same trap.

### Flex Sizing Traps In Fixed-Width Cards
Several distinct mechanisms bit the cohort cards, all looking like "the content does not fit the card". Diagnose by measuring, never by eye: compare each descendant's `getBoundingClientRect().right` against the card's, and read an element's intrinsic floor by temporarily setting `style.width = 'min-content'`.

- **`.nicegui-card` and `.nicegui-column` both set `align-items: flex-start`**, so children are sized by their content instead of stretching. Two consequences: a long label inside a `ui.column` escapes its container and `truncate` has no width to act on (give it `w-full`, plus `truncate` or `break-all`); and a `ui.card_section` comes out narrower than its card (274px in a 320px card), so nothing inside lines up with the card's edges. Every card section needs `w-full`.
- **A flex item defaults to `min-width: auto`, refuses to shrink below its content, and that min-content floor propagates up the whole ancestor chain.** One `whitespace-nowrap` label inside the home card's Status expansion (`100% · 522/522 samples`, 121px) gave the expansion a 299px min-content width, forcing `q-card__section` to 331px inside a 320px card and clipping the phenotype table's last column. Fix: `min-w-0` on every flex item in the chain, `truncate` rather than `whitespace-nowrap`, and a compact rendering that drops the unit noun (`StepRow.detail(compact=True)`).
- **`nicegui-card` carries its own 16px padding on top of the 16px Quasar puts on each `q-card__section`** — 32px a side, leaving 256px of a 320px card. `p-0` on the card, padding only on the sections.
- **A `border-l-*` accent is inside the border box** and eats content width, so equal padding leaves the left inset wider than the right. The cohort cards use `pl-2 pr-3` against a `border-l-4`.

Card geometry for the home page lives in `_CARD_BASE` / `_SECTION_X` in `home.py`, with the reasoning inline. Net effect: 296px of usable content in a 320px card, insets symmetric, no overflow at any container width from 360px up.

Short test data hides all of this: the synthetic fixtures (`24/24`) fitted where the real counts (`522/522 samples`) did not.

### Pedigree Missing Values
The sentinel set `{"", "0", "-9"}` represents unknown/missing in pedigree fields (parent IDs, sex, phenotype). Defined as `_PED_MISSING` in `search.py` and handled in `models.py` via `treat_missing_as_null`.

### TanStack Table
- `DataTable` in `tanstack_table.py` bridges Python to vanilla JS TanStack Table.
- Per-client JS/CSS injection via `client._dt_scripts_injected` attribute (avoids global set memory leak).
- MutationObserver with 30-second timeout for container detection.

### Async
- Use `asyncio.to_thread()` to offload blocking I/O (e.g., `pl.read_csv()`) in NiceGUI page handlers.
- NiceGUI runs a single-threaded event loop; blocking calls freeze the UI for all users.

## Conventions

- Utility code shared across pages lives in `utils/`, not duplicated in page components.
- YAML configs live in `config/`, loaded by corresponding `utils/*.py` modules.
- All new ClinVar terms must be added to `config/clinvar_colors.yaml`.
- All new VEP consequences must be added to `config/vep_consequences.yaml`.
- Column display/group/sort/drop config goes in `config/column_names.yaml`.
