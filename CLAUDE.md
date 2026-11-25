# CLAUDE.md

Project-level instructions for Claude Code.

## Project Overview

NiceGUI-based web application for genetics cohort data visualization. Uses TanStack tables for dataframes, Polars for data processing, and ECharts for plots. Managed with `uv`.

Deliberate exceptions to "ECharts for plots": the cohort PGS violins in `cohort_pgs_tab.py` and the cohort ROH violins in `cohort_roh_tab.py` use Plotly. ECharts has no violin or box series, and both live as sibling tabs on the cohort page. The family report draws static SVG built in Python (`utils/report_svg.py`, and the pedigree through pypedigree in `utils/pedigree_svg.py`): it is printed, and a JS chart can miss the print.

## Commands

- **Run app**: `uv run genetics-viz /path/to/config.yaml` (add `--show` to open a browser tab)
- **Run tests**: `uv run pytest` (`-m "not slow"` skips the subprocess boot test)
- **Lint**: `uv run ruff check .`
- **Format**: `uv run ruff format .`
- **Compile check**: `uv run python -m py_compile src/genetics_viz/path/to/file.py`
- **Always use `uv run`** to execute Python commands, never bare `python`.
- **Always use `uv add`** to add dependencies, never bare `pip install`.

## Tests and CI

- `tests/datagen.py` builds a synthetic data directory and reference bundle in `tmp_path`. CI has no access to the data mount, so every test runs on it. When a new file format is read, extend `datagen` in the format the pipeline actually writes.
- `tests/conftest.py` resets all module-level state before and after each test: the config singleton, the data store registry, the stat-keyed caches, the GeneTrek snapshot, the preferences cache, the per-user gene scorers, the login throttle, the gene index manifests, the Vedette table cache, the UPD tables cache, the cohort SV calls cache, the family statuses cache and the pedigree consistency cache. A new module-level cache must be added to `_reset_globals`, or tests will depend on their order.
- An autouse fixture makes any `urlopen` fail, so no test can reach GitHub by accident.
- **UI tests** (`tests/ui/`) use NiceGUI's user simulation: no browser, timers run for real.
  - `tests/ui/app_main.py` purges `genetics_viz` from `sys.modules` and boots the app again for each test. NiceGUI's per-test reset removes every route, so without the purge nothing would be registered from the second test on.
  - The `user` fixture points `Storage.path` at `tmp_path`. The reset otherwise deletes the real `.nicegui/` sessions.
  - It fails a test on any ERROR log, which is how NiceGUI reports an exception in a lazy tab.
- **UI test gotchas:**
  - **Tabs:** switch them with `ui_helpers.open_tab`, which sends `update:modelValue`. NiceGUI stores listeners under the camelCase name, and clicking a `ui.tab` does nothing in simulation. Tab bars carry `.mark("family-tabs")` and `"cohort-tabs"` for this.
  - **Table and chart contents** render in JS and are invisible to `should_see`. Assert on what the Python side passes: `ui_helpers.capture_tables` records `DataTable` rows, and `ui.echart(...).options` holds the chart data.
- **CI:**
  - `ci.yml` runs on PRs and on main: ruff, pytest on 3.12/3.13 with the lock, and a `package` job that tests the built wheel installed without the lock.
  - `upstream.yml` runs that unlocked install weekly and opens an issue when it fails. Production (`uvx`) ignores `uv.lock`, so an upstream release can break it with no commit here.
  - `publish.yml` runs `ci.yml` first, checks the tag against the version and the CHANGELOG, publishes to PyPI, then creates the GitHub release from the CHANGELOG section (`.github/scripts/release_notes.py`; title = the annotated tag's message). The release steps are in `DEPLOY.md`.
- **Known issues are `xfail(strict=True)`:** once one is fixed, its test starts passing and fails the run until the marker is removed.

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
    filters.py                    # Reusable filter menus (validation statuses, consequences, genesets), unknown-term warning
    sample_counts.py              # Stats dialogs: variants per member + cohort comparison
    gene_prefs_editor.py          # Profile card: GeneTrek tiers per geneset, live preview
    header.py                     # Page header: menus, data dir dropdown, gene lookup, Help, user menu
    validation_loader.py          # Validation TSV loading and badge logic
    diagnostic_loader.py          # Diagnostic TSV loading, saving (role-guarded), ignore toggle, family diagnostics
    notes_loader.py               # Family notes: load, save, delete (role-guarded)
    waves_loader.py               # WAVES validations: load, save (role-guarded), score colours
    icons.py                      # Validation and diagnostic status icons
    variant_dialog.py             # IGV variant inspection dialog
    sample_dialog.py              # Sample dialog: IGV.js bedgraph and CRAM tracks for one sample
    search_stats.py               # Search page Stats dialog: chromosomes, consequences, validations, ideogram
    cell_line_note.py             # The caveat on a cell-line sample's UPD, CNV, de novo and diagnostic reviews
    sv_dialog.py                  # SV inspection dialog: orchestrates the panels below (SVViewState)
    sv_details.py                 # SV dialog: details card, suggested curated coordinates
    sv_igv_views.py               # SV dialog: igv.js views (bedgraph, VAF, read level), SVViewState
    sv_validation_panel.py        # SV dialog: validation form (curated position) and history
    sv_cn_panel.py                # SV dialog: copy-number panel (sample ÷ control panel)
    sv_read_panel.py              # SV dialog: samplot-like read-evidence panel
    vedette_details.py            # SV dialog for a Vedette call: details card, own call, SVSource
    offthread_panel.py            # Container loaded off-thread, superseded by newer refreshes
    download_button.py            # Tab Download button: one TSV, or a menu + "All tables (.zip)"
    tab_groups.py                 # Tab rows in captioned groups, no icons (cohort and family pages)
    family_tabs.py                # Family pages' tab bar: groups, link value, label, tooltip per tab, older link values
    denovo_evidence.py            # Variant dialog: the denovo step's evidence card for one candidate
    tab_badge.py                  # Counts on tabs: the family page's review badges (TabBadge), the cohort page's (TabPills)
    page_link.py                  # Shareable links: tab and filter state in the URL query (PageLink)
    diagnostic_dialog.py          # Diagnostic review dialog
    diagnostics_panel.py          # Family Diagnostics panel (cohort + standalone family pages)
    upd_dialog.py                 # UPD review: evidence, chromosome view, validation, diagnostic (UPD, ROH, Sample QC tabs)
    aneuploidy_dialog.py          # Aneuploidy review: SVs tab card, dialog (evidence, validation, diagnostic), also from a KARYOTYPE flag
    pedigree_tree.py              # Member card's Tree view: pypedigree drawing, selection faded, QC outline, Enlarge
    pedigree_editor.py            # Pedigree editor dialog (Sample QC tab, curators): trees now/after, rows, changes, history; status note; cohort chooser when cohorts differ
    family_report_dialog.py       # Family report: "See report"/"PDF" buttons, dialog, conclusion, print
    family_status.py              # Family status chip + curator dialog (family pages), status bar and yield
    pedigree_consistency.py       # Pedigrees across cohorts: family banner, cohort/data-directory lists, explain dialog
    run_request.py                # Run requests: dialog, request button (or the user's open one), Workflow list
  config/                         # YAML config files (loaded at module level)
    ancestry.yaml                 # Ancestry region colors + PGS z-score color scale
    genetrek.yaml                 # GeneTrek refresh period, tier colours, default profile
    gene_index.yaml               # Gene index refresh period, Parquet row group size
    inheritance.yaml              # MOI flag: GeneTrek mode strings -> codes, verdict colours, PAR
    roh.yaml                      # ROH grades/wording, pattern + class colours, KING degrees
    denovo.yaml                   # denovo step: flag wording/colours/tooltips, strict, origin, QC verdicts, literature references
    upd.yaml                      # upd step: confidence/origin/mode wording, member statuses (role), window colours
    sample_qc.yaml                # Sample QC flag wording (curator/reader), categories, plot lines; also genomeqc's flags
    genome_qc.yaml                # genomeqc step: metric names/formats, seq_class wording, mirrored thresholds
    str.yaml                      # STR class/status colours, support codes, REViewer filename
    report.yaml                   # Family report: wording, PGS traits, aneuploidy cut-off, methods
    family_status.yaml            # Family statuses (order, colours, yield role), solving verdicts, reasons
    pedigree_consistency.yaml     # Cross-cohort pedigree findings: severity, label, short badge, text
    clinvar_colors.yaml           # ClinVar significance term -> color
    cnv.yaml                      # CN mask/display, single-vs-split read view, read-evidence caps/colours
    vedette.yaml                  # Vedette: impact levels, inheritance codes, colours, presets
    column_names.yaml             # Column display names, groups, sorting, drop flags
    cytobands_hg38.tsv            # Cytoband data for ideogram rendering
    continuous_scores.yaml        # Continuous score color ranges
    pipeline_steps.yaml           # ghfc-ngs step labels/order + run status chips
    vep_consequences.yaml         # VEP consequence terms, impacts, colors
    wisecondorx_thresholds.yaml   # WisecondorX call thresholds (robust, permissive; loss and gain)
    validation_badges.yaml        # Validation status and inheritance badge colours and labels
    diagnostic_badges.yaml        # Diagnostic verdict badge colours and labels
    view_presets.yaml             # Column visibility presets
  pages/
    login.py                      # Login page (no auth required)
    cohort/pedigree_consistency.py # /pedigree-consistency: the data directory's cross-cohort findings (curators)
    profile.py                    # User profile (role, password, gene badge tiers)
    search.py                     # Cohort-wide variant search: route and search steps
    search_form.py                # Search page form: source cards, Individuals panel (SearchForm)
    search_results.py             # Search page results: table, presets, Save/Download, row actions
    gene.py                       # /gene/<symbol>, /variant/<variant>, /goto/<n> (data dir switch)
    admin/
      directories.py              # Admin: manage data directories
      users.py                    # Admin: manage users
    cohort/
      home.py                     # Home page with cohort cards
      cohort.py                   # Cohort overview page
      family.py                   # Family detail page with tabs
      standalone_family.py        # /family/<FID>: the family page outside any cohort
      components/
        ancestry_tab.py           # Ancestry PCA tab (per-family)
        pgs_tab.py                # Polygenic scores tab (per-family)
        cohort_ancestry_tab.py    # Ancestry PCA tab (cohort-wide)
        cohort_pgs_tab.py         # Polygenic score violins (cohort-wide)
        roh_tab.py                # Runs of homozygosity tab (per-family)
        upd_tab.py                # UPD tab (per-family): candidates from every source, genome overview, members
        str_tab.py                # STR tab (per-family): EH sub-tab + locus dialog
        cohort_roh_tab.py         # Runs of homozygosity tab + islands dialog (cohort-wide)
        cohort_upd_tab.py         # UPD tab (cohort-wide) + the individuals table's UPD column
        cohort_str_tab.py         # STR tab (cohort-wide): loci counts + flagged genotypes
        cohort_denovo_tab.py      # De novo QC tab (cohort-wide): literature, rate/shares charts, children, flags, curation check
        cohort_snvs_tab.py        # SNVs/Indels tab (cohort-wide): de novo candidates, genes hit in 2+ families
        cohort_svs_tab.py         # SVs tab (cohort-wide): WisecondorX calls, whole-chromosome events, recurrent regions
        sample_qc_tab.py          # Sample QC tab (per-family): issues, kinship, sex, quality, DNA and library
    genome_qc_section.py      # Sample QC tabs' DNA and library section (genomeqc), cell-line badge and callout
        cohort_sample_qc_tab.py   # Sample QC tab (cohort-wide) + cross-family network dialog
        docs_panel.py             # ghfc-ngs params reference, fetched at runtime
        params_panel.py           # Workflow params file, shown verbatim
        status_panel.py           # Pipeline run status (Workflow tab + home card)
        wombat_tab.py             # SNVs/Indels tab (per-family): De novo sub-tab, then the Wombat configs
        wombat_denovo.py          # SNVs/Indels tab: de novo flag filter, per-child summary, curation banner, header QC badge
        wombat_rows.py            # Wombat tab data: result files, table rows, columns
        wombat_filters.py         # Wombat tab filters: WombatFilters, panel, apply_filters
        wombat_stats.py           # Wombat tab Stats dialog
        svs_tab.py                # SV analysis tab (per-family): WisecondorX + Vedette sub-tabs
        vedette_tab.py            # SVs tab, Vedette sub-tab: carrier rows, presets, filters
    validation/
      file.py                     # Per-file validation page
      all.py                      # Aggregated validation page
      statistics.py               # Validation statistics page
      waves.py                    # WAVES validation page (one dialog per sample)
    diagnostic/
      all.py                      # Aggregated diagnostic page
      statistics.py               # Diagnostic statistics page
  utils/
    ancestry.py                   # Ancestry/PGS file discovery, reference panel, colors
    ancestry_plots.py             # Shared ECharts builders for the PCA scatters
    roh.py                        # ROH discovery, reference tier, grades, kinship checks
    roh_plots.py                  # ECharts/SVG builders for the ROH tabs
    str_eh.py                     # ExpansionHunter discovery, reading, cohort cache, locus pivot
    str_plots.py                  # ECharts builders for the STR views
    str_validation.py             # validations/strs.tsv writer (role-guarded) and reader
    sample_qc.py                  # sampleqc (somalier) discovery, reading, flag wording, shaping
    genome_qc.py                  # genomeqc step: discovery (by #step), reading, flags as FlagRows, metric wording
    genome_qc_plots.py            # ECharts builder for the genome QC metrics, one row per seq_class
    sample_qc_acks.py             # validations/sampleqc.tsv acknowledgements (role-guarded) + UPD-derived ones (Acks)
    sample_qc_export.py           # Sample QC exports: flags (with acks), samples, pairs; role-worded
    export.py                     # TSV/zip helpers, export filenames, person columns (FID, IID, sex, phenotype)
    sample_qc_plots.py            # ECharts builders for the Sample QC views
    validation_files.py           # The validation kinds (SNV/SV/STR/UPD) and their files
    validation_status.py          # Validation statuses: consensus of a row's validations, filter, icons/colours
    variant_rows.py               # Wombat rows: VEP/ClinVar/score badges, LCR/gnomAD/MNV/consequence filters
    diagnostic_files.py           # The diagnostic kinds (SNV/SV/STR/UPD) and their files
    diagnostic_labels.py          # Diagnostics panel wording: STR repeat counts, UPD origin
    denovo.py                     # denovo step: discovery, reading, evidence on Wombat rows, strict/flag filter, curation check, cohort summaries
    denovo_plots.py               # ECharts builders for the cohort De novo tab
    upd.py                        # validations/upds.tsv (role-guarded), UPD key, evidence, candidates merged from every source
    aneuploidy.py                 # Whole-chromosome events (WisecondorX), karyotypes, validations/aneuploidies.tsv (role-guarded), candidates
    upd_calls.py                  # upd step tables: discovery (newest label), reading, windows, wording, exports
    upd_plots.py                  # UPD chromosome view (ECharts), genome overview (SVG), cohort chart
    review_counts.py              # What the family and cohort tabs' badges count (QC, de novo, STR, UPD, filter, run)
    family_rows.py                # Cohort Families tab: structure, facts read off-thread, rows, columns, filter; arrows' order
    cohort_snvs.py                # Cohort SNVs/Indels tab: candidate rows, filter, recurrent genes
    cohort_svs.py                 # Cohort SVs tab: lazy WisecondorX scan, call rows, aneuploidies, recurrent regions
    family_report.py              # Family report sections (pure, off-thread): pedigree … diagnostics, methods
    report_html.py                # Family report as one self-contained HTML document (Jinja2)
    report_svg.py                 # Static SVG plots for the report: PGS normal curves, PCA scatter
    pedigree_svg.py               # Pedigree drawing (pypedigree): members -> PED text -> static SVG
    pedigree_history.py           # A cohort's pedigree edits as one HTML document: trees before/after per family
    pedigree_edit.py              # Curator pedigree edits: draft, checks, cohort + family files, history, steps awaiting a rerun
    pedigree_qc.py                # Sample QC against a pedigree draft (mirrors sampleqc_classify), editor suggestions, merge candidates
    pedigree_merge.py             # Merging families under a new or existing FID; curation re-keyed; orphaned directories
    report_summaries.py           # reports/family_summaries.tsv: curator conclusions (role-guarded)
    family_status.py              # reports/family_status.tsv (role-guarded), reasons to look again, yield
    pedigree_consistency.py       # Cohort pedigrees compared (same FID, aliases by barcode), acks (role-guarded)
    run_requests.py               # Requests for a ghfc-ngs run (cohorts/<C>/.ghfc-ngs.requests/), role-guarded; outcome from the run record
    ideogram.py                   # Chromosome backbone SVG (bands, Mb grid)
    cnv_coverage.py               # Coverage bigWigs, CNV control panel, fit table, copy number
    cnv_plots.py                  # ECharts builder for the copy-number panel
    sv_view.py                    # One read window or one per breakpoint, by SV type and size
    vedette.py                    # Vedette SV calls: reading, carrier rows, filters/presets, own call
    sv_reads.py                   # pysam: read pairs/split reads classed by orientation
    sv_breakpoints.py             # Suggested breakpoints: clips/splits, pair interval, CN check
    sv_read_plots.py              # ECharts builder for the read-evidence panel
    auth.py                       # Auth helpers: check_auth, is_authenticated, can_write, can_download; role and grants read live
    access.py                     # What each role may fetch from /data-N (viewer files for readers and students)
    login_throttle.py             # Login backoff per client address and per username (in memory)
    clinvar.py                    # ClinVar color/display utilities (from YAML)
    cohort_state.py               # .ghfc-ngs.state.json parsing, drift, step rows
    cytobands.py                  # Cytoband, chromosome, ideogram constants
    wisecondorx.py                # WisecondorX BED parsing, CNV type and class, thresholds (past_permissive), gene overlap
    locus.py                      # Locus query parsing (chr:start-end, genes) and frame filters
    change_monitor.py             # Background polling of the data directories for changes (live refresh)
    data_availability.py          # Which files a sample or family has (CRAM, bedgraph, ...)
    data.py                       # Multi-store registry, per-user data dir selection
    gene_scoring.py               # Per-user tier scorer + the shared gene badge builders
    genesets.py                   # Geneset filter sources (tiers, GeneTrek, local TSVs) + matchers
    genetrek.py                   # GeneTrek snapshot: API fetch, disk cache, index, links
    inheritance.py                # Trio transmission + MOI consistency verdicts
    user_preferences.py           # <config>.preferences.yaml: per-user prefs, locked atomic writes
    ghfc_ngs_docs.py              # Upstream params.md fetch + link rewriting
    pipeline_params.py            # <NAME>.params.yml discovery and parsing
    score_colors.py               # Continuous score color ranges
    vep.py                        # VEP consequence utilities (from YAML)
    pedigree_labels.py            # Sex/phenotype code normalisation and labels
    pedigree.py                   # Per-family pedigree TSV parser (standalone family page)
    validation_badges.py          # Validation badge colours and labels (from YAML)
    diagnostic_badges.py          # Diagnostic badge colours and labels (from YAML)
    sharding.py                   # Two-level sharded directory path resolution
    tsv.py                        # read_tsv_or_none (string_columns), clean_tsv_field, find_labels (step), inputs_step, append_tsv_row(s), rewrite_tsv_rows, set_ignore, collect_streaming
    wombat.py                     # aggregate_by_variant (Wombat tab, search); family result files, DENOVO_CONFIG
    gene_index.py                 # Cross-cohort gene index: cohort tables -> local Parquet, lookups
    gene_rows.py                  # Gene/variant page rows: carriers, curation joined, /goto links
    search_rows.py                # Search page logic without UI: sources, filters, result rows, columns
    sample_counts.py              # Cohort Wombat/WisecondorX files, per-sample counts, outlier flags
    view_presets.py               # View preset loading with reload support
    column_names.py               # Column name/group/sorting/schema utilities
  templates/
    family_report.html            # The family report document (inline CSS, A4 print rules)
    pedigree_history.html         # The pedigree edit history document (same print rules)
  documentation/
    snvs_validation_guide.md      # The variant dialog's ⓘ guide (inside the package, so the wheel ships it)
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
- **Roles** (`config_model.ROLES`, least rights first): `student` (a reader without downloads), `reader`, `curator`, `administrator`.
- `get_current_user()` reads the session. `get_current_role()` and `granted_data_dirs()` read the **config's** entry for that user (`current_account()`), so an admin's change applies at the next request; `is_authenticated()` is false once the user is removed. Without a loaded config (unit tests), the session alone decides.
- `can_write()` returns `True` if role is `"curator"` or `"administrator"`; `can_download()` for every role but `student` (and an unknown one); `is_admin()` for `"administrator"`.
- **Students get no Download button**: `download_button` returns `None` and draws nothing for them, its senders check again, and Search's Save menu leaves out Download TSV.
- **Passwords are salted scrypt** (`scrypt$N$r$p$salt$hash`, `config_model.hash_password`). Legacy unsalted SHA-512 hex digests still verify; the login rehashes them (`needs_rehash`) and saves the config, so no reset is needed. The test config writes SHA-512 on purpose, so every UI login exercises the migration.
- **The login backs off** (`utils/login_throttle.py`): 5 free failures per username, 20 per client address, then a doubling wait capped at 15 min. Unknown usernames are counted and hashed like real ones. The address is what uvicorn sees: behind a proxy it is the proxy's unless uvicorn trusts its forwarded headers (`FORWARDED_ALLOW_IPS`). The check, the hash and the recorded failure run on the event loop without an await, so parallel attempts cannot all pass the check.
- **Saves record `get_current_user()`**, never the readonly User box: the browser can still send any value for it.
- Write operations are gated with `can_write()` both at UI level (hide or disable the control) and in the writer itself, which raises `PermissionError` for the caller to notify: validations and their Ignore switch (`update_validation_ignore_status`), diagnostics, notes, wave scores, STR/UPD validations, acknowledgements, report summaries, and the trash of a file to validate.

### Multi-Data-Directory
- **A user may be granted some directories only** (`UserConfig.data_directories`, set on Admin → Users; empty means all). `user_data_dirs()` lists the user's registered ones; outside a request (poller, index) it lists them all. Use it, not `get_all_data_stores()`, for anything shown to a user; a thread gets the list taken on the event loop (the gene pages pass it to `data_dir_infos(paths)` and the index lookups' `data_dirs`).
- `get_data_store()` and `get_static_prefix()` follow `current_data_dir()`: the session's `data_dir` if still granted, else the default directory if granted, else the user's first. `check_auth()` logs out a user left with none, and the login refuses them.
- Data directories are served centrally through `data.mount_data_files()` (from `app.py` and `add_data_store`) — NOT in individual pages, and never with `nicegui_app.add_static_files()`, whose route checks nothing. The route answers 401 without a logged-in session, 404 once the directory is removed (or for a `..` path), 403 for a directory not granted, supports Range requests (igv.js) and sends `Cache-Control: private`.
- **By role, the route serves readers and students only the files the pages load** (`utils/access.py`: CRAM/CRAI, bedgraphs and their `.tbi`, coverage bigWigs, REViewer images under `samples/`); curators and administrators get every file. A new browser-loaded file type must be added to `_VIEWER_FILES`, or readers get 403. Links to raw files (the somalier report) are shown to curators only.

### Sharded Directories
- `samples/` and `families/` directories may use two-level sharding: `samples/<shard1>/<shard2>/<id>/`.
- Shard keys: strip `-`, `.`, `_` from entity ID, then shard1 = last char (uppercased), shard2 = second-to-last char (uppercased).
- Auto-detected per data directory (if all children of `samples/` are single-char dirs → sharded). Cached.
- Use `get_sample_path(data_dir, id)` / `get_family_path(data_dir, id)` for filesystem paths.
- Use `get_sample_url(data_dir, id)` / `get_family_url(data_dir, id)` for URL segments.
- All functions in `utils/sharding.py`. Never construct `samples/{id}` or `families/{id}` paths manually.

### Ancestry / Polygenic Scores
- Per-family results live in `<family>/ancestry/`, named `<family_id>.apgs_b<bundle>_<filters>.<kind>.tsv` (kinds: `pcs`, `ancestry`, `pgs_zscore`, and `Q`, the admixture proportions: `IID`, `snp_coverage`, one column per reference component; only the export reads `Q`). The bundle version contains dots, so discovery regexes need a lazy `(.+?)` capture, not `([^.]+)`.
- The PCs are projections onto a **versioned** reference panel, so the panel a family is plotted against is the one named in its own filename tag. Config supplies only the root (`ancestry_reference_dir`); the `v<bundle>/bundle/labels/reference_pcs.tsv.gz` path below it is derived per family. Never read the reference root from the `.qc.json` files — they record the compute-cluster path.
- A missing reference panel is a degraded state, not an error: the Ancestry tab still plots the family's own samples.
- `utils/ancestry.py` holds discovery, probes, the cached reference panel, region colours, and the shared `NULL_VALUES` — the PGS z-score file writes the literal string `nan`, which must be nulled or polars types the column as a string and sorting breaks.
- **Cohort-level results** live in `cohorts/<name>/ancestry/<name>.apgs_b<bundle>_<tag>.<kind>.tsv` — same naming, cohort name as the prefix, plus an extra `family_id` column. `find_ancestry_files_in(directory, prefix)` serves both layouts; `find_cohort_ancestry_files(cohort)` and `find_ancestry_files(data_dir, family_id)` are the wrappers. Not every cohort has them.
- **Several filter-tag variants of one bundle coexist** (`dp6gq15`, `dp10gq20` — minimum depth and GQ). `_tag_strictness` ranks them so the strictest wins by default; a plain string sort picks `dp6gq15` because "6" > "1" at the third character. Pages with more than one variant show a selector.
- Plot builders shared by the family and cohort tabs live in `utils/ancestry_plots.py`. They are pure functions returning ECharts option/series dicts, so neither page component has to import the other's privates.
- Normalise pedigree sex/phenotype through `utils/pedigree_labels.py` before grouping — some cohorts write the codes float-stringified (`-9.0`, `1.0`), which a bare comparison against `"-9"` misses.

### Runs of Homozygosity (ROH)
- The ghfc-ngs `roh` step writes `<family>/roh/<FID>.apgs_b<bundle>_<tag>.{froh,roh,roh_genes,kinship}.tsv` and, per cohort, `cohorts/<name>/roh/<name>.…{froh,roh,roh_genes,kinship,roh_windows}.tsv` (plus `family_id`). Same naming as the ancestry files, so discovery is `find_ancestry_files_in(..., kinds=ROH_KINDS)`; `roh_genes` is optional (needs `annotation_gencode`), `kinship` is absent for a singleton.
- **Read through `read_roh_table`.** It forces the ID columns to strings. Left to inference, `family_id` comes out as an integer from numeric first rows and the read then fails on `AU-ANG-FER-MPX3`. It also strips the `chr` prefix, because the reference hotspots and the cytobands carry none.
- **Reference first, fixed thresholds second.** Each person is compared with unrelated reference genomes of their `ref_group` (a region, or `ALL` = pooled), read from `v<bundle>/bundle/roh/reference_froh.tsv.gz` under `ancestry_reference_dir`, next to `hotspots.tsv.gz`. The fixed 0.03 line is crossed by 150 of 682 unrelated CSA reference genomes, so a threshold-only display would mislabel whole populations. A missing tier is a degraded state, as for the PCA panel.
- **Thresholds and wording live only in `config/roh.yaml`.** They mirror ancestry-pgs constants (`roh.py`), which assign `froh_category` and `homozygosity_pattern` themselves; the app only words and draws them.
- **Wording is role-dependent (`can_write()`).** Viewers get neutral grades and only pedigree-declared kinship pairs. Curators also get the relationship reading, the parents' kinship / expected F, and every pair. The raw files are not served to readers either (`utils/access.py`).
- The cohort tab reads only `froh` (and `roh_windows` on demand, off-thread). The cohort `roh` and `roh_genes` tables are 31 and 93 MB.
- Some pedigrees repeat a row, and a sample can sit in two families: key cohort rows by `(family_id, IID)`.
- `DataTable` does not reorder its initial rows to match `initial_sorting` — sort the rows in Python as well.

### Short Tandem Repeats (STR)
- The ghfc-ngs `expansionhunter` step writes under `str/expansionhunter/` at each level, every file named by the catalog label (`expansionhunter_catalog_name`, e.g. `jmsa50`):
  - sample: `<id>.<L>.{vcf,json}` and `_realigned.bam`;
  - family: `<FID>.<L>.genotypes.tsv`, one row per locus × member;
  - cohort: `<NAME>.<L>.{genotypes,loci,flags}.tsv`.

  Discovery is `find_eh_labels(dir, prefix)`. It is not the ancestry naming: there is no `apgs_b`, and the whole middle of the name is the label.
- **Read through `read_eh_table`.** `.` is the null everywhere. The ID columns and `genotype` are forced to strings: `family_id` is numeric-looking (`88181200007`), and a male's chrX `genotype` is `0`/`1`.
- `allele1` is the shorter allele. A male's chrX locus has one allele, and `allele2` is null.
- **Founders are members with no *genotyped* parent** (`father` and `mother` both null). The pipeline fills the parent columns only for parents genotyped in the family. Cohort statistics and `locus_distribution` count founders only.
- **Classes, statuses, confidences and severities are the pipeline's.** The app only colours and words them, from `config/str.yaml`. `str_plots.allele_class` mirrors the pipeline's order (pathogenic, premutation, intermediate, normal, above_normal) only to colour histogram bars.
- **Flagged means `severity` high or medium, not a status.** `eh_family_table` ranks each row: high (pathogenic), medium (premutation), info (intermediate, carrier, a `report: risk` locus), one rank lower when an allele past normal has `allele*_confidence` `low` for the member's `read_length`. `flags.tsv`, REViewer, the EH tab's count, the cohort panel and "common in cohort" all follow it. Go through `is_flagged` / `flagged_expr` / `loci_flagged`: a table written before severities has no `severity` column (nor `n_high` in `loci.tsv`) and falls back to the old statuses (pathogenic, carrier, intermediate) until the pipeline rebuilds its family. A low-confidence allele's badge ends in `?`.
- **Member column ids are `m<index>`, not sample IDs.** TanStack reads a dot in an accessor key as a path, and the index into the family's member list keeps the ids stable when the member selection changes.
- **Cohort context is optional and loads off-thread** (`load_cohort_context`, cached on `(path, mtime, size)`). The EH table is built only once the context is in. Built earlier, it would be replaced before the browser mounted it, and its container observer would time out.
- **The cohort panel (`cohort_str_tab.py`) counts only the filtered individuals.** Rows are restricted by `restrict_to`, a semi-join on `(family_id, sample)`, since a sample can sit in two families. The founders' statistics and the "common" flag stay whole-cohort, from `loci.tsv`. Flagged genotypes come from the genotype table already in memory, not `flags.tsv`.
- **STR validations live in `validations/strs.tsv`**, keyed `<chrom>:<pos>-<end>:<LOCUS>`. The catalog label is not in the key, so a validation survives a new catalog. The file has the SNV/SV header plus `CuratedRepeats`, so `load_validation_map` and the badges work unchanged. Writes go through `save_str_validation`, which checks `can_write()` itself.
- **STR diagnostics live in `diagnostics/strs.tsv`**, with the SNV/SV diagnostic header and the same key as STR validations. The locus dialog opens the shared `show_diagnostic_dialog` from an anchor *outside* its refreshable section, because an ignore toggle refreshes that section and would otherwise delete the dialog that sent it. The family table is redrawn when the locus dialog closes, not under it.
- **REViewer images** are written by the pipeline's `EXPANSIONHUNTER_REVIEWER` step to the sample's `str/expansionhunter/reviewer/` (`reviewer_subdir` and `reviewer_filename` in `str.yaml`). Things to know:
  - Only flagged genotypes (pathogenic, carrier, intermediate) are drawn.
  - `<sample>.<label>.metrics.tsv` is written for every reviewed sample and marks it reviewed.
  - Images are named by the catalog's **LocusId**, while the table names VariantIds: `FXN_A` is drawn as `FXN`. The catalog is not in the data directory, so `reviewer_path` falls back to the longest drawn locus the variant extends with `_`.
  - The older `str/reviewer/` directory comes from a separate STRchive run with other locus IDs, and is deliberately not shown.
  - Images are shown with `loading="lazy"`, because each one runs to a megabyte.

### Sample QC (somalier)
- The ghfc-ngs `sampleqc` step writes `<family>/qc/<FID>.<label>.{samples,pairs,flags,inputs}.tsv` and `cohorts/<NAME>/qc/<NAME>.<label>.{samples,pairs,flags,inputs}.tsv` plus `.html` (somalier's relate report, linked from the cohort tab). `<label>` is `sampleqc_label` (`somalier_hg38`): the STR naming, so discovery is `tsv.find_labels(dir, prefix, "samples")`. The contract is `modules/sampleqc/scripts/sampleqc_classify` in ghfc-ngs.
- **Read through `read_qc_table`.** `.` is null; IDs are forced to strings. A stub run writes cut-down headers, so missing columns are added as nulls.
- **`genomeqc` shares the `qc/` directories and the file names** (`<prefix>.gqc_v1.{samples,flags,inputs}.tsv`). `find_qc_labels` passes `step="sampleqc"` to `tsv.find_labels`, which drops a label whose `inputs.tsv` names another `#step` (one with no snapshot is kept). Without it the newer `gqc_v1` was read as Sample QC: no flags, no pairs. `datagen.add_genomeqc_tables` (opt-in) writes the step's tables in the real format, newer than sampleqc's.
- **The family files never hold cross-family findings.** The cohort tables are the family ones plus the cross-family pairs *observed related* (unrelated cross-family pairs are not written) and their `level: cohort` flags (`CRYPTIC_RELATED`, `DUPLICATE_OR_MZ`, `SAMPLE_SWAP_CANDIDATE`). The family tab reads the cohort tables off-thread for them; `family_flags_with_cohort` adds only the `cohort`-level rows, since the cohort file repeats the family flags.
- **A cohort pair flag's `family_id` is `"FA,FB"`**: `sample_a` is in FA, `sample_b` in FB. Key samples by `FlagRow.keys` (`(family_id, sample)`), never by sample alone.
- **The pipeline classifies everything** (relationships, karyotypes, flags, severities); the app only words and draws them from `config/sample_qc.yaml`. `KARYOTYPE_<k>` is matched by prefix; an unknown code is shown raw. The plot lines mirror the `sampleqc_classify` constants (CHARR: medium from `charr_warn` 0.02, high from `charr_high` 0.05), and the KING cut-offs come from `roh.yaml`'s `kinship_degrees`.
- **`contamination_charr` is the pipeline's own CHARR (Lu et al. 2023)**, from the extracts; somalier's estimate is `somalier_charr`, kept for reference (several times higher, tracks depth).
- **A `mismatch` pair is not always a contradiction.** The pipeline marks a pair `mismatch` whenever it has a flag, and several are `info`: `CONSANGUINITY`, `RELATED_BEYOND_PEDIGREE` (closer than a pedigree with missing parents says), `MZ_TWIN_CANDIDATE`, and a cross-family `CRYPTIC_RELATED` at the third degree or without the IBS0 of relatives. `pair_severities(flags)` (keyed by `PairRow.flag_key`) lets the scatter draw them as grey diamonds, the network as dashed edges, and the pair tables as "Mismatch (info)".
- **Wording is role-dependent, as for ROH**: viewers get `reader_label`, no evidence, an "Atypical" karyotype, and only pedigree-declared pairs (no parents' pair, no cross-family pairs in the tables). Presentation only.
- **Hubs**: a sample with `HUB_MINIMUM` (5) or more cross-family relatives is outlined in the network. Third-degree pairs across families are mostly shared ancestry or admixture (KING reaches 0.08 for admixed unrelated pairs). The pipeline already leaves out the cross-family pairs of a sample with CHARR ≥ 0.05 and flags it `RELATEDNESS_UNRELIABLE` instead.
- Elsewhere: the family header badge and member icons (family files only, read at page build), the cohort table's `Sample QC` column and filter, the home card's QC figure (off-thread, individuals by worst severity), and a warning above Wombat configs matching `de_?novo` when `denovo_risk` flags hit the family.
- **A wrong relationship holds the de novo table back** (`denovo_gate`, a subset of `denovo_risk`: `NOT_PARENT`, `SAMPLE_UNRELATED_TO_FAMILY`, `SAMPLE_SWAP_CANDIDATE`, `MENDELIAN_INCONSISTENT`). Inherited variants then read as de novo by the ten thousand. The sub-tab is only counted (`count_rows`, no parse), shows the warning and **Load anyway** (`denovo-load-anyway`), which reads and builds the rows in `asyncio.to_thread` (`_prepare` with the scorer taken on the loop). An explained flag does not count; a confirmed one does.
- **Curators acknowledge recurring flags** (MZ twins, declared consanguinity, adoption, tolerated contamination) in `validations/sampleqc.tsv`, through `utils/sample_qc_acks.py`:
  - Keyed by `FlagRow.ack_key`: code, the `family_id` cell (`"FA,FB"` for a cross-family pair), `sample_a`, `sample_b` (empty for one sample). No sites label, so a rerun keeps it. The family and cohort files write a family flag alike, so one key covers both pages. A UPD validation also acknowledges a `UPD_CANDIDATE` flag (see Uniparental Disomy); `load_acks` returns an `Acks` object, looked up by flag.
  - Append-only (`tsv.append_tsv_row`, under `flock`): status `explained` or `confirmed`, a note, user, timestamp. `withdrawn` takes it back; the last row of a key wins. Not a variant validation, so not in `VALIDATION_KINDS`. `save_flag_ack` / `withdraw_flag_ack` check `can_write()` themselves.
  - Acknowledged flags leave the family header badge, the member icons (grey `task_alt` when nothing else is left), the cohort column (one grey "Acknowledged" badge and filter option, not "Any issue"), the cohort panel counts, the home figure, and the de novo warning. The warning keeps a `confirmed` flag: the calls stay unreliable.
  - The family Issues table lists them greyed and last, with the review and note ("Show acknowledged" hides them). The cohort panel leaves them out unless "Include acknowledged" is on.
  - A review redraws the header badge and member icons through `render_sample_qc_tab(on_acks_changed=...)`. The dialog hangs from an anchor outside the refreshable section, as the STR locus dialog does.
  - **Readers see the status, and the note only for categories with `reader_note: true`** in `sample_qc.yaml` (sex, mendelian, sample quality). Elsewhere a note such as "known cousins" would name what `reader_label` hides.

### Genome QC (ghfc-ngs `genomeqc`)
- **What kind of DNA and library each sample is**, which `sampleqc` cannot say: a genome from an EBV-transformed cell line has the person's genotypes (every relatedness, sex and Mendelian check passes) but its clonal somatic mutations read as de novo calls, and culture can add CNVs and isodisomies. The contract is `modules/genomeqc/scripts/genomeqc_tables` in ghfc-ngs; the vault notes are `ghfc-ngs/Genome QC step.md` and `Cell-line and legacy genomes.md`.
- **Files** (label `genomeqc_label`, `gqc_v1`):
  - sample: `samples/…/<ID>/qc/genomeqc/<ID>.<label>.{tsv,stats.txt}` (EBV and mtDNA copies per cell, read length, chr20 depth and alignment metrics; `stats.txt` is samtools' own, curators only);
  - family: `families/…/qc/<FID>.<label>.{variants,inputs}.tsv` (per member SNV/indel profile). **No flags**;
  - cohort: `cohorts/<NAME>/qc/<NAME>.<label>.{samples,flags,inputs}.tsv`: everything, plus `seq_class`, `library`, `dna_source`, `flags`; one flag row per flag (`sample, family_id, seq_class, flag, severity, evidence`).
- **Same `qc/` directories and file names as `sampleqc`.** Each side keeps only its own labels by the `#step` of `inputs.tsv`: Sample QC leniently (`find_labels(step="sampleqc")`: no snapshot is kept), genome QC strictly (`_step_labels`: the step always writes one).
- **Flags exist at cohort level only**: outliers are judged within a `seq_class` (`rl100`/`rl150`, `_pcr` from the sample information table only), 20 samples at least. A family page reads the cohort tables (`load_family_gqc(…, cohort)`, rows whose `family_id` is the family); without them it joins the family's variants table to its members' sample tables, same columns, no class nor flag (`GenomeQc.flagged` false).
- **One row and one `family_id` per sample**, even for a sample in two families.
- **The flags are Sample QC `FlagRow`s** (`genome_qc.flag_rows`): level `genome`, `sample_b` empty, the pipeline's `low` severity (`EBV_TRACE`) as `info`. Their wording is in `sample_qc.yaml` (categories `dna_source` and `sequencing`), so `flag_info`, `ack_key` and the acknowledgements work unchanged. `genome_qc.yaml` holds the metrics (label, format, test side, log scale, threshold lines) and the `seq_class` wording; its thresholds mirror `genomeqc_tables` for the plots and never flag.
- **`CELL_LINE` is a caveat on the de novo, UPD and CNV results, not only an issue to review** (`GenomeQc.cell_lines()`). Readers see "DNA from a cell line"; `EBV_TRACE` is "Atypical DNA signal" for them (EBV in the person is health information), and `dna_source` notes are curator-only.
- **Shown in the Sample QC tabs** (`pages/cohort/components/genome_qc_section.py`, charts in `utils/genome_qc_plots.py`):
  - **The flags are added to each Sample QC flag list where it is used**, not inside `load_*_qc`, whose tables are per sites label:
    - `family_qc_summary(store, fid, cohort)` (header, icons);
    - `qc_family_to_review(…, cohort)`;
    - the family tab's `state["genome_flags"]` (Issues, downloads);
    - the cohort column and pills (`cohort.py`);
    - the cohort tab's `_all_flags`;
    - `cohort_headline(…, extra=)` (home);
    - `family_report`.

    `genome_qc.cohort_flags` reads the flags table only (a few KB, cached), so the header reads it on the event loop, as the de novo QC table.
  - **The "DNA and library" section** loads off-thread (`load_family_gqc`, plus the cohort's rows as the background) and draws one chart per metric, one row per `seq_class`. Each class gets its median and 5-MAD fence (`class_fences`, the pipeline's rule, drawn only) and the fixed lines. EBV is on a log axis from 1e-4 that always reaches 1.
  - **The table passes floats as read**: rounding EBV to its digits made 0.0005 read 0.
  - **The samtools stats links are curator-only and sit outside the table**: a `link` cell opens in the same tab. The data route does not follow a symlink out of the data directory, so a scratch directory of symlinks gives 404 there.
  - **`render_cell_line_badge`** (family header) and **`render_cell_line_callout`** show `CELL_LINE` whether or not it is acknowledged.
- **A cell line's results carry a caveat** (`components/cell_line_note.py`, wording in `genome_qc.yaml` `caveats`: `upd`, `cnv`, `diagnostic`; the de novo banner has its own):
  - UPD tab and dialog, SVs tab (a refreshable that follows the selection) and SV dialog, diagnostic dialog;
  - the family report: the first footnote of the aneuploidy, UPD and diagnostics sections (`CELL_LINE_CAVEATS`), and the `genome_qc` Methods paragraph when the cohort has the tables.

  `genome_qc.cell_lines(cohorts, fid)` reads the flags tables. A dialog opened without a cohort (Search, validation pages) uses `family_cohorts`: every cohort of the data store holding the family.
- Test data: `datagen.add_genomeqc_tables(data_dir, cohort="DEMO")` (opt-in): S3 (1001) `CELL_LINE`, S5 `EBV_TRACE`, S6 a 100 bp genome under 20x (`LOW_DEPTH`), S1 a PCR library (`rl150_pcr`); `cohort=None` writes the sample and family tables only.

### Uniparental Disomy (UPD)
- **Three steps point at a UPD and share no identifier**:
  - the ghfc-ngs `upd` step's events (`families/…/upd/<FID>.<label>.{events,segments,members,inputs}.tsv`, and `windows.tsv` from `upd_v4`; cohort tables concatenated, with `recurrent_families`/`recurrent` on the events and no windows). A span, an origin, a composition (iso/het/mixed), an extent and a confidence per member, analysed as trio, duo or solo (from the ROH runs);
  - the ROH `froh` column `upd_candidate_chrom` (one chromosome per person, no span);
  - the Sample QC pair flag `UPD_CANDIDATE` (`sample_a` the parent whose allele is missing, `sample_b` the child; the chromosome and span only in the evidence, `"… lie on chr2, between 0.1 and 39.7 Mb"`, parsed by `upd.chrom_from_evidence` / `span_from_evidence`).
- **Read the step's tables through `utils/upd_calls.py`**: `find_upd_labels` = `tsv.find_labels(…, "events")`, the newest-written label first (each label is a version of the analysis written beside the older ones; EAGER has `upd_v2` and `upd_v3`). IDs and `chrom` are forced to strings (`family_id` is numeric-looking), `chrom` loses its `chr`, and older labels' missing columns are nulls (`upd_v1` has no `support_fraction`, `shared_with`), so a family's events read `recurrent` as null too. Cached on `(path, mtime, size)` (`clear_upd_cache`).
- **`windows.tsv` counts raw Mendelian classes per 500 kb window**, with the parents the step *used*: `n_mat` (father homozygous, child het), `n_pat`, `n_lacks_father`/`n_lacks_mother` (opposite homozygotes), the child's `n_het`/`n_hom_alt` and `depth_ratio`. A duo has only the present parent's pair (the others `.`), solo none. The app only draws them (`upd_plots.chromosome_option`, mirrored: maternal up, paternal down); it never reads `norm.bcf`, so the classification stays in ghfc-ngs. Labels before `upd_v4` have no windows: the chromosome view draws the lanes only.
- **A UPD is keyed by family, child and chromosome**: `Variant` is `chr<N>:UPD`, in `validations/upds.tsv` and `diagnostics/upds.tsv`. The validation file has the SNV/SV header plus `Start End Composition Extent UpdLabel` (what the curator validated, prefilled from the event; appended, so an older file gains them on its next save); `Inheritance` holds the parental origin (`maternal`, `paternal`, `unknown`), so `load_validation_map` works unchanged. Writes go through `save_upd_validation`, which checks `can_write()` itself. A validation's span and type win over the event's in the report.
- **`upd.upd_candidates` merges the sources into one row per child and chromosome**: an event of the shown confidences (`default_confidences` in `upd.yaml`: high, medium; GHFC-NT `upd_v2` has 4,814 low duo events), the ROH chromosome, a flag naming it, or a standing validation (a validated chromosome stays when a new label drops its event). Once listed, a row carries all its events, whatever their confidence. `family_candidates` and `count_to_review` (not validated present or absent) feed the tab's badge; `probe_family_upd` enables the tab when the family has upd tables, a ROH/QC hint or a validation.
- **Every review opens `show_upd_dialog`** (`components/upd_dialog.py`) from an anchor outside its refreshable section: the UPD tab's candidates, the ROH tab's "UPD candidates" strip and the Sample QC Issues action on a `UPD_CANDIDATE` flag. **A candidate no curator settled is reviewed in the UPD tab**: the ROH strip and the QC action call `on_open_upd(child, chrom)` (the page's `show_upd_tab`, which switches tab and opens the review; a review asked before the tab is drawn waits in `upd_review["pending"]`). A settled one opens in place. The cohort Sample QC panel's family link of such a sample (`hrefField: "_href"`) and the cohort UPD panel's links open `?tab=upd&upd=<sample>:<chrom>`.
- **After a UPD save anywhere, every tab showing a UPD's status follows**: the pages' `on_acks_changed` redraws the header badge, member icons, Diagnostics panel and the UPD tab's count, then runs `upd_listeners` (the UPD, ROH and Sample QC tabs append their inner refresh; never an outer tab content, whose re-render would append its refreshers again). The diagnostic's `Impact` is the evidence source (`UPD step + ROH + Sample QC` or part of it).
- **A standing `present`/`absent` UPD validation acknowledges the flag** (`confirmed`/`explained`; the last non-ignored row wins). It cannot be an `AckKey` (it does not know the parent), so `load_acks` returns an `Acks` object that matches it to the flag by family, child and chromosome. **Look flags up through it** — `acks.get(flag)`, `flag in acks` — never by `flag.ack_key`. Its cache key covers both files. **It also settles the ROH candidate**: `roh.flag_badges(row, validation)` shows `UPD chrN · <origin>` for `present` and drops the badge for `absent`; `roh.is_flagged` ignores an `absent` one (tab, export, report).
- **Wording is role-dependent for the members table only**: a member's `status` (`parentage: S1 excluded`, `not_isolated: parents' kinship …`) is worded through `member_status_label` (readers get `reader_label`, null parts left out) and readers get no `father_used`/`mother_used` (`table_export`): which declared parent was used says when one was excluded as not the parent.
- **Imprinting is the step's `imprinting` column, as written** (per region in trios and duos, per chromosome in solo mode). GeneTrek's imprinting regions are on the roadmap.
- Test data: `datagen.add_upd_case(data_dir)` gives 1001's S3 a chr2 UPD (ROH and flag, family and cohort files). `datagen.add_upd_calls(data_dir)` writes the step's tables for 1001 and DEMO under `upd_v3` (no windows, dated a day earlier) and `upd_v4` (windows): S3 trio paternal iso chr2 0.1–39.7 Mb `high` (the flag's span), S3 chr16 `low` (recurrent at cohort level), S2 solo chr15 `medium` with the Prader-Willi/Angelman text. Both are opt-in: they change flag counts and report sections.

### Aneuploidy
Whole-chromosome gains and losses, and atypical sex chromosomes, validated and diagnosed like any variant. `utils/aneuploidy.py` is the model, `components/aneuploidy_dialog.py` the card and dialog.
- **Two sources, no pipeline step**: WisecondorX calls past the permissive thresholds summed per sample and chromosome (`family_calls` → `chromosome_events`, `whole_chromosome_events` its text; the report and the cohort SVs tab read them there), and the Sample QC karyotype (atypical, or not the declared sex's) with its `KARYOTYPE_<k>` flag. `aneuploidy_min_fraction` comes from `report.yaml` (`DEFAULT_MIN_FRACTION`).
- **One record per sample and chromosome**: `Variant` is `chr21:ANEUPLOIDY`; X and Y are one record, `chrXY:ANEUPLOIDY` (`record_chrom`), since XXY is about both. `validations/aneuploidies.tsv` and `diagnostics/aneuploidies.tsv` (kind `Aneuploidy` in both registries). The validation file has the SNV/SV header plus `Direction Mosaic Karyotype` (ISCN, prefilled by `suggested_karyotype`: `47,XY,+21`, `47,XXY`, `45,X`); `Inheritance` is the origin (`de novo`, `maternal`, `paternal`, `unknown`). `save_aneuploidy_validation` checks `can_write()` itself.
- **`family_candidates`** gives one row per member and chromosome from both sources, plus a standing validation's chromosome; `count_to_review` = not validated present or absent. `load_family_review` adds the calls and QC rows (`segments` for the dialog).
- **Where**: a card above the SVs tab's sub-tabs (`aneuploidy-card`, rows `aneuploidy-row-<sample>-<chrom>`, `aneuploidy-review-<sample>-<chrom>`, `aneuploidy-to-review` for curators), drawn only when the family has a candidate and following the member selection; the Sample QC Issues action on a `KARYOTYPE_<k>` flag opens the same dialog (`aneuploidy-dialog`, `aneuploidy-validation`, `-direction`, `-mosaic`, `-karyotype`, `-origin`, `-validation-save`, `-diagnostic-open`). The SVs tab is enabled only with SV files, so a family with a karyotype flag alone is reviewed from Sample QC. A save runs the page's `on_acks_changed`; the card's reload is in `upd_listeners`.
- **A standing `present`/`absent` `chrXY` validation acknowledges the sample's `KARYOTYPE_<k>` flag** (`confirmed`/`explained`, `Acks.aneuploidy`, `read_aneuploidy_acks`; the ack cache key covers the file).
- **Readers never see a sex-chromosome karyotype**: the card, dialog and Diagnostics panel say `Sex chromosome anomaly` (the flag's `reader_label`) and leave the karyotype out (`describe_aneuploidy(show_detail=)`), as Sample QC does.
- Elsewhere: the Diagnostics panel (`Aneuploidy` tag, the validated karyotype), the report's diagnostics, the Diagnostic pages, the validation pages (a row opens the family's SVs tab), and the cohort SVs tab's whole-chromosome table and download (`validation` per event, `cohort_svs.aneuploidy_validations`).
- Test data: AU-X-1's S5 reads XXY with `KARYOTYPE_XXY` (default data); `datagen.add_report_case` gives 1001's S3 a chr21 gain.

### De Novo Step (SNVs and indels)
- **The ghfc-ngs `denovo` step** reads the Wombat de novo table (the `#wombat_config` line of its inputs, `de_novo_mutations` by default) and writes `families/…/denovo/<FID>.<label>.{dnm,qc,inputs}.tsv` and `cohorts/<NAME>/denovo/<NAME>.<label>.{dnm,qc,inputs}.tsv`. The labels are `tsv.find_labels(…, "dnm")`, newest first. The contract is `modules/denovo/scripts/denovo_tables` in ghfc-ngs; the vault note is `ghfc-ngs/De novo step.md`.
- **One dnm row per Wombat row**: the family table joins the aggregated Wombat de novo rows 1:1 on `(chrom:pos:ref:alt, child)` (checked on CANDY_mpx_ghent: 510 of 510). `denovo.add_evidence` sets the `dnm_*` fields on the Wombat rows and keeps the dnm row as `_dnm` for the variant dialog. A Wombat row with no dnm row gets the app's own `NOT_ASSESSED` flag and stays in Strict view.
- **Read through `utils/denovo.py`.** `.` is null and the IDs are forced to strings (`family_id` is numeric-looking). Tables are cached on `(path, mtime, size)`; `clear_denovo_cache` is in `conftest._reset_globals`.
- **The cohort table is preferred for a family when it is as new** (`load_family_denovo(…, cohort)`). Only it carries `recurrence_families`, the `RECURRENT` flag (which can turn a strict candidate flagged) and the per-child QC verdict (`qc`: ok, high, low, no_callable, with `qc_low`/`qc_high`/`qc_reference`, and `sampleqc_flags`).
- **Strict = no flag but `SIB_SHARED`** (`strict_ignored` in `config/denovo.yaml` mirrors the pipeline's `STRICT_IGNORED`). The QC counts are autosomal. The cohort tab recounts the strict SNVs and indels, transitions and origins from the candidates (`child_rows`), so every figure follows one strict set; the QC table's split is the family pass's, from before `RECURRENT`.
- **Flags, statuses and verdicts are the pipeline's**; the app words them (`config/denovo.yaml`, a tooltip per flag) and filters on them. An unknown flag is shown raw, grey.
- **The De novo sub-tab is first in SNVs/Indels** (`wombat_tab.DENOVO_TAB` = `denovo`, link `wcfg=denovo`; `wcfg=<the source config>` from older links maps to it). The source config gets no sub-tab of its own (`denovo.denovo_source`). Without the step's tables it is the plain Wombat table, with a note.
  - With the step's tables: `DenovoFilter` (link `wdm` strict|all, `wfh` flags hidden, `wfo` flags kept), Wombat exclusions off by default (`WombatFilters(exclusions=False)`, the link defaults follow), the "De novo View" preset.
  - The curation banner lists the rows the filters hide that are curated de novo (`is_curated_denovo`: validated present without an inherited mode, i.e. not paternal/maternal/either/homozygous, or with a diagnostic), and the family's `de novo` validations with no row at all.
  - The variant dialog gets the Python row (`view.row_for`), which holds `_dnm`.
- **With a genome QC, the cohort QC table adds the germline checks** (`seq_class`, `genomeqc_flags`, `qc_reference` `cohort:<class>`, `paternal_share`/`paternal_p` of the strict SNVs, the spectrum and `spectrum_cosine`, the indel profile, `expected_germline_share`, `qc_flags`: `CELL_LINE`, `NOT_GERMLINE_LIKE`, `SPECTRUM_OFF`, worded in `denovo.yaml`).
  - `child_rows` keeps the app's own `paternal_share` (all strict, recounted) and adds the pipeline's as `paternal_share_snv`.
  - `germline_text` words the checks, the expected germline share only for a high count or a flagged child.
  - `is_cell_line` reads both flag columns.
  - No gate: a cell-line child keeps its table, with `render_germline_banner`.
  - `qc_alerts` includes a flagged child whose count is in range.
  - The cohort tab's literature comparison leaves cell lines out (none left: no comparison). The rate chart draws QC lines only for one shared range.
  - Test data: `datagen.add_denovo_calls(…, genome_qc=True)` makes S3 a high, not-germline-like cell line.
- **Cohort De novo tab** (`cohort_denovo_tab.py`): reads the cohort tables and the curation files off-thread and follows the individuals filter. The literature table is `references` in `config/denovo.yaml`, each `metric` a key of `denovo.cohort_summary`. Rates are per child per Gb of haploid callable unmasked sequence, the pipeline's `rate_per_gb`, so 2 × μ × 1e9. The values were read from the papers in Zotero on 2026-11-14; keep the sources with the numbers.
- **The curation check** (`curation_check`): `hidden` = a non-strict candidate curated as de novo; `missing` = a `de novo` validation (present) of an assessed family with no candidate. It shows on the cohort tab, in its pill and in the family banner, for curators only.
- **Badges**: the family SNVs/Indels badge counts strict (or not assessed) HIGH/MODERATE unvalidated candidates when the step ran (`denovo_to_review(…, cohort)`). The cohort De novo pills are QC high, QC low and curation issues (`denovo_cohort_status` off-thread, `denovo_cohort_pills`). The family header's "De novo QC" badge reads only the QC tables (`load_child_qc`).
- Test data: `datagen.add_denovo_calls(data_dir)` gives 1001 (`rare_denovo`) the step's tables. SHANK3 is strict and paternal; PLAIN is flagged `TR_HOMOPOLYMER,PARENT_ALT_POSSIBLE` and validated de novo, so the banner and the check find it. The DEMO cohort tables are written one second later and give S3 the verdict `low`. Opt-in: it changes 1001's de novo table and badge.

### Family Report
Curators send it to the referring clinician of each centre, who decides whether to recontact the family; it is not addressed to patients. "See report" and "PDF" sit at the right end of the title row of both family pages (`render_report_buttons`), for `can_write()` only, so it always uses the curator wording.
- **One HTML string, two uses.** `build_family_report` (pure, blocking, `asyncio.to_thread`) → `render_report_html`. The dialog shows it in an `<iframe srcdoc>`, which keeps Quasar and Tailwind out, so the screen looks like the print. PDF writes it into a hidden iframe and calls the browser's print. The document `<title>` (`<cohort>_<family>_report_<date>`) is the PDF's default name. No server-side PDF library: the VM installs with `uvx`.
- **Plots are SVG strings from `report_svg.py`**, escaped there and marked safe only in `render_report_html` (`svg_of`). `@page` margin boxes give the footer and page numbers (Chrome 131+).
- **The pedigree is drawn by pypedigree** (`utils/pedigree_svg.py`, Alex's library, `bourgeron-lab/pypedigree`; plan in the vault's `Pedigree drawing` note). `render_svg(id_prefix=, responsive=True)`: no script, ids prefixed, CSS scoped under `.pypedigree`. Three workarounds of 0.5.2 (drop the first and third once pypedigree#40 and #41 are released):
  - **the text parser, not `PedigreeBuilder`**: the builder leaves a child with one declared parent unconnected;
  - **a parent named but absent is written as a `not_sampled` row** before parsing: the parser would warn, and `warnings.catch_warnings` is not thread-safe (the report runs in `asyncio.to_thread`);
  - **spacing from the longest label** (`_style`): the layout spaces symbols only, so long dummy IDs overlapped.

  The report sizes it in mm from the `viewBox` (`pedigree_drawing` in `report.yaml`) through `Plot.width`, and writes the observed karyotype under an ID when atypical or not the declared sex's. Tests check structure (`data-id`, `Not sampled`, prefixes), never pixels: figures change between pypedigree minors.
- **Order is the one asked for:** pedigree, ancestry, Sample QC, ROH, PGS, aneuploidy, UPD, diagnostics, with the curator conclusion and "Results at a glance" (`Section.glance`, else `note`) under the pedigree. Each `Section` is available with rows, or missing with a reason.
- **Aneuploidy** is two sources: the QC karyotype for sex chromosomes; for autosomes, same-direction WisecondorX segments covering `aneuploidy_min_fraction` of a chromosome's assessable span (the q arm of an acrocentric one). The pipeline has no aneuploidy step. Both are read through `utils/aneuploidy.py` (see Aneuploidy), and the section adds each candidate's curator validation; a chromosome validated absent is not reported as found.
  - **Only calls past the permissive thresholds count** (`wisecondorx.past_permissive()`: log2 ratio *and* z-score). `whole_chromosome_events` sums whatever it is given, and a chromosome-long segment reaches |z| ≈ 5 with a log2 ratio of ±0.01: before v0.41.2, 88021200005 read +14 (C0XYE8K) and −21 (C0XYE8L) from two such segments. Filter before calling it.
- **PGS percentile = Φ(z)**: the pipeline's z is already ancestry-matched (k-NN SD in PC space). The traits shown are `pgs_traits` in `report.yaml`.
- **Diagnostics** keep `diagnostic_verdicts` (pathogenic, uncertain, conflicting) and add the standing validation and inheritance per sample.
- **Diagnostic details are a snapshot, not a lookup** (`utils/diagnostic_snapshot.py`): the table row's HGVS, zygosity, segregation, carriers, ClinVar, gnomAD, MOI (worded, from the badge tooltip), SV gene overlap and copy number are written with the diagnostic (`SNAPSHOT_COLUMNS`, after the ten old columns), and the report never re-reads Wombat or SV tables. `Source` is always set on a snapshot, so an empty `Source` means a record saved before v0.28.0 ("Not recorded"). Each call site passes `details=` to `show_diagnostic_dialog`, built from the Python row (the row sent back by the table may lack hidden fields).
- **The MANE HGVS is picked before Wombat rows are merged** (`wombat.PICKED_HGVS_COLUMNS`: `_hgvs_transcript`, `_hgvs_mane`, `_hgvsc`, `_hgvsp`): after `aggregate_by_variant` each column joins its values on its own, and which `VEP_HGVSc` belongs to the MANE transcript is lost. Columns starting with `_` are kept out of the Wombat and Search column lists.
- **Methods text and tool versions are in `report.yaml`, hand-copied from ghfc-ngs** (`methods_source` names the commit). The run's own version and commit come from the run record. Update both when the pipeline changes a tool.
- **Curator conclusions** live in `reports/family_summaries.tsv` (`FID Summary User Timestamp`), append-only, last row wins, empty clears. Line breaks are written `\n`. `save_family_summary` checks `can_write()` itself.
- Test data: `datagen.add_report_case(data_dir)` adds the four report traits, a chr21 gain for S3 and diagnostics. It is opt-in because it changes the PGS tables other tests compare whole.

### Family Status
Where each family stands in its review (vault note: Roadmap, "Structured family status"). `utils/family_status.py` is the model, `components/family_status.py` the chip, dialog and bar.
- **Statuses are `config/family_status.yaml`**, in review order: `not_reviewed` (no record, the default), `in_review` (a claim: the record names who), `candidate`, `partially_solved`, `solved`, `unsolved`, `not_analysable`. Each has a `yield` role: `solved`, `reviewed`, `open` or `excluded`.
- **`reports/family_status.tsv`** (`FID Status Comment User Timestamp`), beside the conclusions: append-only through `tsv.append_tsv_row`, the last row of a family wins, keyed by FID only (one status across cohorts). `save_family_status` checks `can_write()`; `save_family_status_async` checks it on the loop and writes in a thread. The file is in `pedigree_merge.curation_files`, so a merge re-keys it. Cached on its stat (`clear_status_cache` in `conftest`, which also clears the per-family diagnostics cache).
- **The status never changes by itself.** `family_state` works out *reasons* to look again, shown as an amber "Reopened" (a pedigree edit) or "Check" badge with the text in the tooltip:
  - `pedigree_waiting` / `pedigree_ready`: the family's last pedigree edit (`edit_statuses`, against the cohort's last run) is after the status. Phenotype edits count: they rerun nothing in ghfc-ngs but change who is affected;
  - `diagnostic_after`: a standing diagnostic's `Timestamp` is after the status;
  - `solved_without_diagnostic` / `diagnostic_not_solved`: a `solving_verdicts` diagnostic stands or not, against the status, whatever the dates (an Ignore toggle keeps its timestamp, so it is caught here). No reason for `in_review` or `not_analysable`;
  - `qc_after`: an unacknowledged high flag in QC tables (Sample QC, genome QC) written after the status (file mtime: flags carry no date);
  - `results_after`: a `*.inputs.tsv` of the family written after the status (`step_outputs`, which walks the family directory). **Family pages only**: across a cohort that is one walk per family over SMB.

  Reasons about dates need a status record; saving a status again clears them.
- **The yield** (`count_statuses`): solved = `solved` + `partially_solved`, over the reviewed (`yield` solved or reviewed), and the analysable families (every counted one). A family is not counted when `not_analysable`, or with no affected sequenced member and no record.
- **Every role sees the status**; the comment is curators only (`status_badges(state, curator)`), and the report states the status without it.
- **Where**: the family pages' chip (`family-status`, `family-status-reason`; dialog `family-status-dialog`, `-select`, `-comment`, `-save`, `-history`), drawn from the status file at once, its reasons read off-thread after (`FamilyStatusChip.update`, step outputs walked once per page) and again after each curation save (it is in `data_table_refreshers` and `on_acks_changed`); the Families tab's first column (`status_fields`, filter link `fstatus`, options the labels plus `Reopened` and `Needs re-review`), read in `load_family_facts`, empty until then; the bar and yield above the table over the shown families (`families-status-summary`), and the curators' "to recheck" pill (`status_pills`); the home cards (`home-status-<cohort>`, `cohort_status_counts` in the off-thread pass); the report's "Review status" line.
- Test data: `datagen.add_family_statuses(data_dir)`: 1001 in review then solved, 1002 not analysable, AU-X-1 none. Opt-in.

### Pedigrees Across Cohorts
The cohorts of a data directory share one `families/` tree, and ghfc-ngs builds a family from the rows of the cohort it runs (`main.nf`, `FamilyPedigree.onDisk`): a family two cohorts list differently is rebuilt, and its results overwritten, by each cohort's run. Barcodes have no aliases, so a barcode under two FIDs is one family known by two codes. Measured on GHFC-GRCh38 (2026-11-18): four families differ between cohorts (C0733-011-021, AU-RD-LAR-016, 420498, C0733-011-126), 59 are the same people under two FIDs (CANDY_mpx ↔ dataset_paper_20251010), one is part of another (820910 ⊂ 88081200008).
- **`utils/pedigree_consistency.py`** reads every cohort pedigree as the pipeline does (`rows_from_lines`: positional columns, first row per member, `pedigree_edit.norm_*`), and `scan(store)` returns `Finding`s (cached on the files' stats; `clear_consistency_cache` in `conftest`). Codes and wording are `config/pedigree_consistency.yaml`:
  - one FID in several cohorts: `SAME_FID_ROWS_DIFFER` (parents, sex, dummies), `SAME_FID_MEMBERS_DIFFER`, `FID_COLLISION` (no member in common), all high; `SAME_FID_PHENOTYPE_DIFFERS` (info: a cohort per phenotype explains it);
  - FIDs sharing a sequenced barcode (dummies never link: their IDs embed a barcode): `ALIAS_SAME_MEMBERS`, `ALIAS_SUBSET` (medium), `ALIAS_OVERLAP`, and `BARCODE_ROWS_DIFFER` (a shared sample with another sex, or other parents both families hold), high.
- **Acknowledgements**: `validations/pedigree_consistency.tsv` (`Key Code Families Status Note User Timestamp`), append-only, last row wins, `explained` or `withdrawn`. `Finding.key` holds a digest of the rows involved, so a change of rows raises a finding again. Info findings are never acknowledged nor counted. `save_ack` checks `can_write()`.
- **Where**:
  - every role: the family pages' banner (`family-consistency`: the other places, linked);
  - curators also get, in the banner, the details and the explain dialog (`consistency-explain`, `-dialog`, `-note`, `-save`);
  - curators: the home warning (`home-consistency`), `/pedigree-consistency` (every finding), the cohort's findings on the Workflow tab (`workflow-consistency`), and the Families tab's Pedigree column badges (`consistency_badges`: `Differs: <cohort>`, `Alias: <FID>`);
  - an open high finding is a family status reason (`pedigree_differs`, through `family_status.disagreements_by_family`); aliases are not.
- **The pedigree editor** no longer refuses a family whose cohorts differ. `load_family(…, reference=)` keeps each file's rows (`FamilySource.drafts`) and edits the reference's (`original`). A save writes the reference rows, as edited, into every file (`FamilySource.target`: each file keeps its own phenotype unless the edit changed the reference's), and the record lists each other cohort's changes (`file_changes`, `NOHEAD: S3: father 0 -> S1`) with the strongest rerun (`rerun`). The dialog asks for the reference first (`ped-reference`, `ped-reference-<cohort>`) and says what the save writes (`ped-cohorts-differ`).
- **"Seen in N families" counts a family once** (`dir_families`, see Gene Index), whether two cohorts list it under one FID or two.
- **Families ghfc-ngs held** (ghfc-ngs#66, `pedigree_conflicts: hold` by default): before scheduling, a run reads the other cohorts' pedigrees and leaves a family they list with other rows (phenotype aside) as it was, recording it in `last_run.pedigree_conflicts` (`[{family, cohorts, changes: {cohort: [clause]}}]`; absent = not recorded, `[]` = none). `cohort_state.RunRecord.held_families` (`HeldFamily`, `describe()` → `differs from NOHEAD: S3: father S1 -> 0`). Shown on the Workflow tab (`workflow-held-families`, every role; the changes for curators), on the family pages (`family-held`, every cohort whose last run held it) and as a `Held` badge in the Families tab's Pedigree column. Test data: `datagen.add_held_families(data_dir)` (DEMO's last run held 1001; opt-in).
- Test data: `datagen.add_cohort_disagreements(data_dir)`: NOHEAD lists 1001 with S3's parents at 0; cohort ALIAS knows AU-X-1 as AU-77 and S3 alone as K3. Opt-in (a third cohort).

### Gene Index (gene and variant pages)
`/gene/<symbol>` and `/variant/<chr:pos:ref:alt>` (header box `gene-lookup`) gather a gene or variant from every cohort of **every data directory**, through `utils/gene_index.py`.
- **Only cohort-level tables are indexed**: Wombat results (`cohort_wombat_sources`), WisecondorX calls (not `Below threshold`, 87% of a table), the newest STR catalog's genotypes, the `roh_genes` table. A family without a cohort table is missing; that is said in the CHANGELOG, not worked around.
- **Variant counts** (`variant_counts.parquet` per directory): the snv files summed up by variant and cohort (carriers with an alternate allele, families), rebuilt in the build process when an snv source changed, and held in memory (~30 MB for 900k variant × cohort pairs). The family Wombat tables read them on the event loop from memory only (`load=False`; the refresh loads them, the held de novo path loads off-thread): `cohort_carriers`, `cohort_freq` (over the cohort's pedigree), `seen_families` (every cohort the user may open, this family counted even when its cohort has no table; a link to `/variant/`) and `seen_cohorts` (`add_variant_counts`). A lookup by scanning the snv files took 4-6 s over SMB. **A family counts once per data directory** (`dir_families`, the same on each cohort row of a variant): one family listed by two cohorts, or under two FIDs, shares its carriers' barcodes, so each carrier takes the smallest `(cohort, FID)` it is listed under and the distinct ones are counted. Overlapping families stay two. The gene and variant pages count alike (`gene_rows.distinct_families`). A counts file without the column (`_counts_current`) is served with nulls (the rows then sum the cohorts, as before) and rebuilt at the next refresh, alone: no `INDEX_VERSION` bump. On GHFC-GRCh38 it fixed 25,876 variants (21,959 seen in at most three families), mostly families in GHFC-NT and dataset_paper_20251010 under one ID; peak memory of the build 0.24 → 0.40 GB (SPARK WES 1.5 → 1.9 GB).
- **One Parquet file per source, sorted by `gene_key`** (upper-cased symbol), under `gene_index_dir` (default `<config dir>/.gene_index/<sha1(data dir)>/`), with a `manifest.json` of each source's `(mtime, size)`, saved after every table. A refresh re-reads only changed sources, drops vanished ones, and runs off-thread at startup and every `refresh_minutes`. The manifest on disk is published before the check, so a restart serves lookups at once. The default sits next to the config file, which on the VM is on the share.
- **Each table is built in a process of its own** (`_build`: spawn, `nice 10`, `oom_score_adj 1000`, half the cores through `POLARS_MAX_THREADS`). v0.26.0 built in the app: CANDY_mpx's WisecondorX table took 13.4 GB, the VM's OOM killer restarted the app every ~45 s, and each restart began the index again. Now a killed build is that table's error (`KILLED`), kept until its `(mtime, size)` changes; the app is untouched. Tests build in-process (`GENETICS_VIZ_GENE_INDEX_IN_PROCESS`, set by `conftest`); one slow test covers the child.
- **Memory**: Wombat tables are scanned lazily (`scan_cohort_wombat`) and grouped per variant, sample and gene in the streaming engine (SPARK iWES-v3's 2 GB table: 14.7 GB eager, 4.6 GB now); a gene column is exploded apart from the other columns (`_row` join). Measure a new source's peak memory on the real tables before releasing: the share is mounted on the Mac.
- **Every kind's files have exactly `SCHEMAS[kind]`** (`_conform`), so all cohorts scan as one table. Changing a schema means bumping `INDEX_VERSION`, which rebuilds everything.
- **Only the dozen Wombat columns used are parsed.** A headerless table goes through `read_cohort_wombat(columns=...)`, which selects by index, since polars 1.x and 2.x disagree on how `new_columns` combines with a selection.
- **Validations and diagnostics are not indexed**: read live per data directory (`gene_rows.load_curation`) and joined by the keys the family tabs save (`chr:pos:ref:alt`, `chr:start-end:del|dup`, `str_variant_key`). Diagnostics are also found by their `Gene` column, whose SV genes carry an overlap suffix (`SCN2A:exonic`); the variant page reads the validation history straight from the files, so a variant in no cohort table still shows its curation.
- **Look genes up through `gene_rows.gene_keys`**: the current symbol plus previous/alias names (`matching_identifiers`), since VEP symbols lag HGNC.
- **Rows of another data directory link through `/goto/<n>?next=<path>`** (`_dir`, `_next`; `n` from the `/data-N` prefix, never a path), which sets `app.storage.user["data_dir"]` and follows `safe_next`.
- The pages only read. Their diagnostic column is a `badge_list` (`DiagnosticBadges`), not the `diagnostic` cell, which would offer a "+" that opens nothing.
- Test data: `datagen.add_gene_case(data_dir)` adds cohort OTHER sharing DEMO's SHANK3 de novo, with validations and diagnostics; `datagen.build_second_data_dir(root)` gives a second data directory. Both opt-in.

### Family Member Card (Table | Tree)
- **The member card has a Table | Tree switch** (`member-view`, `member_selector.py`); the view is `app.storage.user["member_view"]`, a preference, not a link key. The tree is built on first showing (`components/pedigree_tree.py`).
- **Selection stays in the table.** The tree only shows it: every change (checkbox, only, All `members-all`, None `members-none`) and every QC review (`update_qc_icons`) redraws the whole SVG (about 1 ms). Nothing is looked up by `data-id` in JavaScript, and `ui.html` would run no script anyway.
- **Classes on the nodes, styled by the tree's own `<style>`**: `deselected` (faded; `dimmed` is pypedigree's own), `qc-high` / `qc-medium` (the worst open flag, the table icon's colour). They go through `PedigreePerson.classes`, the PED `css_class` column.
- The page's tree and the Enlarge dialog's (`member-tree-enlarge`, `member-tree-big`) are two drawings, id prefixes `member-tree-` and `member-tree-big-`. Each box is sized from the `viewBox`: the card's width, `PX_PER_UNIT` per layout unit and `MAX_HEIGHT_PX` high, whichever is smallest.
- In tests, `user.find` sees visible elements only: assert a hidden tree with `should_not_see(marker="member-tree")`, and match `node deselected`, not `deselected`, which the CSS holds too.

### Pedigree Edits (curators)
`utils/pedigree_edit.py` is the model and writer; `components/pedigree_editor.py` the dialog (vault note `Pedigree editor`).
- **ghfc-ngs reads only the cohort pedigree** and writes `<FID>.pedigree.tsv` from it (`FAMILIAL_PEDIGREE`: the cohort's header line, then the family's rows, first per member). A family file edited alone is overwritten at the next run. So an edit rewrites the family's lines in **every cohort pedigree listing the family** (`family_cohort_files`; when they differ, the reference cohort's rows are written into all, see Pedigrees Across Cohorts), leaves every other line byte for byte, puts new rows after the family's last line, then writes the family file as the pipeline would (`family_file_text`).
- **Reruns**: each pedigree-reading step keeps `<prefix>.inputs.tsv`, the family's rows normalised (`fid sample father mother sex`: unknown parent `''`/`0`/`-9`/`.` → `0`, sex `1`/`2` else `0`, no phenotype) and reruns when they differ. So a parent, sex or dummy edit needs a normal run; a phenotype edit reruns nothing; a change of sequenced members needs `--clean-stale-families` and is refused for now (`rerun_needed`). `pending_steps` compares every `*.inputs.tsv` of the family directory with the current rows, by member (Wombat's are sorted; the STR snapshot has no dummies). `0.0`/`-9.0` are missing for the app (`models.MISSING_PARENT`) but an ID for the pipeline: `problems` warns.
- **`save_pedigree_edit`** checks `can_write()` itself, needs a comment, refuses an error from `problems` (father not male, cycle, a new dummy not `DUMMY_<child>_<F|M…>` or already in another family…), a live run (`live_run`: `.ghfc-ngs.progress.json` saying `running` and updated in the last `LIVE_RUN_MINUTES`; the pipeline never removes the file, a finished run leaves `success`/`failed` in it and a killed one a `running` that stops moving) and a cohort file changed since `load_family` (sha256). Under `validations/.pedigree_edit.lock` it writes every file beside its target first (a permission error stops it before any change), renames them over (mode kept), then keeps `families/…/<FID>/pedigree_edits/<FID>.<EditId>.{before,after}.tsv` and appends the family log and `validations/pedigree_edits.tsv` (the cohort page's index). An empty phenotype is written `0`: ghfc-ngs drops trailing empty fields and stops on a five-column row.
- **The dialog** opens from **Edit pedigree** (`pedigree-edit`) in the family's Sample QC tab, curators only, hung from its own anchor outside the tab's refreshable section. `load_family` runs off-thread; the draft is edited in place (`apply` keeps an undo stack) and redrawn whole: the trees (`ped-edit-now`, `ped-edit-new`, changed members in a `changed` class), one row per member and implied dummy (`ped-father-<iid>` … `ped-pheno-<iid>`, `ped-swap-<iid>`, `ped-remove-<iid>`; the parent select's `__new_dummy__` adds a dummy for that member), the full-sibs control (`ped-shared-*`), changes (`ped-changes`), problems (`ped-problems`). Save (`ped-save`) needs a comment and no error; `save_pedigree_edit_async` checks the role on the loop and writes in a thread. After a save the tab reloads the cohorts (`store.reload`, off-thread) and the page. The history (`ped-history`) shows each edit's trees and restores its rows before (`ped-restore-<id>`) into the draft, saved as a new edit.
- **The Sample QC tab's status note** (`render_pedigree_status`, curators, off-thread): the steps whose `inputs.tsv` still holds older rows (`pedigree-pending`, in the tabs' order) and the last edit (`pedigree-last-edit`).
- **The Families tab's Pedigree column** (curators): the last edit's date from `validations/pedigree_edits.tsv`, amber while no ghfc-ngs run started after it (`edit_statuses` against the run record's `last_run.started_at`; a phenotype-only edit never waits). The checkbox filter "To rerun" (link `fped`) keeps those. Its header and filter set an 85 px floor: a curator's table no longer fits a 1,280 px window (49 px of scroll; it fits from about 1,330 px).
- **Suggestions from Sample QC** (`utils/pedigree_qc.py`, the dialog's `ped-suggestions`): `Kinship` and `compare` mirror ghfc-ngs `sampleqc_classify` (`163ce59`: KING cut-offs, `Pedigree.kinship`/`relationship`, `compare`'s flags), so `check_pairs` gives the flags the pipeline would write for a draft, against the observed `pairs.tsv` as written; a `UPD_CANDIDATE` pair (acknowledged or not) is no `NOT_PARENT`. `Fit` = (contradictions, pairs related beyond the pedigree); the summary's `ped-fit` line words it now and with the changes. `suggest` tries candidate edits on a copy and keeps those that improve the fit: a dummy for a parent who is not one, a dummy of their own for half sibs declared full, the parent full sibs should share, swaps (the best `MAX_SWAPS`). Sibship ones (a shared dummy for full sibs, one each or of their own for half sibs, with a parent written 0: the pipeline flags none of it) are kept when the fit is not worse. Sex from an XX/XY karyotype, with a caution when the sample's pairs contradict the pedigree (a swap). Applying one (`ped-apply-<key>`) records its flags with the edit. Keep the mirror in step when `sampleqc_classify` changes.
- **Across families** (phase 4): the dialog's "Another family" (`ped-other-family`, `ped-other-open`; families related in the cohort's Sample QC pairs first) opens a second family alongside. Only families listed in exactly the same cohort pedigrees can be edited together (`companion_families`; the writer refuses otherwise). A member row's menu then offers `swap_across` (`ped-swap-<iid>-<other>`: each barcode takes the other's row and references) and `move_member` (`ped-move-<iid>`: its children keep one shared dummy; it arrives a founder). `cross_problems` checks the families together (no sequenced member lost or doubled, none emptied); `save_cross_family_edit(_async)` writes them under one lock and one edit ID, each cohort file rendered family by family (`PedigreeText.from_text` between), one record per family (`rerun` `clean_stale_families` where members changed, the comment naming the other family). `QcContext.with_family` reads both families' pairs and the cohort pairs between them as one family; `suggest_across` offers moves (with no parents, or as the relative's full sib, parent, child) and swaps that improve the combined fit.
- **The Workflow tab lists the families waiting for a run** (curators, `render_waiting_families`, `workflow-pedigree-waiting`): edits (not phenotype-only) since the run record's last start, the rerun each needs, and the joint call's samples against the pedigree's (`bcf_samples`: the `vcfs/<FID>.norm.bcf` header, pysam; ghfc-ngs's own drift check). A family whose members differ from its call needs `--clean-stale-families`, and the note names them.
- **Merging families** (phase 5, `utils/pedigree_merge.py`): `load_merge` reads two or more families in the same cohort pedigrees; the target is a new FID (suggested: the IDs joined by `_`, `suggested_fid`) or one of theirs, checked by `target_problems` (no space, comma or slash; a new one free in the cohort pedigrees and under `families/`). The editor opens the merged draft (`MergeSource.as_source`) in a merge mode of the same dialog (`ped-merge-fid`; from the editor's `ped-merge-with`/`ped-merge-start`, or the cohort Sample QC tab's merge candidates), its QC the families' as one (`pedigree_qc.merged_context`, cohort pairs between them included). `save_merge(_async)` replaces the families' lines where the first stood (`render_merge_text`), writes the target's family file, keeps every family's rows before and the merged rows after (the old families' histories too), then **re-keys the curation** (`rekey_curation` over `curation_files`: validations, diagnostics, `validations/sampleqc.tsv`, `notes/notes.tsv`, `reports/family_summaries.tsv`; not a `"FA,FB"` acknowledgement) after copying each file it changes into `pedigree_edits/<target>.<EditId>.curation/`. Records: `merged` for the old families, the target's (`normal` for a new FID, `clean_stale_families` for an existing one) last. The old directories stay (ghfc-ngs leaves them): `orphaned_families` lists them on the Workflow tab (`workflow-orphaned-families`). A new curation file keyed by FID must be added to `curation_files`.
- **Renaming a family** is a merge of that one family (`load_merge(..., rename=True)`, `MergeSource.is_rename`; the editor's `ped-rename-start`, the same dialog and `ped-merge-fid` input). The target starts as the family's own ID, which `target_problems` refuses ("Type the family's new ID"). Records say `renamed to <FID>` (`RENAMED_TO`; `orphaned_families` reads both prefixes) and the old family's `Rerun` is `merged`. **A former family's ID may be taken again** (`MergeSource.former`: FIDs with a `merged` record in the index; `reuses_directory`): its leftover directory is allowed and the rerun is `clean_stale_families`. A directory of any other family is still refused. **Saving a merge under `suggested_fid` asks first** (`ped-confirm-fid`, `ped-change-fid`, `ped-keep-fid`): curators kept the joined IDs by mistake.
- **Merge candidates** (curators, the cohort Sample QC tab's `merge-candidates`): `pedigree_qc.merge_candidates` groups families joined by cross-family pairs of `MERGE_CLASSES` (first and second degree; a duplicate is not one). A group whose families are in different cohort pedigrees has its button disabled. In a merged draft, a pair observed parent-offspring and unlinked gets "Make D K's father" (both ways, by sex): for families entered one member per FID.
- **The cohort Sample QC tab's Download menu** adds the cohort pedigree as written (`download-sampleqc-pedigree`, every role that may download) and, for curators, the edit history (`download-sampleqc-pedigree-history`, `utils/pedigree_history.py` + `templates/pedigree_history.html`, `Download(media_type="text/html")`). The history lists the index's records whose `Cohorts` hold the cohort, each taken from its family's log (the index has no `Flags`), grouped by family, last edited first. Its trees are drawn with pypedigree (`pedigree_svg`) from the kept rows: the first edit's `before`, the last edit's `after`, members changed between them classed `changed`. A merged-away family has no rows under its ID after (its record's clause is shown), a merge target none before. `changed_ids` and the plain `draft_people` live there; the editor's `draft_people` adds the QC classes. Printed from the browser, as the report.
- **Changes are worded as the pipeline's run summary** (`FamilyPedigree.describe`): `S3: father S1 -> DUMMY_S3_F`, `DUMMY_S3_F: added`.
- Test data: `datagen.snapshot_rows(members, *extra)` writes a step's rows by the pipeline's rule; the sampleqc, genomeqc and denovo family snapshots now hold them, so nothing waits for a rerun in the default data.

### Family Page Tabs
- **Both family pages draw their analysis tabs through `components/family_tabs.family_tab_bar`.** `FAMILY_TAB_GROUPS` holds each group's caption and its tabs' link value (`tab=`), label and tooltip; `FAMILY_TABS`, `TAB_LABELS` and `TAB_SLUGS` derive from it, and a tab mapped to false in the availability dict is disabled. A new tab is one entry there plus its panel and lazy loader on each page.
- **The Wombat tab is SNVs/Indels** (`snvs`). A renamed tab's old link value goes in `TAB_ALIASES` (`wombat` → `snvs`), and the pages read `tab=` through `canonical_tab`.
- **The row matches the cohort page's** (`components/tab_groups.py`): no icons, a small caption over each group (`tab_group`), mixed case, dense tabs from the left (`GROUPED_TABS_PROPS`). The groups and their order are the cohort page's: Checks (Sample QC), Variants (SNVs/Indels, SVs, STR, UPD), Background (Ancestry, ROH, Polygenic Scores). Sample QC is first, but the page still opens on SNVs/Indels (`default_tab_name`). The icons were dropped in v0.39.0: an icon for ROH or STR said nothing.
- **A `ui.tab` inside a group is still a tab of the bar**: it registers with the nearest `q-tabs` through the group's `div`s, so the bar's value, `update:modelValue` and the panels work unchanged.
- **A tooltip is text `should_see` matches.** The ROH tab's "Runs of homozygosity" is also the panel's title, so a test waiting for the panel waits for something only the panel draws.

### Family Page Review Badges
- **The family tabs carry a count of what waits for a curator** (`components/tab_badge.py`, counts in `utils/review_counts.py`), on both family pages, through `family_review_badges`.
  - **Sample QC**: unacknowledged high/medium flags of the family files, `UPD_CANDIDATE` left to the UPD tab when it exists. Shown to every role, as the header badge.
  - **SNVs/Indels** (mark `snvs-tab-count`): de novo calls aggregated per variant and carrier, HIGH or MODERATE impact (`VEP_IMPACT`, else the consequences' impact), never validated (`UNVALIDATED`). With the denovo step's tables: its config's strict (or not assessed) candidates. Without them: every `wombat.DENOVO_CONFIG` config after the LCR/gnomAD/MNV exclusions. A family whose de novo table Sample QC holds back is not read: `!`.
  - **STR**: flagged genotypes (`flagged_expr`) of the newest catalog, no validation.
  - **UPD**: `upd.count_to_review`.
  - Wombat, STR and UPD are for `can_write()` only, resolved on the event loop before the counts run in `asyncio.to_thread`.
- **Counts cover the whole family**, not the member selection, and are computed once the page is drawn (`ui.timer(0.2)`), so a page never waits on them.
- **Refresh points**:
  - `on_acks_changed` (QC, UPD);
  - `render_wombat_tab(on_review_changed=…)`, called from `refresh_curation` after a variant dialog save;
  - `render_str_tab(on_review_changed=…)`, called from the locus dialog's `on_saved`.

  A badge refresh's timer hangs from the tab, which outlives the dialog that saved.
- Shared with the Wombat tab: `utils/wombat.py` holds the family result files' discovery and reading (`find_wombat_files`, `load_config_frame`, `variant_key`, `DENOVO_CONFIG`), re-exported by `wombat_rows`. The de novo gate flags come from `sample_qc_acks.unexplained_denovo_risk`, and the STR curated statuses from `str_validation.curated_statuses`.

### Family Diagnostics Panel
- One renderer, `components/diagnostics_panel.py`, for the cohort and standalone family pages; each page wraps it in its own refreshable.
- STR and UPD keys are unreadable, so `utils/diagnostic_labels.py` words them: an STR is its gene badge, an `STR` tag and the count (`CuratedRepeats` of the member's standing validation, else `~` + the longer allele from the newest catalog's family genotype table); a UPD is `UPD` + `chr<N> · <origin>`. The raw key is in the tooltip.

### Vedette SV Calls
Vedette is Alex's SV caller (github.com/darkxander/vedette), under active development; the tab is there to visualise its calls and record visual validations that feed back into it. Requests to Alex were issues #5–#11 there (all closed by v1.9.0); `vedette labels --from-genetics-viz validations/vedette.tsv` reads our validations back as training labels.
- **Files** (from Vedette v1.4.0, format in its `docs/EXPORT.md`): family `svs/vedette/<FID>.vedette.tsv` (one row per row of Vedette's events page), `<FID>.vedette.qc.tsv` (one row per member: wave score, carrier load, de novo candidates, karyotype, `flags`), `<FID>.vedette.json` (`schema_version`, `vedette_version`, `model_version`, cohort, `pon`, `members` order, `parquet` name) and, from v1.9.0 where pyarrow is installed, `<FID>.vedette.parquet`; sample `svs/vedette/<ID>.filtered.vcf.gz` (the sample's own calls, every FILTER kept). No cohort level yet, so Vedette is not in the gene index or Search.
- **The Parquet file is preferred** (`read_family_rows`): one line per row × member, typed, folded back by `rows_from_parquet` into the TSV's row dicts (`GT_<id>`-style keys, lists `;`-joined) plus `OwnDetails` (`{(row_id, member): own call fields}`: id, FILTER, bounds, CIPOS/CIEND, GQ, PE, SR, CG, RD, CNQ, model). Both paths must give the same rows (`test_parquet_and_tsv_give_the_same_rows`). Parquet values are typed: an integer 0 is falsy, so never display a cell with `x or "-"` (`_shown`).
- **Rows are events, not sites**: cohort sites carried by the same members are collapsed (`row_id`, `site_ids`, `n_sites`), plus `kind` = `subthreshold` (NDD genes, `filter` LowSQ), `mosaic` and `chromosomal` rows (no GT, `mosaic_cn`). From v1.7.0 the export also holds the rare non-PASS cohort sites (carried by ≤ 20 % of the cohort, `export_nonpass_max_carrier_frac`): `is_nonpass` (a `site` row whose `filter` is not PASS) rows are hidden unless the "Show non-PASS sites" switch is on (`show_nonpass`, link `vnp`, independent of the presets).
- **Coordinates**: `start`/`end` are the first and last affected bases, 1-based inclusive; `vcf_pos = start - 1`. VCFs with `##vedette_coords=padded` (v1.4.0+) have `POS` on the base before and `END` on the last base; older ones `POS` = first affected base, `END` = last + 1. `find_own_call` converts either. The reads give 0-based positions: the suggested start is the last base kept (= `vcf_pos`). `SVSource.start_offset` (1 for Vedette, 0 for WisecondorX, whose curated starts were saved that way) adds to the start shown and saved, so Vedette curated starts are first affected bases like its `start`.
- **Read by column name, degrading** (`utils/vedette.py`): pre-v1.4 exports (one row per PASS site, `site_id`, French `impact` label, `genes` truncated with `(+N)`, an `inheritance` column, no JSON: caller from the family VCF's `##source=`) still open. Prefer `impact_score/level/reason`, `exonic_genes`/`genic_genes` (`;` lists, exonic badges bordered), `INH_<id>`/`DNQ_<id>`, `carriers`. Codes are worded in `config/vedette.yaml` (`benign` impact, `karyotype_<K>` inheritance, kinds, own call, support, QC flags); an unknown one is shown raw, grey. `supported_schema_version` warns above the table when the JSON's is newer.
- **One row per carrier** (the `carriers` column: GT ≥ 0/1 and GQ ≥ Vedette's minimum; before v1.4, alt GT), so validation, diagnostics and the dialog work per sample as for WisecondorX. From v1.9.0 the carrier's own call is per member (`OWN_/SQ_/SRC_/PRECISE_<id>`): the Own call column and filter are the carrier's (`_own_call_is_member`), else the event's `own_call`. Rows and gene badges are built off-thread (`_load`), the scorer taken on the event loop first. The table is cached on `(path, mtime, size)` (`clear_vedette_cache`).
- **Filters never drop a row for a missing value** (no GQ on a mosaic row, no carrier fraction on a sub-threshold one). The frequency filter is the cohort carrier fraction (`max_freq`, link `vfreq`).
- **Row ids are renumbered when the cohort is merged again**: validations and diagnostics are keyed by `site_key` (v1.9.0: `chr:start-end:del|dup` of the row's cohort site, its own bounds, stable while a merge keeps them), else the row's locus (`_locus_key`). `SVSource.variant_key` passes it to the dialog. Across a re-merge, Vedette's `site-map` gives old → new keys (`same`, `moved`, `split`, `merged`…); Alex publishes it at the next re-merge, and nothing reads it yet. The row id is only a note in the validation row.
- **Validations in `validations/vedette.tsv`** (`VEDETTE_VALIDATION_HEADER`: the SV columns + `SiteId CallerVersion ModelVersion`; caller `vedette-<version>` and model from the JSON), a `VALIDATION_KINDS` entry of its own so Alex gets a clean, versioned set. The family report and the dialog's suggestions merge it with `svs.tsv` (`merge_validation_maps`). Diagnostics stay in the shared `diagnostics/svs.tsv`.
- **`show_sv_dialog(source=SVSource(...))`** carries what differs by caller: validation kind and extra columns, details card, breakpoint precision (narrow read windows), title. The default is WisecondorX's. `resolve_own_call` takes the carrier's own call from the export (Parquet: whole; TSV `OWN_` = `none`: nothing to look up) and only for older exports matches it in the sample's VCF (one tabix query, by type and ≥ 50 % reciprocal overlap). `row_is_precise`: `PRECISE_<id>`, else the own call's, else rows under 1 kb.
- **pysam raises on an INFO key the header does not declare**: read optional keys through `_info`.
- **Read evidence draws CIGAR `D`/`I` of ≥ `min_indel_bp`** (`InRead`, grouped by position and length): small events have no pair, split or clip. For a deletion they count as split reads at each end in the suggested breakpoints.
- Test data: `datagen.VEDETTE_ROWS` (family 1001 in the v1.9 format, one row of each kind plus a non-PASS site, written by default as TSV and Parquet), its QC file and JSON (S3 flagged), S3's own-call VCF (padded), and `CRAM_INREAD_DELS` reads holding `VEDETTE_DEL` in their CIGAR. `VEDETTE_SITE_KEY` is wider than site 1's row on purpose.

### CNV Coverage and SV Read Evidence
- The ghfc-ngs `cnv_coverage` step writes `samples/.../<ID>/coverage/<ID>.cov_b<bin>q<mapq>.{bw,json}`: every bin holds `2 × depth / autosomal median` (about 2 on a normal autosome). The `cnv_panel` step of a controls cohort writes a versioned, self-contained `cohorts/<controls>/cnv_panel/v<version>/` (per-bin median and unscaled MAD for `autosomes`, `XX`, `XY`, plus `panel.json` with PARs and expected ploidy). A cohort fitted against it gets `cohorts/<NAME>/coverage/<NAME>.<tag>.panel_v<version>.fit.tsv`. The contract is documented in ghfc-ngs `documentation/params.md`.
- **The division happens here, at view time**: CN = expected ploidy × sample / panel median, PAR-aware, in `utils/cnv_coverage.py`. Per-sample files never depend on a panel, so a new panel version rewrites nothing.
- **The panel is found through the cohort's params `cnv_panel_dir`, a compute-cluster path.** Its tail from the last `cohorts/` component is resolved under the cohort's data directory, then every other one. The cluster path itself is never opened (same trap as the run record and the ancestry `.qc.json`). A controls cohort has no `cnv_panel_dir`; it uses the panel it built itself (`cnv_panel_version`, else the newest `cnv_panel/v*/`).
- **Degraded, never an error**: no panel, a panel of another tag (`cov_b<bin>q<mapq>`, refused, as in the pipeline), or a karyotype with no panel group (`X0`, `XXY`, unknown) → the sample's own normalised value, unmasked, with the reason shown. Karyotype is the fit table's `coverage_sex`, else the pedigree sex.
- Bins are masked (drawn grey) where the panel median is under half the expected ploidy or its MAD/median is high (`config/cnv.yaml`). Regions over `max_exact_bins` are read through the bigWig zoom levels.
- **The SV dialog** opens on Read Evidence, which carries the mode toggle. When the carrier has a coverage bigWig, the Copy Number View below it replaces the IGV by1000 bedgraph view. It is folded, and drawn on first opening: an ECharts chart built in a hidden panel has no width. The bedgraph view stays as the fallback, open. Neither chart zooms on the mouse wheel, so the wheel scrolls the dialog. Drag pans, and the rulers zoom.
- **Read windows** (`utils/sv_view.py`): one window when span × 1.6 fits igv.js's alignment `visibilityWindow` (also set on the CRAM tracks from the same config value), one per breakpoint otherwise; INS is always one, BND always two. An Auto / Single / Split toggle overrides it, for the IGV reads and the read-evidence panel alike.
- **Read evidence** (`utils/sv_reads.py`) decodes the CRAMs with pysam and needs `reference_fasta` in the app config (the exact alignment FASTA). Without it nothing is opened: htslib would otherwise try to download the reference by MD5.
  - Only evidence is drawn: discordant pairs, split reads and soft-clip clusters (several reads clipped at one position).
  - Pairs are classed by orientation against the window's own proper-pair insert sizes: `+-` with a large insert = deletion, `-+` duplication, `++`/`--` inversion, mate on another chromosome. Concordant pairs are only counted; their 2.5-97.5 percentile band is shaded.
  - A split read is its two aligned parts joined by an arc above them, so the link never hides a pair's line.
  - Each sample × window cell is one grid with two y-axes: insert size (log) on the left and, on the right, the CN from the coverage bigWig divided by the panel, the same computation as the Copy Number View.
  - A sample with coverage but no CRAM still gets its CN row. A window wider than the visibility window draws CN only.
  - A thin read-depth line at 25 bp (`2 × depth / autosomal_median_depth` from the coverage JSON) shows the exact drop point. It is not divided by the panel.
  - "Show MAPQ 0 reads" keeps the reads under `min_mapq` and draws their evidence faded. They never count toward the insert band, the depth or the suggestions.
  - **Suggested breakpoints** (`utils/sv_breakpoints.py`), from the carrier, per side: evidence left of the SV's middle speaks for its start, right of it for its end.
    - **Exact position**: clip clusters and split-read junctions, scored clips + 2 × splits. A split read's junction is the clipped end of each part, from its own CIGAR and the SA CIGAR. That is right for every type: a deletion joins `a2`→`b1`, a tandem duplication `b2`→`a1`.
    - **Orientation**: evidence that doesn't fit the type counts for `mismatch_weight`. A deletion's start is clipped after it; a duplication is the reverse. Inversions and others accept either side.
    - **Interval**: from the type's discordant pairs. A `+` read lies before the breakpoint within `insert_hi`, a `-` read after it.
    - **Copy number is used only for types that change it (`del`, `dup`)**, never for inversions or other balanced events. It checks the exact position, and when the reads give none it supplies an approximate coverage step (dashed, "≈", "Use approximate position"). The coverage is the 25 bp depth without panel-masked bins, so repeats losing depth in everyone don't read as steps.
    - **Confidence**: strong, moderate or weak. A coverage contradiction, or an exact position outside the pairs' interval, makes it weak.
    - "Use as curated position" fills the curated inputs.
  - **The read views follow the curated bounds.** `read_bounds` and `_update_roi`'s `roi_coords` are the same dict, so every curated change (suggestion chips, IGV pickers, the refresh button, the reads' suggestion) re-centres IGV, Read Evidence and the Copy Number View.
  - WisecondorX calls get ±5 kb breakpoint windows (`flank_bp_imprecise`, `precise=False`), since their breakpoints sit on 1 kb bins. A split read's far part can suggest an end outside the window; taking it re-centres the windows there.
- Tests: `datagen` writes the bigWigs with pybigtools and a **reference-free CRAM** (`no_ref=1`) with a hand-written `.fai`, since a real chr2 reference would be 240 MB.

### Cohort Discovery
- **A cohort is a directory under `cohorts/`.** The directory name is the cohort name. Hidden directories (a NAS `.snapshot`) are skipped.
- `DataStore.cohorts` holds only cohorts with a readable `<NAME>.pedigree.tsv`, so the header dropdown, validation pages and search keep their existing meaning. Directories with no pedigree — or one that will not parse — land in `DataStore.incomplete_cohorts` as `CohortStub`s, carded on the home page but not navigable. `CohortStub.reason` is `"missing"` or `"unreadable"`; the two need different fixes and do not share a message.
- Discovery lives in **one** place, `_iter_cohort_dirs()`. It used to be hand-copied into `load()`, `take_snapshot()` and `reload()`, which drifted the moment the rules changed.
- `take_snapshot()` watches the pedigree **and** `<NAME>.params.yml` (the newer of the two mtimes, `0.0` when absent) and keys every cohort directory, so a directory appearing or gaining a pedigree registers as a change.
- `create_header(cohort_name)` must never be handed a name outside `store.cohorts` as the select's value — Quasar raises on it, which used to turn every unknown cohort in the URL into a 500 and left the pages' own "not found" branches unreachable.

### Cohort Run State
Two files may sit in a cohort directory, both written by the ghfc-ngs workflow: `<NAME>.params.yml` (parameters, including the `steps:` array) and `.ghfc-ngs.state.json` (the run record — schema documented in `COHORT_STATE.md` of bourgeron-lab/ghfc-ngs). The cohort page's Workflow tab (run status, parameters, params reference) exists only when the params file does.

- **genetics-viz never measures completion.** The record already carries per-step `{done, total, pct}` with denominators from the pedigree. Re-deriving it would duplicate the pipeline's planning logic and drift from it, and it is unaffordable anyway: a cold existence check costs ~5.6 ms over the CIFS mount, so a 2700-sample cohort is about a minute of stat calls.
- **A null step is unmeasured, not 0%.** `ancestry` is null unless requested, `extractor` always is. Render "not measured", never an empty bar.
- **Check `schema_version` before trusting any field.** Anything but `1` is refused whole rather than half-interpreted.
- **A missing state file means "not run since run recording shipped"**, not "never run". The Status view then lists the params file's steps, all unmeasured.
- **`status: running` is not a liveness signal.** There is no heartbeat; a killed run leaves the record behind until the next run closes it out. `RunRecord.is_stale_running` flags one older than 24 h.
- **The `path` fields in the record are compute-cluster paths** (`/pasteur/helix/...`) and do not resolve locally, so drift detection hashes `cohort.pedigree_file` / `cohort.params_file`. Same trap as the ancestry `.qc.json` files.
- The record is deliberately **not** in the change-monitor snapshot: the pipeline rewrites it on every run and `reload()` re-parses every pedigree in the data directory. `read_cohort_state` caches on the file's own `(path, mtime, size)` and picks up rewrites by itself.
- `config/pipeline_steps.yaml` holds step labels, display order and status chip colours only. A step the pipeline reports but the file does not know is appended under its raw name.
- **A step's completion, drift and start come from its own last runs**, not `last_run`: a `--steps sampleqc,genomeqc` run becomes `last_run` and says nothing of the other steps. `CohortState.by_step` is ghfc-ngs's `last_run_by_step` (bourgeron-lab/ghfc-ngs#71: per step, `last` and `last_success`, finished runs that requested and planned it), or `derive_by_step` over the last success and the history for an older file (`by_step_recorded` False). A live run counts for the steps it was asked for (`live_run_for`).
  - `build_step_rows`: each row from its step's run; `latest` False shows "as of <date>". `step_drift`: per step, against its last success.
  - **A pedigree edit waits for every pedigree-reading step** (`PEDIGREE_STEPS`, the steps with an `InputsSnapshot` in ghfc-ngs, among the params' steps): `pedigree_run_started(cohort)` is the earliest of their starts, `None` when one never ran. `family_status.cohort_edits` and the Workflow tab's waiting families pass it to `edit_statuses` / `waiting_families`; never `last_run.started_at`.

### Run Requests
Curators ask for a ghfc-ngs run; **genetics-viz never starts one**. The VM has no cluster account, and a password or a key per user would put a credential in the app. The request is a file the ghfc-ngs-swiftbar dashboard reads over its own SSH poll, as the person who runs the pipeline (vault: `Workflow run requests`).
- **Files**, in `cohorts/<C>/.ghfc-ngs.requests/`. Each writer creates only its own and never edits another's:
  - `<ID>.json`, the request (`schema_version` 1, `steps` or `null` for a normal run, `clean_stale_families`, `reason`, `families`, `origin`, `requested_by`, `requested_at`), and `<ID>.withdrawn.json`, written here;
  - `<ID>.handled.json`, `started` (`job_id`, `steps_sent`, `clean_sent`, `covered_by`) or `dismissed` (`note`), written by swiftbar.

  `<ID>` is `YYYYMMDDTHHMMSS-<6 of [a-z0-9]>` (`ID_PATTERN`), checked before it goes in a path. Writes go to a temporary file in the folder, then a rename. The folder's `.lock` (`flock`) keeps two clicks from both passing the one-open-request check.
- **`utils/run_requests.py`**:
  - `write_request(_async)` checks `can_write()` and refuses a missing reason, a step outside the params file's (`allowed_steps`, which splits `"a,b"` as ghfc-ngs does) and a second open request by the same user on the cohort;
  - `withdraw_request(_async)`: the author's, or any for `is_admin()`;
  - `read_requests` gives the state (`open`, `withdrawn`, `dismissed`, `started`) and, for a started one, the outcome: the run record whose `launch.slurm_job_id` (`RunRecord.slurm_job_id`) is the handled `job_id`. A file that is not a request (other schema, bad ID, not JSON) is skipped. Free text is capped.
- **Opt-in per data directory**: `DataDirectoryConfig.run_requests` (Admin → Directories switch `admin-run-requests-<i>`, saved only when on), `requests_enabled(data_dir)`. Off without a loaded config. Only for a directory swiftbar polls: nothing else reads requests.
- **Curators and administrators only** (`components/run_request.requests_shown`):
  - Workflow tab `workflow-requests` (the button and the last 30 days, `recent`);
  - the Workflow button's count (`request_pills`);
  - cohort Sample QC `sampleqc-request-run` (sampleqc/genomeqc among the params);
  - the family Sample QC note's `pedigree-request-run-<cohort>` (one per cohort listing the family; `_request_steps`: the pending steps among the params, `wombat` with `wombat_cohort`, else a normal run);
  - `workflow-pedigree-request` in the waiting families (a normal run, clean when one needs it).
- **Markers**: the dialog is `run-request-dialog` (`-mode`, `-step-<step>`, `-clean`, `-reason`, `-command`, `-send`). `render_request_control` shows the button, or the user's open request (`<mark>-open`) with `run-request-withdraw-<id>`; the list's rows are `run-request-<id>`, with `run-request-list-withdraw-<id>`.
- Test data: `datagen.add_run_requests(data_dir)` (opt-in): DEMO gets an open request by admin, one started with fewer steps (job `RUN_REQUEST_JOB`, also written in the run record) and one dismissed.

### Cohort Page Tabs
One row of tabs in captioned groups, without icons (`_TAB_GROUPS` in `cohort.py`: caption, then each tab's label, `tab=` link value, tooltip; drawn through `components/tab_groups.py`): Browse (`Families`, the default, and `Individuals`), Checks (`Sample QC`, `De novo QC`), Variants (`SNVs/Indels`, `SVs`, `STR`, `UPD`), Background (`Ancestry`, `ROH`, `Polygenic Scores`). Up to v0.34.0 the analyses sat in a column beside the table; once the table had nine columns (1,366 px) that column was 138 px wide in a 1,600 px window. The Statistics tab (carriers per geneset) was removed in v0.39.0; `tab=stats` falls back to Individuals.
- **Workflow is a button on the title line** (`cohort-workflow`, `_WORKFLOW`), not a tab: it is about the run, not the cohort. Its panel is `ui.tab_panel("Workflow")` with no `ui.tab`, opened by `show_tab("Workflow")`, so `tab=workflow` links and `open_tab(..., "Workflow")` still work; the row then shows no tab selected, and the button turns `unelevated` (`mark_workflow`).
- **The Families tab** (`utils/family_rows.py`) has one row per family. Pedigree, Sample QC, UPD and diagnostics columns are built at page build; ancestry (founders' regions), the strongest FROH grade, de novo QC verdicts, the de novo review count and the conclusions come from `load_family_facts`, run in `asyncio.to_thread` by `load_families` after the page is drawn, then `build_family_rows` + `recompute`. STR flags arrive with `load_badges`, which rebuilds the rows the same way.
  - **One scope for both tables** (`recompute`): individuals = those passing the individuals filters, in families passing the families filters; families = those passing their filters with, when individuals filters are on, a member passing them. `filtered_state["individuals"]` is the scope every analysis tab follows. Family filters are in the link under `f*` keys (`_FAMILY_FILTER_LINK`); "Edit filters" opens `state["filter_tab"]`, the table last filtered.
  - **The family Sample QC cell is one badge** (`compact_qc_cell`: worst severity, members flagged, categories in the tooltip; `_qc_cats` keeps the categories for the filter) plus "Cell line". One badge per category, as for individuals, made the column 590 px wide.
  - **The de novo review count is the cohort `dnm` table's** strict rows with a HIGH or MODERATE consequence and no non-ignored validation (`dn_open_counts`). The family tab's badge also counts Wombat rows the step did not assess; reading every family's Wombat table here would be too slow.
  - **The family page's arrows walk the Families tab's order** (`family_rows.save_family_order` / `family_order`, in `app.storage.user["family_order"][<cohort>]`). `remember_order` saves the table's rows (filtered, and sorted by `DataTable.sort_rows`, which `recompute` now applies before `update_data`) after each filter, sort or page change, and forgets them when neither filter nor sort is on, including on a bare page load. The family page uses the list only when it holds the family, else every family by ID; `family-position` shows "N of M".
  - **Width**: it fitted a 1,280 px window for curators until the Pedigree column (see Pedigree Edits) because the text filters are `size`d (`filter.size`, an opt-in in `data_table.js`: an input is about 20 characters wide by default, which set the column floor) and there is no Members column (Affected says `1 / 3`).
- **Every variant kind has a cohort view under Variants.** `De novo QC` (`tab=denovo`, `cohort_denovo_tab.py`) is the step's QC only; the candidates are `SNVs/Indels` (`tab=snvs`, `utils/cohort_snvs.py`): the cohort `dnm` rows of the children shown, strict HIGH/MODERATE by default, with `denovo.load_curation`'s status and diagnostic, and `recurrent_genes` (two or more families). Gene badges are built in the loading thread with the scorer taken on the loop.
  - **`SVs`** (`tab=svs`, `utils/cohort_svs.py`) reads `sample_counts.cohort_wcx_file` through `scan_calls`: `pl.scan_csv` of the needed columns only, filtered to the permissive thresholds (`past_permissive()`, the report's filter too) and 100 kb before a streaming collect (`tsv.collect_streaming`: polars 1.17, which the project allows, refuses `engine="streaming"` and takes `streaming=True`; CI's polars 1.17.1 job caught it), cached on `(path, mtime, size)` (`clear_cache` in `conftest`). Never read a cohort WisecondorX table eagerly: CANDY_mpx's took 13.4 GB. A call row is per call and family of its sample (a sample can sit in two); `Variant` is the SV validation key `chr:start-end:del|dup`. Whole-chromosome events reuse `family_report.whole_chromosome_events`, so the tab and the report agree; recurrent regions merge overlapping same-direction calls of two or more families.
  - **Polygenic Scores lists its outliers** (`ancestry.pgs_outlier_rows`): the individuals of the filter and regions at |z| ≥ the cut-off (`PGS_OUTLIER_THRESHOLDS`) for the plotted trait or every trait, percentile Φ(z) through `report_svg.percentile`/`ordinal`, so the tab and the report word a score alike. Nothing is read twice: it pairs the tab's own z-score rows with the individuals (`scoped()`).
  - Eleven tabs fit a 1,280 px window (the last ends at 1,216 of 1,240 px); a twelfth will need the labels shortened or two rows.
- **The tab row and the filter bar are one sticky block** under the fixed page header. The header's height changes when its menus wrap, so `_HEADER_HEIGHT_JS` keeps `--gv-header-height` up to date (ResizeObserver) and the block's `top` reads it.
- **The filter bar** shows the count and one removable chip per filter (`_filter_chips`, worded from the column definitions). A chip's cross, or "Clear all", calls `apply_filters(..., to_table=True)`, which sets the table's filter row through `DataTable.set_filters` (no filter event follows, so the page applies them itself). The filters stay in the table header; "Edit filters" switches to Individuals. Workflow hides the bar: it is about the cohort, not the filter.
- **Every tab but Families, Individuals and Workflow loads on first view** through `lazy_content(label, render, ...)`, which registers an opener in `openers`; `Workflow`'s opener is its own, `async` one, because it warms the run-record and document caches off-thread. Tabs without data are disabled; `Workflow` is *created only when the cohort has a params file* — absent, not disabled.
- **The tab change is `tabs.on_value_change`**, so `tabs.set_value(label)` from Python (Edit filters) runs the same handler as a click. Tab values are labels (`ui.tab_panels(tabs, value=<label>)`), never `ui.tab` objects.
- **Links**: `tab=<slug>`. Older links named the tab `panel=<slug>` (the right-hand panel) or `view=workflow`; `_linked_tab` still reads them and the page rewrites the address to `tab=`.
- **Badges** (`TabPills`, counts in `utils/review_counts.py`) follow the filters, so they are rebuilt in memory by `update_pills` at each filter change; only the STR flags and the run record are read, once, off-thread (`load_badges`). Families and Individuals: how many the filters leave (`filter_pills`, either noun). Sample QC: individuals by worst unacknowledged flag (`qc_worst_by_individual`, `UPD_CANDIDATE` left to UPD when the cohort has that tab), high and medium apart. UPD: High/Medium candidate chromosomes no curator validated present or absent (the rows' `_upd_open`, from `upd_column_fields`). STR: open flagged genotypes of the newest catalog, from the cohort `flags.tsv` (else the genotype table). Workflow (on its button): the last run's status icon and colour (`workflow_pills`; nothing without a record or for an unknown status). UPD and STR are for `can_write()` only. Colours are the Sample QC severities', as on the family page.
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

### GeneTrek (gene badges, genesets, modes of inheritance)
GeneTrek (`genetrek_url` in the app config, default the dev instance) is the lab's gene database. The app reads it **as a snapshot**: the whole table is one ~4 MB parquet (`/api/download/database`), downloaded once per GeneTrek `release` (`/api/health`) and cached under `genetrek_cache_dir` (default `<config dir>/.genetrek/`). `utils/genetrek.py` holds it.
- **Never call the API per gene.** Tables render thousands of badges whose tooltips are computed in Python before the table is sent. The snapshot also means GeneTrek never learns which genes a user looks at; only a click reveals one.
- **Startup loads the disk cache, then checks the release, both off-thread** (`start_refresh`), and re-checks every `refresh_hours`. Pages render with neutral badges until the snapshot is in. A failed refresh keeps the previous snapshot and records the error; it never raises (`except Exception` is deliberate there).
- **Resolve identifiers through the snapshot** (`resolve`): symbol, then ENSG, then previous symbol, then alias — the last two only when they name exactly one gene. VEP symbols lag HGNC (`KIAA1109` is `BLTP1`). Filters use `matching_identifiers` to expand a set with those old names instead of resolving every cell.
- **Geneset versions coexist** (`hc_ndd_v1` … `hc_ndd_v7`). Their alias (`hc_ndd`) names the current one, so preferences are keyed by `GenesetInfo.key` (alias first). A key GeneTrek drops is kept in the file and ignored.
- **Badges are coloured by tiers, not a sum.** Each user assigns sets to Tier 1-4, "show" (tooltip only) or hidden; a gene takes the best tier's colour. A sum would double-count versions and roll-ups. `GeneScorer` is per user: call `get_gene_scorer()` on the event loop (it reads `app.storage.user`), never inside `asyncio.to_thread`.
- **Build badges with `build_gene_badges` / `gene_badge_html`** (`utils/gene_scoring.py`), never by hand: they carry `href` (GeneTrek `?search=SYMBOL`), `textColor` and the multi-line tooltip, and cap a cell at six with a "+N genes" badge linking to GeneTrek's table (`search_url`, a `#search=` fragment so patient gene lists and loci stay out of server logs). SV regions (`chr:start-end`, the `curated_locus` cell) link the same way through a row's `_genetrek_url`; the ROH genes table and the STR `Locus` columns (`str_eh.add_locus_badges`, the locus coloured as its gene) use gene badges too.
- **Multi-line table tooltips** use `data-gene-tip`, shown by the floating tooltip in `data_table.js`; the CSS `data-tooltip` is one clipped line. `badge_list` cells switch automatically when a tooltip has a newline.
- **Preferences live in `<config stem>.preferences.yaml`, not the config file.** `save_config` rewrites the config whole from memory and truncates before locking; profile saves there would race with admin saves. `user_preferences.py` locks a `.lock` file, re-reads, changes one entry and `os.replace`s. Saving checks the user itself; the default is admin-only.
- **MOI flag** (`utils/inheritance.py`, `config/inheritance.yaml`): a gene's modes are pooled from every GeneTrek field in category `Inheritance`, because the combined `inheritance` holds one mode. Recessive second hits are looked for only among the rows of the same table (`annotate_moi`), and the tooltip says so. A male's chrX call outside PAR is hemizygous whatever the GT says.

### Variants Per Sample (family Stats dialogs)
The Wombat and WisecondorX Stats dialogs chart each selected member's count in the current table view, then place each member on the cohort's per-sample distribution (`utils/sample_counts.py`, drawn by `components/sample_counts.py`).
- **The cohort comparison is unfiltered on both sides.** It counts the cohort-level file (`cohorts/<C>/wombat/<C>.rare.<vep>[.annotated].<config>.results.tsv`, or `<C>_aberrations[.annotated].bed`), never the table's filters. Older cohorts have only the `.annotated.` name, so the lookup tries both.
- **Several cohort de novo tables have no header** (CANDY_mpx, CANDY_mpx_ghent, PMS): a blank line, then rows. Read them through `read_cohort_wombat`: a family file of the same config lends its header when the widths match; otherwise `MissingHeaderError` is shown, never a cohort of zeros or a first row taken for column names. The Stats dialog borrows from the family's own file; the Search page tries the cohort's families in turn (`family_header_sources`). The default `datagen` data has one, DEMO's `strict` table (`build_headerless_cohort_wombat`).
- **Search lists cohort Wombat tables through `cohort_wombat_sources`**, sorted, so its default source does not change when a file is added. The `annotated` part of the name is not part of the config.
- **Founders and non-founders are compared within their own group** (`Sample.is_founder`), with every pedigree sample zero-filled. A de novo config has rows only for children, and pooling them with the parents' zeros would flag any child with two calls.
- An outlier is outside Q1-Q3 by more than 1.5 × IQR (`OUTLIER_IQR_FACTOR`) and on the far side of the median.
- Counts load off-thread and are cached on `(path, mtime, size)`. The comparison sits outside the dialog's refreshable section, so the SNV/Indel toggle does not recount.

### Tab Downloads (TSV exports)
Every cohort tab but Individuals and Workflow, and every family tab, has a `Download` button (`components/download_button.py`, marked `download-<tab>`; menu entries `download-<tab>-<key>`, zip `download-<tab>-zip`).
- **Scope is what is shown**: the cohort table's filtered individuals plus the tab's own selectors (PGS regions, ROH reference groups, bundle/tag, catalog), or the family's selected members. Columns are all of them, not just the plotted one.
- **Built from explicit column lists, never by dumping row dicts.** Badge columns hold `""`, and the row dicts carry fields readers must not see. The builders sit next to their data (`ancestry_export`/`pgs_export`, `froh_export`/`kinship_export`/`segments_export`, `genotypes_export`, `sample_qc_export`) and serve both levels, so a family file and a cohort file of one kind have the same columns.
- **Role-aware, as the tabs**: readers get `reader_label`, the neutral ROH grade, no evidence, no observed relationship, no `froh_category`/parents' kinship, an "Atypical" karyotype, declared pairs only, and an ack note only for `reader_note` categories.
- Per-sample files start `FID, IID, sex, phenotype` (pedigree codes via `export.person`); the PGS export is samples × traits at both levels, not the family table's transpose. ROH percentiles are headed `_vs_ref_group`: they are against the reference group, not the cohort.
- Builders run in `asyncio.to_thread` at click time: resolve `can_write()` at render and close over it. Content reaches `ui.download` as **bytes** — the test `UserDownload` fetches a `str` as a URL. A header-only result shows "Nothing to download".
- Wombat and SV exports are the table's current `DataTable.rows` (filtered and sorted in Python) with the shown columns, through `is_export_column`, as the Search page's Download TSV.

### Shareable Links (view state in the URL)
The family pages (both), the cohort page and the Search page keep their view in the URL query, so the address bar is a link a curator can send (`components/page_link.py`).
- **Read at build, written as it changes.** A page takes `request: Request`, builds `PageLink.from_request(request)` and opens with its values; every change calls `link.set(...)` then `link.push()` (`history.replaceState`: no reload, no history entry per click). Only non-default values are written, so a page opened bare is in its usual state.
- **Keys:** family pages `tab`, `m` (members), the Wombat tab's shown config `wcfg` and its `wlcr`/`wgnomad`/`wmnv`/`wgs`/`wimp`/`wval`/`wsort`/`wpage`, the SV tab `svtab` and `svcalls`, the Vedette sub-tab `vpre` (preset, or `custom` with `vtype`/`vkind`/`vown`/`vsize`/`vfreq`/`vgq`/`vimp`/`vinh`), `vgs`/`vval` and `vnp` (non-PASS shown), the UPD tab's `updl` (label, when not the newest), `updc` (confidences) and `upd=<sample>:<chrom>` (the review open); cohort page `tab` (older links: `panel`, `view=workflow`), the individuals table's column filters (`fid`, `sid`, `sex`, `pheno`, `diag`, `qc`, `upd`, `father`, `mother`) and the families table's (`ffid`, `fstruct`, `ftrio`, `fmulti`, `fanc`, `froh`, `fqc`, `fdn`, `frev`, `fdiag`, `fconc`); Search `q`, the whole form as compact JSON, written when Search runs and run on open.
- Lists are comma-joined, `-` is an emptied list (not the same as absent). Consequences shorten to impact levels (`wimp=HIGH`). Search names sources by display name, never by directory: a link carries no server path.
- Only the shown Wombat config is in the link; the others keep their state while hidden. A linked value the page cannot honour (a tab without data, an unknown config or member) falls back to the default.
- **The query is in the server's access log**, as the path already is. The fragment would not be, but Python could only read it after the page was built, through a JavaScript round trip.
- **Login keeps the link**: `check_auth` redirects to `/login?next=<path?query>`, and the login goes there after success. `safe_next` accepts only a path on this site (not `//host`, `/\host`, a full URL or `/login`).
- Tables: `DataTable(initial_filters=..., on_state_change=...)` opens with column filter values and reports sort/page changes.

### Pedigree Missing Values
The sentinel set `{"", "0", "-9"}` represents unknown/missing in pedigree fields (parent IDs, sex, phenotype). Defined as `_PED_MISSING` in `search.py` and handled in `models.py` via `treat_missing_as_null`.

### Dummy Pedigree Members
A dummy is an individual who was never sequenced and exists only to link real members: the shared parents of MZ twins, the grandparents joining an aunt to a cousin. ghfc-ngs uses the same convention (pedigree-only members, never scheduled as samples).
- **The ID prefix `DUMMY_` is the only test**: `models.is_dummy_id`. The rest of the name is the barcode of the real child it hangs from, then `F`/`M` steps (`DUMMY_B00IDL0_F`, `DUMMY_B00INOK_MF` = mother's father). A dummy may have its own row (phenotype `0`), dummy parents of its own, or no row at all, named only in a child's father/mother column: `Family.add_implied_dummies` makes it a member (father male, mother female).
- **Both parsers split them out**: `Family.samples` holds the sequenced members only, `Family.dummies` the dummies, `Family.members` both. Everything that reads data or counts iterates `samples`, so it skips dummies without asking; `get_family_members` / `load_family_pedigree` return sequenced members unless `include_dummies=True`. Only pedigree views ask for them: the member selector (via `sequenced_ids` for the selection), the report's family table.
- **A dummy is never selected** in the member selector (no checkbox), so no tab looks it up. Parent IDs are the remaining leak: a lookup by a father/mother ID (the variant dialog's parent CRAMs, `parents_to_validate`) must go through `is_sequenced_parent`, which is false for a missing, `-` or dummy parent.
- **The pedigree still resolves through dummies**: twins sharing dummy parents are siblings wherever parent IDs are compared. `Sample.is_founder` means "no sequenced parent", so twins with dummy parents are founders (they have no de novo calls). A dummy's phenotype is never read: dummies are not in `samples`, so `family_pedigree` sees no phenotype for a dummy parent.
- **Marked with one wording**, `pedigree_labels.DUMMY_*`: a muted row and "dummy" badge in the selector, `dummy_badges` (a grey `badge_list` badge) for parent cells of TanStack tables, a greyed italic row with a dashed rule and a footnote in the report.
- Test data: `datagen.add_dummy_case(data_dir)` adds families 2001 (MZ twins S7/S8 with dummy parents) and 2002 (three generations; S9's dummy father has no row) to DEMO. Opt-in, because it changes DEMO's counts. `tests/ui/test_dummy_pages.py` spies on `sharding.get_entity_path` to catch any sample lookup of a dummy.

### TanStack Table
- `DataTable` in `tanstack_table.py` bridges Python to vanilla JS TanStack Table.
- Per-client JS/CSS injection via `client._dt_scripts_injected` attribute (avoids global set memory leak).
- Tables are created when their container enters the document, with no timeout (`static/js/data_table_mount.js`, one shared MutationObserver). Quasar mounts only the shown tab panel, and a left keep-alive panel sits detached, so a table built for a hidden tab gets its container whenever the tab is opened. `update_data` / `set_column_visibility` calls made before that are queued.
- A tab bar whose panels each hold a table builds each panel on first view (the Wombat sub-tabs, as the family page's tabs).

### Async
- Use `asyncio.to_thread()` to offload blocking I/O (e.g., `pl.read_csv()`) in NiceGUI page handlers.
- NiceGUI runs a single-threaded event loop; blocking calls freeze the UI for all users.

## Conventions

- **Validation statuses are defined once**, in `utils/validation_status.py`: `UNVALIDATED` (`"TODO"`, a row nobody validated), `FILTER_STATUSES` (menu order), their icons and chart colours, `consensus()` (several curators' validations → one status and inheritance) and `status_selected` / `filter_by_status` (an empty selection filters nothing; `in phase MNV` counts as `present`). Never spell `"TODO"` or the status list in a page.
- **Wombat rows are badged and filtered through `utils/variant_rows.py`** (Wombat tab, Search, per-file validation page), with `create_consequence_filter_menu` for the Impacts menu. An empty consequence selection filters nothing, like the validation filter.
- **Validation kinds are listed once**, in `utils/validation_files.py`. The All validations and Statistics pages loop over it, so a new kind is one entry there plus its writer.
- **Diagnostic kinds are listed once**, in `utils/diagnostic_files.py`. The following loop over it: the Diagnostic All and Statistics pages, the cohort and home summaries, and `load_family_diagnostics`. The per-kind table callers (wombat, svs, search, validation/file) name their own file. `save_diagnostic_entry` and `update_diagnostic_ignore_status` check `can_write()` themselves.
- Utility code shared across pages lives in `utils/`, not duplicated in page components.
- **Columns are added to a curation file at the end only.** `append_tsv_rows` extends a header that is a strict prefix of the writer's, under its lock, so an older file gains the columns on its next save; older records keep their shorter lines.
- **Curation files (validations, diagnostics, notes, waves, acks, summaries) are written only through `utils/tsv.py`**: `append_tsv_row(s)` to add records, `rewrite_tsv_rows` (or `set_ignore`) to edit one in place. Both check the header and write under one `flock`, and flush before unlocking: a buffered write reaches the file only at close, and the next writer would find it empty and write the header again. Never `csv.DictWriter`, which quotes a field holding `"`.
- **Read curation files with `quoting=csv.QUOTE_NONE`.** They are written unquoted, and the default reader takes a comment's leading `"` for quoting and drops it. Rows an older `DictWriter` rewrite quoted (a comment holding `"`, in a file whose ignore toggle or note deletion ran before this) now show their quotes doubled.
- YAML configs live in `config/`, loaded by corresponding `utils/*.py` modules.
- All new ClinVar terms must be added to `config/clinvar_colors.yaml`.
- All new VEP consequences must be added to `config/vep_consequences.yaml`.
- Column display/group/sort/drop config goes in `config/column_names.yaml`.
