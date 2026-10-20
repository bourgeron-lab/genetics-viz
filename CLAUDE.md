# CLAUDE.md

Project-level instructions for Claude Code.

## Project Overview

NiceGUI-based web application for genetics cohort data visualization. Uses TanStack tables for dataframes, Polars for data processing, and ECharts for plots. Managed with `uv`.

Deliberate exceptions to "ECharts for plots": `stats_panel.py`, the cohort PGS violins in `cohort_pgs_tab.py` and the cohort ROH violins in `cohort_roh_tab.py` use Plotly. ECharts has no violin or box series, and all three live as sibling tabs on the cohort page.

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
- `tests/conftest.py` resets all module-level state before and after each test: the config singleton, the data store registry, the stat-keyed caches, the GeneTrek snapshot, the preferences cache and the per-user gene scorers. A new module-level cache must be added to `_reset_globals`, or tests will depend on their order.
- An autouse fixture makes any `urlopen` fail, so no test can reach GitHub by accident.
- **UI tests** (`tests/ui/`) use NiceGUI's user simulation: no browser, timers run for real.
  - `tests/ui/app_main.py` purges `genetics_viz` from `sys.modules` and boots the app again for each test. NiceGUI's per-test reset removes every route, so without the purge nothing would be registered from the second test on.
  - The `user` fixture points `Storage.path` at `tmp_path`. The reset otherwise deletes the real `.nicegui/` sessions.
  - It fails a test on any ERROR log, which is how NiceGUI reports an exception in a lazy tab.
- **UI test gotchas:**
  - **Tabs:** switch them with `ui_helpers.open_tab`, which sends `update:modelValue`. NiceGUI stores listeners under the camelCase name, and clicking a `ui.tab` does nothing in simulation. Tab bars carry `.mark("family-tabs")`, `"cohort-panel-tabs"` and `"cohort-page-tabs"` for this.
  - **Table and chart contents** render in JS and are invisible to `should_see`. Assert on what the Python side passes: `ui_helpers.capture_tables` records `DataTable` rows, and `ui.echart(...).options` holds the chart data.
- **CI:**
  - `ci.yml` runs on PRs and on main: ruff, pytest on 3.12/3.13 with the lock, and a `package` job that tests the built wheel installed without the lock.
  - `upstream.yml` runs that unlocked install weekly and opens an issue when it fails. Production (`uvx`) ignores `uv.lock`, so an upstream release can break it with no commit here.
  - `publish.yml` runs `ci.yml` first and checks the tag against the version.
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
    filters.py                    # Reusable filter menus (validation statuses, genesets)
    sample_counts.py              # Stats dialogs: variants per member + cohort comparison
    gene_prefs_editor.py          # Profile card: GeneTrek tiers per geneset, live preview
    header.py                     # Shared page header with data dir dropdown + user menu
    validation_loader.py          # Validation TSV loading and badge logic
    variant_dialog.py             # IGV variant inspection dialog
    sv_dialog.py                  # SV inspection dialog (copy number, IGV, read evidence)
    sv_cn_panel.py                # SV dialog: copy-number panel (sample ÷ control panel)
    sv_read_panel.py              # SV dialog: samplot-like read-evidence panel
    offthread_panel.py            # Container loaded off-thread, superseded by newer refreshes
    diagnostic_dialog.py          # Diagnostic review dialog
  config/                         # YAML config files (loaded at module level)
    ancestry.yaml                 # Ancestry region colors + PGS z-score color scale
    genetrek.yaml                 # GeneTrek refresh period, tier colours, default profile
    inheritance.yaml              # MOI flag: GeneTrek mode strings -> codes, verdict colours, PAR
    roh.yaml                      # ROH grades/wording, pattern + class colours, KING degrees
    sample_qc.yaml                # Sample QC flag wording (curator/reader), categories, plot lines
    str.yaml                      # STR class/status colours, support codes, REViewer filename
    clinvar_colors.yaml           # ClinVar significance term -> color
    cnv.yaml                      # CN mask/display, single-vs-split read view, read-evidence caps/colours
    column_names.yaml             # Column display names, groups, sorting, drop flags
    cytobands_hg38.tsv            # Cytoband data for ideogram rendering
    continuous_scores.yaml        # Continuous score color ranges
    pipeline_steps.yaml           # ghfc-ngs step labels/order + run status chips
    vep_consequences.yaml         # VEP consequence terms, impacts, colors
    view_presets.yaml             # Column visibility presets
  pages/
    login.py                      # Login page (no auth required)
    profile.py                    # User profile (role, password, gene badge tiers)
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
        roh_tab.py                # Runs of homozygosity tab (per-family)
        str_tab.py                # STR tab (per-family): EH sub-tab + locus dialog
        cohort_roh_tab.py         # Runs of homozygosity sub-tab + islands dialog (cohort-wide)
        cohort_str_tab.py         # STR sub-tab (cohort-wide): loci counts + flagged genotypes
        sample_qc_tab.py          # Sample QC tab (per-family): issues, kinship, sex, quality
        cohort_sample_qc_tab.py   # Sample QC sub-tab (cohort-wide) + cross-family network dialog
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
    roh.py                        # ROH discovery, reference tier, grades, kinship checks
    roh_plots.py                  # ECharts/SVG builders for the ROH tabs
    str_eh.py                     # ExpansionHunter discovery, reading, cohort cache, locus pivot
    str_plots.py                  # ECharts builders for the STR views
    str_validation.py             # validations/strs.tsv writer (role-guarded) and reader
    sample_qc.py                  # sampleqc (somalier) discovery, reading, flag wording, shaping
    sample_qc_acks.py             # validations/sampleqc.tsv: curator flag acknowledgements (role-guarded)
    sample_qc_plots.py            # ECharts builders for the Sample QC views
    validation_files.py           # The validation kinds (SNV/SV/STR) and their files
    diagnostic_files.py           # The diagnostic kinds (SNV/SV/STR) and their files
    ideogram.py                   # Chromosome backbone SVG (bands, Mb grid)
    cnv_coverage.py               # Coverage bigWigs, CNV control panel, fit table, copy number
    cnv_plots.py                  # ECharts builder for the copy-number panel
    sv_view.py                    # One read window or one per breakpoint, by SV type and size
    sv_reads.py                   # pysam: read pairs/split reads classed by orientation
    sv_breakpoints.py             # Suggested breakpoints: clips/splits, pair interval, CN check
    sv_read_plots.py              # ECharts builder for the read-evidence panel
    auth.py                       # Auth helpers: check_auth, can_write, get_current_user
    clinvar.py                    # ClinVar color/display utilities (from YAML)
    cohort_state.py               # .ghfc-ngs.state.json parsing, drift, step rows
    cytobands.py                  # Cytoband, chromosome, ideogram constants
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
    sharding.py                   # Two-level sharded directory path resolution
    tsv.py                        # read_tsv_or_none (string_columns), clean_tsv_field, find_labels, append_tsv_row
    wombat.py                     # aggregate_by_variant, shared by the Wombat tab and search
    sample_counts.py              # Cohort Wombat/WisecondorX files, per-sample counts, outlier flags
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

### Runs of Homozygosity (ROH)
- The ghfc-ngs `roh` step writes `<family>/roh/<FID>.apgs_b<bundle>_<tag>.{froh,roh,roh_genes,kinship}.tsv` and, per cohort, `cohorts/<name>/roh/<name>.…{froh,roh,roh_genes,kinship,roh_windows}.tsv` (plus `family_id`). Same naming as the ancestry files, so discovery is `find_ancestry_files_in(..., kinds=ROH_KINDS)`; `roh_genes` is optional (needs `annotation_gencode`), `kinship` is absent for a singleton.
- **Read through `read_roh_table`.** It forces the ID columns to strings. Left to inference, `family_id` comes out as an integer from numeric first rows and the read then fails on `AU-ANG-FER-MPX3`. It also strips the `chr` prefix, because the reference hotspots and the cytobands carry none.
- **Reference first, fixed thresholds second.** Each person is compared with unrelated reference genomes of their `ref_group` (a region, or `ALL` = pooled), read from `v<bundle>/bundle/roh/reference_froh.tsv.gz` under `ancestry_reference_dir`, next to `hotspots.tsv.gz`. The fixed 0.03 line is crossed by 150 of 682 unrelated CSA reference genomes, so a threshold-only display would mislabel whole populations. A missing tier is a degraded state, as for the PCA panel.
- **Thresholds and wording live only in `config/roh.yaml`.** They mirror ancestry-pgs constants (`roh.py`), which assign `froh_category` and `homozygosity_pattern` themselves; the app only words and draws them.
- **Wording is role-dependent (`can_write()`).** Viewers get neutral grades and only pedigree-declared kinship pairs. Curators also get the relationship reading, the parents' kinship / expected F, and every pair. This is presentation only: the data directory is served statically without an auth check.
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
- **Classes and statuses are the pipeline's.** The app only colours and words them, from `config/str.yaml`. `str_plots.allele_class` mirrors the pipeline's order (pathogenic, intermediate, normal, above_normal) only to colour histogram bars.
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
- **The family files never hold cross-family findings.** The cohort tables are the family ones plus the cross-family pairs *observed related* (unrelated cross-family pairs are not written) and their `level: cohort` flags (`CRYPTIC_RELATED`, `DUPLICATE_OR_MZ`, `SAMPLE_SWAP_CANDIDATE`). The family tab reads the cohort tables off-thread for them; `family_flags_with_cohort` adds only the `cohort`-level rows, since the cohort file repeats the family flags.
- **A cohort pair flag's `family_id` is `"FA,FB"`**: `sample_a` is in FA, `sample_b` in FB. Key samples by `FlagRow.keys` (`(family_id, sample)`), never by sample alone.
- **The pipeline classifies everything** (relationships, karyotypes, flags, severities); the app only words and draws them from `config/sample_qc.yaml`. `KARYOTYPE_<k>` is matched by prefix; an unknown code is shown raw. The plot lines mirror the `sampleqc_classify` constants, and the KING cut-offs come from `roh.yaml`'s `kinship_degrees`.
- **Wording is role-dependent, as for ROH**: viewers get `reader_label`, no evidence, an "Atypical" karyotype, and only pedigree-declared pairs (no parents' pair, no cross-family pairs in the tables). Presentation only.
- **Hubs**: a sample with `HUB_MINIMUM` (5) or more cross-family relatives is outlined in the network. One EAGER founder is "third-degree" to 37 samples of other families; that is contamination or population structure, not kinship.
- Elsewhere: the family header badge and member icons (family files only, read at page build), the cohort table's `Sample QC` column and filter, the home card's QC figure (off-thread, individuals by worst severity), and a warning above Wombat configs matching `de_?novo` when `denovo_risk` flags hit the family.
- **Curators acknowledge recurring flags** (MZ twins, declared consanguinity, adoption, tolerated contamination) in `validations/sampleqc.tsv`, through `utils/sample_qc_acks.py`:
  - Keyed by `FlagRow.ack_key`: code, the `family_id` cell (`"FA,FB"` for a cross-family pair), `sample_a`, `sample_b` (empty for one sample). No sites label, so a rerun keeps it. The family and cohort files write a family flag alike, so one key covers both pages.
  - Append-only (`tsv.append_tsv_row`, under `flock`): status `explained` or `confirmed`, a note, user, timestamp. `withdrawn` takes it back; the last row of a key wins. Not a variant validation, so not in `VALIDATION_KINDS`. `save_flag_ack` / `withdraw_flag_ack` check `can_write()` themselves.
  - Acknowledged flags leave the family header badge, the member icons (grey `task_alt` when nothing else is left), the cohort column (one grey "Acknowledged" badge and filter option, not "Any issue"), the cohort panel counts, the home figure, and the de novo warning. The warning keeps a `confirmed` flag: the calls stay unreliable.
  - The family Issues table lists them greyed and last, with the review and note ("Show acknowledged" hides them). The cohort panel leaves them out unless "Include acknowledged" is on.
  - A review redraws the header badge and member icons through `render_sample_qc_tab(on_acks_changed=...)`. The dialog hangs from an anchor outside the refreshable section, as the STR locus dialog does.
  - **Readers see the status, and the note only for categories with `reader_note: true`** in `sample_qc.yaml` (sex, mendelian, sample quality). Elsewhere a note such as "known cousins" would name what `reader_label` hides.

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
- **Right-hand panel, inside `General`**: `Statistics` / `Ancestry` / `Polygenic Scores` / `ROH` / `STR`, nested inside the General tab panel. Anything describing the pipeline run belongs on `Workflow` instead: it is a property of the cohort, not of the current filter, and does not fit a 400px column.
- The right-hand panel also has `Sample QC`.
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
- **MOI flag** (`utils/inheritance.py`, `config/inheritance.yaml`): a gene's modes are pooled from every GeneTrek field in category `Inheritance`, because the combined `inheritance` holds one mode. Recessive second hits are looked for only among the rows of the same table (`annotate_moi`), and the tooltip says so. A male's chrX call outside PAR is hemizygous whatever the GT says. The DNM tab is not wired: `render_dnm_tab` has no caller.

### Variants Per Sample (family Stats dialogs)
The Wombat and WisecondorX Stats dialogs chart each selected member's count in the current table view, then place each member on the cohort's per-sample distribution (`utils/sample_counts.py`, drawn by `components/sample_counts.py`).
- **The cohort comparison is unfiltered on both sides.** It counts the cohort-level file (`cohorts/<C>/wombat/<C>.rare.<vep>[.annotated].<config>.results.tsv`, or `<C>_aberrations[.annotated].bed`), never the table's filters. Older cohorts have only the `.annotated.` name, so the lookup tries both.
- **Several cohort de novo tables have no header** (CANDY_mpx, CANDY_mpx_ghent, PMS): a blank line, then rows. Read them through `read_cohort_wombat`: a family file of the same config lends its header when the widths match; otherwise `MissingHeaderError` is shown, never a cohort of zeros or a first row taken for column names. The Stats dialog borrows from the family's own file; the Search page tries the cohort's families in turn (`family_header_sources`). The default `datagen` data has one, DEMO's `strict` table (`build_headerless_cohort_wombat`).
- **Search lists cohort Wombat tables through `cohort_wombat_sources`**, sorted, so its default source does not change when a file is added. The `annotated` part of the name is not part of the config.
- **Founders and non-founders are compared within their own group** (`Sample.is_founder`), with every pedigree sample zero-filled. A de novo config has rows only for children, and pooling them with the parents' zeros would flag any child with two calls.
- An outlier is outside Q1-Q3 by more than 1.5 × IQR (`OUTLIER_IQR_FACTOR`) and on the far side of the median.
- Counts load off-thread and are cached on `(path, mtime, size)`. The comparison sits outside the dialog's refreshable section, so the SNV/Indel toggle does not recount.

### Pedigree Missing Values
The sentinel set `{"", "0", "-9"}` represents unknown/missing in pedigree fields (parent IDs, sex, phenotype). Defined as `_PED_MISSING` in `search.py` and handled in `models.py` via `treat_missing_as_null`.

### TanStack Table
- `DataTable` in `tanstack_table.py` bridges Python to vanilla JS TanStack Table.
- Per-client JS/CSS injection via `client._dt_scripts_injected` attribute (avoids global set memory leak).
- Tables are created when their container enters the document, with no timeout (`static/js/data_table_mount.js`, one shared MutationObserver). Quasar mounts only the shown tab panel, and a left keep-alive panel sits detached, so a table built for a hidden tab gets its container whenever the tab is opened. `update_data` / `set_column_visibility` calls made before that are queued.
- A tab bar whose panels each hold a table builds each panel on first view (the Wombat sub-tabs, as the family page's tabs).

### Async
- Use `asyncio.to_thread()` to offload blocking I/O (e.g., `pl.read_csv()`) in NiceGUI page handlers.
- NiceGUI runs a single-threaded event loop; blocking calls freeze the UI for all users.

## Conventions

- **Validation kinds are listed once**, in `utils/validation_files.py`. The All validations and Statistics pages loop over it, so a new kind is one entry there plus its writer.
- **Diagnostic kinds are listed once**, in `utils/diagnostic_files.py`. The following loop over it: the Diagnostic All and Statistics pages, the cohort and home summaries, and `load_family_diagnostics`. The per-kind table callers (wombat, svs, search, validation/file) name their own file. `save_diagnostic_entry` and `update_diagnostic_ignore_status` check `can_write()` themselves.
- Utility code shared across pages lives in `utils/`, not duplicated in page components.
- YAML configs live in `config/`, loaded by corresponding `utils/*.py` modules.
- All new ClinVar terms must be added to `config/clinvar_colors.yaml`.
- All new VEP consequences must be added to `config/vep_consequences.yaml`.
- Column display/group/sort/drop config goes in `config/column_names.yaml`.
