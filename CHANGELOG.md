# Changelog

All notable changes to genetics-viz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.12.0] - 2026-09-22

### Added
- **A cohort is now any directory under `cohorts/`.** Directories whose `<NAME>.pedigree.tsv` is missing — or present but unparseable — are listed on the home page as *Incomplete Cohorts*: amber, non-clickable cards naming the file they need, instead of being silently skipped. Six such directories existed on the GHFC data directory and were invisible. `DataStore.cohorts` still holds only explorable cohorts, so the header dropdown, validation pages and search are unchanged; the rest arrive in `DataStore.incomplete_cohorts`, which distinguishes a missing pedigree from an unreadable one because the two need different fixes.
- **Full-width `General` / `Workflow` tabs on the cohort page.** `General` is the individuals table and the right-hand panels that follow its filters, unchanged. `Workflow` gathers everything about the pipeline run and exists only for a cohort with a `<NAME>.params.yml`.
- **Pipeline status, read from the workflow's own run record.** The ghfc-ngs workflow writes `.ghfc-ngs.state.json` into each cohort directory with per-step `{done, total, pct}` counts (schema: `COHORT_STATE.md` in bourgeron-lab/ghfc-ngs), and genetics-viz never recomputes them — re-deriving would duplicate the pipeline's planning logic, and a cold existence check costs ~5.6 ms over the CIFS mount, so a 2700-sample cohort would be a minute of stat calls. Home page cards gain a foldable **Status** item: a status chip and an "N/M steps complete" summary collapsed, per-step progress bars when opened. The Workflow tab shows the full record — run timings and duration, last successful run, pipeline version and commit, run history — with warnings for a stale `running` record (there is no heartbeat, so a killed run leaves one behind), counts measured before the run, and outputs a failed run may not have published. A step the record reports as `null` renders "not measured", never 0%: unmeasured and zero are different facts.
- **Drift detection.** The Workflow tab re-hashes the local pedigree and params file against the checksums in the record, flagging "the pedigree has been edited since this run". The paths inside the record are compute-cluster paths, so the local files are hashed instead.
- **The cohort's `params.yml` beside the pipeline's own reference for it**, two side-by-side panels under the status, each scrolling internally so the page never grows a long scrollbar. The parameters file is shown verbatim — these files carry a great deal of their meaning in comments. The reference is `documentation/params.md`, fetched at runtime from the public ghfc-ngs repository (stdlib `urllib`, 10s timeout, cached an hour — no new dependency) rather than vendored, since it changes often. Read outside its own repository, its repo-relative links are rewritten to absolute GitHub URLs opening in a new tab, with `/tree/main/` for the one directory link; its 73 in-page anchors resolve inside the panel via markdown2's `header-ids`, with clicks intercepted so the page neither jumps nor collects a URL fragment. A failed fetch degrades to a note naming the error plus a link to the original, leaving the rest of the tab intact.
- **`config/pipeline_steps.yaml`** — step labels, display order and run-status chip colours. Display metadata only: completion is measured by the pipeline, never recomputed here. A step the pipeline reports but this file does not know is shown under its raw name.

### Fixed
- **An unknown cohort name in the URL returned a 500.** `create_header` passed the name straight into a `ui.select` value; Quasar raises when the value is not among the options, which fired before the pages' own "cohort not found" branches could run, making them dead code. Those branches now work, and a pedigree-less cohort directory explains what it is missing instead of erroring.
- **The cohort cards lost most of their width to layout defects**, which the new Status panel made visible: `nicegui-card`'s own 16px padding doubled with the 16px on each `q-card__section` (32px a side), and because `nicegui-card` is a flex column with `align-items: flex-start` the sections did not stretch, coming out 274px wide in a 320px card — so nothing inside lined up with the card's edges. Sections now carry `w-full`, the card carries `p-0`, and padding is `pl-2 pr-3` to offset the 4px accent border. Usable content went from 250px to 296px with symmetric insets, and no element overflows at any container width from 360px up.
- **Two flex-sizing overflows on the new panels.** Long file paths escaped their card because `.nicegui-column` sets `align-items: flex-start`, so a child is sized by its content and `truncate` has no width to act on; those labels now carry `w-full`. Separately, the step rows' `whitespace-nowrap` counts gave the home card's Status expansion a 299px min-content width — a flex item will not shrink below its content, and that floor propagates up the ancestor chain — which forced the card section to 331px inside a fixed 320px card and clipped the phenotype table's last column. Every flex item in the chain now carries `min-w-0` and truncates, and the card's compact rendering drops the unit noun (`100% · 522/522` rather than `100% · 522/522 samples`); the Workflow tab, which has the room, keeps it.

### Changed
- Cohort discovery was hand-copied into `DataStore.load()`, `take_snapshot()` and `reload()`; all three now share one `_iter_cohort_dirs()` generator. `take_snapshot()` additionally watches `<NAME>.params.yml` and keys every cohort directory rather than only pedigree-bearing ones, so a directory appearing or gaining a pedigree is detected. The run record is deliberately *not* watched: the pipeline rewrites it on every run, and a reload re-parses every pedigree in the data directory — its reader caches on the file's own stat instead.

## [0.11.0] - 2026-09-21

### Added
- **Ancestry and Polygenic Scores tabs on the cohort page** — the right-hand Statistics panel becomes three tabs. Both new tabs read the cohort-level pipeline output and follow the filters applied to the individuals table.
- **Cohort Ancestry tab** — every cohort member on PC1×PC2 over the reference panel, which a checkbox (ticked by default) can hide. Points are coloured by predicted region and outlined red when flagged as outliers; the tooltip names the sample, region, population and confidence. Beside the plot, the region distribution as a pie plus a count/percent table. For EAGER that is 586 samples across 7 regions.
- **Cohort Polygenic Scores tab** — violins of the cohort's z-scores for a selected trait. The grouping axis is Phenotype, Sex or predicted Region; an optional split draws each group as two sex halves; and a region multiselect restricts which samples are included. Group sizes are printed on the axis and the individual points are overlaid, because filtering by region can leave groups of a handful of samples.
- **Bundle/filter-tag selector** on the cohort and family ancestry tabs, shown when a directory holds more than one variant.

### Fixed
- **The strictest filter-tag variant is now preferred.** Ancestry directories now hold two variants of the same bundle — `dp6gq15` and `dp10gq20` — and the tags were sorted as plain strings, so `dp6gq15` won because "6" sorts above "1" at the third character. The family tabs released in 0.10.0 therefore switched variant silently as soon as the second one appeared. Tags are now ranked by their depth and GQ thresholds.
- Pedigree sex and phenotype codes written float-stringified (`-9.0`, `1.0`, as the PMS cohort does) were not recognised as their integer equivalents, so they grouped as unknown.

### Changed
- The shared PCA plot builders moved from the family Ancestry tab into `utils/ancestry_plots.py`, and the phenotype label helpers from `home.py` into `utils/pedigree_labels.py`, so the cohort tabs reuse them rather than importing another page's private functions.
- The cohort violins are Plotly, the second deliberate exception to the project's ECharts default — ECharts has no violin series, and the sibling Statistics tab is already Plotly.

## [0.10.0] - 2026-09-17

### Added
- **Ancestry tab on the family page** — the family's samples projected onto the versioned `ancestry-pgs` reference panel. Four scatter plots: PC1×PC2 and PC3×PC4 over the full panel, then the same pair zoomed to a window covering the family plus every reference sample of the regions predicted for it (~13× magnification for a single-region family). The reference cloud is one series per region, so the legend toggles regions, and it is `silent` so it never steals a tooltip from a family point. Family samples are labelled diamonds filled with their predicted region colour, red-bordered when flagged as outliers. Above the plots, a table of predicted region, population, both confidences, kNN distances, local reference density and the outlier flag.
- **Polygenic Scores tab** — the per-family `pgs_zscore.tsv`, which the pipeline writes wide (one row per sample, ~100 score columns), transposed to one row per trait and one column per member, with a diverging colour scale on the z-scores and a text filter on the trait name.
- **`ancestry_reference_dir` config key** — root of the references tree. The bundle version subdirectory is derived per family from its own filename tag (`apgs_b1.0.0` → `v1.0.0`), so each family is plotted against the panel it was actually projected onto rather than a globally pinned one. The reference panel is cached per bundle and read once, not per page load.
- **`config/ancestry.yaml`** — colourblind-safe reference region colours and the PGS z-score colour scale.
- Both tabs carry existence probes and are disabled for families with no `ancestry/` directory, which is currently the majority. They are wired into the cohort-scoped and the standalone family page.

### Fixed
- **`save_config` no longer drops unknown top-level keys.** It rewrites the whole YAML from known fields, so `poll_interval` was already being silently deleted the next time an administrator added a user or a data directory; `ancestry_reference_dir` would have inherited the same fate. Both are now emitted when they carry information.
- A missing or unconfigured reference panel is treated as a degraded state rather than a failure: the Ancestry tab still plots the family's own samples and names the missing config key, instead of erroring.

### Changed
- CLAUDE.md structure listing corrected — the continuous-score config file is `config/continuous_scores.yaml`, not `config/score_colors.yaml`.

## [0.9.2] - 2026-09-16

### Added
- **Tabbed family information panel** — Notes and Diagnostics moved into a **General** tab, joined by **Ancestry** and **Polygenic Scores** placeholders. General opens by default, and a half-typed note survives a visit to another tab.
- **`/health` endpoint** — unauthenticated liveness probe reporting the running version from `importlib.metadata`, so a deployment can assert that the release it asked for is the one answering. Left unauthenticated on purpose so a pool monitor can use it instead of a bare TCP connect.

### Changed
- **Pedigree table is much shorter** — row pitch dropped from 73px to 29px and the table from 640px to 260px for an eight-member family, by removing the inter-row flex gap and using a dense checkbox.
- The member table became a shared `render_member_selector` component, so the standalone family page — which carried a near-verbatim copy with the same bugs — is fixed by the same change.

### Fixed
- **Pedigree table alignment** — the grid style was applied to the header row and each member row separately, making every row an independent grid sized to its own content. Six of seven columns drifted by up to 75px. A single grid container now holds every cell.
- **Diagnostics merged across curators** — a variant recorded by two curators appeared twice; it is now one row, with a tooltip listing each curator's own verdict and date.
- **Diagnostics merged across individuals** — a variant reaching the same verdict in several family members is one row listing every barcode in pedigree order. Members that reached *different* verdicts stay on separate rows.
- **Empty wombat files no longer log as failures** — wombat leaves a zero-byte file when a step had nothing to report (a family without two sequenced parents has no de novo mutations) and polars raises `NoDataError` on it, which was being caught by the handler meant for real failures and logging a traceback on every visit. A new `read_tsv_or_none` helper returns `None` for such a file, leaving `FileNotFoundError` to propagate.

## [0.9.1] - 2026-08-25

### Fixed
- **Pinned TanStack table-core to 8.21.3 on the CDN** — the jsDelivr URL was unpinned, so it resolved to whatever npm tagged `latest`. TanStack Table v9.1.2 became `latest` and dropped the `getCoreRowModel` export, so every table threw `TanStack.getCoreRowModel is not a function` at load time and never rendered.

## [0.9.0] - 2026-05-20

### Added
- **Top-level tab disabling** — the Wombat and SVs tabs on the family page are greyed out and unclickable when no data files exist for a family.
- **Sub-tab disabling** — within Wombat, sub-tabs whose DataFrame is empty are disabled rather than showing blank content, distinguishing "step ran but found nothing" from "step was not run".
- **Variant count badges** — enabled sub-tabs show a count, e.g. `config_name (42)`, so the volume is visible without opening the tab.
- **Smart default selection** — the first tab with actual data is auto-selected, skipping disabled ones.
- Lightweight `probe_wombat_data()` / `probe_svs_data()` helpers check for data existence without a full DataFrame load, so tab state can be decided before rendering.

## [0.8.1] - 2026-05-19

### Fixed
- **Deferred config save until the password dialog is closed** — in admin user management, `save_config()` ran before showing the generated password, so the app could reload before the administrator had copied it. The YAML write now happens only when the dialog is dismissed.

## [0.8.0] - 2026-05-19

### Added
- **Live data refresh** — filesystem changes are detected automatically; no restart when cohorts are added, removed or updated on disk.
- **Background polling** — one server-wide task checks all data directories every 30 seconds for new, removed or modified cohort pedigree files, using `asyncio.to_thread()` for the stat calls so the event loop stays responsive.
- **Toast notifications** — connected browsers are told what changed.
- **Manual refresh button** in the header for an immediate check.
- **Atomic reload** — `DataStore` rebuilds its cohort dict in a local variable before swapping, so concurrent page loads always see a consistent state.
- **`poll_interval: <seconds>`** config key to adjust the polling interval (default 30).

## [0.7.9] - 2026-04-08

### Added
- **Help menu in header bar** — small `help_outline` icon button (right side, before the user menu) with three options: Email (mailto to the Teams channel address with "Bug report" subject), Teams channel (opens the genetics-viz channel in a new tab), and Report bug (opens a new GitHub issue in a new tab).

## [0.7.8] - 2026-04-08

### Fixed
- **Family page member table** — rebuilt the member selection table using native NiceGUI CSS-grid layout instead of a static HTML table + JS DOM moves. The previous approach was fragile because Vue/NiceGUI reconciliation would revert the manual moves. Checkboxes and "only" buttons now render reliably inside their cells. Same fix applied to standalone family page.

## [0.7.7] - 2026-04-08

### Fixed
- **Family page member table checkboxes** — fixed timing bug where checkboxes and "only" buttons sometimes rendered below the table instead of inside the Select column. The DOM-move JavaScript now polls for both source elements and target cells with bounded retries.
- **DataTable container race condition** — replaced fragile MutationObserver with a bounded polling loop in `tanstack_table.py`. Fixes intermittent issue where the Wombat/SVs dataframe wouldn't render and required a page reload.

## [0.7.6] - 2026-04-07

### Changed
- **Cohort card phenotype table styling** — shortened phenotype labels (`2 (aff)`, `1 (unaff)`, `-9 (unk)`) with `whitespace-nowrap` to prevent wrapping; wrapped the table in a rounded border with a separator line between the header row and data rows.

## [0.7.5] - 2026-04-07

### Changed
- **Cohort card phenotype breakdown** — replaced the single Affected/Diagnosed columns with a per-phenotype table showing N samples, Pat/Unc counts, and diagnostic yield % for each phenotype category (2=affected, 1=unaffected, -9=unknown). Only categories with samples are displayed.

## [0.7.4] - 2026-04-07

### Changed
- **Cohort card diagnosed count** — now shows two numbers: pathogenic (red) / uncertain-only (amber), where uncertain-only counts samples with at least one uncertain diagnostic but no pathogenic. Percentage still based on pathogenic / affected.

## [0.7.3] - 2026-04-07

### Added
- **Family navigation arrows** — prev/next arrow buttons on the cohort family page header to navigate between families in the current project. Hidden for single-family cohorts; previous hidden on first family, next hidden on last.

## [0.7.2] - 2026-04-07

### Added
- **Family Notes panel** — new panel above Diagnostics on both cohort and standalone family pages for adding free-text notes to families or specific samples. Notes stored in `notes/notes.tsv` with hard delete support. Includes add form (with optional sample selector), delete button, info tooltip, and write protection.

## [0.7.1] - 2026-04-03

### Added
- **Version number in header** — displays version below the app name in the header bar
- **Cohort card statistics** — home page cohort cards now show affected count, pathogenic diagnostic count, and diagnostic yield percentage alongside families and samples
- **Diagnostic column in cohort table** — new Diagnostic column with colored badges (pathogenic/uncertain/benign) and multiselect filter (including NA for undiagnosed) on the cohort detail page

## [0.7.0] - 2026-04-03

### Added
- **Quick Search on home page** — search samples by barcode or families by FID with async O(1) sharded filesystem lookup; results show data availability badges and action buttons
- **Sample visualization dialog** — fullscreen IGV.js viewer with bedgraph (CNV), CRAM (alignments), and VAF tracks for any sample, with locus navigation input
- **Standalone family page** (`/family/{fid}`) — browse family data independently of any cohort/project, with data availability panel, member table, diagnostics, and analysis tabs (Wombat, SVs)
- **Standalone pedigree parser** (`utils/pedigree.py`) — parses per-family pedigree files with extended column name support (FatherBarcode, MotherBarcode, Pheno_*)
- **Data availability checker** (`utils/data_availability.py`) — checks existence of data files for samples and families

### Changed
- `show_sv_dialog()` and `show_variant_dialog()` accept optional `family_members_override` and `sample_parents_override` parameters for standalone pages without cohort context (backward compatible)
- **Sharding detection** — replaced "all children must be single-char" heuristic with try-sharded-first strategy to support hybrid directories with both shard buckets and direct entity folders

## [0.6.9] - 2026-03-24

### Added
- **CRAM split view re-centering** — when ROI is updated (via suggestion click, curated position refresh, or New Start/End buttons), the Read-Level Split View automatically re-centers both panes on the new boundary positions (±1500bp windows)

## [0.6.8] - 2026-03-24

### Fixed
- **IGV.js ROI display** — reverted from v3.8.0 to v2.15.13 (latest 2.x) due to ROI rendering incompatibility in v3. ROI overlays now display correctly again.

## [0.6.7] - 2026-03-24

### Fixed
- **IGV.js v3 ROI display** — ROI regions were not rendering because `loadROI()` became async in v3. Moved ROI loading from `createBrowser` config to explicit `await browser.loadROI()` calls after browser creation. Also wrapped dynamic `_update_roi()` in an async IIFE.

## [0.6.6] - 2026-03-24

### Added
- **SV suggestion inheritance preselection** — clicking a coordinate suggestion from a parent now auto-selects "paternal"/"maternal" in the inheritance dropdown
- **Dual overlap percentages** — SV suggestions now show both "ours" (overlap as % of current SV) and "theirs" (overlap as % of suggested SV) for at-a-glance size comparison

### Changed
- **IGV.js upgraded from v2.15.11 to v3.8.0** across all 6 page files (7 script tags)

## [0.6.5] - 2026-03-23

### Fixed
- **SV ideogram colors in stats dialog** — SVs from to_validate files with non-standard `call` column values (e.g. "primary") now correctly show GAIN/LOSS colors on the ideogram. Uses shared `infer_sv_type()` for consistent SV type detection across all code paths.

## [0.6.4] - 2026-03-23

### Added
- **SV coordinate suggestions** — panel between SV Details and CNV Coverage View suggesting curated coordinates from overlapping validated SVs; prioritizes parents, then family, then cohort; click to auto-fill curated position fields and update ROI

## [0.6.3] - 2026-03-23

### Fixed
- **SV validation display for GAIN variants** — validations on "gain" type SVs now display correctly in the validation file table. Root cause: inconsistent SV type inference between save (sv_dialog.py, checked only `call` column) and lookup (file.py, checked `wisecondorX` first). Extracted shared `infer_sv_type()` utility in `utils/wisecondorx.py` used by all callers.

## [0.6.2] - 2026-03-16

### Added
- **VAF View panel** in SV dialog — displays VAF bedgraph tracks (`*.vaf.bedgraph.gz`) as a collapsible panel between CNV Coverage View and Read-Level Split View, rendered as scatter points with fixed 0–1 Y-axis scale

## [0.6.1] - 2026-03-16

### Added
- **Sharded directory support** — two-level sharding for `samples/` and `families/` directories (`utils/sharding.py`), with auto-detection and backward compatibility for flat layouts

### Changed
- All sample/family path constructions now use `get_sample_path()`/`get_family_path()` and `get_sample_url()`/`get_family_url()` from the sharding module

## [0.6.0] - 2026-03-09

### Added
- **Authentication & authorization** — YAML-configured user accounts with SHA-512 password hashing and role-based access control (reader, curator, administrator)
- **Login page** — username/password login with session persistence
- **Profile page** — view current role and change password
- **Admin: user management** — add, remove, change role, and reset password from the web interface
- **Admin: data directory management** — add, remove, and set default data directory from the web interface
- **Multi-data-directory support** — per-user data directory selection via header dropdown; YAML config lists multiple directories with descriptions and default flag
- **YAML config file** — single config file replaces CLI data_dir argument; holds data directories, user list, and auto-generated storage secret
- **Write protection** — save operations (validations, diagnostics) gated by `can_write()` at both UI and backend level
- **Config model** (`config_model.py`) — dataclasses for config loading/saving with file locking and password helpers

### Changed
- **CLI** — accepts `config_file` (path to YAML) instead of `data_dir`
- **App init** — loads config, initializes multi-store registry, registers per-directory static files, sets storage secret from config
- **Header** — data directory dropdown replaces reload button; user menu with profile/logout added
- **Static file URLs** — migrated from hardcoded `/data/` to `get_static_prefix()` for per-directory routing
- **All page handlers** — auth guard added (`check_auth()` redirect at top of every page)
- **Validation/diagnostic saves** — username read from session instead of OS user; `can_write()` guards added

## [0.5.0] - 2026-02-24

### Added
- **Exclude samples filter** — search Individuals panel: enter sample IDs (space/comma/blur) as deletable chips to exclude from results
- **SV support in statistics dialog** — SV deduplication, gain/loss classification in consequence chart, SV checkbox filter, variant type breakdown in subtitle
- **SV validation file support** — auto-detects SV format in validation files, loads `svs.tsv`, opens SV dialog instead of IGV
- **SV type inference** — extracts SV type from wisecondorX call, generic call, type column, or ratio sign
- **WisecondorX tooltips** — CNV call cells show ratio/zscore on hover in search results
- **Numeric column formatting** — integer columns display with thin-space (U+202F) thousands separator and proper Unicode minus sign
- **Column type config** — `type: int/float` in `column_names.yaml` auto-infers numerical sorting and number cell rendering

### Changed
- Sex and Phenotype dropdowns now on same row in search Individuals panel
- Validation file table preserves sorting and pagination across refreshes
- Statistics dialog subtitle shows variant type breakdown (SNVs, Indels, SVs)

### Fixed
- `KeyError: 'error'` on ratio Min/Max input fields — `props(remove="error")` was deleting the key from internal props dict

## [0.4.1] - 2026-02-18

### Added
- **Inheritance inference** — variant dialog auto-infers inheritance mode from VCF genotypes
- **Search TSV export** — download search results as TSV or save as validation file
- **Validation file stats** — stats dialog (chromosome distribution, consequence/status charts, ideogram) now available on validation file pages
- **Validation file badges** — gene, consequence, ClinVar, and score badges in validation file tables
- **Validation file column selector** — column visibility dialog with presets for validation pages
- **CRAM version check script** — `utils/check_cram_version.sh` checks IGV.js compatibility via magic bytes
- **FID column** added to all view presets

### Refactored
- Extracted `utils/locus.py` — centralized locus query parsing and filtering (from search.py)
- Extracted `utils/genesets.py` — shared geneset loading (was duplicated in 3 places)
- Extracted `utils/wisecondorx.py` — WisecondorX parsing, CNV classification, and color utilities (from svs_tab.py)
- Extracted `components/search_stats.py` — reusable stats dialog (from search.py)

### Changed
- Validation all/statistics pages now async with offloaded I/O

## [0.3.0] - 2026-02-16

### Added
- **TanStack Table** - New high-performance DataTable component replacing NiceGUI ag-grid, with virtual scrolling, column sorting, and custom cell renderers
- **Cohort-wide Search** - Tabbed search parameters panel with Variants and Individuals tabs
  - Individual filters: sex, phenotype, and "only samples with both parents"
  - Filters use pedigree data from the Cohort object (no duplicate parsing)
- **Variant Statistics dialog** - Stats button on search results showing:
  - Stacked bar chart of variants per chromosome by validation status
  - Consequence distribution pie chart (highest priority per variant)
  - Validation status distribution pie chart
  - Interactive SVG ideogram with cytoband rendering and variant positions
  - SNV/Indel type filter checkboxes with live refresh
- **ClinVar `Likely_risk_allele`** term added to clinvar_colors.yaml
- **Column names config** (`column_names.yaml`) for display labels, groups, sorting, width, and drop flags
- **View presets** for quick column visibility switching
- **Gene scoring** utilities and score color coding
- **Stats panel** component for carrier frequency box/bar plots
- **CLAUDE.md** project-level instructions for Claude Code

### Changed
- **Pedigree parsing** now handles `-9` as missing/unknown (in addition to `0` and empty)
- **Pedigree header detection** strips leading `#` (supports `#FID` headers)
- **Search panel** uses dense Quasar props for a more compact layout
- Pedigree data loaded from the already-parsed `Cohort` object via `DataStore` instead of re-reading the file

### Refactored
- Extracted `utils/vep.py` - consolidated VEP consequence utilities from 3 files
- Extracted `utils/clinvar.py` - consolidated ClinVar utilities from 2 files
- Extracted `utils/view_presets.py` - view preset logic from wombat_tab
- Extracted `utils/cytobands.py` - cytoband/ideogram constants shared between wombat_tab and search
- Fixed `header.py` coupling - no longer imports from page components

### Fixed
- **TanStack table memory leak** - replaced global `_injected_pages` set with per-client attribute
- **MutationObserver timeout** - 30-second timeout prevents indefinite observation
- `Sample.is_founder` now correctly handles `-9` parent IDs

### Removed
- Dead code: `variant.py` (unused variant page)
- Dead code: `waves_backup.py` (obsolete backup)

## [0.2.0] - 2026-01-20

### Changed
- **Breaking Change**: New validation TSV format with `Comment` and `Ignore` columns
  - Use migration script: `utils/snvs_validations_migration_0.1.1_to_0.2.0.sh`
- Validation form now defaults to "present" status
- Validation/all page aggregates by Variant/Sample with unique user lists
- "in phase MNV" normalized to "present" for conflict detection but displayed distinctly

### Added
- Enhanced inheritance options: "not maternal", "not paternal", "homozygous"
- New validation status: "in phase MNV" for multi-nucleotide variants
- Comment field for free-text notes on validations
- Ignore functionality to exclude specific validations from statistics and conflict detection
- Interactive validation guide accessible via info button in dialog
- Validation history shows ignore toggle for each entry
- Ignored validations displayed with reduced opacity
- Statistics page excludes ignored validations and shows separate ignored count

### Fixed
- Fixed context issues when toggling ignore status from within dialog
- Fixed table refresh after validation changes

## [0.1.1] - 2026-01-15

### Added
- Initial stable release
- Core validation functionality
- IGV.js integration
- WAVES validation support
- Multi-cohort management
- Family structure visualization
- DNM and WOMBAT analysis tables
- Variant validation tracking with inheritance patterns
