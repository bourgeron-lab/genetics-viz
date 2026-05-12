# Changelog

All notable changes to genetics-viz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
