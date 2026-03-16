# Changelog

All notable changes to genetics-viz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
