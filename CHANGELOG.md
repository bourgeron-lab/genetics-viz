# Changelog

All notable changes to genetics-viz will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
