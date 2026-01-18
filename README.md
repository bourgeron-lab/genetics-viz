# Genetics-Viz 🧬

A web-based visualization tool for genetics cohort data.

## Features

- Browse multiple cohorts from a single data directory
- View family structures and pedigree information
- Interactive data tables with sorting and filtering
- Clean, modern web interface built with NiceGUI

## Installation

### From PyPI (when published)

```bash
pip install genetics-viz
```

### From source

```bash
# Clone the repository
git clone https://github.com/fcliquet/genetics-viz.git
cd genetics-viz

# Install with uv
uv sync

# Or install with pip
pip install -e .
```

## Usage

```bash
# Start the web application
genetics-viz /path/to/data/directory

# With custom host and port
genetics-viz /path/to/data/directory --host 0.0.0.0 --port 8000

# With auto-reload for development
genetics-viz /path/to/data/directory --reload
```

## Data Directory Structure

The data directory should have the following structure:

```
data_directory/
└── cohorts/
    ├── cohort1/
    │   └── cohort1.pedigree.tsv
    ├── cohort2/
    │   └── cohort2.pedigree.tsv
    └── ...
```

### Pedigree File Format

Pedigree files should be tab-separated with the following columns. The header is optional - if present, it must start with "FID":

**With header:**

```
FID IID PAT MAT SEX PHENOTYPE
FAM001 IND001 0 0 M affected
```

**Without header (positional columns):**

```
FAM001 IND001 0 0 M affected
```

| Column | Possible Names (case-insensitive) |
|--------|-----------------------------------|
| Family ID | `FID`, `family_id`, `familyid`, `family` |
| Individual ID | `IID`, `individual_id`, `sample_id`, `sample` |
| Father ID | `PAT`, `father_id`, `fatherid`, `father` |
| Mother ID | `MAT`, `mother_id`, `motherid`, `mother` |
| Sex | `SEX`, `gender` |
| Phenotype | `PHENOTYPE`, `affected`, `status` |

## Development

```bash
# Install development dependencies
uv sync --dev

# Run tests
uv run pytest

# Run linter
uv run ruff check .

# Format code
uv run ruff format .
```

## License

MIT License
