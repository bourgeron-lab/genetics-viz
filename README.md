# Genetics-Viz 🧬

A web-based visualization tool for genetics cohort data, providing interactive analysis and validation of genetic variants.

## Features

### Core Features

- 📊 **Multi-Cohort Management** - Browse and analyze multiple cohorts from a single data directory
- 👨‍👩‍👧‍👦 **Family Structure Visualization** - View pedigree information and family relationships
- 🧬 **Variant Analysis** - Interactive tables for DNM (de novo mutations) and WOMBAT analysis
- ✅ **Variant Validation** - Track and validate genetic variants with inheritance patterns
- 🔬 **IGV Integration** - Built-in IGV.js browser for sequence visualization (CRAM files)
- 🌊 **WAVES Validation** - Specialized validation workflow for bedGraph/coverage analysis
- 🎨 **Modern UI** - Clean, responsive interface built with NiceGUI

### Validation Features

- Save validation status (present/absent/uncertain/different/in phase MNV)
- Track inheritance patterns (de novo/paternal/maternal/not paternal/not maternal/either/homozygous)
- Add optional comments to validations
- Mark validations as ignored (excluded from statistics and conflict detection)
- View validation history with timestamps and ignore status
- Interactive validation guide accessible via info button
- Filter variants by validation status
- Automatic conflict detection (ignoring validations marked as ignored)
- Export validation data

## Installation

### Quick Start with uvx (Recommended)

The easiest way to run genetics-viz without installation:

```bash
uvx genetics-viz /path/to/data/directory
```

### From PyPI

```bash
pip install genetics-viz
```

### From Source

```bash
# Clone the repository
git clone https://github.com/bourgeron-lab/genetics-viz.git
cd genetics-viz

# Install with uv (recommended)
uv sync
uv run genetics-viz /path/to/data/directory

# Or install with pip
pip install -e .
genetics-viz /path/to/data/directory
```

### Alternative: Local Python/Virtualenv

```bash
# Create a virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install genetics-viz
pip install genetics-viz

# Run the application
genetics-viz /path/to/data/directory
```

## Usage

### Command Line Options

```bash
# Basic usage
genetics-viz /path/to/data/directory

# With custom host and port
genetics-viz /path/to/data/directory --host 0.0.0.0 --port 8080

# Full help
genetics-viz --help
```

### Web Interface

Once started, open your browser to `http://localhost:8000` (or the specified port).

The interface provides:

- **Home Page** - List of available cohorts
- **Cohort View** - Family list and overview
- **Family View** - DNM and WOMBAT analysis tabs
- **Validation Pages** - Track variant validations (file-specific and all validations)
- **WAVES Validation** - Specialized coverage/bedGraph validation workflow

## Data Directory Structure

The tool expects the following directory structure:

```
data_directory/
├── cohorts/
│   ├── cohort1/
│   │   ├── cohort1.pedigree.tsv
│   │   └── families/
│   │       ├── FAM001/
│   │       │   ├── FAM001.wombat.*.tsv (WOMBAT analysis files)
│   │       │   └── FAM001.dnm.*.tsv (DNM analysis files)
│   │       └── FAM002/
│   │           └── ...
│   └── cohort2/
│       └── ...
├── samples/
│   ├── SAMPLE001/
│   │   └── sequences/
│   │       ├── SAMPLE001.GRCh38_GIABv3.cram
│   │       ├── SAMPLE001.GRCh38_GIABv3.cram.crai
│   │       └── SAMPLE001.GRCh38.bedGraph.gz (for WAVES)
│   └── SAMPLE002/
│       └── ...
└── validations/
    ├── snvs.tsv (variant validations)
    └── waves.tsv (WAVES validations)
```

### Required Files

#### Pedigree File Format

Pedigree files (`cohort_name.pedigree.tsv`) should be tab-separated. The header is optional - if present, it must start with "FID":

**With header:**

```
FID IID PAT MAT SEX PHENOTYPE
FAM001 SAMPLE001 0 0 M affected
FAM001 SAMPLE002 0 0 F unaffected
```

**Without header:**

```
FAM001 SAMPLE001 0 0 M affected
FAM001 SAMPLE002 0 0 F unaffected
```

**Column Mapping** (case-insensitive):

| Column | Possible Names |
|--------|----------------|
| Family ID | `FID`, `family_id`, `familyid`, `family` |
| Individual ID | `IID`, `individual_id`, `sample_id`, `sample` |
| Father ID | `PAT`, `father_id`, `fatherid`, `father` |
| Mother ID | `MAT`, `mother_id`, `motherid`, `mother` |
| Sex | `SEX`, `gender` |
| Phenotype | `PHENOTYPE`, `affected`, `status` |

#### CRAM Files (for IGV visualization)

- Format: `SAMPLE_ID.GRCh38_GIABv3.cram`
- Index: `SAMPLE_ID.GRCh38_GIABv3.cram.crai`
- Location: `samples/SAMPLE_ID/sequences/`

#### BedGraph Files (for WAVES validation)

- Format: `SAMPLE_ID.GRCh38.bedGraph.gz`
- Location: `samples/SAMPLE_ID/sequences/`

#### Analysis Files

- **DNM files**: `FAMILY_ID.dnm.*.tsv` (must contain `chr:pos:ref:alt` column)
- **WOMBAT files**: `FAMILY_ID.wombat.*.tsv` (must contain `#CHROM`, `POS`, `REF`, `ALT` columns)

## GHFC Lab Usage

### Prerequisites

You need to either:

- Be on the Institut Pasteur network, OR
- Be connected via VPN

### Mounting ghfc_wgs from Helix

#### On macOS

```bash
# Mount the network drive
# In Finder: Go > Connect to Server (⌘K)
# Enter: smb://helix.pasteur.fr/ghfc_wgs
# Or via command line:
open 'smb://helix.pasteur.fr/projects/ghfc_wgs'
```

The drive will be mounted at `/Volumes/ghfc_wgs`

#### On Linux

```bash
# Create mount point
sudo mkdir -p /mnt/ghfc_wgs

# Mount via CIFS
sudo mount -t cifs //helix.pasteur.fr/projects/ghfc_wgs /mnt/ghfc_wgs -o username=YOUR_USERNAME,domain=PASTEUR

# Or add to /etc/fstab for automatic mounting:
# //helix.pasteur.fr/projects/ghfc_wgs /mnt/ghfc_wgs cifs username=YOUR_USERNAME,password=YOUR_PASSWORD,domain=PASTEUR,uid=1000,gid=1000 0 0
```

### Running genetics-viz for GHFC Data

#### Method 1: Using uvx (Recommended - No Installation)

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Run directly with uvx
uvx genetics-viz /Volumes/ghfc_wgs/WGS/GHFC-GRCh38

# On Linux (adjust mount point):
uvx genetics-viz /mnt/ghfc_wgs/WGS/GHFC-GRCh38
```

#### Method 2: Using uv with Local Installation

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install genetics-viz
uv pip install genetics-viz

# Run the application
genetics-viz /Volumes/ghfc_wgs/WGS/GHFC-GRCh38
```

#### Method 3: Traditional Python/pip

```bash
# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate

# Install genetics-viz
pip install genetics-viz

# Run the application
genetics-viz /Volumes/ghfc_wgs/WGS/GHFC-GRCh38
```

### Access the Application

Once started, open your browser to:

```
http://localhost:8000
```

To access from other machines on the network:

```bash
genetics-viz /Volumes/ghfc_wgs/WGS/GHFC-GRCh38 --host 0.0.0.0 --port 8000
```

Then access via: `http://YOUR_MACHINE_IP:8000`

## Validation Workflow

### SNV Validation

1. Navigate to a variant table (DNM or WOMBAT tabs, or Validation pages)
2. Click "View in IGV" button for a variant
3. In the dialog:
   - Review variant details with collapsible sections
   - Add additional samples (parents, siblings, or by barcode)
   - Examine CRAM tracks in IGV viewer
   - Click the info button (ℹ️) to view validation guidelines
   - Set validation status (default: present) and inheritance pattern
   - Add an optional comment
   - Click "Save Validation"
4. The validation is saved to `validations/snvs.tsv`
5. View validation history below the form
   - Toggle the "Ignore" switch to exclude validations from statistics
   - Ignored validations appear with reduced opacity
6. Validation/all page aggregates multiple validations per variant/sample
   - Shows unique list of users who validated each variant
   - Computes final status from non-ignored validations

### WAVES Validation

1. Go to "Validation" > "Waves" in the menu
2. Select a cohort and pedigree
3. Select a sample from the pedigree
4. Click "View on IGV" for the sample
5. In the dialog:
   - Review bedGraph coverage tracks for the sample
   - Add additional samples for comparison
   - Set validation status
   - Click "Save Validation"
6. The validation is saved to `validations/waves.tsv`

## Development

```bash
# Clone repository
git clone https://github.com/bourgeron-lab/genetics-viz.git
cd genetics-viz

# Install with development dependencies
uv sync --dev

# Run tests
uv run pytest

# Run linter
uv run ruff check .

# Format code
uv run ruff format .

# Run with auto-reload for development
uv run genetics-viz --reload /path/to/data
```

## Validation File Formats

### SNV Validations (`validations/snvs.tsv`)

**Version 0.2.0+ format:**

```
FID Variant Sample User Inheritance Validation Comment Ignore Timestamp
FAM001 chr1:12345:A:T SAMPLE001 username de novo present Initial validation 0 2026-01-18T10:30:00
FAM001 chr1:12345:A:T SAMPLE001 reviewer homozygous present Confirmed 0 2026-01-19T14:20:00
FAM002 chr2:67890:G:C SAMPLE002 username unknown uncertain Low coverage 1 2026-01-18T11:00:00
```

**Columns:**

- **FID**: Family ID
- **Variant**: chr:pos:ref:alt format
- **Sample**: Sample ID
- **User**: Username who performed validation
- **Inheritance**: de novo, paternal, maternal, not paternal, not maternal, either, homozygous, or unknown
- **Validation**: present, absent, uncertain, different, or "in phase MNV"
- **Comment**: Optional free-text comment
- **Ignore**: 0 (included) or 1 (excluded from statistics and conflict detection)
- **Timestamp**: ISO format timestamp

**Migration from v0.1.1:**

If upgrading from v0.1.1, use the provided migration script:

```bash
./utils/snvs_validations_migration_0.1.1_to_0.2.0.sh /path/to/data/validations/snvs.tsv
```

This adds the `Comment` and `Ignore` columns with default values.

### WAVES Validations (`validations/waves.tsv`)

```
Cohort Pedigree Sample User Validation Timestamp
cohort1 FAM001 SAMPLE001 username present 2026-01-18T10:30:00
```

## Troubleshooting

### Cannot access GHFC data

- Verify VPN connection or Pasteur network access
- Check that ghfc_wgs is properly mounted
- Verify mount path (`/Volumes/ghfc_wgs` on macOS, `/mnt/ghfc_wgs` on Linux)

### IGV not displaying

- Ensure CRAM files and indices (.crai) exist
- Check that files follow naming convention: `SAMPLE_ID.GRCh38_GIABv3.cram`
- Verify IGV.js is loading (check browser console)

### Pedigree file not recognized

- Ensure tab-separated format
- Verify required columns are present
- Check file naming: `cohort_name.pedigree.tsv`

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

For detailed changes between versions, see [CHANGELOG.md](CHANGELOG.md).

## License

MIT License - See LICENSE file for details

## Citation

If you use this tool in your research, please cite:

```
Genetics-Viz: A web-based visualization tool for genetics cohort data
GitHub: https://github.com/bourgeron-lab/genetics-viz
```

## Support

For issues, questions, or feature requests, please open an issue on GitHub:
<https://github.com/bourgeron-lab/genetics-viz/issues>
