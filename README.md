# Genetics-Viz 🧬

A web-based visualization tool for genetics cohort data, providing interactive analysis and validation of genetic variants.

## Features

### Core Features

- 📊 **Multi-Cohort Management** - Browse and analyze multiple cohorts from a single data directory, with a Families tab, an Individuals tab and a cohort view per analysis
- 👨‍👩‍👧‍👦 **Family Structure Visualization** - Pedigree table or drawn tree (pypedigree), dummy members, and a pedigree editor for curators
- 🧬 **SNVs and Indels** - Interactive TanStack tables of the Wombat results, with the ghfc-ngs de novo step's flags, evidence and per-child QC
- 🧱 **Structural Variants** - WisecondorX and Vedette calls, copy number against a control panel, read evidence and suggested breakpoints
- 🔁 **STRs and UPD** - ExpansionHunter genotypes with REViewer images, and uniparental disomy candidates from the upd step, ROH and Sample QC
- 🧪 **Sample and Genome QC** - somalier relatedness, sex and contamination checks, DNA source and library (cell lines, read length, PCR), with curator acknowledgements
- 🌍 **Ancestry, ROH and Polygenic Scores** - PCA against a versioned reference panel, runs of homozygosity against reference groups, ancestry-matched PGS z-scores
- 🔍 **Cohort-Wide Search** - Search variants across all samples with filters on locus, genesets, impact, individuals (sex, phenotype, parental status), and validation status
- 🧭 **Gene and Variant Pages** - "Have we seen this before?": every SNV, SV, STR genotype, ROH segment, validation and diagnostic of a gene or a variant, across all cohorts and data directories
- 🏷️ **GeneTrek Gene Badges** - Genes coloured by the user's geneset tiers, with constraint, expression and modes of inheritance in the tooltip
- ✅ **Curation** - Validations and diagnostics for SNVs, SVs, STRs and UPDs, family statuses and conclusions, and a printable family report for the referring clinician
- 🔬 **IGV Integration** - Built-in IGV.js browser for sequence visualization (CRAM files)
- 🌊 **WAVES Validation** - Specialized validation workflow for bedGraph/coverage analysis
- ⚙️ **Workflow** - The ghfc-ngs run status, parameters and reference, and run requests from curators
- 🔐 **Authentication & Authorization** - YAML-configured user accounts with role-based access (student, reader, curator, administrator) and data directories granted per user
- 📂 **Multi-Data-Directory** - Switch between multiple data directories per user session via YAML config
- 🛠️ **Admin Pages** - Manage users and data directories from the web interface

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
uvx genetics-viz /path/to/config.yaml
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
uv run genetics-viz /path/to/config.yaml

# Or install with pip
pip install -e .
genetics-viz /path/to/config.yaml
```

### Alternative: Local Python/Virtualenv

```bash
# Create a virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install genetics-viz
pip install genetics-viz

# Run the application
genetics-viz /path/to/config.yaml
```

## Configuration

### YAML Config File

The application requires a YAML configuration file listing data directories and users:

```yaml
data_directories:
  - path: /path/to/data/directory1
    description: "Primary WGS data"
    default: true
  - path: /path/to/data/directory2
    description: "Secondary dataset"

user_list:
  - username: admin
    password: "<scrypt hash>"
    role: administrator
  - username: curator1
    password: "<scrypt hash>"
    role: curator
  - username: viewer
    password: "<scrypt hash>"
    role: reader
  - username: student1
    password: "<scrypt hash>"
    role: student
    # Optional: the data directories this user may open (default: all)
    data_directories:
      - /path/to/data/directory2

# Optional. Seconds between checks of the data directories for changes.
# poll_interval: 30

# Optional. The FASTA the CRAMs were aligned to (GRCh38_GIABv3), for the SV
# dialog's read evidence. Without it the CRAMs are not decoded in Python.
# reference_fasta: /path/to/GRCh38_GIABv3.fa

# Optional. Root of the ancestry/PGS references tree, holding one
# v<bundle_version>/ subdirectory per reference bundle. The version to use is
# derived from the per-family data filenames, so it is not configured here.
# Without this key, the Ancestry tab still plots each family's own samples but
# cannot draw the reference panel behind them.
ancestry_reference_dir: /path/to/references/ancestry-pgs

# Optional. The GeneTrek instance gene badges link to and take their genesets
# from (default: the dev instance). Its data is cached in genetrek_cache_dir
# (default: .genetrek/ next to this file) and refreshed when GeneTrek
# publishes a new release, so the app still works offline on the last one.
genetrek_url: https://genetrek-py.dev.pasteur.cloud/
# genetrek_cache_dir: /path/to/cache

# Optional. Where the gene index behind the gene and variant pages is kept
# (default: .gene_index/ next to this file). It copies, by gene, every
# cohort's Wombat, WisecondorX, ExpansionHunter and ROH genes tables, so keep
# it as private as the data. It is built at startup, off the request path,
# and only tables that changed are read again.
# gene_index_dir: /path/to/gene_index

# Auto-generated on first run — do not edit manually
storage_secret: "<hex-string>"
```

Each user's gene badge settings (which GeneTrek genesets colour badges, and
how) are saved in `<config name>.preferences.yaml` next to the config file,
from the Profile page. An administrator can make their own settings the default
for everyone.

### Generating Password Hashes

Passwords are stored as salted scrypt. Add users from **Admin → Users**, which
generates a password; for the first administrator, hash one with:

```bash
uv run python -c "from genetics_viz.config_model import hash_password; print(hash_password('your_password'))"
```

A legacy unsalted SHA-512 hex digest (`sha512sum`) still works, and is rehashed
to scrypt at that user's first login.

### User Roles

| Role | Permissions |
|------|-------------|
| **student** | Browse, reader wording, no downloads (no Download buttons, no TSV from Search) |
| **reader** | Browse and download the tables, cannot save validations or diagnostics |
| **curator** | Reader + save validations and diagnostics, family reports, every raw data file |
| **administrator** | All curator permissions + manage users and data directories |

Readers and students fetch from a data directory only the files the pages load
themselves (IGV tracks, REViewer images); curators and administrators may
fetch any file. A user with `data_directories` sees only those directories.

### Initial Setup

1. Create a YAML config file with at least one data directory and one administrator user
2. Generate the admin password hash (see above)
3. Run: `genetics-viz /path/to/config.yaml`
4. Log in with admin credentials at `http://localhost:8080/login`
5. Use the admin pages to add more users and data directories

## Usage

### Command Line Options

```bash
# Basic usage
genetics-viz /path/to/config.yaml

# With custom host and port
genetics-viz /path/to/config.yaml --host 0.0.0.0 --port 8080

# Full help
genetics-viz --help
```

### Web Interface

Once started, open your browser to `http://localhost:8080` (or the specified port). You will be redirected to the login page.

The interface provides:

- **Home Page** - Cohort cards with counts, QC and review status
- **Cohort View** - Families and Individuals tabs, then a tab per check, variant type and background analysis, and the Workflow panel
- **Family View** - Sample QC, SNVs/Indels, SVs, STR, UPD, Ancestry, ROH and Polygenic Scores tabs, the family status, diagnostics and report
- **Search Page** - Cohort-wide variant search with tabbed filters (Variants and Individuals)
- **Gene and Variant Pages** - Type a gene symbol or a variant (`chr22:50674000:G:A`) in the header box: `/gene/<symbol>` and `/variant/<variant>` gather it from every cohort of every data directory
- **Variant Statistics** - Charts and ideogram views for search results
- **Validation Pages** - Track variant validations (file-specific and all validations)
- **Diagnostic Pages** - Track variant diagnostic conclusions
- **WAVES Validation** - Specialized coverage/bedGraph validation workflow
- **Profile Page** - View role and change password
- **Admin Pages** - Manage data directories and users (administrators only)

## Data Directory Structure

The tool reads the layout the [ghfc-ngs](https://github.com/bourgeron-lab/ghfc-ngs) workflow
writes. Every part but the cohort pedigree is optional: a tab without its files is disabled.

```
data_directory/
├── cohorts/
│   └── COHORT/
│       ├── COHORT.pedigree.tsv          (required)
│       ├── COHORT.params.yml            (ghfc-ngs parameters: Workflow panel)
│       ├── .ghfc-ngs.state.json         (ghfc-ngs run record)
│       ├── .ghfc-ngs.requests/          (run requests, when enabled)
│       ├── wombat/   COHORT.rare.<vep>[.annotated].<config>.results.tsv
│       ├── denovo/   COHORT.<label>.{dnm,qc,inputs}.tsv
│       ├── qc/       COHORT.<label>.{samples,pairs,flags,inputs}.tsv (sampleqc, genomeqc)
│       ├── ancestry/ COHORT.apgs_b<bundle>_<filters>.{pcs,ancestry,pgs_zscore,Q}.tsv
│       ├── roh/      COHORT.apgs_b<bundle>_<filters>.{froh,roh,roh_genes,kinship,roh_windows}.tsv
│       ├── str/expansionhunter/  COHORT.<catalog>.{genotypes,loci,flags}.tsv
│       ├── upd/      COHORT.<label>.{events,segments,members,inputs}.tsv
│       ├── svs/wisecondorx/      COHORT_aberrations[.annotated].bed
│       ├── coverage/ COHORT.<tag>.panel_v<version>.fit.tsv
│       └── cnv_panel/v<version>/ (a controls cohort's copy-number panel)
├── families/<shards>/FID/               (same subdirectories, per family, plus:)
│   ├── FID.pedigree.tsv
│   ├── wombat/   FID.rare.<vep>[.annotated].<config>.tsv
│   ├── svs/vedette/      FID.vedette.{tsv,parquet,qc.tsv,json}
│   ├── svs/wisecondorx/  FID_aberrations.annotated.bed
│   └── pedigree_edits/   (pedigree edits made in the app)
├── samples/<shards>/ID/
│   ├── sequences/  ID.GRCh38_GIABv3.cram(.crai), ID.by1000.bedgraph.gz(.tbi) (WAVES)
│   ├── coverage/   ID.cov_b<bin>q<mapq>.{bw,json}
│   ├── str/expansionhunter/  ID.<catalog>.{vcf,json}, reviewer/
│   ├── qc/genomeqc/          ID.<label>.{tsv,stats.txt}
│   └── svs/vedette/          ID.filtered.vcf.gz
├── params/genesets/*.tsv     (local gene sets for the geneset filters)
├── validations/              (snvs, svs, strs, upds, vedette, waves, sampleqc, pedigree edits)
├── diagnostics/              (snvs, svs, strs, upds)
├── notes/notes.tsv
└── reports/                  (family conclusions and statuses)
```

`samples/` and `families/` may be sharded two levels deep (`samples/1/S/S1/`): the last and
second-to-last characters of the ID, with `-`, `.` and `_` removed. This is detected per data
directory.

### Required Files

A cohort is any directory under `cohorts/`, and the directory name is the cohort name. A
directory whose `<name>.pedigree.tsv` is missing or unreadable is listed on the home page as an
incomplete cohort rather than being skipped, and is not explorable until the pedigree is in
place.

#### Optional Workflow Files

`COHORT.params.yml` is the parameter file the [ghfc-ngs](https://github.com/bourgeron-lab/ghfc-ngs)
workflow was launched with. When it is present, the cohort page gains a **Workflow** button:
the run status per step, the parameters verbatim and the workflow's parameter reference. The home
page card gains a foldable **Status** item with the same per-step breakdown.

`.ghfc-ngs.state.json` is the workflow's own run record (see `COHORT_STATE.md` in the ghfc-ngs
repository) and is written only by the pipeline. It supplies the run status, timings and the
per-step completion counts; genetics-viz reads it and never measures completion itself. Its
absence is reported as "no run recorded", which is not the same as the cohort never having been
processed.

#### Pedigree File Format

Pedigree files (`cohort_name.pedigree.tsv`) should be tab-separated. The header is optional - if present, it must start with "FID" (a leading `#` is stripped automatically):

**With header:**

```
#FID	IID	PAT	MAT	SEX	PHENOTYPE
FAM001	SAMPLE001	SAMPLE003	SAMPLE004	1	2
FAM001	SAMPLE002	0	0	2	1
```

**Without header (positional columns):**

```
FAM001	SAMPLE001	SAMPLE003	SAMPLE004	1	2
FAM001	SAMPLE002	0	0	2	1
```

Missing/unknown values for parent IDs are `0`, `-9`, or empty. These values are also treated as unknown for sex and phenotype when building filter options.

**Column Mapping** (case-insensitive, `#` prefix stripped):

| Column | Possible Names |
|--------|----------------|
| Family ID | `FID`, `family_id`, `familyid`, `family` |
| Individual ID | `IID`, `individual_id`, `sample_id`, `sampleid`, `sample` |
| Father ID | `PAT`, `father_id`, `fatherid`, `father`, `paternal_id` |
| Mother ID | `MAT`, `mother_id`, `motherid`, `mother`, `maternal_id` |
| Sex | `SEX`, `gender` |
| Phenotype | `PHENOTYPE`, `affected`, `status`, `affection` |

#### CRAM Files (for IGV visualization)

- Format: `SAMPLE_ID.GRCh38_GIABv3.cram`
- Index: `SAMPLE_ID.GRCh38_GIABv3.cram.crai`
- Location: `samples/SAMPLE_ID/sequences/`

#### BedGraph Files (for WAVES validation)

- Format: `SAMPLE_ID.by1000.bedgraph.gz`, with its tabix index (`.tbi`)
- Location: `samples/SAMPLE_ID/sequences/`

#### Analysis Files

- **Wombat files**: `families/.../FAMILY_ID/wombat/FAMILY_ID.rare.<vep>[.annotated].<config>.tsv` (must contain `#CHROM`, `POS`, `REF`, `ALT` columns). De novo calls are a Wombat config (`de_novo_mutations` by default), checked by the ghfc-ngs `denovo` step when its tables exist
- **Locations**: `families/` and `samples/` may be sharded; see the tree above

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
uvx genetics-viz /path/to/genetics_viz.yaml

# On Linux (adjust path):
uvx genetics-viz /path/to/genetics_viz.yaml
```

#### Method 2: Using uv with Local Installation

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install genetics-viz
uv pip install genetics-viz

# Run the application
genetics-viz /path/to/genetics_viz.yaml
```

### Access the Application

Once started, open your browser to:

```
http://localhost:8080
```

To access from other machines on the network:

```bash
genetics-viz /path/to/genetics_viz.yaml --host 0.0.0.0 --port 8080
```

Then access via: `http://YOUR_MACHINE_IP:8080`

## Validation Workflow

### SNV Validation

1. Navigate to a variant table (the family SNVs/Indels tab, Search, or the Validation pages)
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
uv run genetics-viz --reload /path/to/config.yaml
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
