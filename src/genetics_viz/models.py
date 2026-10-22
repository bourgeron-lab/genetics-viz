"""
Data models for genetics-viz.

This module provides classes for loading and managing cohort data,
including pedigree information, families, and samples.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

logger = logging.getLogger(__name__)

# Standard pedigree column names (positional, no header)
PEDIGREE_COLUMNS = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]

#: Parent IDs meaning "no parent". Some cohorts write the codes
#: float-stringified (``0.0``, ``-9.0``), as they do for sex and phenotype.
MISSING_PARENT = frozenset({"", "0", "-9", "0.0", "-9.0"})

#: Header names for each pedigree field, in order of preference. Matched
#: case-insensitively; FatherBarcode/MotherBarcode are the per-family format.
PEDIGREE_COLUMN_CANDIDATES: dict[str, tuple[str, ...]] = {
    "family_id": ("FID", "FAMILY_ID", "FAMILYID", "FAMILY", "#FAMILY_ID"),
    "sample_id": (
        "IID",
        "INDIVIDUAL_ID",
        "SAMPLE_ID",
        "SAMPLEID",
        "SAMPLE",
        "INDIVIDUAL",
        "ID",
    ),
    "father_id": (
        "PAT",
        "FATHER_ID",
        "FATHERID",
        "FATHER",
        "PATERNAL_ID",
        "FATHERBARCODE",
    ),
    "mother_id": (
        "MAT",
        "MOTHER_ID",
        "MOTHERID",
        "MOTHER",
        "MATERNAL_ID",
        "MOTHERBARCODE",
    ),
    "sex": ("SEX", "GENDER"),
    "phenotype": ("PHENOTYPE", "AFFECTED", "STATUS", "AFFECTION"),
}

_HEADER_NAMES = frozenset(
    name.lstrip("#") for names in PEDIGREE_COLUMN_CANDIDATES.values() for name in names
)


def identify_pedigree_columns(columns: list[str]) -> dict[str, str | None]:
    """Map each pedigree field to the column holding it, by header name.

    Any column starting with ``PHENO_`` stands in for a missing phenotype.
    """
    by_upper = {c.upper(): c for c in columns}
    mapping: dict[str, str | None] = {}
    for key, candidates in PEDIGREE_COLUMN_CANDIDATES.items():
        for name in candidates:
            if name in by_upper:
                mapping[key] = by_upper[name]
                break
    if "phenotype" not in mapping:
        for upper_name, original_name in by_upper.items():
            if upper_name.startswith("PHENO_"):
                mapping["phenotype"] = original_name
                break
    return mapping


def pedigree_has_header(pedigree_file: Path) -> bool:
    """Whether the pedigree's first line is a header rather than a member.

    A header names at least one known pedigree column. A prefix test on
    ``FID`` would also take a headerless file whose first family is
    ``FID001`` for a header, and would miss ``sample_id``-first headers.
    """
    with open(pedigree_file) as f:
        first_line = f.readline()
    fields = (field.strip().lstrip("#").upper() for field in first_line.split("\t"))
    return any(f in _HEADER_NAMES or f.startswith("PHENO_") for f in fields)


def read_positional_pedigree(pedigree_file: Path, skip_header: bool) -> pl.DataFrame:
    """Read a pedigree as strings, naming the columns by position.

    The first six are ``PEDIGREE_COLUMNS`` (FID, IID, PAT, MAT, SEX,
    PHENOTYPE), any further ones ``COL7``, ``COL8``...
    """
    df = pl.read_csv(
        pedigree_file,
        separator="\t",
        has_header=False,
        skip_rows=1 if skip_header else 0,
        infer_schema_length=0,  # Read all as strings
    )
    num_cols = len(df.columns)
    col_names = PEDIGREE_COLUMNS[:num_cols]
    # Pad with generic names if more columns than expected
    while len(col_names) < num_cols:
        col_names.append(f"COL{len(col_names) + 1}")
    # Rename by position, not by the generated names: polars 1.x calls
    # headerless columns column_1..column_N, polars 2.x column_0..column_N-1.
    return df.rename(dict(zip(df.columns, col_names)))


@dataclass
class Sample:
    """Represents an individual sample in a pedigree."""

    sample_id: str
    family_id: str
    father_id: str | None = None
    mother_id: str | None = None
    sex: str | None = None
    phenotype: str | None = None

    @property
    def is_founder(self) -> bool:
        """Check if this sample is a founder (no parents in pedigree)."""
        return (self.father_id is None or self.father_id in MISSING_PARENT) and (
            self.mother_id is None or self.mother_id in MISSING_PARENT
        )


@dataclass
class Family:
    """Represents a family in a cohort."""

    family_id: str
    samples: list[Sample] = field(default_factory=list)

    @property
    def num_samples(self) -> int:
        """Return the number of samples in this family."""
        return len(self.samples)

    @property
    def num_founders(self) -> int:
        """Return the number of founders in this family."""
        return sum(1 for s in self.samples if s.is_founder)

    def get_sample(self, sample_id: str) -> Sample | None:
        """Get a sample by ID."""
        for sample in self.samples:
            if sample.sample_id == sample_id:
                return sample
        return None


@dataclass
class Cohort:
    """Represents a cohort/project containing multiple families."""

    name: str
    path: Path
    pedigree_file: Path
    #: The ghfc-ngs workflow parameters file, when the directory has one.
    #: Resolved during discovery so a page render never stats for it.
    params_file: Path | None = None
    families: dict[str, Family] = field(default_factory=dict)
    _dataframe: pl.DataFrame | None = field(default=None, repr=False)

    @property
    def num_families(self) -> int:
        """Return the number of families in this cohort."""
        return len(self.families)

    @property
    def num_samples(self) -> int:
        """Return the total number of samples across all families."""
        return sum(f.num_samples for f in self.families.values())

    @property
    def dataframe(self) -> pl.DataFrame:
        """Return the raw pedigree dataframe."""
        if self._dataframe is None:
            self._dataframe = self._load_pedigree()
        return self._dataframe

    def _load_pedigree(self) -> pl.DataFrame:
        """
        Load the pedigree file into a DataFrame.

        Handles both files with and without headers. A header is skipped:
        the columns are named by position either way.
        """
        return read_positional_pedigree(
            self.pedigree_file, skip_header=pedigree_has_header(self.pedigree_file)
        )

    @classmethod
    def from_directory(cls, path: Path, params_file: Path | None = None) -> "Cohort":
        """
        Create a Cohort from a directory containing a pedigree file.

        The pedigree file should be named {cohort_name}.pedigree.tsv

        ``params_file`` is passed in by the directory scan, which has already
        looked for it; it is resolved here when called directly.
        """
        name = path.name
        pedigree_file = path / f"{name}.pedigree.tsv"

        if not pedigree_file.exists():
            raise FileNotFoundError(f"Pedigree file not found: {pedigree_file}")

        if params_file is None:
            # Deferred for the same reason as in _iter_cohort_dirs.
            from genetics_viz.utils.pipeline_params import find_params_file

            params_file = find_params_file(path, name)

        cohort = cls(
            name=name,
            path=path,
            pedigree_file=pedigree_file,
            params_file=params_file,
        )
        cohort._parse_pedigree()
        return cohort

    def _parse_pedigree(self) -> None:
        """Parse the pedigree file and populate families and samples."""
        df = self.dataframe

        # Map columns to standard names
        col_mapping = self._identify_columns(df)

        family_col = col_mapping.get("family_id")
        sample_col = col_mapping.get("sample_id")
        father_col = col_mapping.get("father_id")
        mother_col = col_mapping.get("mother_id")
        sex_col = col_mapping.get("sex")
        phenotype_col = col_mapping.get("phenotype")

        if family_col is None or sample_col is None:
            raise ValueError(
                f"Could not identify required columns (family_id, sample_id) "
                f"in pedigree file. Found columns: {df.columns}"
            )

        self.families = {}

        for row in df.iter_rows(named=True):
            family_id = str(row[family_col])
            sample_id = str(row[sample_col])

            def get_value(
                col: str | None, treat_missing_as_null: bool = False
            ) -> str | None:
                if col is None:
                    return None
                val = row.get(col)
                # Check for None or empty string
                if val is None or val == "":
                    return None
                # For father/mother fields, "0" and "-9" mean no parent
                if treat_missing_as_null and val in MISSING_PARENT:
                    return None
                return str(val)

            sample = Sample(
                sample_id=sample_id,
                family_id=family_id,
                father_id=get_value(father_col, treat_missing_as_null=True),
                mother_id=get_value(mother_col, treat_missing_as_null=True),
                sex=get_value(sex_col),
                phenotype=get_value(phenotype_col),
            )

            if family_id not in self.families:
                self.families[family_id] = Family(family_id=family_id)

            self.families[family_id].samples.append(sample)

    def _identify_columns(self, df: pl.DataFrame) -> dict[str, str | None]:
        """Identify column names from various naming conventions."""
        return identify_pedigree_columns(df.columns)

    def get_families_summary(self) -> list[dict]:
        """Get a summary of all families as a list of dicts (for NiceGUI tables)."""
        data = []
        for family in self.families.values():
            data.append(
                {
                    "Family ID": family.family_id,
                    "Members": family.num_samples,
                    "Founders": family.num_founders,
                }
            )
        return data

    def get_family_members(self, family_id: str) -> list[dict]:
        """Get members of a specific family as a list of dicts (for NiceGUI tables)."""
        family = self.families.get(family_id)
        if family is None:
            return []

        data = []
        for sample in family.samples:
            data.append(
                {
                    "Sample ID": sample.sample_id,
                    "Father": sample.father_id or "-",
                    "Mother": sample.mother_id or "-",
                    "Sex": sample.sex or "-",
                    "Phenotype": sample.phenotype or "-",
                }
            )
        return data


@dataclass(frozen=True)
class CohortStub:
    """A cohorts/ directory with no pedigree file.

    A cohort is a directory under ``cohorts/``; the pedigree is what makes it
    explorable. A directory missing one is reported rather than skipped, so a
    cohort that exists but is not ready to browse is visible as such.
    """

    name: str
    path: Path
    #: cohorts/<name>/<name>.pedigree.tsv.
    expected_pedigree: Path
    params_file: Path | None = None
    #: Why the cohort is not explorable: "missing" when there is no pedigree
    #: file, "unreadable" when there is one but it could not be parsed. The two
    #: call for different fixes, so they are not collapsed into one message.
    reason: str = "missing"
    #: The parse error, when reason is "unreadable".
    error: str | None = None


def _iter_cohort_dirs(
    cohorts_dir: Path,
) -> Iterator[tuple[Path, Path, Path | None]]:
    """Yield ``(directory, pedigree_file, params_file)`` for every cohort dir.

    The single definition of cohort discovery: every directory under
    ``cohorts/`` is a cohort, except hidden ones (a NAS ``.snapshot``). ``pedigree_file`` is the expected path and may not
    exist -- callers decide what that means. ``params_file`` is resolved, so it
    is None when the directory has none.

    This exists because the same scan is needed by ``load``, ``take_snapshot``
    and ``reload``; three hand-copied versions of it drifted the moment the
    rules changed.
    """
    # Deferred: genetics_viz.utils.__init__ imports utils.data, which imports
    # this module, so any utils import at module level here is a cycle.
    from genetics_viz.utils.pipeline_params import find_params_file

    if not cohorts_dir.exists():
        return
    for cohort_path in sorted(cohorts_dir.iterdir()):
        if cohort_path.name.startswith(".") or not cohort_path.is_dir():
            continue
        name = cohort_path.name
        yield (
            cohort_path,
            cohort_path / f"{name}.pedigree.tsv",
            find_params_file(cohort_path, name),
        )


def _mtime_or_zero(path: Path | None) -> float:
    """Modification time of a watched file, 0.0 when it is absent."""
    if path is None:
        return 0.0
    try:
        return path.stat().st_mtime
    except OSError:
        return 0.0


@dataclass
class DirectorySnapshot:
    """Point-in-time snapshot of a DataStore's cohort directory state."""

    #: Cohort directory name -> the newest mtime among the files that change
    #: what the UI shows: the pedigree and the workflow parameters file. Every
    #: cohort directory is keyed, including those with no pedigree, so one
    #: appearing or gaining a pedigree registers as a change.
    cohort_mtimes: dict[str, float] = field(default_factory=dict)
    cohort_names: frozenset[str] = field(default_factory=frozenset)


@dataclass
class ChangeReport:
    """Describes what changed between two directory snapshots."""

    data_dir: str
    added_cohorts: list[str] = field(default_factory=list)
    removed_cohorts: list[str] = field(default_factory=list)
    modified_cohorts: list[str] = field(default_factory=list)

    @property
    def has_changes(self) -> bool:
        return bool(self.added_cohorts or self.removed_cohorts or self.modified_cohorts)

    def summary_lines(self) -> list[str]:
        lines: list[str] = []
        if self.added_cohorts:
            lines.append(f"New cohorts: {', '.join(self.added_cohorts)}")
        if self.removed_cohorts:
            lines.append(f"Removed cohorts: {', '.join(self.removed_cohorts)}")
        if self.modified_cohorts:
            lines.append(f"Updated cohorts: {', '.join(self.modified_cohorts)}")
        return lines


@dataclass
class DataStore:
    """
    Central data store for all cohorts in a data directory.

    This class manages loading and caching of cohort data.
    """

    data_dir: Path
    #: Pedigree-bearing cohorts only, so every existing consumer of this dict
    #: -- the header dropdown, the validation pages, search -- keeps its meaning.
    cohorts: dict[str, Cohort] = field(default_factory=dict)
    #: Cohort directories with no usable pedigree. Carded on the home page, but
    #: not explorable.
    incomplete_cohorts: list[CohortStub] = field(default_factory=list)
    _loaded: bool = field(default=False, repr=False)
    _snapshot: DirectorySnapshot | None = field(default=None, repr=False)

    @property
    def cohorts_dir(self) -> Path:
        """Return the path to the cohorts directory."""
        return self.data_dir / "cohorts"

    def load(self) -> None:
        """Load all cohorts from the data directory."""
        if self._loaded:
            return

        if not self.cohorts_dir.exists():
            raise FileNotFoundError(f"Cohorts directory not found: {self.cohorts_dir}")

        self.cohorts, self.incomplete_cohorts = self._scan()
        self._loaded = True
        self._snapshot = self.take_snapshot()

    def _scan(self) -> tuple[dict[str, Cohort], list[CohortStub]]:
        """Scan the cohorts directory into loaded cohorts and pedigree-less stubs."""
        cohorts: dict[str, Cohort] = {}
        stubs: list[CohortStub] = []

        for cohort_path, pedigree_file, params_file in _iter_cohort_dirs(
            self.cohorts_dir
        ):
            if not pedigree_file.exists():
                stubs.append(
                    CohortStub(
                        name=cohort_path.name,
                        path=cohort_path,
                        expected_pedigree=pedigree_file,
                        params_file=params_file,
                        reason="missing",
                    )
                )
                continue

            try:
                cohort = Cohort.from_directory(cohort_path, params_file=params_file)
                cohorts[cohort.name] = cohort
            except Exception as e:
                # The pedigree is there but unusable, which leaves the cohort
                # just as unexplorable -- and needs a different fix than an
                # absent file, so it is reported separately.
                logger.warning("Failed to load cohort %s: %s", cohort_path.name, e)
                stubs.append(
                    CohortStub(
                        name=cohort_path.name,
                        path=cohort_path,
                        expected_pedigree=pedigree_file,
                        params_file=params_file,
                        reason="unreadable",
                        error=str(e),
                    )
                )

        return cohorts, stubs

    def take_snapshot(self) -> DirectorySnapshot:
        """Stat the cohort directory and return a lightweight snapshot of mtimes.

        Watches the pedigree and the workflow parameters file: the first decides
        whether a cohort is explorable, the second whether its Parameters and
        Status tabs exist.

        ``.ghfc-ngs.state.json`` is deliberately *not* watched. The pipeline
        rewrites it on every run, and a change here triggers ``reload``, which
        re-parses every pedigree in the data directory; a status write must not
        cost that. The state reader caches on the file's own stat instead, so it
        picks up a rewrite without help from here.
        """
        mtimes = {
            cohort_path.name: max(
                _mtime_or_zero(pedigree_file), _mtime_or_zero(params_file)
            )
            for cohort_path, pedigree_file, params_file in _iter_cohort_dirs(
                self.cohorts_dir
            )
        }
        return DirectorySnapshot(
            cohort_mtimes=mtimes,
            cohort_names=frozenset(mtimes),
        )

    @staticmethod
    def compare_snapshot(
        old: DirectorySnapshot, new: DirectorySnapshot
    ) -> ChangeReport:
        """Compare two snapshots and return a report of what changed."""
        added = sorted(new.cohort_names - old.cohort_names)
        removed = sorted(old.cohort_names - new.cohort_names)
        modified = sorted(
            name
            for name in old.cohort_names & new.cohort_names
            if old.cohort_mtimes[name] != new.cohort_mtimes[name]
        )
        return ChangeReport(
            data_dir="",
            added_cohorts=added,
            removed_cohorts=removed,
            modified_cohorts=modified,
        )

    def reload(self) -> None:
        """Re-scan the data directory and atomically replace the cohorts dict."""
        if not self.cohorts_dir.exists():
            return

        new_cohorts, new_stubs = self._scan()
        self.cohorts = new_cohorts
        self.incomplete_cohorts = new_stubs
        self._loaded = True
        self._snapshot = self.take_snapshot()

    def get_cohort(self, name: str) -> Cohort | None:
        """Get a cohort by name."""
        self.load()
        return self.cohorts.get(name)

    def get_cohorts_summary(self) -> list[dict]:
        """Get a summary of all cohorts as a list of dicts."""
        self.load()

        data = []
        for cohort in self.cohorts.values():
            data.append(
                {
                    "Cohort": cohort.name,
                    "Families": cohort.num_families,
                    "Samples": cohort.num_samples,
                }
            )

        return data
