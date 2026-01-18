"""
Data models for genetics-viz.

This module provides classes for loading and managing cohort data,
including pedigree information, families, and samples.
"""

from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

# Standard pedigree column names (positional, no header)
PEDIGREE_COLUMNS = ["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"]


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
        return (
            self.father_id is None or self.father_id == "0" or self.father_id == ""
        ) and (self.mother_id is None or self.mother_id == "0" or self.mother_id == "")


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

        Handles both files with and without headers.
        If header is present (starts with "FID"), it is skipped.
        """
        # Read first line to check for header
        with open(self.pedigree_file) as f:
            first_line = f.readline().strip()

        has_header = first_line.upper().startswith("FID")

        # Always read without header, skip first line if it's a header
        df = pl.read_csv(
            self.pedigree_file,
            separator="\t",
            has_header=False,
            skip_rows=1 if has_header else 0,
            infer_schema_length=0,  # Read all as strings
        )

        # Assign standard column names based on position
        num_cols = len(df.columns)
        col_names = PEDIGREE_COLUMNS[:num_cols]
        # Pad with generic names if more columns than expected
        while len(col_names) < num_cols:
            col_names.append(f"COL{len(col_names) + 1}")
        df = df.rename({f"column_{i + 1}": name for i, name in enumerate(col_names)})

        return df

    @classmethod
    def from_directory(cls, path: Path) -> "Cohort":
        """
        Create a Cohort from a directory containing a pedigree file.

        The pedigree file should be named {cohort_name}.pedigree.tsv
        """
        name = path.name
        pedigree_file = path / f"{name}.pedigree.tsv"

        if not pedigree_file.exists():
            raise FileNotFoundError(f"Pedigree file not found: {pedigree_file}")

        cohort = cls(name=name, path=path, pedigree_file=pedigree_file)
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
                col: str | None, treat_zero_as_null: bool = False
            ) -> str | None:
                if col is None:
                    return None
                val = row.get(col)
                # Check for None or empty string
                if val is None or val == "":
                    return None
                # For father/mother fields, "0" means no parent (founder)
                if treat_zero_as_null and val == "0":
                    return None
                return str(val)

            sample = Sample(
                sample_id=sample_id,
                family_id=family_id,
                father_id=get_value(father_col, treat_zero_as_null=True),
                mother_id=get_value(mother_col, treat_zero_as_null=True),
                sex=get_value(sex_col),
                phenotype=get_value(phenotype_col),
            )

            if family_id not in self.families:
                self.families[family_id] = Family(family_id=family_id)

            self.families[family_id].samples.append(sample)

    def _identify_columns(self, df: pl.DataFrame) -> dict[str, str | None]:
        """Identify column names from various naming conventions."""
        columns = {c.upper(): c for c in df.columns}

        mapping: dict[str, str | None] = {}

        # Family ID column
        for name in ["FID", "FAMILY_ID", "FAMILYID", "FAMILY", "#FAMILY_ID"]:
            if name in columns:
                mapping["family_id"] = columns[name]
                break

        # Sample/Individual ID column
        for name in [
            "IID",
            "INDIVIDUAL_ID",
            "SAMPLE_ID",
            "SAMPLEID",
            "SAMPLE",
            "INDIVIDUAL",
            "ID",
        ]:
            if name in columns:
                mapping["sample_id"] = columns[name]
                break

        # Father ID column
        for name in ["PAT", "FATHER_ID", "FATHERID", "FATHER", "PATERNAL_ID"]:
            if name in columns:
                mapping["father_id"] = columns[name]
                break

        # Mother ID column
        for name in ["MAT", "MOTHER_ID", "MOTHERID", "MOTHER", "MATERNAL_ID"]:
            if name in columns:
                mapping["mother_id"] = columns[name]
                break

        # Sex column
        for name in ["SEX", "GENDER"]:
            if name in columns:
                mapping["sex"] = columns[name]
                break

        # Phenotype column
        for name in ["PHENOTYPE", "AFFECTED", "STATUS", "AFFECTION"]:
            if name in columns:
                mapping["phenotype"] = columns[name]
                break

        return mapping

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


@dataclass
class DataStore:
    """
    Central data store for all cohorts in a data directory.

    This class manages loading and caching of cohort data.
    """

    data_dir: Path
    cohorts: dict[str, Cohort] = field(default_factory=dict)
    _loaded: bool = field(default=False, repr=False)

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

        self.cohorts = {}

        for cohort_path in sorted(self.cohorts_dir.iterdir()):
            if not cohort_path.is_dir():
                continue

            pedigree_file = cohort_path / f"{cohort_path.name}.pedigree.tsv"
            if not pedigree_file.exists():
                # Skip directories without pedigree files
                continue

            try:
                cohort = Cohort.from_directory(cohort_path)
                self.cohorts[cohort.name] = cohort
            except Exception as e:
                print(f"Warning: Failed to load cohort {cohort_path.name}: {e}")

        self._loaded = True

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
