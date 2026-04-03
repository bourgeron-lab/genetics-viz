"""Standalone pedigree file parser for per-family pedigree TSV files.

Parses pedigree files independently of the Cohort model, supporting
the per-family format (FatherBarcode/MotherBarcode/sample_id columns)
as well as the standard cohort format (PAT/MAT/IID columns).
"""

from pathlib import Path

import polars as pl

from genetics_viz.models import Family, Sample


def _identify_columns(df: pl.DataFrame) -> dict[str, str | None]:
    """Identify pedigree column names from various naming conventions.

    Extends the Cohort._identify_columns() mappings with per-family
    pedigree column names (FatherBarcode, MotherBarcode, sample_id).
    """
    columns = {c.upper(): c for c in df.columns}
    mapping: dict[str, str | None] = {}

    for key, candidates in {
        "family_id": ["FID", "FAMILY_ID", "FAMILYID", "FAMILY", "#FAMILY_ID"],
        "sample_id": [
            "IID",
            "INDIVIDUAL_ID",
            "SAMPLE_ID",
            "SAMPLEID",
            "SAMPLE",
            "INDIVIDUAL",
            "ID",
        ],
        "father_id": [
            "PAT",
            "FATHER_ID",
            "FATHERID",
            "FATHER",
            "PATERNAL_ID",
            "FATHERBARCODE",
        ],
        "mother_id": [
            "MAT",
            "MOTHER_ID",
            "MOTHERID",
            "MOTHER",
            "MATERNAL_ID",
            "MOTHERBARCODE",
        ],
        "sex": ["SEX", "GENDER"],
        "phenotype": ["PHENOTYPE", "AFFECTED", "STATUS", "AFFECTION"],
    }.items():
        for name in candidates:
            if name in columns:
                mapping[key] = columns[name]
                break

    # Fallback: any column starting with PHENO_ is a phenotype column
    if "phenotype" not in mapping:
        for upper_name, original_name in columns.items():
            if upper_name.startswith("PHENO_"):
                mapping["phenotype"] = original_name
                break

    return mapping


def load_family_pedigree(pedigree_path: Path) -> list[dict]:
    """Parse a per-family pedigree TSV into a list of member dicts.

    Returns the same format as Cohort.get_family_members():
        [{"Sample ID": ..., "Father": ..., "Mother": ..., "Sex": ..., "Phenotype": ...}]
    """
    df = pl.read_csv(pedigree_path, separator="\t", infer_schema_length=0)
    col_mapping = _identify_columns(df)

    sample_col = col_mapping.get("sample_id")
    father_col = col_mapping.get("father_id")
    mother_col = col_mapping.get("mother_id")
    sex_col = col_mapping.get("sex")
    phenotype_col = col_mapping.get("phenotype")

    if sample_col is None:
        return []

    _MISSING = {None, "", "0", "-9"}
    members = []
    for row in df.iter_rows(named=True):
        sample_id = str(row[sample_col])
        father = str(row[father_col]) if father_col and row.get(father_col) else "-"
        mother = str(row[mother_col]) if mother_col and row.get(mother_col) else "-"
        sex = str(row[sex_col]) if sex_col and row.get(sex_col) else "-"
        phenotype = (
            str(row[phenotype_col]) if phenotype_col and row.get(phenotype_col) else "-"
        )

        if father in _MISSING:
            father = "-"
        if mother in _MISSING:
            mother = "-"

        members.append(
            {
                "Sample ID": sample_id,
                "Father": father,
                "Mother": mother,
                "Sex": sex,
                "Phenotype": phenotype,
            }
        )
    return members


def load_family_object(pedigree_path: Path, family_id: str) -> Family:
    """Parse a per-family pedigree TSV into a Family model object.

    Returns a Family with Sample children, supporting .num_samples
    and .num_founders properties.
    """
    df = pl.read_csv(pedigree_path, separator="\t", infer_schema_length=0)
    col_mapping = _identify_columns(df)

    sample_col = col_mapping.get("sample_id")
    father_col = col_mapping.get("father_id")
    mother_col = col_mapping.get("mother_id")
    sex_col = col_mapping.get("sex")
    phenotype_col = col_mapping.get("phenotype")

    family = Family(family_id=family_id)

    if sample_col is None:
        return family

    _MISSING = {"", "0", "-9"}

    for row in df.iter_rows(named=True):
        sample_id = str(row[sample_col])

        def _val(col: str | None, treat_missing: bool = False) -> str | None:
            if col is None:
                return None
            val = row.get(col)
            if val is None or val == "":
                return None
            if treat_missing and str(val) in _MISSING:
                return None
            return str(val)

        sample = Sample(
            sample_id=sample_id,
            family_id=family_id,
            father_id=_val(father_col, treat_missing=True),
            mother_id=_val(mother_col, treat_missing=True),
            sex=_val(sex_col),
            phenotype=_val(phenotype_col),
        )
        family.samples.append(sample)

    return family
