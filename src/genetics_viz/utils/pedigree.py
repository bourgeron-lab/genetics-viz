"""Standalone pedigree file parser for per-family pedigree TSV files.

Parses pedigree files independently of the Cohort model. A file with a
header is read by column name, which covers the per-family format
(FatherBarcode/MotherBarcode/sample_id columns) as well as the standard
cohort names (PAT/MAT/IID). A file without one - as the ghfc-ngs pipeline
writes ``<FID>.pedigree.tsv`` - is read by position: FID, IID, father,
mother, sex, phenotype. Header detection is the cohort pedigree's own.
"""

from pathlib import Path

import polars as pl

from genetics_viz.models import (
    MISSING_PARENT,
    Family,
    Sample,
    identify_pedigree_columns,
    pedigree_has_header,
    read_positional_pedigree,
)


def _read_family_samples(pedigree_path: Path, family_id: str) -> list[Sample]:
    """Read a per-family pedigree into Samples, headed or headerless.

    Parent codes in ``MISSING_PARENT`` become ``None``. Sex and phenotype are
    kept as written, as in the cohort model; ``utils/pedigree_labels``
    normalises them for display.
    """
    if pedigree_has_header(pedigree_path):
        df = pl.read_csv(pedigree_path, separator="\t", infer_schema_length=0)
    else:
        df = read_positional_pedigree(pedigree_path, skip_header=False)
    col_mapping = identify_pedigree_columns(df.columns)

    sample_col = col_mapping.get("sample_id")
    if sample_col is None:
        return []
    father_col = col_mapping.get("father_id")
    mother_col = col_mapping.get("mother_id")
    sex_col = col_mapping.get("sex")
    phenotype_col = col_mapping.get("phenotype")

    samples = []
    for row in df.iter_rows(named=True):

        def _val(col: str | None, treat_missing: bool = False) -> str | None:
            if col is None:
                return None
            val = row.get(col)
            if val is None or val == "":
                return None
            if treat_missing and val in MISSING_PARENT:
                return None
            return str(val)

        samples.append(
            Sample(
                sample_id=str(row[sample_col]),
                family_id=family_id,
                father_id=_val(father_col, treat_missing=True),
                mother_id=_val(mother_col, treat_missing=True),
                sex=_val(sex_col),
                phenotype=_val(phenotype_col),
            )
        )
    return samples


def load_family_pedigree(pedigree_path: Path) -> list[dict]:
    """Parse a per-family pedigree TSV into a list of member dicts.

    Returns the same format as Cohort.get_family_members():
        [{"Sample ID": ..., "Father": ..., "Mother": ..., "Sex": ..., "Phenotype": ...}]
    """
    return [
        {
            "Sample ID": sample.sample_id,
            "Father": sample.father_id or "-",
            "Mother": sample.mother_id or "-",
            "Sex": sample.sex or "-",
            "Phenotype": sample.phenotype or "-",
        }
        for sample in _read_family_samples(pedigree_path, family_id="")
    ]


def load_family_object(pedigree_path: Path, family_id: str) -> Family:
    """Parse a per-family pedigree TSV into a Family model object.

    Returns a Family with Sample children, supporting .num_samples
    and .num_founders properties.
    """
    return Family(
        family_id=family_id,
        samples=_read_family_samples(pedigree_path, family_id),
    )
