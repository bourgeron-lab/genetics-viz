"""TSV exports of what a page shows: cell formatting, filenames, zip bundles.

Every tab of the cohort and family pages offers a download of its data. The
tables they draw are not exported as they are: badge columns hold ``""`` (the
text is in ``<column>_badges``), and the row dicts carry fields hidden from
readers (evidence, observed relationships, parents' kinship). Each export is
built from an explicit column list instead, by a builder next to the data it
shapes (``ancestry``, ``roh``, ``str_eh``, ``sample_qc_export``), so a family
file and a cohort file of one kind have the same columns.
"""

import csv
import io
import re
import zipfile
from datetime import date
from typing import Any, Dict, Iterable, List, Optional, Sequence

import polars as pl

from .pedigree_labels import normalize_pheno, normalize_sex

#: Keys that are internal/display-only and never exported from a table's rows.
INTERNAL_KEYS = {
    "_source_type",
    "_cohort_name",
    "_is_snv",
    "_original_locus",
    "_curated_tooltip",
    "_wisecondorx_tooltip",
    "IsCurated",
    "ConsequenceBadges",
    "ClinVarBadges",
    "GeneBadges",
    "VEP_Gene_badges",
    "MoiBadges",
    "_moi_rank",
    "_genetrek_url",
    "Diagnostic_badge",
}


def is_export_column(key: str) -> bool:
    """Whether a table column holds data rather than display state."""
    if key in INTERNAL_KEYS or key == "actions":
        return False
    if key.endswith("_badge") or key.endswith("_badges"):
        return False
    return True


def cell(value: Any) -> Any:
    """One TSV cell: empty for missing, a list joined, anything else as is."""
    if value is None:
        return ""
    if isinstance(value, float) and value != value:
        return ""
    if isinstance(value, (list, tuple, set)):
        return ", ".join(str(v) for v in value if v is not None)
    return value


def rows_to_tsv(
    rows: Iterable[Dict[str, Any]],
    columns: Sequence[str],
    headers: Optional[Sequence[str]] = None,
) -> str:
    """TSV text of ``rows``, only ``columns``, under ``headers`` (default: the ids)."""
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(list(headers or columns))
    for row in rows:
        writer.writerow([cell(row.get(col)) for col in columns])
    return output.getvalue()


def frame_to_tsv(df: pl.DataFrame) -> str:
    """TSV text of a frame; nulls as empty cells, lists joined."""
    return rows_to_tsv(df.to_dicts(), df.columns)


_UNSAFE = re.compile(r"[^A-Za-z0-9._+-]+")


def _safe(part: str) -> str:
    return _UNSAFE.sub("_", str(part)).strip("_")


def export_filename(
    prefix: str,
    tab: str,
    table: Optional[str] = None,
    tag: Optional[str] = None,
    ext: str = "tsv",
) -> str:
    """``<prefix>.<tab>[.<table>][.<tag>].<YYYY-MM-DD>.<ext>``."""
    parts = [prefix, tab, table, tag, date.today().isoformat()]
    return ".".join(_safe(p) for p in parts if p) + f".{ext}"


def zip_texts(files: Dict[str, str]) -> bytes:
    """A zip archive of named text files, UTF-8."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, text in files.items():
            archive.writestr(name, text.encode("utf-8"))
    return buffer.getvalue()


def person(fid: Any, iid: Any, sex: Any, phenotype: Any) -> Dict[str, Any]:
    """The identity columns every per-sample export starts with."""
    return {
        "FID": fid or "",
        "IID": iid or "",
        "sex": normalize_sex(None if sex == "-" else sex),
        "phenotype": normalize_pheno(None if phenotype == "-" else phenotype),
    }


PERSON_COLUMNS = ["FID", "IID", "sex", "phenotype"]


def people_from_individuals(individuals: Iterable[Dict[str, Any]]) -> List[Dict]:
    """The cohort table's filtered rows as :func:`person` dicts, each once."""
    seen: Dict[tuple, Dict[str, Any]] = {}
    for ind in individuals:
        if not ind.get("Sample ID"):
            continue
        key = (ind.get("FID") or "", ind["Sample ID"])
        if key not in seen:
            seen[key] = person(key[0], key[1], ind.get("Sex"), ind.get("Phenotype"))
    return list(seen.values())


def people_from_family(
    store: Any, cohort_name: str, family_id: str, members: Iterable[str]
) -> List[Dict]:
    """The selected members of a family as :func:`person` dicts, pedigree order.

    A member the pedigree does not list (or a family opened outside its
    cohort) keeps an unknown sex and phenotype.
    """
    wanted = list(dict.fromkeys(members))
    cohort = store.cohorts.get(cohort_name) if cohort_name else None
    family = cohort.families.get(family_id) if cohort is not None else None
    known = {s.sample_id: s for s in family.samples} if family is not None else {}
    out = []
    for sample_id in wanted:
        s = known.get(sample_id)
        out.append(
            person(
                family_id,
                sample_id,
                s.sex if s else None,
                s.phenotype if s else None,
            )
        )
    return out


def people_frame(people: List[Dict[str, Any]]) -> pl.DataFrame:
    """:func:`person` dicts as a frame of string columns."""
    return pl.DataFrame(
        {c: [p[c] for p in people] for c in PERSON_COLUMNS},
        schema={c: pl.String for c in PERSON_COLUMNS},
    )
