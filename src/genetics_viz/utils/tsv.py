"""Tab-separated file reading helpers."""

import csv
import fcntl
import gzip
import logging
import re
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

import polars as pl

logger = logging.getLogger(__name__)


def clean_tsv_field(text: str) -> str:
    """Make free text safe as one TSV field: no tabs, no line breaks.

    The curation files (notes, diagnostics, validations) are appended to line by
    line, so a tab or newline typed into a comment would split the record.
    """
    return (
        text.replace("\t", " ")
        .replace("\r\n", " ")
        .replace("\n", " ")
        .replace("\r", " ")
    )


def _tsv_line(fields: Iterable[Any]) -> str:
    return "\t".join(clean_tsv_field(str(v)) for v in fields) + "\n"


def append_tsv_row(path: Path, header: List[str], fields: Iterable[Any]) -> None:
    """Append one record to a curation file, writing ``header`` first if it is new.

    Every field goes through :func:`clean_tsv_field`. The write holds an
    exclusive ``flock``, and the header check is made under it, at the end of
    the file as it is now: two first saves must not both write the header.
    """
    append_tsv_rows(path, header, [fields])


def append_tsv_rows(
    path: Path, header: List[str], rows: Iterable[Iterable[Any]]
) -> None:
    """Append several records under one lock, as :func:`append_tsv_row`.

    With no rows, this only creates the file with its header.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            f.seek(0, 2)
            if f.tell() == 0:
                f.write("\t".join(header) + "\n")
            f.writelines(_tsv_line(fields) for fields in rows)
            # Before unlocking: buffered lines would reach the file only at
            # close, and the next writer would find it empty.
            f.flush()
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def rewrite_tsv_rows(
    path: Path, change: Callable[[List[Dict[str, str]]], bool]
) -> bool:
    """Edit a curation file in place: read, ``change`` the rows, write back.

    ``change`` edits the list of row dicts in place and returns whether it
    changed anything; the file is rewritten only then. Read and write happen
    under one exclusive ``flock``, so a save in between is neither lost nor
    truncated away. Returns ``change``'s answer, ``False`` for a missing file.

    Rows are read with ``QUOTE_NONE`` and written back unquoted, as
    :func:`append_tsv_row` writes them: a ``csv.DictWriter`` would quote a
    comment holding ``"``, and the next read would keep the quotes.
    """
    if not path.exists():
        return False
    with open(path, "r+") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
            rows = list(reader)
            fieldnames = list(reader.fieldnames or [])
            if not change(rows):
                return False
            f.seek(0)
            f.write("\t".join(fieldnames) + "\n")
            f.writelines(
                _tsv_line(row.get(name) or "" for name in fieldnames) for row in rows
            )
            f.truncate()
            f.flush()
            return True
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def set_ignore(
    path: Path,
    family_id: str,
    variant_key: str,
    sample: str,
    timestamp: str,
    ignore_value: str,
) -> bool:
    """Set ``Ignore`` on one validation or diagnostic record, in place.

    The record is the one of this family, variant, sample and timestamp: both
    kinds of file key their rows alike. Returns False when none matches.
    """

    def change(rows: List[Dict[str, str]]) -> bool:
        for row in rows:
            if (
                row.get("FID") == family_id
                and row.get("Variant") == variant_key
                and row.get("Sample") == sample
                and row.get("Timestamp") == timestamp
            ):
                row["Ignore"] = ignore_value
                return True
        return False

    return rewrite_tsv_rows(path, change)


def _header(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def read_tsv_or_none(
    path: Path, *, string_columns: Optional[Iterable[str]] = None, **kwargs: Any
) -> pl.DataFrame | None:
    """Read a tab-separated file, returning ``None`` when it holds no rows at all.

    Wombat writes a zero-byte file when a step produced nothing - a family
    without two sequenced parents has no de novo mutations to report, for
    instance - and :func:`polars.read_csv` raises ``NoDataError`` on such a
    file rather than returning an empty frame. That is an expected, ordinary
    outcome here, not a failure, so it comes back as ``None`` and callers can
    skip the file without wrapping every read in a try/except that would also
    swallow real errors.

    A file carrying only a header still returns an empty DataFrame, which a
    ``len(df) == 0`` check already covers, and a *missing* file still raises
    ``FileNotFoundError`` - that one is not expected.

    ``string_columns`` names identifier columns to read as strings when the
    file has them. Left to inference, a column of family IDs whose first rows
    happen to be numeric is typed as an integer, and the read then fails on the
    first ID that is not (``AU-ANG-FER-MPX3``).

    Other keyword arguments are passed through to :func:`polars.read_csv`;
    ``separator`` defaults to a tab.
    """
    kwargs.setdefault("separator", "\t")
    if string_columns:
        wanted = set(string_columns)
        overrides = {c: pl.Utf8 for c in _header(Path(path)) if c in wanted}
        if overrides:
            kwargs["schema_overrides"] = {
                **overrides,
                **kwargs.get("schema_overrides", {}),
            }
    try:
        return pl.read_csv(path, **kwargs)
    except pl.exceptions.NoDataError:
        return None


def find_labels(directory: Path, prefix: str, kind: str) -> List[str]:
    """Labels of the ``<prefix>.<label>.<kind>.tsv`` files in ``directory``.

    Newest first, then by name. The label is the whole middle of the name: a
    lazy capture would cut a dotted label short. Used where the pipeline names
    a step's files by a free label (an STR catalog, a sample QC sites file)
    rather than by a bundle version and filter tag.
    """
    if not directory.is_dir():
        return []
    pattern = re.compile(rf"^{re.escape(prefix)}\.(.+)\.{re.escape(kind)}\.tsv$")
    found: List[Tuple[float, str]] = []
    try:
        for entry in directory.iterdir():
            m = pattern.match(entry.name)
            if m:
                found.append((entry.stat().st_mtime, m.group(1)))
    except OSError:
        logger.warning("Cannot list %s", directory, exc_info=True)
        return []
    return [label for _, label in sorted(found, key=lambda t: (-t[0], t[1]))]
