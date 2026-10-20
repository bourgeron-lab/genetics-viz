"""Tab-separated file reading helpers."""

import fcntl
import gzip
import logging
import re
from pathlib import Path
from typing import Any, Iterable, List, Optional, Tuple

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


def append_tsv_row(path: Path, header: List[str], fields: Iterable[Any]) -> None:
    """Append one record to a curation file, writing ``header`` first if it is new.

    Every field goes through :func:`clean_tsv_field`. The write holds an
    exclusive ``flock``, and the header check is made under it, at the end of
    the file as it is now: two first saves must not both write the header.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            f.seek(0, 2)
            if f.tell() == 0:
                f.write("\t".join(header) + "\n")
            f.write("\t".join(clean_tsv_field(str(v)) for v in fields) + "\n")
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


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
