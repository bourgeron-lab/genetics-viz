"""Tab-separated file reading helpers."""

import gzip
from pathlib import Path
from typing import Any, Iterable, Optional

import polars as pl


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
