"""Tab-separated file reading helpers."""

from pathlib import Path
from typing import Any

import polars as pl


def read_tsv_or_none(path: Path, **kwargs: Any) -> pl.DataFrame | None:
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

    Keyword arguments are passed through to :func:`polars.read_csv`;
    ``separator`` defaults to a tab.
    """
    kwargs.setdefault("separator", "\t")
    try:
        return pl.read_csv(path, **kwargs)
    except pl.exceptions.NoDataError:
        return None
