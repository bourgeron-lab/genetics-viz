"""Shared column display names and groups from YAML config."""

import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import polars as pl
import yaml


def _load_column_names() -> Dict[str, Dict[str, str]]:
    """Load column name/group mappings from YAML."""
    config_path = Path(__file__).parent.parent / "config" / "column_names.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f) or {}


COLUMN_NAMES: Dict[str, Dict[str, str]] = _load_column_names()


def get_display_label(col: str) -> str:
    """Return the display name for a column. Falls back to raw col ID."""
    entry = COLUMN_NAMES.get(col)
    if entry:
        return entry.get("name", col)
    return col


def get_column_group(col: str) -> str:
    """Return the group name for a column. Empty string means ungrouped."""
    entry = COLUMN_NAMES.get(col)
    if entry:
        return entry.get("group", "")
    return ""


def get_dropped_columns() -> set:
    """Return the set of column IDs marked with ``drop: true`` in YAML."""
    return {
        col for col, entry in COLUMN_NAMES.items()
        if entry.get("drop") is True
    }


_TYPE_MAP: Dict[str, pl.DataType] = {
    "string": pl.Utf8,
    "int": pl.Int64,
    "float": pl.Float64,
    "bool": pl.Boolean,
}


def get_schema_overrides() -> Dict[str, pl.DataType]:
    """Return Polars schema_overrides for columns with an explicit type in YAML."""
    overrides: Dict[str, pl.DataType] = {}
    for col, entry in COLUMN_NAMES.items():
        col_type = entry.get("type")
        if col_type and col_type in _TYPE_MAP:
            overrides[col] = _TYPE_MAP[col_type]
    return overrides


def get_column_sorting(col: str) -> str:
    """Return the sorting strategy for a column. Empty string means default."""
    entry = COLUMN_NAMES.get(col)
    if entry:
        return entry.get("sorting", "")
    return ""


def get_column_min_width(col: str) -> int | None:
    """Return the min_width (px) for a column, or None."""
    entry = COLUMN_NAMES.get(col)
    if entry:
        v = entry.get("min_width")
        return int(v) if v is not None else None
    return None


def get_column_max_width(col: str) -> int | None:
    """Return the max_width (px) for a column, or None."""
    entry = COLUMN_NAMES.get(col)
    if entry:
        v = entry.get("max_width")
        return int(v) if v is not None else None
    return None


def apply_width_constraints(col_def: Dict[str, Any], col: str) -> None:
    """Add minWidth / maxWidth from YAML to a column definition dict."""
    min_w = get_column_min_width(col)
    if min_w is not None:
        col_def["minWidth"] = min_w
    max_w = get_column_max_width(col)
    if max_w is not None:
        col_def["maxWidth"] = max_w


# ---- Genomic sort key ----

_GENOMIC_RE = re.compile(r"^(?:chr)?(.+?):(\d+)(.*)", re.IGNORECASE)
_CHROM_ORDER: Dict[str, int] = {str(i): i for i in range(1, 23)}
_CHROM_ORDER.update({"X": 23, "Y": 24, "M": 25, "MT": 25})


def genomic_sort_key(value: Any) -> Tuple[int, int, str]:
    """Return a sort key tuple for a genomic coordinate string.

    Handles formats like ``chr1:123456:A:C``, ``10:234567:G:T``,
    ``chrX:123-567``.  Chromosomes are ordered 1-22, X, Y, M/MT.
    """
    if value is None or value == "":
        return (999, 0, "")
    s = str(value)
    m = _GENOMIC_RE.match(s)
    if not m:
        return (999, 0, s)
    chrom = m.group(1).upper()
    pos = int(m.group(2))
    rest = m.group(3)
    rank = _CHROM_ORDER.get(chrom, 26)
    return (rank, pos, rest)


def reorder_columns_by_group(columns: List[str]) -> List[str]:
    """Reorder columns so that all columns sharing a group are consecutive.

    When the first column of a group is encountered, all other columns
    belonging to that same group are pulled forward to follow it immediately
    (preserving their relative order).  Ungrouped columns stay in place.
    """
    result: List[str] = []
    added: set = set()
    seen_groups: set = set()

    for col in columns:
        if col in added:
            continue
        group = get_column_group(col)
        if not group:
            result.append(col)
            added.add(col)
        elif group not in seen_groups:
            seen_groups.add(group)
            for c in columns:
                if c not in added and get_column_group(c) == group:
                    result.append(c)
                    added.add(c)
        # else: already added when the group was first seen

    return result
