"""Shared column display names and groups from YAML config."""

from pathlib import Path
from typing import Dict, List

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
