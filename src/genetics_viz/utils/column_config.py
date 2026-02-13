"""Column display names and grouping configuration."""

from pathlib import Path
from typing import Dict, Optional, Tuple

import yaml


def _load_column_config() -> Dict[str, dict]:
    """Load column names config from YAML."""
    config_path = Path(__file__).parent.parent / "config" / "column_names.yaml"
    if not config_path.exists():
        return {}
    with open(config_path, "r") as f:
        data = yaml.safe_load(f)
    return data or {}


COLUMN_CONFIG = _load_column_config()


def get_column_display_name(col: str) -> str:
    """Get display name for a column from config, falling back to column name."""
    entry = COLUMN_CONFIG.get(col)
    if entry and entry.get("name"):
        return entry["name"]
    return col


def get_column_group(col: str) -> Optional[str]:
    """Get the group name for a column, or None if not grouped."""
    entry = COLUMN_CONFIG.get(col)
    if entry and entry.get("group"):
        return entry["group"]
    return None


def get_column_info(col: str) -> Tuple[str, Optional[str]]:
    """Get (display_name, group_name) for a column."""
    return get_column_display_name(col), get_column_group(col)
