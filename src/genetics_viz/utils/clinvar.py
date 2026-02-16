"""ClinVar significance colour coding utilities.

Centralised loading and lookup functions for ClinVar significance
colours.  Previously duplicated across wombat_tab.py and search.py.
"""

from pathlib import Path
from typing import Dict

import yaml

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "clinvar_colors.yaml"


def _load_clinvar_colors() -> Dict[str, str]:
    """Load ClinVar colors from YAML config file."""
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f)


CLINVAR_COLORS: Dict[str, str] = _load_clinvar_colors()


def get_clinvar_color(significance: str) -> str:
    """Get color for a ClinVar significance term (case-insensitive)."""
    if not significance:
        return "#757575"
    sig_lower = significance.lower()
    for key, color in CLINVAR_COLORS.items():
        if key.lower() == sig_lower:
            return color
    return "#757575"


def format_clinvar_display(significance: str) -> str:
    """Format ClinVar significance for display: replace _ with space."""
    return significance.replace("_", " ")


def reload_clinvar_config() -> None:
    """Reload ClinVar colors from YAML file."""
    global CLINVAR_COLORS
    CLINVAR_COLORS = _load_clinvar_colors()
