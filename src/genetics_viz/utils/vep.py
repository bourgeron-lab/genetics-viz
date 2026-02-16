"""VEP (Variant Effect Predictor) consequence utilities.

Centralised loading and lookup functions for VEP consequence data
(impact levels, colours, priority ordering).  Previously duplicated
across wombat_tab.py, search.py and validation/file.py.
"""

from pathlib import Path
from typing import Dict

import yaml

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "vep_consequences.yaml"


def _load_vep_consequences() -> Dict[str, tuple]:
    """Load VEP consequences from YAML config file."""
    with open(_CONFIG_PATH, "r") as f:
        data = yaml.safe_load(f)
    return {term: (info["impact"], info["color"]) for term, info in data.items()}


VEP_CONSEQUENCES: Dict[str, tuple] = _load_vep_consequences()

# Lower number = higher priority (insertion-order from YAML)
VEP_CONSEQUENCE_PRIORITY: Dict[str, int] = {
    term: idx for idx, term in enumerate(VEP_CONSEQUENCES.keys())
}


def get_consequence_color(consequence: str) -> str:
    """Get color for a consequence term."""
    return VEP_CONSEQUENCES.get(consequence, ("MODIFIER", "#636363"))[1]


def get_consequence_impact(consequence: str) -> str:
    """Get impact level for a consequence term."""
    return VEP_CONSEQUENCES.get(consequence, ("MODIFIER", "#636363"))[0]


def get_consequence_priority(consequence: str) -> int:
    """Get priority for a consequence term (lower = higher priority)."""
    return VEP_CONSEQUENCE_PRIORITY.get(consequence, 9999)


def get_highest_priority_consequence(consequence_str: str) -> int:
    """Get the highest priority (lowest number) from a comma/ampersand-separated consequence string."""
    if not consequence_str:
        return 9999

    consequences = []
    for part in str(consequence_str).split(","):
        for cons in part.split("&"):
            cons = cons.strip()
            if cons:
                consequences.append(cons)

    if not consequences:
        return 9999

    return min(get_consequence_priority(cons) for cons in consequences)


def get_highest_consequence_term(consequence_str: str) -> str:
    """Get the highest-priority consequence term from a comma/ampersand-separated string."""
    if not consequence_str:
        return "Unknown"
    consequences = []
    for part in str(consequence_str).split(","):
        for cons in part.split("&"):
            cons = cons.strip()
            if cons:
                consequences.append(cons)
    if not consequences:
        return "Unknown"
    return min(consequences, key=lambda c: VEP_CONSEQUENCE_PRIORITY.get(c, 9999))


def format_consequence_display(consequence: str) -> str:
    """Format consequence for display: remove _variant suffix and replace _ with space."""
    return consequence.replace("_variant", "").replace("_", " ")


def reload_vep_config() -> None:
    """Reload VEP consequences from YAML file."""
    global VEP_CONSEQUENCES, VEP_CONSEQUENCE_PRIORITY
    VEP_CONSEQUENCES = _load_vep_consequences()
    VEP_CONSEQUENCE_PRIORITY = {
        term: idx for idx, term in enumerate(VEP_CONSEQUENCES.keys())
    }
