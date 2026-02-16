"""View preset configuration for column visibility.

Previously lived in wombat_tab.py but is shared by search.py and
the header reload mechanism.
"""

from pathlib import Path
from typing import Any, Dict, List

import yaml

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "view_presets.yaml"


def _load_view_presets() -> List[Dict[str, Any]]:
    """Load view presets from YAML config file."""
    with open(_CONFIG_PATH, "r") as f:
        data = yaml.safe_load(f)
    return data.get("presets", [])


VIEW_PRESETS: List[Dict[str, Any]] = _load_view_presets()


def select_preset_for_config(config_name: str, presets: List[Dict]) -> Dict:
    """Select the first preset whose keywords contain the config_name,
    or return the first preset if none match."""
    config_lower = config_name.lower()

    for preset in presets:
        keywords = preset.get("keywords", [])
        if any(keyword.lower() in config_lower for keyword in keywords):
            return preset

    return presets[0] if presets else {"name": "Default", "columns": []}


def reload_view_presets() -> None:
    """Reload view presets from YAML file."""
    global VIEW_PRESETS
    VIEW_PRESETS = _load_view_presets()
