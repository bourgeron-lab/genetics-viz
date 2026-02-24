"""Validation badge configuration from YAML."""

from pathlib import Path
from typing import Any, Dict, List, Tuple

import yaml


def _load_config() -> Dict[str, Any]:
    config_path = Path(__file__).parent.parent / "config" / "validation_badges.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f) or {}


_CFG: Dict[str, Any] = _load_config()

INHERITANCE_CFG: Dict[str, Dict[str, str]] = _CFG.get("inheritance", {})
STATUS_CFG: Dict[str, Dict[str, str]] = _CFG.get("status", {})
DEFAULT_CFG: Dict[str, str] = _CFG.get("default", {})


def build_validation_badge(
    status: str,
    inheritance: str,
    validations: List[Tuple[str, str, str, str, str, str, str]],
) -> Dict[str, Any]:
    """Build badge and detail data for a validation cell.

    Parameters
    ----------
    status : aggregated validation status (e.g. "present", "absent")
    inheritance : aggregated inheritance (e.g. "de novo", "")
    validations : raw non-ignored validation tuples from validation_map

    Returns
    -------
    dict with ``badge`` (rendering info) and ``details`` (tooltip rows)
    """
    s_cfg = STATUS_CFG.get(status, {})
    i_cfg = INHERITANCE_CFG.get(inheritance, {})

    badge: Dict[str, str] = {
        "symbol": s_cfg.get("symbol", ""),
        "sc": s_cfg.get("color", "#9e9e9e"),
        "bg": i_cfg.get("bg", DEFAULT_CFG.get("bg", "#ffffff")),
        "text": i_cfg.get("text", ""),
        "tc": i_cfg.get("text_color", DEFAULT_CFG.get("text_color", "#374151")),
    }
    label = s_cfg.get("label", "")
    if label:
        badge["label"] = label

    details: List[Dict[str, str]] = []
    for v in validations:
        d_cfg = STATUS_CFG.get(v[0], {})
        # Extract date-only from timestamp (drop time portion)
        ts = v[6] if len(v) > 6 else ""
        date_only = ts.split(" ")[0].split("T")[0] if ts else ""
        details.append(
            {
                "s": v[0],
                "sy": d_cfg.get("symbol", ""),
                "sc": d_cfg.get("color", "#9e9e9e"),
                "i": v[1],
                "c": v[2],
                "t": date_only,
                "u": v[7] if len(v) > 7 else "",
            }
        )

    return {"badge": badge, "details": details}
