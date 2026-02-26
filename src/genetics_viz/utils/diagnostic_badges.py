"""Diagnostic badge configuration from YAML."""

from pathlib import Path
from typing import Any, Dict, List, Tuple

import yaml


def _load_config() -> Dict[str, Any]:
    config_path = Path(__file__).parent.parent / "config" / "diagnostic_badges.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f) or {}


_CFG: Dict[str, Any] = _load_config()

STATUS_CFG: Dict[str, Dict[str, str]] = _CFG.get("status", {})
DEFAULT_CFG: Dict[str, str] = _CFG.get("default", {})


def reload_diagnostic_config() -> None:
    """Hot-reload diagnostic badge config."""
    global _CFG, STATUS_CFG, DEFAULT_CFG
    _CFG = _load_config()
    STATUS_CFG = _CFG.get("status", {})
    DEFAULT_CFG = _CFG.get("default", {})


def build_diagnostic_badge(
    status: str,
    diagnostics: List[Tuple[str, str, str, str]],
) -> Dict[str, Any]:
    """Build badge and detail data for a diagnostic cell.

    Parameters
    ----------
    status : aggregated diagnostic status (e.g. "pathogenic", "conflicting")
    diagnostics : raw non-ignored diagnostic tuples
        Each tuple: (diagnostic_value, user, timestamp, comment)

    Returns
    -------
    dict with ``badge`` (rendering info) and ``details`` (tooltip rows)
    """
    s_cfg = STATUS_CFG.get(status, {})

    badge: Dict[str, str] = {
        "symbol": s_cfg.get("symbol", ""),
        "sc": s_cfg.get("color", "#9e9e9e"),
        "bg": s_cfg.get("bg", DEFAULT_CFG.get("bg", "#ffffff")),
        "tc": s_cfg.get("text_color", DEFAULT_CFG.get("text_color", "#374151")),
        "text": status[:3].upper() if status else "",
    }

    details: List[Dict[str, str]] = []
    for d in diagnostics:
        d_cfg = STATUS_CFG.get(d[0], {})
        ts = d[2] if len(d) > 2 else ""
        date_only = ts.split(" ")[0].split("T")[0] if ts else ""
        details.append(
            {
                "s": d[0],
                "sy": d_cfg.get("symbol", ""),
                "sc": d_cfg.get("color", "#9e9e9e"),
                "u": d[1] if len(d) > 1 else "",
                "t": date_only,
                "c": d[3] if len(d) > 3 else "",
            }
        )

    return {"badge": badge, "details": details}
