"""Diagnostic verdicts and their badges, from ``config/diagnostic_badges.yaml``.

The verdicts are the ACMG/AMP five tiers, strongest first (:data:`VERDICTS`),
plus :data:`CONFLICTING`, computed when curators disagree. Pages never spell
them: they use the helpers here, so a tier is added in the YAML alone.
"""

from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import yaml


def _load_config() -> Dict[str, Any]:
    config_path = Path(__file__).parent.parent / "config" / "diagnostic_badges.yaml"
    with open(config_path, "r") as f:
        return yaml.safe_load(f) or {}


_CFG: Dict[str, Any] = _load_config()

STATUS_CFG: Dict[str, Dict[str, str]] = _CFG.get("status", {})
DEFAULT_CFG: Dict[str, str] = _CFG.get("default", {})

#: Verdict computed when curators disagree on one variant and sample.
CONFLICTING = "conflicting"
PATHOGENIC = "pathogenic"
LIKELY_PATHOGENIC = "likely pathogenic"
UNCERTAIN = "uncertain"
LIKELY_BENIGN = "likely benign"
BENIGN = "benign"

#: The verdicts a curator records, strongest first.
VERDICTS: Tuple[str, ...] = tuple(k for k in STATUS_CFG if k != CONFLICTING)
#: The same, then :data:`CONFLICTING` (filters, charts).
ALL_STATUSES: Tuple[str, ...] = (*VERDICTS, CONFLICTING)


def reload_diagnostic_config() -> None:
    """Hot-reload diagnostic badge config."""
    global _CFG, STATUS_CFG, DEFAULT_CFG
    _CFG = _load_config()
    STATUS_CFG = _CFG.get("status", {})
    DEFAULT_CFG = _CFG.get("default", {})


def verdict_label(verdict: str) -> str:
    """``likely pathogenic`` → ``Likely pathogenic``; an unknown value as is."""
    return STATUS_CFG.get(verdict, {}).get("label", verdict)


def verdict_symbol(verdict: str) -> str:
    """The badge text: ``P``, ``LP``, ``VUS``, ``LB``, ``B`` or ``!``."""
    return STATUS_CFG.get(verdict, {}).get("symbol", "")


def verdict_icon(verdict: str) -> Tuple[str, str]:
    """``(Material icon, Quasar colour)``; ``("", "grey")`` for an unknown value."""
    cfg = STATUS_CFG.get(verdict)
    if not cfg:
        return "", "grey"
    return cfg.get("icon", "help_outline"), cfg.get("quasar", "grey")


def verdict_color(verdict: str, default: str = "#9e9e9e") -> str:
    """The verdict's hex colour (charts, table cells)."""
    return STATUS_CFG.get(verdict, {}).get("hex", default)


def verdict_rank(verdict: str) -> int:
    """Higher for a stronger pathogenic verdict: P > LP > VUS > LB > B > other.

    Used where one verdict must stand for several records (a sample's, a
    family's). :data:`CONFLICTING` and unknown values rank 0.
    """
    if verdict not in VERDICTS:
        return 0
    return len(VERDICTS) - VERDICTS.index(verdict)


def strongest(verdicts: Iterable[str]) -> str:
    """The highest-ranked of ``verdicts`` (``""`` when none is a verdict)."""
    best = max(verdicts, key=verdict_rank, default="")
    return best if verdict_rank(best) else ""


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
        # A known verdict is its symbol alone; an unknown one is shown raw.
        "text": "" if s_cfg else (status[:3].upper() if status else ""),
    }

    details: List[Dict[str, str]] = []
    for d in diagnostics:
        d_cfg = STATUS_CFG.get(d[0], {})
        ts = d[2] if len(d) > 2 else ""
        date_only = ts.split(" ")[0].split("T")[0] if ts else ""
        details.append(
            {
                "s": verdict_label(d[0]),
                "sy": d_cfg.get("symbol", ""),
                "sc": d_cfg.get("color", "#9e9e9e"),
                "u": d[1] if len(d) > 1 else "",
                "t": date_only,
                "c": d[3] if len(d) > 3 else "",
            }
        )

    return {"badge": badge, "details": details}
