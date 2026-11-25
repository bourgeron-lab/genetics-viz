"""Icon definitions and utilities for genetics-viz."""

from typing import Dict, Tuple

from genetics_viz.utils.diagnostic_badges import ALL_STATUSES, verdict_icon
from genetics_viz.utils.validation_status import STATUS_ICONS

# Map validation status to (icon_name, color)
VALIDATION_STATUS_ICONS: Dict[str, Tuple[str, str]] = STATUS_ICONS


def get_validation_icon(status: str) -> Tuple[str, str]:
    """Get the icon name and color for a validation status.

    Args:
        status: The validation status string

    Returns:
        Tuple of (icon_name, color)
    """
    return VALIDATION_STATUS_ICONS.get(status, ("", "grey"))


# Map diagnostic status to (icon_name, color), from config/diagnostic_badges.yaml
DIAGNOSTIC_STATUS_ICONS: Dict[str, Tuple[str, str]] = {
    status: verdict_icon(status) for status in ALL_STATUSES
}


def get_diagnostic_icon(status: str) -> Tuple[str, str]:
    """Get the icon name and color for a diagnostic status.

    Args:
        status: The diagnostic status string

    Returns:
        Tuple of (icon_name, color)
    """
    return DIAGNOSTIC_STATUS_ICONS.get(status, ("", "grey"))
