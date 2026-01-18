"""Icon definitions and utilities for genetics-viz."""

from typing import Dict, Tuple

# Map validation status to (icon_name, color)
VALIDATION_STATUS_ICONS: Dict[str, Tuple[str, str]] = {
    "present": ("check_circle", "green"),
    "absent": ("cancel", "red"),
    "uncertain": ("help", "orange"),
    "conflicting": ("bolt", "amber-9"),
    "TODO": ("assignment", "grey"),
}


def get_validation_icon(status: str) -> Tuple[str, str]:
    """Get the icon name and color for a validation status.

    Args:
        status: The validation status string

    Returns:
        Tuple of (icon_name, color)
    """
    return VALIDATION_STATUS_ICONS.get(status, ("", "grey"))
