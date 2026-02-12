"""Utility functions for calculating gradient colors for continuous scores."""

import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml


def _load_score_configs() -> Dict:
    """Load continuous score configurations from YAML."""
    config_path = Path(__file__).parent.parent / "config" / "continuous_scores.yaml"
    with open(config_path, "r") as f:
        data = yaml.safe_load(f)
    return data.get("scores", {})


SCORE_CONFIGS = _load_score_configs()


def get_score_columns() -> List[str]:
    """Get list of all configured score column names."""
    return list(SCORE_CONFIGS.keys())


def reload_score_configs():
    """Reload score configurations from YAML file."""
    global SCORE_CONFIGS
    SCORE_CONFIGS = _load_score_configs()


def hex_to_rgb(hex_color: str) -> Tuple[int, int, int]:
    """Convert hex color to RGB tuple."""
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))


def rgb_to_hex(r: int, g: int, b: int) -> str:
    """Convert RGB tuple to hex color."""
    return f"#{r:02x}{g:02x}{b:02x}"


def interpolate_color(color1: str, color2: str, ratio: float) -> str:
    """Interpolate between two hex colors.

    Args:
        color1: Start color (hex)
        color2: End color (hex)
        ratio: Value between 0 (color1) and 1 (color2)

    Returns:
        Interpolated hex color
    """
    r1, g1, b1 = hex_to_rgb(color1)
    r2, g2, b2 = hex_to_rgb(color2)

    r = int(r1 + (r2 - r1) * ratio)
    g = int(g1 + (g2 - g1) * ratio)
    b = int(b1 + (b2 - b1) * ratio)

    return rgb_to_hex(r, g, b)


def get_score_color(score_name: str, value: float) -> Optional[Dict[str, str]]:
    """Get gradient color for a continuous score value.

    Args:
        score_name: Name of the score (case-insensitive)
        value: Score value

    Returns:
        Dict with 'color' (hex) and 'label' (category), or None if score not configured
    """
    # Case-insensitive lookup
    score_name_lower = score_name.lower()
    config = None
    for key in SCORE_CONFIGS:
        if key.lower() == score_name_lower:
            config = SCORE_CONFIGS[key]
            break

    if not config:
        return None

    # Apply log transformation if configured
    if config.get("log", False) and value > 0:
        value = math.log10(value)

    points = config["points"]

    # Handle edge cases
    if value <= points[0]["threshold"]:
        return {
            "color": points[0]["color"],
            "label": points[0]["label"]
        }
    if value >= points[-1]["threshold"]:
        return {
            "color": points[-1]["color"],
            "label": points[-1]["label"]
        }

    # Find the two points that frame the value
    for i in range(len(points) - 1):
        if points[i]["threshold"] <= value <= points[i + 1]["threshold"]:
            lower = points[i]
            upper = points[i + 1]

            # Calculate interpolation ratio
            threshold_range = upper["threshold"] - lower["threshold"]
            if threshold_range == 0:
                ratio = 0
            else:
                ratio = (value - lower["threshold"]) / threshold_range

            # Interpolate color
            color = interpolate_color(lower["color"], upper["color"], ratio)

            # Use label of closer point
            label = lower["label"] if ratio < 0.5 else upper["label"]

            return {
                "color": color,
                "label": label
            }

    # Fallback (shouldn't reach here)
    return {
        "color": "#757575",  # Gray
        "label": "Unknown"
    }
