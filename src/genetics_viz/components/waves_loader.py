"""Shared waves validation data loading utilities."""

import csv
from pathlib import Path
from typing import Dict, List

from genetics_viz.utils.tsv import append_tsv_row

WAVES_HEADER = ["sample", "user", "timestamp", "wave"]


def load_waves_validations(validation_file_path: Path) -> Dict[str, List[int]]:
    """Load waves validation data from validations/waves.tsv.

    Args:
        validation_file_path: Path to the validations/waves.tsv file

    Returns:
        Dictionary mapping sample_id to list of wave scores
    """
    waves_map: Dict[str, List[int]] = {}

    if not validation_file_path.exists():
        return waves_map

    with open(validation_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            sample_id = row.get("sample")
            wave_str = row.get("wave")

            if sample_id and wave_str:
                try:
                    wave = int(wave_str)
                    if sample_id not in waves_map:
                        waves_map[sample_id] = []
                    waves_map[sample_id].append(wave)
                except ValueError:
                    pass

    return waves_map


def load_waves_validations_full(validation_file_path: Path) -> Dict[str, List[Dict]]:
    """Load full waves validation records from validations/waves.tsv.

    Args:
        validation_file_path: Path to the validations/waves.tsv file

    Returns:
        Dictionary mapping sample_id to list of validation records with keys:
        - wave: int score
        - user: str username
        - timestamp: str ISO timestamp
    """
    waves_map: Dict[str, List[Dict]] = {}

    if not validation_file_path.exists():
        return waves_map

    with open(validation_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            sample_id = row.get("sample")
            wave_str = row.get("wave")
            user = row.get("user", "")
            timestamp = row.get("timestamp", "")

            if sample_id and wave_str:
                try:
                    wave = int(wave_str)
                    if sample_id not in waves_map:
                        waves_map[sample_id] = []
                    waves_map[sample_id].append(
                        {
                            "wave": wave,
                            "user": user,
                            "timestamp": timestamp,
                        }
                    )
                except ValueError:
                    pass

    return waves_map


def save_wave_validation(
    validation_file_path: Path,
    sample_id: str,
    user: str,
    wave: int,
    timestamp: str,
) -> None:
    """Save a wave validation to the validations/waves.tsv file.

    Args:
        validation_file_path: Path to the validations/waves.tsv file
        sample_id: Sample ID
        user: User name
        wave: Wave score (0, 1, 2, or 3)
        timestamp: ISO timestamp string
    """
    append_tsv_row(
        validation_file_path, WAVES_HEADER, [sample_id, user, timestamp, wave]
    )


def get_wave_category(waves: List[int] | None) -> str:
    """Get the category for a sample based on its wave scores.

    Args:
        waves: List of wave scores for the sample, or None if no validations

    Returns:
        Category string: "TODO", "good", "low wave", "medium wave", or "high wave"
    """
    if not waves:
        return "TODO"

    # Return the highest (worst) wave score
    max_wave = max(waves)
    return {
        0: "good",
        1: "low wave",
        2: "medium wave",
        3: "high wave",
    }.get(max_wave, "TODO")


def get_wave_color(category: str) -> str:
    """Get the color for a wave category.

    Args:
        category: Wave category

    Returns:
        CSS color string
    """
    return {
        "good": "green",
        "low wave": "yellow",
        "medium wave": "orange",
        "high wave": "red",
        "TODO": "gray",
    }.get(category, "gray")


def get_wave_score_color(wave: int) -> str:
    """Get the color for a wave score.

    Args:
        wave: Wave score (0, 1, 2, or 3)

    Returns:
        CSS color string
    """
    return {
        0: "green",
        1: "yellow",
        2: "orange",
        3: "red",
    }.get(wave, "gray")
