"""Shared waves validation data loading utilities."""

import csv
import fcntl
from pathlib import Path
from typing import Dict, List


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
        reader = csv.DictReader(f, delimiter="\t")
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
        reader = csv.DictReader(f, delimiter="\t")
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
    # Ensure the validations directory exists
    validation_file_path.parent.mkdir(parents=True, exist_ok=True)

    # Check if file exists and has header
    file_exists = validation_file_path.exists()
    needs_header = not file_exists

    if file_exists:
        # Check if file is empty or has no header
        with open(validation_file_path, "r") as f:
            first_line = f.readline().strip()
            if not first_line or first_line != "sample\tuser\ttimestamp\twave":
                needs_header = True

    # Open file for appending with exclusive lock
    with open(validation_file_path, "a") as f:
        # Acquire exclusive lock
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)

        try:
            # Write header if needed
            if needs_header:
                f.write("sample\tuser\ttimestamp\twave\n")

            # Write the validation record
            f.write(f"{sample_id}\t{user}\t{timestamp}\t{wave}\n")
        finally:
            # Release lock
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


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
