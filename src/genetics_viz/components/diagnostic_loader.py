"""Shared diagnostic data loading utilities."""

import csv
import fcntl
from pathlib import Path
from typing import Any, Dict, List, Tuple

from genetics_viz.utils.diagnostic_badges import build_diagnostic_badge

# TSV header for diagnostics files
DIAGNOSTIC_HEADER = (
    "FID\tVariant\tGene\tImpact\tSample\tUser\tTimestamp\tComment\tIgnore\tDiagnostic\n"
)


def ensure_diagnostic_file(file_path: Path) -> None:
    """Create the diagnostics file with header if it doesn't exist."""
    file_path.parent.mkdir(parents=True, exist_ok=True)
    if not file_path.exists():
        with open(file_path, "w") as f:
            f.write(DIAGNOSTIC_HEADER)


def load_diagnostic_map(
    diagnostic_file_path: Path, family_id: str | None = None
) -> Dict[Tuple[str, str], List[Tuple[str, str, str, str, str]]]:
    """Load diagnostic data from snvs.tsv or svs.tsv into a lookup map.

    Args:
        diagnostic_file_path: Path to the diagnostics file
        family_id: Optional family ID to filter by

    Returns:
        Dictionary mapping (variant_key, sample_id) to list of
        (diagnostic_value, user, timestamp, comment, ignore)
    """
    diagnostic_map: Dict[Tuple[str, str], List[Tuple[str, str, str, str, str]]] = {}

    if not diagnostic_file_path.exists():
        return diagnostic_map

    with open(diagnostic_file_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            fid = row.get("FID")
            variant_key = row.get("Variant")
            sample_id = row.get("Sample")
            diagnostic = row.get("Diagnostic", "")
            user = row.get("User", "")
            timestamp = row.get("Timestamp", "")
            comment = row.get("Comment", "")
            ignore = row.get("Ignore", "0")

            if family_id is not None and fid != family_id:
                continue

            if variant_key and sample_id:
                map_key = (variant_key, sample_id)
                if map_key not in diagnostic_map:
                    diagnostic_map[map_key] = []
                diagnostic_map[map_key].append(
                    (diagnostic, user, timestamp, comment, ignore)
                )

    return diagnostic_map


def add_diagnostic_status_to_row(
    row: Dict[str, Any],
    diagnostic_map: Dict[Tuple[str, str], List[Tuple[str, str, str, str, str]]],
    variant_key: str,
    sample_id: str,
) -> None:
    """Add Diagnostic and Diagnostic_badge fields to a row.

    Args:
        row: The row dict to modify
        diagnostic_map: Mapping from (variant_key, sample_id) to diagnostics
        variant_key: The variant key
        sample_id: The sample ID
    """
    map_key = (variant_key, sample_id)

    if map_key in diagnostic_map:
        all_diagnostics = diagnostic_map[map_key]
        # Filter out ignored entries
        diagnostics = [d for d in all_diagnostics if d[4] != "1"]

        if not diagnostics:
            row["Diagnostic"] = ""
            row["Diagnostic_badge"] = None
            return

        diagnostic_values = [d[0] for d in diagnostics]
        unique_diagnostics = set(diagnostic_values)

        if len(unique_diagnostics) > 1:
            row["Diagnostic"] = "conflicting"
        else:
            row["Diagnostic"] = diagnostic_values[0]

        # Build badge with non-ignored entries for tooltip
        badge_data = [(d[0], d[1], d[2], d[3]) for d in diagnostics]
        row["Diagnostic_badge"] = build_diagnostic_badge(row["Diagnostic"], badge_data)
    else:
        row["Diagnostic"] = ""
        row["Diagnostic_badge"] = None


def save_diagnostic_entry(
    diagnostic_file: Path,
    family_id: str,
    variant_key: str,
    gene: str,
    impact: str,
    sample: str,
    user: str,
    timestamp: str,
    comment: str,
    diagnostic: str,
) -> None:
    """Append a diagnostic entry to the TSV file with file locking."""
    ensure_diagnostic_file(diagnostic_file)
    with open(diagnostic_file, "a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            f.write(
                f"{family_id}\t{variant_key}\t{gene}\t{impact}\t{sample}\t"
                f"{user}\t{timestamp}\t{comment}\t0\t{diagnostic}\n"
            )
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def update_diagnostic_ignore_status(
    diagnostic_file: Path,
    family_id: str,
    variant_key: str,
    sample: str,
    timestamp: str,
    ignore_value: str,
) -> bool:
    """Update the Ignore status for a specific diagnostic row.

    Uses read-all, modify, write-back pattern with file locking.

    Returns:
        True if update was successful, False otherwise
    """
    if not diagnostic_file.exists():
        return False

    rows = []
    fieldnames: list[str] = []
    with open(diagnostic_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows.append(row)

    updated = False
    for row in rows:
        if (
            row.get("FID") == family_id
            and row.get("Variant") == variant_key
            and row.get("Sample") == sample
            and row.get("Timestamp") == timestamp
        ):
            row["Ignore"] = ignore_value
            updated = True
            break

    if not updated:
        return False

    with open(diagnostic_file, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)

    return True


def load_family_diagnostics(
    snv_file: Path,
    sv_file: Path,
    family_id: str,
    sample_ids: List[str],
) -> List[Dict[str, str]]:
    """Load all diagnostic entries for a family and set of samples.

    Returns a list of raw row dicts (one per diagnostic entry).
    """
    entries: List[Dict[str, str]] = []
    sample_set = set(sample_ids)

    for diag_file in [snv_file, sv_file]:
        if not diag_file.exists():
            continue
        source = "snv" if diag_file == snv_file else "sv"
        with open(diag_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if (
                    row.get("FID") == family_id
                    and row.get("Sample") in sample_set
                    and row.get("Ignore", "0") != "1"
                ):
                    entries.append({**row, "_source": source})

    return entries
