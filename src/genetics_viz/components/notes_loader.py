"""Shared notes data loading utilities.

Notes are stored in notes/notes.tsv in the data directory.
Each note has a FID, optional sample ID, message, user, and timestamp.
"""

import csv
from pathlib import Path
from typing import Dict, List

from genetics_viz.utils.tsv import append_tsv_row, append_tsv_rows, rewrite_tsv_rows

NOTES_HEADER = ["FID", "Sample", "Message", "User", "Timestamp"]


def ensure_notes_file(file_path: Path) -> None:
    """Create the notes file with header if it doesn't exist."""
    append_tsv_rows(file_path, NOTES_HEADER, [])


def load_family_notes(notes_file: Path, family_id: str) -> List[Dict[str, str]]:
    """Load all notes for a family.

    Returns a list of row dicts sorted by timestamp descending (newest first).
    """
    entries: List[Dict[str, str]] = []

    if not notes_file.exists():
        return entries

    with open(notes_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
        for row in reader:
            if row.get("FID") == family_id:
                entries.append(dict(row))

    entries.sort(key=lambda r: r.get("Timestamp", ""), reverse=True)
    return entries


def save_note(
    notes_file: Path,
    family_id: str,
    sample: str,
    message: str,
    user: str,
    timestamp: str,
) -> None:
    """Append a note entry to the TSV file with file locking."""
    append_tsv_row(
        notes_file, NOTES_HEADER, [family_id, sample, message, user, timestamp]
    )


def delete_note(
    notes_file: Path,
    family_id: str,
    timestamp: str,
    user: str,
) -> bool:
    """Hard-delete a note from the TSV file.

    Matches on FID + Timestamp + User as composite key.

    Returns True if a note was deleted, False otherwise.
    """

    def remove(rows: List[Dict[str, str]]) -> bool:
        kept = [
            r
            for r in rows
            if not (
                r.get("FID") == family_id
                and r.get("Timestamp") == timestamp
                and r.get("User") == user
            )
        ]
        if len(kept) == len(rows):
            return False
        rows[:] = kept
        return True

    return rewrite_tsv_rows(notes_file, remove)
