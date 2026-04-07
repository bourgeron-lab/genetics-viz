"""Shared notes data loading utilities.

Notes are stored in notes/notes.tsv in the data directory.
Each note has a FID, optional sample ID, message, user, and timestamp.
"""

import csv
import fcntl
from pathlib import Path
from typing import Dict, List

NOTES_HEADER = "FID\tSample\tMessage\tUser\tTimestamp\n"


def ensure_notes_file(file_path: Path) -> None:
    """Create the notes file with header if it doesn't exist."""
    file_path.parent.mkdir(parents=True, exist_ok=True)
    if not file_path.exists():
        with open(file_path, "w") as f:
            f.write(NOTES_HEADER)


def load_family_notes(notes_file: Path, family_id: str) -> List[Dict[str, str]]:
    """Load all notes for a family.

    Returns a list of row dicts sorted by timestamp descending (newest first).
    """
    entries: List[Dict[str, str]] = []

    if not notes_file.exists():
        return entries

    with open(notes_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
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
    ensure_notes_file(notes_file)
    # Sanitize message: replace tabs and newlines to keep TSV intact
    clean_msg = message.replace("\t", " ").replace("\n", " ").replace("\r", "")
    with open(notes_file, "a") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            f.write(f"{family_id}\t{sample}\t{clean_msg}\t{user}\t{timestamp}\n")
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def delete_note(
    notes_file: Path,
    family_id: str,
    timestamp: str,
    user: str,
) -> bool:
    """Hard-delete a note from the TSV file.

    Matches on FID + Timestamp + User as composite key.
    Uses read-all, remove, write-back with file locking.

    Returns True if a note was deleted, False otherwise.
    """
    if not notes_file.exists():
        return False

    rows: List[Dict[str, str]] = []
    fieldnames: list[str] = []
    with open(notes_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for row in reader:
            rows.append(dict(row))

    original_len = len(rows)
    rows = [
        r
        for r in rows
        if not (
            r.get("FID") == family_id
            and r.get("Timestamp") == timestamp
            and r.get("User") == user
        )
    ]

    if len(rows) == original_len:
        return False

    with open(notes_file, "w") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)

    return True
