"""Curator conclusions of the family report: ``<data_dir>/reports/family_summaries.tsv``.

The report a curator sends to the referring clinician opens with their own
words: what was found, and whether recontacting the family is advised. One row
per save, keyed by family; the last row of a family wins, and an empty summary
clears it. Append-only under an exclusive ``flock``, like the curation files.

A conclusion keeps its paragraphs: a line break is written ``\\n`` (and a
backslash ``\\\\``), since :func:`tsv.append_tsv_row` turns real ones into
spaces.
"""

import csv
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Optional, Tuple

from .auth import can_write
from .tsv import append_tsv_row

logger = logging.getLogger(__name__)

HEADER = ["FID", "Summary", "User", "Timestamp"]


def summaries_file(data_dir: Path) -> Path:
    return data_dir / "reports" / "family_summaries.tsv"


@dataclass(frozen=True)
class FamilySummary:
    """The standing conclusion of one family."""

    text: str
    user: str
    timestamp: str

    @property
    def date(self) -> str:
        return self.timestamp[:10]


def _encode(text: str) -> str:
    return text.replace("\\", "\\\\").replace("\r\n", "\n").replace("\n", "\\n")


def _decode(text: str) -> str:
    out = []
    chars = iter(text)
    for ch in chars:
        if ch == "\\":
            nxt = next(chars, "")
            out.append("\n" if nxt == "n" else nxt)
        else:
            out.append(ch)
    return "".join(out)


def save_family_summary(
    data_dir: Path,
    family_id: str,
    text: str,
    *,
    user: str,
    timestamp: Optional[str] = None,
) -> None:
    """Record ``text`` as the family's conclusion. Viewers are refused."""
    # The role check sits here as well as on the form: the editor is only
    # hidden from viewers, which a crafted event would not respect.
    if not can_write():
        raise PermissionError("Only curators and administrators can write a report")
    if not user:
        raise ValueError("No user to record the conclusion under")
    append_tsv_row(
        summaries_file(data_dir),
        HEADER,
        [
            family_id,
            _encode(text.strip()),
            user,
            timestamp or datetime.now().isoformat(timespec="seconds"),
        ],
    )


def read_family_summaries(data_dir: Path) -> Dict[str, FamilySummary]:
    """``FID -> FamilySummary``: the last row of each family, unless empty."""
    path = summaries_file(data_dir)
    if not path.exists():
        return {}
    latest: Dict[str, FamilySummary] = {}
    with open(path, "r", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
            fid = row.get("FID") or ""
            if fid:
                latest[fid] = FamilySummary(
                    text=_decode(row.get("Summary") or ""),
                    user=row.get("User") or "",
                    timestamp=row.get("Timestamp") or "",
                )
    return {fid: s for fid, s in latest.items() if s.text}


#: The file's ``(path, mtime, size)`` -> its summaries. A save appends, so it
#: changes the size and the next read picks it up.
_summary_cache: Dict[Tuple[str, float, int], Dict[str, FamilySummary]] = {}


def latest_family_summary(data_dir: Path, family_id: str) -> Optional[FamilySummary]:
    """The family's standing conclusion, cached on the file's stat. Never raises."""
    path = summaries_file(data_dir)
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (str(path), stat.st_mtime, stat.st_size)
    if key not in _summary_cache:
        for old in [k for k in _summary_cache if k[0] == key[0]]:
            del _summary_cache[old]
        try:
            _summary_cache[key] = read_family_summaries(data_dir)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            return None
    return _summary_cache[key].get(family_id)


def clear_summary_cache() -> None:
    _summary_cache.clear()
