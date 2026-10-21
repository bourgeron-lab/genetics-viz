"""Curator acknowledgements of Sample QC flags: ``<data_dir>/validations/sampleqc.tsv``.

Real cohorts raise the same flags at every run: monozygotic twins, a declared
consanguineous couple, an adopted child, a contamination judged tolerable. A
curator reviews such a flag once and records why, and it stops counting in the
family badge, the member icons, the cohort column and panel, the home card and
the de novo warning. The Sample QC tabs still list it, greyed.

- One row per decision, keyed by :attr:`FlagRow.ack_key`: the flag code, the
  flag's ``family_id`` cell (``"FA,FB"`` for a cross-family pair), ``sample_a``
  and ``sample_b`` (empty for a one-sample flag). The sites label is not in the
  key, so an acknowledgement survives a rerun or a new label.
- ``Status`` is ``explained`` or ``confirmed`` (a confirmed issue);
  ``withdrawn`` takes an acknowledgement back. The last row of a key wins.
- ``Note`` is the curator's free text. Viewers see the status, and the note
  only for the categories whose reader wording names no relationship
  (``reader_note`` in ``config/sample_qc.yaml``).

Not a variant validation, so it is not one of ``VALIDATION_KINDS``. Append-only
under an exclusive ``flock``, like the other curation files.

A ``UPD_CANDIDATE`` flag is also acknowledged by a UPD validation
(``validations/upds.tsv``, :mod:`genetics_viz.utils.upd`): ``present`` as
``confirmed``, ``absent`` as ``explained``. That validation knows the child
and the chromosome but not the parent the flag names, so it cannot be an
:data:`AckKey`; :class:`Acks` matches it to the flag instead. Look flags up
through :class:`Acks` (``acks.get(flag)``, ``flag in acks``), never by key.
"""

import csv
import logging
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from . import sample_qc, upd
from .auth import can_write
from .sample_qc import FlagRow, flag_badges, flag_title
from .tsv import append_tsv_row

logger = logging.getLogger(__name__)

HEADER = ["Flag", "FID", "SampleA", "SampleB", "Status", "Note", "User", "Timestamp"]

STATUSES = ("explained", "confirmed")
WITHDRAWN = "withdrawn"

AckKey = Tuple[str, str, str, str]


def ack_file(data_dir: Path) -> Path:
    return data_dir / "validations" / "sampleqc.tsv"


@dataclass(frozen=True)
class Ack:
    """The standing acknowledgement of one flag."""

    status: str
    note: str
    user: str
    timestamp: str
    #: ``ack`` for a row of ``sampleqc.tsv``, ``upd`` for a UPD validation.
    source: str = "ack"

    @property
    def label(self) -> str:
        return ack_status_info(self.status)["label"]

    @property
    def color(self) -> str:
        return ack_status_info(self.status)["color"]

    @property
    def date(self) -> str:
        return self.timestamp[:10]


def ack_status_info(status: str) -> Dict[str, str]:
    """Label, colour and help text of an acknowledgement status."""
    cfg = sample_qc.ACKNOWLEDGEMENT.get("statuses", {}).get(status, {})
    return {
        "label": cfg.get("label", status),
        "color": cfg.get("color", "#757575"),
        "help": cfg.get("help", ""),
    }


def muted_color() -> str:
    return sample_qc.ACKNOWLEDGEMENT.get("muted_color", "#bdbdbd")


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def _append(
    data_dir: Path, flag: FlagRow, status: str, note: str, user: str, timestamp: str
) -> None:
    # The role check sits here as well as on the form: the button is only
    # hidden from viewers, which a crafted event would not respect.
    if not can_write():
        raise PermissionError(
            "Only curators and administrators can acknowledge QC flags"
        )
    if not user:
        raise ValueError("No user to record the acknowledgement under")
    code, family_id, sample_a, sample_b = flag.ack_key
    fields = [code, family_id, sample_a, sample_b, status, note.strip(), user]
    append_tsv_row(
        ack_file(data_dir),
        HEADER,
        fields + [timestamp or datetime.now().isoformat(timespec="seconds")],
    )


def save_flag_ack(
    data_dir: Path,
    flag: FlagRow,
    *,
    user: str,
    status: str,
    note: str = "",
    timestamp: Optional[str] = None,
) -> None:
    """Acknowledge ``flag`` (or change its acknowledgement). Viewers are refused."""
    if status not in STATUSES:
        raise ValueError(f"Unknown acknowledgement status: {status!r}")
    _append(data_dir, flag, status, note, user, timestamp or "")


def withdraw_flag_ack(
    data_dir: Path, flag: FlagRow, *, user: str, timestamp: Optional[str] = None
) -> None:
    """Take back the acknowledgement of ``flag``: it counts again."""
    _append(data_dir, flag, WITHDRAWN, "", user, timestamp or "")


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_acks(data_dir: Path) -> Dict[AckKey, Ack]:
    """Every standing acknowledgement: the last row of each key, unless withdrawn."""
    path = ack_file(data_dir)
    if not path.exists():
        return {}
    out: Dict[AckKey, Ack] = {}
    with open(path, "r", newline="") as f:
        # QUOTE_NONE: a note starting with a quote is text, not a quoted field.
        for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
            key = (
                row.get("Flag") or "",
                row.get("FID") or "",
                row.get("SampleA") or "",
                row.get("SampleB") or "",
            )
            if not key[0] or not key[2]:
                continue
            status = row.get("Status") or ""
            if status == WITHDRAWN:
                out.pop(key, None)
            elif status in STATUSES:
                out[key] = Ack(
                    status=status,
                    note=row.get("Note") or "",
                    user=row.get("User") or "",
                    timestamp=row.get("Timestamp") or "",
                )
    return out


#: UPD validation -> acknowledgement status; ``uncertain`` leaves the flag counting.
_UPD_STATUS = {"present": "confirmed", "absent": "explained"}

UpdKey = Tuple[str, str]


def read_upd_acks(data_dir: Path) -> Dict[UpdKey, Dict[str, Ack]]:
    """``(FID, child) -> chromosome -> Ack`` from each standing UPD validation."""
    rows = upd.read_upd_validations(data_dir)
    latest: Dict[Tuple[str, str, str], Dict[str, str]] = {}
    for row in rows:
        if row.get("Ignore", "0") == "1":
            continue
        key = (row.get("FID") or "", row.get("Sample") or "", row.get("Variant") or "")
        if all(key):
            latest[key] = row
    out: Dict[UpdKey, Dict[str, Ack]] = {}
    for (fid, sample, variant), row in latest.items():
        validation = row.get("Validation") or ""
        status = _UPD_STATUS.get(validation)
        if status is None:
            continue
        origin = upd.origin_label(row.get("Inheritance") or "")
        note = f"UPD validated {validation}" + (f" ({origin})" if origin else "")
        if row.get("Comment"):
            note += f": {row['Comment']}"
        out.setdefault((fid, sample), {})[upd.chrom_of_key(variant)] = Ack(
            status=status,
            note=note,
            user=row.get("User") or "",
            timestamp=row.get("Timestamp") or "",
            source="upd",
        )
    return out


class Acks:
    """The standing acknowledgements, looked up by flag.

    ``by_key`` holds the rows of ``sampleqc.tsv``; ``upd`` the acknowledgements
    UPD validations give ``UPD_CANDIDATE`` flags, matched on the family, the
    child (``sample_b``) and the chromosome the evidence names. A flag whose
    evidence names none matches any validated chromosome of the child. When a
    flag has both, the newer one stands.
    """

    def __init__(
        self,
        by_key: Optional[Dict[AckKey, Ack]] = None,
        upd_acks: Optional[Dict[UpdKey, Dict[str, Ack]]] = None,
    ) -> None:
        self.by_key: Dict[AckKey, Ack] = by_key or {}
        self.upd: Dict[UpdKey, Dict[str, Ack]] = upd_acks or {}

    def _from_upd(self, flag: FlagRow) -> Optional[Ack]:
        if not self.upd or not upd.is_upd_flag(flag):
            return None
        by_chrom = self.upd.get((flag.family_id, upd.child_of(flag)))
        if not by_chrom:
            return None
        chrom = upd.chrom_from_evidence(flag)
        if chrom is not None:
            return by_chrom.get(chrom)
        return max(by_chrom.values(), key=lambda a: a.timestamp)

    def get(self, flag: FlagRow) -> Optional[Ack]:
        own = self.by_key.get(flag.ack_key)
        derived = self._from_upd(flag)
        if own is None or derived is None:
            return own or derived
        return derived if derived.timestamp > own.timestamp else own

    def __contains__(self, flag: object) -> bool:
        return isinstance(flag, FlagRow) and self.get(flag) is not None

    def __getitem__(self, flag: FlagRow) -> Ack:
        ack = self.get(flag)
        if ack is None:
            raise KeyError(flag.ack_key)
        return ack

    def __iter__(self) -> Iterator[AckKey]:
        return iter(self.by_key)

    def __len__(self) -> int:
        return len(self.by_key) + sum(len(v) for v in self.upd.values())


#: ``(sampleqc.tsv stat, upds.tsv stat)`` -> acknowledgements. A save appends,
#: so it changes the size and the next read picks it up.
_ack_cache: Dict[Tuple[Tuple[str, float, int], ...], Acks] = {}


def _stat_key(path: Path) -> Tuple[str, float, int]:
    try:
        stat = path.stat()
    except OSError:
        return (str(path), 0.0, -1)
    return (str(path), stat.st_mtime, stat.st_size)


def load_acks(data_dir: Path) -> Acks:
    """:func:`read_acks` and :func:`read_upd_acks`, cached on the files' stat. Never raises."""
    key = (_stat_key(ack_file(data_dir)), _stat_key(upd.upd_validation_file(data_dir)))
    if key[0][2] < 0 and key[1][2] < 0:
        return Acks()
    if key not in _ack_cache:
        for old in [k for k in _ack_cache if k[0][0] == key[0][0]]:
            del _ack_cache[old]
        try:
            _ack_cache[key] = Acks(read_acks(data_dir), read_upd_acks(data_dir))
        except Exception:
            logger.warning("Failed to read the acknowledgements", exc_info=True)
            return Acks()
    return _ack_cache[key]


def clear_ack_cache() -> None:
    _ack_cache.clear()


# ---------------------------------------------------------------------------
# Shaping
# ---------------------------------------------------------------------------


def split_acknowledged(
    flags: Iterable[FlagRow], acks: Acks
) -> Tuple[List[FlagRow], List[FlagRow]]:
    """``(still counting, acknowledged)``, each in the input order."""
    active: List[FlagRow] = []
    acked: List[FlagRow] = []
    for f in flags:
        (acked if f in acks else active).append(f)
    return active, acked


def active_flags(flags: Iterable[FlagRow], acks: Acks) -> List[FlagRow]:
    """The flags no curator has acknowledged."""
    return split_acknowledged(flags, acks)[0]


def note_visible(flag: FlagRow, show_detail: bool) -> bool:
    """Whether the user may read the note: curators always, viewers by category."""
    if show_detail:
        return True
    return bool(sample_qc.CATEGORIES.get(flag.info.category, {}).get("reader_note"))


def ack_summary(flag: FlagRow, ack: Ack, show_detail: bool) -> str:
    """``Explained by alice, 2026-10-20: known MZ twins`` (the note if allowed)."""
    text = f"{ack.label} by {ack.user or 'a curator'}"
    if ack.date:
        text += f", {ack.date}"
    if ack.note:
        text += f": {ack.note}" if note_visible(flag, show_detail) else " (note)"
    return text


def acknowledged_badges(
    flags: Iterable[FlagRow], acks: Acks, show_detail: bool
) -> List[Dict[str, str]]:
    """:func:`flag_badges` of acknowledged flags, greyed, the review in the tooltip."""
    out = []
    for f in flags:
        ack = acks.get(f)
        if ack is None:
            continue
        for badge in flag_badges([f], show_detail):
            badge["color"] = muted_color()
            badge["tooltip"] += f"\nAcknowledged — {ack_summary(f, ack, show_detail)}"
            out.append(badge)
    return out


def acknowledged_titles(
    flags: Iterable[FlagRow], acks: Acks, show_detail: bool
) -> List[str]:
    """One line per acknowledged flag: its title and the review, for tooltips."""
    return list(
        dict.fromkeys(
            f"{flag_title(f, show_detail)} — {acks[f].label.lower()}"
            for f in flags
            if f in acks
        )
    )
