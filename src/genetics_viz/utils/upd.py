"""Curator validation of uniparental disomy: ``<data_dir>/validations/upds.tsv``.

Three pipeline steps point at a UPD, and they share no identifier:

- the ``upd`` step's events (:mod:`genetics_viz.utils.upd_calls`): a span,
  an origin, a composition and a confidence, per member;
- the ``roh`` step's ``upd_candidate_chrom`` (one chromosome per person, no
  span, no parent);
- the ``sampleqc`` step's ``UPD_CANDIDATE`` pair flag. ``sample_a`` is the
  parent whose allele the child lacks, ``sample_b`` the child; the chromosome
  and span are only in the evidence (``"… lie on chr2, between 0.1 and 39.7
  Mb; …"``).

A UPD is therefore keyed by family, child and chromosome. ``Variant`` is
``chr<N>:UPD``, so a validation made from the UPD, ROH or Sample QC tab is the
same record. :func:`upd_candidates` merges the three into one row per child
and chromosome.

The span, composition and extent the curator validated, and the upd label
checked, are recorded after the shared columns (``Start End Composition
Extent UpdLabel``): the report says what was validated even after a newer
label calls the event differently.

Same layout and vocabulary as the SNV and SV files, so ``validation_loader``,
the badges and the validation pages read it unchanged. ``Inheritance`` holds
the parental origin: ``maternal`` (both copies from the mother), ``paternal``
or ``unknown``.

A ``present`` or ``absent`` validation acknowledges the matching Sample QC
flag (see :mod:`genetics_viz.utils.sample_qc_acks`).

Append-only under an exclusive ``flock``, like the other curation files.
"""

import csv
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

from . import upd_calls
from .auth import can_write
from .roh import find_roh_files, read_roh_table
from .sample_qc import FlagRow, load_family_qc
from .tsv import append_tsv_row
from .upd_calls import normalise_chrom
from .validation_files import validation_path

__all__ = ["normalise_chrom"]

HEADER = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    "Comment", "Ignore", "Timestamp",
    # Appended in v0.32.0: an older file gains them on its next save.
    "Start", "End", "Composition", "Extent", "UpdLabel",
]  # fmt: skip

VALIDATIONS = ["present", "absent", "uncertain"]

ORIGINS = ["unknown", "maternal", "paternal"]

#: Codes a curator records, worded through ``config/upd.yaml``; ``""`` is
#: not recorded.
COMPOSITIONS = ["", *upd_calls.VALIDATION_FORM.get("compositions", [])]
EXTENTS = ["", *upd_calls.VALIDATION_FORM.get("extents", [])]

FLAG_CODE = "UPD_CANDIDATE"

logger = logging.getLogger(__name__)

_EVIDENCE_CHROM = re.compile(r"\blie on (?:chr)?([0-9XYM]+)\b")
_EVIDENCE_SPAN = re.compile(r"\bbetween ([0-9.]+) and ([0-9.]+) Mb")


def upd_validation_file(data_dir: Path) -> Path:
    return validation_path(data_dir, "UPD")


def upd_key(chrom: Any) -> str:
    """``chr8:UPD``."""
    return f"chr{normalise_chrom(chrom)}:UPD"


def chrom_of_key(key: str) -> str:
    """``"8"`` for ``chr8:UPD``."""
    return normalise_chrom(key.split(":", 1)[0])


def chrom_from_evidence(flag: FlagRow) -> Optional[str]:
    """The chromosome a ``UPD_CANDIDATE`` flag's evidence names, else ``None``."""
    match = _EVIDENCE_CHROM.search(flag.evidence or "")
    return match.group(1) if match else None


def span_from_evidence(flag: FlagRow) -> Optional[Tuple[int, int]]:
    """The ``(start, end)`` in bp of a ``UPD_CANDIDATE`` flag's opposite
    homozygotes (``"between 0.1 and 39.7 Mb"``), else ``None``."""
    match = _EVIDENCE_SPAN.search(flag.evidence or "")
    if not match:
        return None
    return int(float(match.group(1)) * 1e6), int(float(match.group(2)) * 1e6)


def child_of(flag: FlagRow) -> str:
    """The child of a ``UPD_CANDIDATE`` pair: ``sample_b``."""
    return flag.sample_b or flag.sample_a


def is_upd_flag(flag: FlagRow) -> bool:
    return flag.code == FLAG_CODE and bool(flag.sample_b)


def default_origin(
    flags: Iterable[FlagRow],
    father: Optional[str],
    mother: Optional[str],
    events: Iterable[Dict[str, Any]] = (),
) -> str:
    """The origin the evidence points to, else ``unknown``.

    The upd step's events (trio and duo) name it; a Sample QC flag names the
    parent whose allele the child lacks, so both copies came from the other
    one. Evidence for both parents leaves it ``unknown``.
    """
    origins = {
        e.get("origin") for e in events if e.get("origin") in ("maternal", "paternal")
    }
    if not origins:
        for f in flags:
            if father and f.sample_a == father:
                origins.add("maternal")
            elif mother and f.sample_a == mother:
                origins.add("paternal")
    return origins.pop() if len(origins) == 1 else "unknown"


def save_upd_validation(
    data_dir: Path,
    *,
    family_id: str,
    sample: str,
    chrom: Any,
    user: str,
    validation: str,
    origin: str = "unknown",
    comment: str = "",
    timestamp: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    composition: str = "",
    extent: str = "",
    upd_label: str = "",
) -> None:
    """Append one validation. Refused for viewers and for unknown values.

    The role check sits here as well as on the form: the save button is only
    hidden from viewers, which a crafted event would not respect. ``start``
    and ``end`` (bp) are the span validated, both or neither; ``upd_label``
    the upd step's label the curator looked at.
    """
    if not can_write():
        raise PermissionError("Only curators and administrators can save validations")
    if validation not in VALIDATIONS:
        raise ValueError(f"Unknown validation: {validation!r}")
    if origin not in ORIGINS:
        raise ValueError(f"Unknown parental origin: {origin!r}")
    if not normalise_chrom(chrom):
        raise ValueError("No chromosome for the UPD")
    if not user:
        raise ValueError("No user to record the validation under")
    if composition not in COMPOSITIONS:
        raise ValueError(f"Unknown composition: {composition!r}")
    if extent not in EXTENTS:
        raise ValueError(f"Unknown extent: {extent!r}")
    if (start is None) != (end is None):
        raise ValueError("A span needs both a start and an end")
    if start is not None and end is not None and not 0 <= int(start) < int(end):
        raise ValueError("The span's start must lie before its end")

    fields = [
        family_id, upd_key(chrom), sample, user, origin, validation,
        comment, "0",
        timestamp or datetime.now().isoformat(),
        "" if start is None else str(int(start)),
        "" if end is None else str(int(end)),
        composition, extent, upd_label,
    ]  # fmt: skip
    append_tsv_row(upd_validation_file(data_dir), HEADER, fields)


def read_upd_validations(
    data_dir: Path, family_id: Optional[str] = None, key: Optional[str] = None
) -> List[Dict[str, str]]:
    """UPD validations, oldest first, optionally of one family and one chromosome."""
    path = upd_validation_file(data_dir)
    if not path.exists():
        return []
    with open(path, "r", newline="") as f:
        return [
            row
            for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
            if (family_id is None or row.get("FID") == family_id)
            and (key is None or row.get("Variant") == key)
        ]


def latest_validation(
    rows: Iterable[Dict[str, str]], sample: str
) -> Optional[Dict[str, str]]:
    """The sample's standing validation: its last row that is not ignored."""
    latest = None
    for row in rows:
        if row.get("Sample") == sample and row.get("Ignore", "0") != "1":
            latest = row
    return latest


def origin_label(origin: str) -> str:
    """``maternal`` → ``maternal``; ``unknown`` and empty → ``""``."""
    return "" if origin in ("", "unknown") else origin


# ---------------------------------------------------------------------------
# Evidence
# ---------------------------------------------------------------------------


@dataclass
class UpdEvidence:
    """What the pipeline says about one child's chromosome."""

    #: The ``froh`` row names this chromosome as ``upd_candidate_chrom``.
    roh: bool = False
    #: The family's ``UPD_CANDIDATE`` flags on this child and chromosome.
    flags: List[FlagRow] = field(default_factory=list)
    father: Optional[str] = None
    mother: Optional[str] = None
    #: The upd step's events on this child and chromosome, every confidence.
    events: List[Dict[str, Any]] = field(default_factory=list)
    #: The upd step's label they come from (``""`` without upd tables).
    upd_label: str = ""

    @property
    def event(self) -> Optional[Dict[str, Any]]:
        """The most confident event."""
        return upd_calls.best_event(self.events)

    @property
    def source(self) -> str:
        """``UPD step + ROH + Sample QC`` or part of it; the diagnostic's Impact."""
        parts = (
            (["UPD step"] if self.events else [])
            + (["ROH"] if self.roh else [])
            + (["Sample QC"] if self.flags else [])
        )
        return " + ".join(parts) or "Curator"


def flags_on(flags: Iterable[FlagRow], sample: str, chrom: Any) -> List[FlagRow]:
    """The ``UPD_CANDIDATE`` flags about ``sample`` on ``chrom``.

    A flag whose evidence names no chromosome is kept: nothing says it is
    about another one.
    """
    wanted = normalise_chrom(chrom)
    return [
        f
        for f in flags
        if is_upd_flag(f)
        and child_of(f) == sample
        and chrom_from_evidence(f) in (None, wanted)
    ]


def gather_evidence(
    data_dir: Path,
    family_id: str,
    sample: str,
    chrom: Any,
    upd_label: Optional[str] = None,
) -> UpdEvidence:
    """The family's upd step, ROH and Sample QC findings on ``sample``'s ``chrom``.

    Small family tables; the upd step's of ``upd_label`` (default: the
    newest). Never raises: a table that will not read is simply no evidence.
    """
    ev = UpdEvidence()
    wanted = normalise_chrom(chrom)
    try:
        calls = upd_calls.load_family_upd(data_dir, family_id, upd_label)
    except Exception:
        logger.warning("Failed to read the UPD calls of %s", family_id, exc_info=True)
        calls = None
    if calls is not None:
        ev.upd_label = calls.label
        ev.events = sorted(
            calls.events_of(sample, wanted), key=lambda e: -upd_calls.confidence_rank(e)
        )
        member = calls.member(sample) or {}
        ev.father = member.get("father_used") or None
        ev.mother = member.get("mother_used") or None
    try:
        files = find_roh_files(data_dir, family_id)
        froh = read_roh_table(files.path("froh")) if files else None
    except Exception:
        logger.warning("Failed to read the ROH of %s", family_id, exc_info=True)
        froh = None
    if froh is not None:
        for row in froh.to_dicts():
            if row.get("IID") == sample:
                ev.roh = normalise_chrom(row.get("upd_candidate_chrom")) == wanted
                ev.father = ev.father or row.get("father") or None
                ev.mother = ev.mother or row.get("mother") or None
    try:
        qc = load_family_qc(data_dir, family_id)
    except Exception:
        logger.warning("Failed to read the sample QC of %s", family_id, exc_info=True)
        qc = None
    if qc is not None:
        ev.flags = flags_on(qc.flags, sample, wanted)
        for row in qc.sample_rows():
            if row.get("sample") == sample:
                ev.father = ev.father or row.get("father") or None
                ev.mother = ev.mother or row.get("mother") or None
    return ev


def roh_runs(data_dir: Path, family_id: str, sample: str, chrom: Any) -> List[Dict]:
    """The ROH step's runs of ``sample`` on ``chrom`` (``start``, ``end``, …).

    Blocking; never raises.
    """
    wanted = normalise_chrom(chrom)
    try:
        files = find_roh_files(data_dir, family_id)
        roh = read_roh_table(files.path("roh")) if files else None
    except Exception:
        logger.warning("Failed to read the ROH of %s", family_id, exc_info=True)
        return []
    if roh is None or "IID" not in roh.columns:
        return []
    return [
        r
        for r in roh.to_dicts()
        if r.get("IID") == sample and normalise_chrom(r.get("chrom")) == wanted
    ]


# ---------------------------------------------------------------------------
# Candidates: one row per child and chromosome
# ---------------------------------------------------------------------------


def standing_validations(
    rows: Iterable[Dict[str, str]],
) -> Dict[Tuple[str, str], Dict[str, str]]:
    """``(sample, chrom) -> last row that is not ignored``, of one family's rows."""
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in rows:
        if row.get("Ignore", "0") == "1" or not row.get("Sample"):
            continue
        chrom = chrom_of_key(row.get("Variant") or "")
        if chrom:
            out[(row["Sample"], chrom)] = row
    return out


def validated_span(row: Optional[Dict[str, str]]) -> Optional[Tuple[int, int]]:
    """The ``(start, end)`` a validation recorded, else ``None``."""
    try:
        return int((row or {}).get("Start") or ""), int((row or {}).get("End") or "")
    except ValueError:
        return None


@dataclass
class UpdCandidate:
    """Everything known about one child's chromosome."""

    sample: str
    chrom: str
    #: The upd step's events on it, every confidence, the most confident first.
    events: List[Dict[str, Any]] = field(default_factory=list)
    roh: bool = False
    flags: List[FlagRow] = field(default_factory=list)
    #: The standing validation, if any.
    validation: Optional[Dict[str, str]] = None

    @property
    def key(self) -> str:
        return upd_key(self.chrom)

    @property
    def event(self) -> Optional[Dict[str, Any]]:
        return upd_calls.best_event(self.events)

    @property
    def confidence(self) -> Optional[str]:
        return (self.event or {}).get("confidence")

    @property
    def status(self) -> str:
        """The standing validation (``present``, ``absent``, ``uncertain``), or ``""``."""
        return (self.validation or {}).get("Validation", "")

    @property
    def reviewed(self) -> bool:
        """A ``present`` or ``absent`` validation; ``uncertain`` still needs a look."""
        return self.status in ("present", "absent")

    @property
    def span(self) -> Optional[Tuple[int, int]]:
        """The validated span, else the best event's."""
        span = validated_span(self.validation)
        if span is not None:
            return span
        event = self.event
        if event and event.get("start") is not None and event.get("end") is not None:
            return int(event["start"]), int(event["end"])
        return None


def upd_candidates(
    events: Iterable[Dict[str, Any]],
    froh_rows: Iterable[Dict[str, Any]],
    flags: Iterable[FlagRow],
    validations: Iterable[Dict[str, str]],
    confidences: Optional[Iterable[str]] = None,
    members: Optional[List[str]] = None,
) -> List[UpdCandidate]:
    """One candidate per (child, chromosome) any source points at.

    A chromosome is listed when an event of ``confidences`` (default: every
    one), the ROH ``upd_candidate_chrom``, a ``UPD_CANDIDATE`` flag naming it
    or a standing validation is about it; it then carries all its events,
    whatever their confidence. ``members`` restricts and orders the children.
    Sorted by member, the most confident first, then chromosome.
    """
    shown = set(confidences) if confidences is not None else None
    wanted = set(members) if members is not None else None
    found: Dict[Tuple[str, str], UpdCandidate] = {}

    def get(sample: str, chrom: Any) -> Optional[UpdCandidate]:
        chrom = normalise_chrom(chrom)
        if not sample or not chrom or (wanted is not None and sample not in wanted):
            return None
        return found.setdefault((sample, chrom), UpdCandidate(sample, chrom))

    all_events = list(events)
    for e in all_events:
        if shown is None or e.get("confidence") in shown:
            get(e.get("sample") or "", e.get("chrom"))
    for row in froh_rows:
        cand = get(row.get("IID") or "", row.get("upd_candidate_chrom"))
        if cand is not None:
            cand.roh = True
    for f in flags:
        if not is_upd_flag(f):
            continue
        chrom = chrom_from_evidence(f)
        cand = get(child_of(f), chrom) if chrom else None
        if cand is not None:
            cand.flags.append(f)
    for (sample, chrom), row in standing_validations(validations).items():
        cand = get(sample, chrom)
        if cand is not None:
            cand.validation = row
    for e in all_events:
        cand = found.get((e.get("sample") or "", normalise_chrom(e.get("chrom"))))
        if cand is not None:
            cand.events.append(e)
    for cand in found.values():
        cand.events.sort(key=lambda e: -upd_calls.confidence_rank(e))

    order = {m: i for i, m in enumerate(members or [])}
    chrom_order = {c: i for i, c in enumerate([str(i) for i in range(1, 23)] + ["X"])}
    return sorted(
        found.values(),
        key=lambda c: (
            order.get(c.sample, len(order)),
            c.sample,
            -upd_calls.confidence_info(c.confidence)["rank"],
            chrom_order.get(c.chrom, len(chrom_order)),
        ),
    )


# ---------------------------------------------------------------------------
# A family's UPD hints, for the tab, its badge and the page
# ---------------------------------------------------------------------------


@dataclass
class UpdHints:
    """The ROH and Sample QC findings that point at a UPD in one family."""

    #: ``froh`` rows naming an ``upd_candidate_chrom``.
    froh_rows: List[Dict[str, Any]] = field(default_factory=list)
    #: ``UPD_CANDIDATE`` flags of the family's own QC files.
    flags: List[FlagRow] = field(default_factory=list)

    @property
    def any(self) -> bool:
        return bool(self.froh_rows or self.flags)


def family_upd_hints(data_dir: Path, family_id: str) -> UpdHints:
    """The family's ROH and Sample QC UPD hints. Small tables; never raises."""
    hints = UpdHints()
    try:
        files = find_roh_files(data_dir, family_id)
        froh = read_roh_table(files.path("froh")) if files else None
    except Exception:
        logger.warning("Failed to read the ROH of %s", family_id, exc_info=True)
        froh = None
    if froh is not None and "upd_candidate_chrom" in froh.columns:
        hints.froh_rows = [
            r for r in froh.to_dicts() if r.get("upd_candidate_chrom") and r.get("IID")
        ]
    try:
        qc = load_family_qc(data_dir, family_id)
    except Exception:
        logger.warning("Failed to read the sample QC of %s", family_id, exc_info=True)
        qc = None
    if qc is not None:
        hints.flags = [f for f in qc.flags if is_upd_flag(f)]
    return hints


def family_candidates(
    data_dir: Path,
    family_id: str,
    members: Optional[List[str]] = None,
    confidences: Optional[Iterable[str]] = None,
) -> List[UpdCandidate]:
    """:func:`upd_candidates` of the family, newest upd label. Blocking, small."""
    calls = upd_calls.load_family_upd(data_dir, family_id)
    hints = family_upd_hints(data_dir, family_id)
    return upd_candidates(
        calls.events if calls is not None else [],
        hints.froh_rows,
        hints.flags,
        read_upd_validations(data_dir, family_id),
        upd_calls.DEFAULT_CONFIDENCES if confidences is None else confidences,
        members,
    )


def probe_family_upd(data_dir: Path, family_id: str) -> bool:
    """Whether the family's UPD tab has anything: upd tables, a hint, a validation."""
    if upd_calls.probe_upd_data(data_dir, family_id):
        return True
    if family_upd_hints(data_dir, family_id).any:
        return True
    return bool(read_upd_validations(data_dir, family_id))


def count_to_review(candidates: Iterable[UpdCandidate]) -> int:
    """Candidates no curator settled (no ``present`` or ``absent`` validation)."""
    return sum(1 for c in candidates if not c.reviewed)


def parse_upd_link(value: Optional[str]) -> Optional[Tuple[str, str]]:
    """``("S3", "2")`` for a link's ``upd=S3:2``; ``None`` when malformed."""
    sample, sep, chrom = (value or "").rpartition(":")
    chrom = normalise_chrom(chrom)
    if not sep or not sample or not chrom:
        return None
    return sample, chrom
