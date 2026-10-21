"""Curator validation of uniparental disomy: ``<data_dir>/validations/upds.tsv``.

Two pipeline steps point at a UPD, and they share no identifier:

- the ``roh`` step's ``upd_candidate_chrom`` (one chromosome per person, no
  span, no parent);
- the ``sampleqc`` step's ``UPD_CANDIDATE`` pair flag. ``sample_a`` is the
  parent whose allele the child lacks, ``sample_b`` the child; the chromosome
  and span are only in the evidence (``"… lie on chr2, between 0.1 and 39.7
  Mb; …"``).

A UPD is therefore keyed by family, child and chromosome. ``Variant`` is
``chr<N>:UPD``, so a validation made from the ROH tab and one made from the
Sample QC tab are the same record.

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
from typing import Any, Dict, Iterable, List, Optional

from .auth import can_write
from .roh import find_roh_files, read_roh_table
from .sample_qc import FlagRow, load_family_qc
from .tsv import append_tsv_row
from .validation_files import validation_path

HEADER = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    "Comment", "Ignore", "Timestamp",
]  # fmt: skip

VALIDATIONS = ["present", "absent", "uncertain"]

ORIGINS = ["unknown", "maternal", "paternal"]

FLAG_CODE = "UPD_CANDIDATE"

logger = logging.getLogger(__name__)

_EVIDENCE_CHROM = re.compile(r"\blie on (?:chr)?([0-9XYM]+)\b")


def upd_validation_file(data_dir: Path) -> Path:
    return validation_path(data_dir, "UPD")


def normalise_chrom(chrom: Any) -> str:
    """``"chr8"``, ``"8"`` and ``8`` are all ``"8"``."""
    text = str(chrom or "").strip()
    return text[3:] if text.lower().startswith("chr") else text


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


def child_of(flag: FlagRow) -> str:
    """The child of a ``UPD_CANDIDATE`` pair: ``sample_b``."""
    return flag.sample_b or flag.sample_a


def is_upd_flag(flag: FlagRow) -> bool:
    return flag.code == FLAG_CODE and bool(flag.sample_b)


def default_origin(
    flags: Iterable[FlagRow], father: Optional[str], mother: Optional[str]
) -> str:
    """The origin the Sample QC flags point to, else ``unknown``.

    The flag names the parent whose allele the child lacks, so both copies came
    from the other one. Flags against both parents leave it ``unknown``.
    """
    origins = set()
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
) -> None:
    """Append one validation. Refused for viewers and for unknown values.

    The role check sits here as well as on the form: the save button is only
    hidden from viewers, which a crafted event would not respect.
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

    fields = [
        family_id, upd_key(chrom), sample, user, origin, validation, comment, "0",
        timestamp or datetime.now().isoformat(),
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

    @property
    def source(self) -> str:
        """``ROH``, ``Sample QC``, ``ROH + Sample QC``; the diagnostic's Impact."""
        parts = (["ROH"] if self.roh else []) + (["Sample QC"] if self.flags else [])
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
    data_dir: Path, family_id: str, sample: str, chrom: Any
) -> UpdEvidence:
    """The family's ROH and Sample QC findings on ``sample``'s ``chrom``.

    Two small family tables. Never raises: a table that will not read is
    simply no evidence.
    """
    ev = UpdEvidence()
    wanted = normalise_chrom(chrom)
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
                ev.father = row.get("father") or None
                ev.mother = row.get("mother") or None
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
