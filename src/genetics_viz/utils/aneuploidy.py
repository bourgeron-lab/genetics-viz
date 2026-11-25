"""Aneuploidies: whole-chromosome events, and their curator validation.

Two pipeline results point at an aneuploidy, and the pipeline has no
aneuploidy step of its own:

- **autosomes** (and X, Y): same-direction WisecondorX segments past the
  permissive thresholds that cover ``aneuploidy_min_fraction`` of a
  chromosome's assessable span (the q arm of an acrocentric one), from the
  family's ``svs/wisecondorx/<FID>_aberrations.annotated.bed``;
- **sex chromosomes**: the ``sampleqc`` karyotype from the X and Y copy
  numbers, atypical (``XXY``, ``X0``) or not the declared sex's, and its
  ``KARYOTYPE_<k>`` flag.

The family report, the cohort SVs tab and the family SVs tab's Aneuploidy
card all read the events through :func:`chromosome_events`, so they agree.

A curator validates one chromosome of one sample. ``Variant`` is
``chr21:ANEUPLOIDY``; X and Y are one record, ``chrXY:ANEUPLOIDY``, since a
karyotype such as XXY is about both. The file
(``<data_dir>/validations/aneuploidies.tsv``) has the SNV/SV header, so
``validation_loader``, the badges and the validation pages read it unchanged,
plus what the curator recorded: ``Direction`` (gain, loss), ``Mosaic`` and
``Karyotype`` (ISCN, ``47,XY,+21``). ``Inheritance`` holds the origin
(``de novo``, ``maternal``, ``paternal``, ``unknown``).

A standing ``present`` or ``absent`` validation of ``chrXY`` acknowledges the
sample's ``KARYOTYPE_<k>`` flag (:mod:`genetics_viz.utils.sample_qc_acks`).

Append-only under an exclusive ``flock``, like the other curation files.
"""

import csv
import logging
import math
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import polars as pl
import yaml

from .auth import can_write
from .cytobands import CHROM_ORDER, CYTOBANDS, norm_chrom
from . import sample_qc
from .sample_qc import FlagRow, is_typical_karyotype
from .sharding import get_family_path
from .tsv import append_tsv_row
from .validation_files import validation_path
from .wisecondorx import parse_wisecondorx_bed, past_permissive

logger = logging.getLogger(__name__)

#: The validation and diagnostic type label (``VALIDATION_KINDS``).
KIND = "Aneuploidy"

HEADER = [
    "FID", "Variant", "Sample", "User", "Inheritance", "Validation",
    "Comment", "Ignore", "Timestamp",
    "Direction", "Mosaic", "Karyotype",
]  # fmt: skip

VALIDATIONS = ["present", "absent", "uncertain"]
ORIGINS = ["unknown", "de novo", "maternal", "paternal"]
#: ``""`` is not recorded.
DIRECTIONS = ["", "gain", "loss"]
MOSAIC = ["", "no", "yes"]

#: The chromosome of the sex chromosomes' record.
SEX_CHROMS = "XY"


def _min_fraction() -> float:
    path = Path(__file__).parent.parent / "config" / "report.yaml"
    try:
        with open(path) as f:
            config = yaml.safe_load(f) or {}
        return float(config.get("aneuploidy_min_fraction", 0.75))
    except (OSError, ValueError, yaml.YAMLError):
        return 0.75


#: The share of a chromosome's assessable span the events must cover: the
#: report's ``aneuploidy_min_fraction`` (``config/report.yaml``).
DEFAULT_MIN_FRACTION = _min_fraction()

_ACROCENTRIC = {"13", "14", "15", "21", "22"}
_SUFFIX = ":ANEUPLOIDY"


# ---------------------------------------------------------------------------
# Keys
# ---------------------------------------------------------------------------


def record_chrom(chrom: Any) -> str:
    """``"21"`` for ``chr21``; ``"XY"`` for X, Y or XY."""
    text = str(chrom or "").strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    text = text.upper()
    if text in ("X", "Y", "XY"):
        return SEX_CHROMS
    return norm_chrom(text) if text else ""


def aneuploidy_key(chrom: Any) -> str:
    """``chr21:ANEUPLOIDY``, ``chrXY:ANEUPLOIDY``."""
    return f"chr{record_chrom(chrom)}{_SUFFIX}"


def chrom_of_key(key: str) -> str:
    """``"21"`` for ``chr21:ANEUPLOIDY``."""
    return record_chrom(key.split(":", 1)[0])


def chrom_label(chrom: Any) -> str:
    """``chr21``, or ``sex chromosomes``."""
    chrom = record_chrom(chrom)
    return "sex chromosomes" if chrom == SEX_CHROMS else f"chr{chrom}"


def _chrom_rank(chrom: str) -> int:
    if chrom == SEX_CHROMS:
        return len(CHROM_ORDER)
    try:
        return CHROM_ORDER.index(chrom)
    except ValueError:
        return len(CHROM_ORDER) + 1


# ---------------------------------------------------------------------------
# Whole-chromosome events from WisecondorX
# ---------------------------------------------------------------------------


def _float(value: Any) -> Optional[float]:
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def assessable_span(chrom: str) -> Optional[Tuple[float, float]]:
    """``(start, end)`` in Mb of what read depth can assess: the q arm of an
    acrocentric chromosome, the whole chromosome otherwise."""
    bands = CYTOBANDS.get(chrom) or []
    if not bands:
        return None
    end = max(b["end"] for b in bands)
    start = 0.0
    if chrom in _ACROCENTRIC:
        start = max((b["end"] for b in bands if b["stain"] == "acen"), default=0.0)
    return start, end


@dataclass(frozen=True)
class ChromEvent:
    """A gain or loss covering most of one chromosome."""

    chrom: str
    #: ``gain`` or ``loss``.
    direction: str
    #: The share of the assessable span the segments cover.
    fraction: float

    @property
    def text(self) -> str:
        """``+21 (97%)``, ``−X (88%)``."""
        sign = "+" if self.direction == "gain" else "−"
        return f"{sign}{self.chrom} ({self.fraction * 100:.0f}%)"


def chromosome_events(
    calls: Optional[pl.DataFrame], min_fraction: float
) -> Dict[str, List[ChromEvent]]:
    """``sample -> events``: gains and losses covering most of a chromosome.

    Same-direction segments (gain for a positive log2 ratio) are summed per
    sample and chromosome, against the chromosome's assessable length. Pass
    the calls past the permissive thresholds only (:func:`family_calls`): a
    segment below them can span a chromosome with no copy-number change.
    """
    out: Dict[str, List[ChromEvent]] = {}
    if calls is None or len(calls) == 0:
        return out
    needed = {"chr", "start", "end", "ratio", "sample"}
    if not needed <= set(calls.columns):
        return out
    covered: Dict[Tuple[str, str, str], float] = {}
    for r in calls.to_dicts():
        ratio, start, end = _float(r["ratio"]), _float(r["start"]), _float(r["end"])
        if ratio is None or start is None or end is None or ratio == 0:
            continue
        chrom = norm_chrom(r["chr"])
        span = assessable_span(chrom)
        if span is None:
            continue
        overlap = min(end / 1e6, span[1]) - max(start / 1e6, span[0])
        if overlap > 0:
            key = (str(r["sample"]), chrom, "gain" if ratio > 0 else "loss")
            covered[key] = covered.get(key, 0.0) + overlap
    order = {c: i for i, c in enumerate(CHROM_ORDER)}
    for (sample, chrom, direction), mb in sorted(
        covered.items(), key=lambda kv: (kv[0][0], order.get(kv[0][1], 99))
    ):
        span = assessable_span(chrom)
        length = span[1] - span[0] if span else 0.0
        if length > 0 and mb / length >= min_fraction:
            out.setdefault(sample, []).append(ChromEvent(chrom, direction, mb / length))
    return out


def whole_chromosome_events(
    calls: Optional[pl.DataFrame], min_fraction: float
) -> Dict[str, List[str]]:
    """``sample -> ["+21 (97%)", "−18 (90%)"]``: :func:`chromosome_events` as text."""
    return {
        sample: [e.text for e in events]
        for sample, events in chromosome_events(calls, min_fraction).items()
    }


def wcx_file(data_dir: Path, family_id: str) -> Path:
    """The family's WisecondorX calls."""
    return (
        get_family_path(data_dir, family_id)
        / "svs"
        / "wisecondorx"
        / f"{family_id}_aberrations.annotated.bed"
    )


def family_calls(data_dir: Path, family_id: str) -> Optional[pl.DataFrame]:
    """The family's WisecondorX calls past the permissive thresholds, ``None``
    without a file. A file that will not parse is logged and gives ``None``."""
    path = wcx_file(data_dir, family_id)
    if not path.exists():
        return None
    try:
        calls = parse_wisecondorx_bed(path)
    except Exception:
        logger.warning("Failed to read %s", path, exc_info=True)
        return None
    if calls is not None and {"ratio", "zscore"} <= set(calls.columns):
        calls = calls.filter(past_permissive())
    return calls


# ---------------------------------------------------------------------------
# Validations
# ---------------------------------------------------------------------------


def aneuploidy_validation_file(data_dir: Path) -> Path:
    return validation_path(data_dir, KIND)


def save_aneuploidy_validation(
    data_dir: Path,
    *,
    family_id: str,
    sample: str,
    chrom: Any,
    user: str,
    validation: str,
    origin: str = "unknown",
    direction: str = "",
    mosaic: str = "",
    karyotype: str = "",
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
        raise ValueError(f"Unknown origin: {origin!r}")
    if direction not in DIRECTIONS:
        raise ValueError(f"Unknown direction: {direction!r}")
    if mosaic not in MOSAIC:
        raise ValueError(f"Unknown mosaic value: {mosaic!r}")
    if not record_chrom(chrom):
        raise ValueError("No chromosome for the aneuploidy")
    if not user:
        raise ValueError("No user to record the validation under")
    fields = [
        family_id, aneuploidy_key(chrom), sample, user, origin, validation,
        comment, "0",
        timestamp or datetime.now().isoformat(),
        direction, mosaic, " ".join(karyotype.split()),
    ]  # fmt: skip
    append_tsv_row(aneuploidy_validation_file(data_dir), HEADER, fields)


def read_aneuploidy_validations(
    data_dir: Path, family_id: Optional[str] = None, key: Optional[str] = None
) -> List[Dict[str, str]]:
    """Aneuploidy validations, oldest first, optionally of one family and one key."""
    path = aneuploidy_validation_file(data_dir)
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


def standing_validations(
    rows: Iterable[Dict[str, str]],
) -> Dict[Tuple[str, str], Dict[str, str]]:
    """``(sample, chrom) -> standing validation`` of one family's rows."""
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    for row in rows:
        if row.get("Ignore", "0") == "1":
            continue
        chrom = chrom_of_key(row.get("Variant") or "")
        if row.get("Sample") and chrom:
            out[(row["Sample"], chrom)] = row
    return out


# ---------------------------------------------------------------------------
# Karyotypes
# ---------------------------------------------------------------------------


def _flag_prefix() -> str:
    # Read through the module: a config reload rebinds KARYOTYPE_FLAG.
    return str(sample_qc.KARYOTYPE_FLAG.get("prefix") or "KARYOTYPE_")


def is_karyotype_flag(flag: FlagRow) -> bool:
    """A ``KARYOTYPE_<k>`` flag, about one sample."""
    return flag.code.startswith(_flag_prefix()) and not flag.sample_b


def flag_karyotype(flag: FlagRow) -> str:
    """``XXY`` for ``KARYOTYPE_XXY``."""
    prefix = _flag_prefix()
    return flag.code[len(prefix) :] if flag.code.startswith(prefix) else ""


def declared_karyotype(sex: str) -> str:
    """``XY`` for a declared male (``1``), ``XX`` for a female, else ``""``."""
    return {"1": "XY", "2": "XX"}.get(str(sex), "")


def sex_complement(karyotype: str) -> str:
    """ISCN sex chromosomes of an observed karyotype: ``X0`` is ``X``."""
    text = (karyotype or "").upper().replace("0", "")
    return text if text and set(text) <= {"X", "Y"} else ""


def suggested_karyotype(
    chrom: str,
    direction: str,
    observed: str = "",
    declared: str = "",
) -> str:
    """An ISCN karyotype to prefill: ``47,XY,+21``, ``45,XX,-18``, ``47,XXY``.

    Autosomes take the sample's observed sex chromosomes (else the declared
    sex's); sex chromosomes are the observed karyotype itself. ``""`` when
    there is not enough to write one.
    """
    chrom = record_chrom(chrom)
    if chrom == SEX_CHROMS:
        sex = sex_complement(observed)
        return f"{44 + len(sex)},{sex}" if sex else ""
    sex = sex_complement(observed) if is_typical_karyotype(observed) else ""
    sex = sex or declared
    if not sex or direction not in ("gain", "loss"):
        return ""
    if direction == "gain":
        return f"{44 + len(sex) + 1},{sex},+{chrom}"
    return f"{44 + len(sex) - 1},{sex},-{chrom}"


# ---------------------------------------------------------------------------
# Candidates
# ---------------------------------------------------------------------------


@dataclass
class AneuploidyCandidate:
    """One sample and chromosome some source, or a curator, points at."""

    sample: str
    #: ``"21"``, or ``"XY"`` for the sex chromosomes.
    chrom: str
    #: WisecondorX's whole-chromosome events on it (X and Y for ``XY``).
    events: List[ChromEvent] = field(default_factory=list)
    #: The Sample QC karyotype, for ``XY`` when atypical or not the declared sex's.
    karyotype: str = ""
    #: The sample's Sample QC karyotype, whatever it is (``""`` when not checked).
    observed: str = ""
    #: ``XY``/``XX`` from the pedigree sex, ``""`` when unknown.
    declared: str = ""
    #: The sample's ``KARYOTYPE_<k>`` flags (``XY`` only).
    flags: List[FlagRow] = field(default_factory=list)
    #: The standing validation, if any.
    validation: Optional[Dict[str, str]] = None

    @property
    def key(self) -> str:
        return aneuploidy_key(self.chrom)

    @property
    def status(self) -> str:
        """The standing validation, ``""`` when none."""
        return (self.validation or {}).get("Validation", "")

    @property
    def direction(self) -> str:
        """The validation's direction, else the events', else the karyotype's."""
        recorded = (self.validation or {}).get("Direction") or ""
        if recorded:
            return recorded
        directions = {e.direction for e in self.events}
        if len(directions) == 1:
            return directions.pop()
        if self.karyotype and self.declared:
            gained = len(sex_complement(self.karyotype)) - len(self.declared)
            return "gain" if gained > 0 else "loss" if gained < 0 else ""
        return ""

    @property
    def label(self) -> str:
        """``+21``, ``−18``, ``XXY``, or the chromosome when nothing says more."""
        if self.chrom == SEX_CHROMS:
            recorded = (self.validation or {}).get("Karyotype") or ""
            return (
                self.karyotype
                or next((flag_karyotype(f) for f in self.flags), "")
                or recorded
                or "sex chromosomes"
            )
        sign = {"gain": "+", "loss": "−"}.get(self.direction, "")
        return f"{sign}{self.chrom}" if sign else f"chr{self.chrom}"

    @property
    def sources(self) -> List[str]:
        out = []
        if self.events:
            out.append("WisecondorX " + ", ".join(e.text for e in self.events))
        if self.karyotype:
            declared = f", declared {self.declared}" if self.declared else ""
            out.append(f"Sample QC {self.karyotype}{declared}")
        elif self.flags:
            out.append("Sample QC flag")
        return out

    @property
    def suggested_karyotype(self) -> str:
        return suggested_karyotype(
            self.chrom, self.direction, self.karyotype or self.observed, self.declared
        )


def family_candidates(
    data_dir: Path,
    family_id: str,
    members: Sequence[Tuple[str, str]],
    qc_rows: Iterable[Dict[str, Any]] = (),
    flags: Iterable[FlagRow] = (),
    calls: Optional[pl.DataFrame] = None,
    min_fraction: float = DEFAULT_MIN_FRACTION,
) -> List[AneuploidyCandidate]:
    """One row per member and chromosome, in member then chromosome order.

    ``members`` are ``(sample, pedigree sex code)``; ``qc_rows`` the family's
    Sample QC sample rows; ``calls`` the WisecondorX calls past the
    permissive thresholds (:func:`family_calls`). A standing validation keeps
    its chromosome listed when the sources no longer point at it.
    """
    member_order = {iid: i for i, (iid, _) in enumerate(members)}
    declared = {iid: declared_karyotype(sex) for iid, sex in members}
    out: Dict[Tuple[str, str], AneuploidyCandidate] = {}
    qc_rows = list(qc_rows)
    observed = {
        str(r.get("sample") or ""): str(r.get("karyotype") or "")
        for r in qc_rows
        if str(r.get("karyotype") or "") not in ("", "unknown", "—")
    }

    def candidate(sample: str, chrom: str) -> AneuploidyCandidate:
        key = (sample, chrom)
        if key not in out:
            out[key] = AneuploidyCandidate(
                sample,
                chrom,
                declared=declared.get(sample, ""),
                observed=observed.get(sample, ""),
            )
        return out[key]

    for sample, events in chromosome_events(calls, min_fraction).items():
        if sample not in member_order:
            continue
        for e in events:
            candidate(sample, record_chrom(e.chrom)).events.append(e)
    for row in qc_rows:
        sample = str(row.get("sample") or "")
        karyotype = str(row.get("karyotype") or "")
        if sample not in member_order or karyotype in ("", "unknown", "—"):
            continue
        expected = declared.get(sample, "")
        if not is_typical_karyotype(karyotype) or (expected and karyotype != expected):
            candidate(sample, SEX_CHROMS).karyotype = karyotype
    for f in flags:
        if is_karyotype_flag(f) and f.sample_a in member_order:
            candidate(f.sample_a, SEX_CHROMS).flags.append(f)
    standing = standing_validations(read_aneuploidy_validations(data_dir, family_id))
    for (sample, chrom), row in standing.items():
        if sample in member_order:
            candidate(sample, chrom).validation = row
    return sorted(
        out.values(),
        key=lambda c: (member_order.get(c.sample, 0), _chrom_rank(c.chrom)),
    )


def count_to_review(candidates: Iterable[AneuploidyCandidate]) -> int:
    """Candidates no curator validated present or absent."""
    return sum(1 for c in candidates if c.status not in ("present", "absent"))


# ---------------------------------------------------------------------------
# One family's review
# ---------------------------------------------------------------------------


def member_sexes(
    store: Any, family_id: str, cohort_name: str = ""
) -> List[Tuple[str, str]]:
    """``(sample, sex code)`` of the family's sequenced members.

    From the named cohort, else the first cohort of the store listing the
    family (the standalone family page has none).
    """
    from .pedigree_labels import normalize_sex

    cohorts = getattr(store, "cohorts", {}) or {}
    candidates = [cohorts[cohort_name]] if cohort_name in cohorts else []
    candidates += [c for c in cohorts.values() if family_id in c.families]
    for cohort in candidates:
        family = cohort.families.get(family_id)
        if family is not None:
            return [(s.sample_id, normalize_sex(s.sex)) for s in family.samples]
    return []


@dataclass
class FamilyReview:
    """What the Aneuploidy card and dialog read for one family."""

    candidates: List[AneuploidyCandidate]
    #: The WisecondorX calls past the permissive thresholds, ``None`` without a file.
    calls: Optional[pl.DataFrame]
    #: Sample QC sample rows by sample.
    qc_rows: Dict[str, Dict[str, Any]]

    def candidate(self, sample: str, chrom: Any) -> Optional[AneuploidyCandidate]:
        chrom = record_chrom(chrom)
        return next(
            (c for c in self.candidates if c.sample == sample and c.chrom == chrom),
            None,
        )

    def segments(self, sample: str, chrom: Any) -> List[Dict[str, Any]]:
        """The sample's calls on ``chrom`` (X and Y for the sex chromosomes)."""
        if self.calls is None or not len(self.calls):
            return []
        chrom = record_chrom(chrom)
        wanted = {"X", "Y"} if chrom == SEX_CHROMS else {chrom}
        return [
            r
            for r in self.calls.to_dicts()
            if str(r.get("sample")) == sample and norm_chrom(r.get("chr", "")) in wanted
        ]


def load_family_review(
    data_dir: Path,
    family_id: str,
    members: Sequence[Tuple[str, str]],
    min_fraction: float = DEFAULT_MIN_FRACTION,
) -> FamilyReview:
    """The family's candidates, calls and QC rows. Blocking: call off the loop
    for a large WisecondorX file."""
    qc = sample_qc.load_family_qc(data_dir, family_id)
    qc_rows = qc.sample_rows() if qc is not None else []
    calls = family_calls(data_dir, family_id)
    return FamilyReview(
        candidates=family_candidates(
            data_dir,
            family_id,
            members,
            qc_rows,
            qc.flags if qc is not None else [],
            calls,
            min_fraction,
        ),
        calls=calls,
        qc_rows={str(r.get("sample")): r for r in qc_rows},
    )
