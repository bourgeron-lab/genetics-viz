"""Family status: where each family stands in its review, set by a curator.

A family is *not reviewed* until a curator says otherwise. The statuses and
their wording are in ``config/family_status.yaml``; the records live in
``<data_dir>/reports/family_status.tsv`` (``FID Status Comment User
Timestamp``), next to the report conclusions. One row per save, append-only
under the curation files' ``flock``; the last row of a family wins. Keyed by
family only, so a family listed in two cohorts has one status.

The status never changes by itself. What happened since it was set is worked
out at read time (:func:`family_state`) and shown beside it as *reasons* to
look again: a pedigree edit (the family is *reopened*), a diagnostic recorded
later, a diagnostic that contradicts the status, a QC flag or step results
written later. Saving a status again clears the reasons about dates.
"""

import csv
import logging
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import yaml

from .auth import can_write
from .diagnostic_files import diagnostic_paths
from .tsv import append_tsv_row

logger = logging.getLogger(__name__)

HEADER = ["FID", "Status", "Comment", "User", "Timestamp"]

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "family_status.yaml"


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


_CONFIG = _load_config()
#: Status key -> ``label``, ``color``, ``yield``, ``description``; in review order.
STATUSES: Dict[str, Dict[str, str]] = dict(_CONFIG.get("statuses") or {})
REASONS: Dict[str, Dict[str, str]] = dict(_CONFIG.get("reasons") or {})
SOLVING_VERDICTS: Tuple[str, ...] = tuple(_CONFIG.get("solving_verdicts") or ())

NOT_REVIEWED = "not_reviewed"
IN_REVIEW = "in_review"


def _with_yield(kind: str) -> frozenset:
    return frozenset(k for k, v in STATUSES.items() if v.get("yield") == kind)


#: Statuses by how they count in the yield.
SOLVED = _with_yield("solved")
REVIEWED = _with_yield("reviewed") | SOLVED
EXCLUDED = _with_yield("excluded")


def status_info(status: Optional[str]) -> Dict[str, str]:
    """A status's wording; an unknown one is shown raw, grey."""
    key = status or NOT_REVIEWED
    return STATUSES.get(
        key, {"label": key, "color": "#9e9e9e", "yield": "open", "description": ""}
    )


def status_label(status: Optional[str]) -> str:
    return status_info(status)["label"]


def status_file(data_dir: Path) -> Path:
    return Path(data_dir) / "reports" / "family_status.tsv"


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StatusRecord:
    """One save of a family's status."""

    family_id: str
    status: str
    comment: str
    user: str
    timestamp: str

    @property
    def date(self) -> str:
        return self.timestamp[:10]

    @property
    def moment(self) -> Optional[datetime]:
        return as_moment(self.timestamp)


def as_moment(timestamp: Any) -> Optional[datetime]:
    """An ISO timestamp or a file's ``st_mtime`` as an aware datetime."""
    if timestamp is None or timestamp == "":
        return None
    if isinstance(timestamp, (int, float)):
        return datetime.fromtimestamp(timestamp).astimezone()
    try:
        moment = datetime.fromisoformat(str(timestamp))
    except ValueError:
        return None
    return moment if moment.tzinfo else moment.astimezone()


def save_family_status(
    data_dir: Path,
    family_id: str,
    status: str,
    comment: str = "",
    *,
    user: str,
    timestamp: Optional[str] = None,
) -> StatusRecord:
    """Record ``status`` for the family. Viewers are refused."""
    _check_role()
    return _write_status(data_dir, family_id, status, comment, user, timestamp)


async def save_family_status_async(
    data_dir: Path, family_id: str, status: str, comment: str = "", *, user: str
) -> StatusRecord:
    """:func:`save_family_status` for a page: the role is checked on the event
    loop (it reads the session), the file is written in a thread."""
    import asyncio

    _check_role()
    return await asyncio.to_thread(
        _write_status, data_dir, family_id, status, comment, user, None
    )


def _check_role() -> None:
    # Checked here as well as on the control, which a crafted event bypasses.
    if not can_write():
        raise PermissionError("Only curators and administrators can set a status")


def _write_status(
    data_dir: Path,
    family_id: str,
    status: str,
    comment: str,
    user: str,
    timestamp: Optional[str],
) -> StatusRecord:
    if status not in STATUSES:
        raise ValueError(f"Unknown family status: {status}")
    if not user:
        raise ValueError("No user to record the status under")
    if not family_id:
        raise ValueError("No family to record the status of")
    record = StatusRecord(
        family_id,
        status,
        " ".join((comment or "").split()),
        user,
        timestamp or datetime.now().isoformat(timespec="seconds"),
    )
    append_tsv_row(
        status_file(data_dir),
        HEADER,
        [family_id, status, record.comment, user, record.timestamp],
    )
    return record


def _read_history(path: Path) -> Dict[str, List[StatusRecord]]:
    out: Dict[str, List[StatusRecord]] = {}
    with open(path, "r", newline="") as f:
        for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
            fid = row.get("FID") or ""
            status = row.get("Status") or ""
            if not fid or not status:
                continue
            out.setdefault(fid, []).append(
                StatusRecord(
                    fid,
                    status,
                    row.get("Comment") or "",
                    row.get("User") or "",
                    row.get("Timestamp") or "",
                )
            )
    return out


#: The file's ``(path, mtime, size)`` -> every family's records, oldest first.
_cache: Dict[Tuple[str, float, int], Dict[str, List[StatusRecord]]] = {}


def status_history(data_dir: Path) -> Dict[str, List[StatusRecord]]:
    """``FID -> records``, oldest first; cached on the file's stat. Never raises."""
    path = status_file(data_dir)
    try:
        stat = path.stat()
    except OSError:
        return {}
    key = (str(path), stat.st_mtime, stat.st_size)
    if key not in _cache:
        for old in [k for k in _cache if k[0] == key[0]]:
            del _cache[old]
        try:
            _cache[key] = _read_history(path)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            return {}
    return _cache[key]


def family_statuses(data_dir: Path) -> Dict[str, StatusRecord]:
    """``FID -> its standing record`` (the last row of each family)."""
    return {fid: rows[-1] for fid, rows in status_history(data_dir).items() if rows}


def family_history(data_dir: Path, family_id: str) -> List[StatusRecord]:
    """The family's records, newest first."""
    return list(reversed(status_history(data_dir).get(family_id, [])))


def clear_status_cache() -> None:
    _cache.clear()
    _diag_cache.clear()


# ---------------------------------------------------------------------------
# Diagnostics per family
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FamilyDiagnostics:
    """A family's standing diagnostics, as far as its status cares."""

    verdicts: frozenset = frozenset()
    #: The newest standing diagnostic's timestamp.
    latest: str = ""

    @property
    def solving(self) -> bool:
        return bool(self.verdicts & set(SOLVING_VERDICTS))


_diag_cache: Dict[Tuple[Tuple[str, float, int], ...], Dict[str, FamilyDiagnostics]] = {}


def diagnostics_by_family(data_dir: Path) -> Dict[str, FamilyDiagnostics]:
    """``FID -> FamilyDiagnostics`` from every diagnostic file, ignored rows left out.

    Cached on the files' stats. Never raises: a file that fails is skipped.
    """
    paths = [p for _, p in diagnostic_paths(data_dir)]
    stats = []
    for path in paths:
        try:
            st = path.stat()
            stats.append((str(path), st.st_mtime, st.st_size))
        except OSError:
            continue
    key = tuple(stats)
    if key in _diag_cache:
        return _diag_cache[key]
    verdicts: Dict[str, set] = {}
    latest: Dict[str, str] = {}
    for name, _mtime, _size in stats:
        try:
            with open(name, "r", newline="") as f:
                for row in csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE):
                    if row.get("Ignore", "0") == "1":
                        continue
                    fid = row.get("FID") or ""
                    verdict = row.get("Diagnostic") or ""
                    if not fid or not verdict:
                        continue
                    verdicts.setdefault(fid, set()).add(verdict)
                    ts = row.get("Timestamp") or ""
                    if _later(ts, latest.get(fid, "")):
                        latest[fid] = ts
        except Exception:
            logger.warning("Failed to read %s", name, exc_info=True)
    out = {
        fid: FamilyDiagnostics(frozenset(v), latest.get(fid, ""))
        for fid, v in verdicts.items()
    }
    _diag_cache.clear()
    _diag_cache[key] = out
    return out


def _later(a: str, b: str) -> bool:
    """Whether timestamp ``a`` is after ``b`` (an unreadable one never is)."""
    ma, mb = as_moment(a), as_moment(b)
    if ma is None:
        return False
    return mb is None or ma > mb


# ---------------------------------------------------------------------------
# Where a family stands
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Reason:
    """Why a family should be looked at again."""

    code: str
    label: str
    text: str

    @property
    def reopens(self) -> bool:
        """A pedigree edit: the family is *reopened*."""
        return self.code.startswith("pedigree_")


def _reason(code: str, **values: str) -> Reason:
    info = REASONS.get(code, {})
    text = str(info.get("text", code))
    if values:
        text = text.format(**values)
    return Reason(code, str(info.get("label", code)), text)


@dataclass(frozen=True)
class FamilyState:
    """A family's status, and the reasons to look at it again."""

    family_id: str
    record: Optional[StatusRecord] = None
    reasons: Tuple[Reason, ...] = ()
    #: Whether a sequenced member is affected; a family without one is left
    #: out of the yield unless a curator gave it a status.
    affected: bool = True

    @property
    def status(self) -> str:
        return self.record.status if self.record else NOT_REVIEWED

    @property
    def reopened(self) -> bool:
        return any(r.reopens for r in self.reasons)

    @property
    def waiting_rerun(self) -> bool:
        return any(r.code == "pedigree_waiting" for r in self.reasons)

    @property
    def counted(self) -> bool:
        """Whether the family counts in the yield at all."""
        if self.status in EXCLUDED:
            return False
        return self.affected or self.record is not None


def family_state(
    family_id: str,
    record: Optional[StatusRecord],
    *,
    edit: Any = None,
    diagnostics: Optional[FamilyDiagnostics] = None,
    qc_high_written: Any = None,
    results: Sequence[Tuple[str, Any]] = (),
    affected: bool = True,
    disagreements: Sequence[str] = (),
) -> FamilyState:
    """Work out the reasons to look at a family again.

    ``edit`` is the family's ``pedigree_edit.EditStatus``; ``qc_high_written``
    when the family has an open high-severity QC flag, the newest write
    (ISO timestamp or ``st_mtime``) of the tables holding the flags;
    ``results`` the ``(step, written)`` of each step output kept for the
    family. ``disagreements`` word the open high-severity pedigree findings
    between cohorts (``pedigree_consistency``). Reasons about dates need a
    status; the diagnostic consistency and the disagreements hold whatever
    the dates.
    """
    reasons: List[Reason] = []
    status = record.status if record else NOT_REVIEWED
    since = record.moment if record and status != NOT_REVIEWED else None
    if since is not None and edit is not None:
        edited = as_moment(edit.last.timestamp)
        if edited is not None and edited > since:
            reasons.append(
                _reason("pedigree_waiting" if edit.awaiting_rerun else "pedigree_ready")
            )
    diag = diagnostics or FamilyDiagnostics()
    if since is not None:
        latest = as_moment(diag.latest)
        if latest is not None and latest > since:
            reasons.append(_reason("diagnostic_after"))
    if status in SOLVED and not diag.solving:
        reasons.append(_reason("solved_without_diagnostic"))
    elif diag.solving and status not in SOLVED | EXCLUDED | {IN_REVIEW}:
        reasons.append(_reason("diagnostic_not_solved"))
    if since is not None:
        written = as_moment(qc_high_written)
        if written is not None and written > since:
            reasons.append(_reason("qc_after"))
        steps = sorted(
            {
                step
                for step, when in results
                if (m := as_moment(when)) is not None and m > since
            }
        )
        if steps:
            reasons.append(_reason("results_after", steps=", ".join(steps)))
    if disagreements:
        reasons.append(_reason("pedigree_differs", what="; ".join(disagreements)))
    return FamilyState(family_id, record, tuple(reasons), affected)


# ---------------------------------------------------------------------------
# Counts and yield
# ---------------------------------------------------------------------------


@dataclass
class StatusCounts:
    """A set of families by status, and the diagnostic yield."""

    #: Status -> families, counted families only, in the config's order.
    by_status: Dict[str, int] = field(default_factory=dict)
    total: int = 0
    excluded: int = 0
    reopened: int = 0
    to_check: int = 0

    @property
    def solved(self) -> int:
        return sum(n for s, n in self.by_status.items() if s in SOLVED)

    @property
    def reviewed(self) -> int:
        return sum(n for s, n in self.by_status.items() if s in REVIEWED)

    @property
    def analysable(self) -> int:
        return sum(self.by_status.values())

    @property
    def not_reviewed(self) -> int:
        return self.by_status.get(NOT_REVIEWED, 0)

    def yield_text(self, compact: bool = False) -> str:
        """``Solved 3 of 10 reviewed (30%) · 40 analysable``; ``compact`` drops
        the analysable families (``Solved 3/10 reviewed (30%)``)."""
        if not self.analysable:
            return "No family to review"
        sep = "/" if compact else " of "
        text = f"Solved {self.solved:,}{sep}{self.reviewed:,} reviewed"
        if self.reviewed:
            text += f" ({100 * self.solved / self.reviewed:.0f}%)"
        if compact:
            return text
        return text + f" · {self.analysable:,} analysable"


def count_statuses(states: Iterable[FamilyState]) -> StatusCounts:
    counts = StatusCounts(by_status={k: 0 for k in STATUSES if k not in EXCLUDED})
    for state in states:
        counts.total += 1
        if not state.counted:
            counts.excluded += 1
            continue
        counts.by_status[state.status] = counts.by_status.get(state.status, 0) + 1
        counts.reopened += state.reopened
        counts.to_check += bool(state.reasons)
    return counts


# ---------------------------------------------------------------------------
# Badges
# ---------------------------------------------------------------------------


def status_tooltip(state: FamilyState, curator: bool) -> str:
    """Status, who and when, the comment (curators) and the reasons."""
    info = status_info(state.status)
    lines = [info["label"]]
    record = state.record
    if record is not None:
        lines[0] += f", {record.user}, {record.date}"
        if curator and record.comment:
            lines.append(record.comment)
    elif not state.affected:
        lines.append("No affected sequenced member: not in the yield")
    if info.get("description"):
        lines.append(str(info["description"]))
    lines.extend(f"⚠ {r.text}" for r in state.reasons)
    return "\n".join(lines)


def status_badges(state: FamilyState, curator: bool) -> List[Dict[str, str]]:
    """The status badge, then one amber badge when there is a reason to look again."""
    info = status_info(state.status)
    tooltip = status_tooltip(state, curator)
    badges = [
        {
            "label": info["label"]
            if state.affected or state.record
            else f"{info['label']} (no affected)",
            # Grey only for a family nobody gave a status and out of the yield.
            "color": info["color"] if state.counted or state.record else "#bdbdbd",
            "tooltip": tooltip,
        }
    ]
    if state.reasons:
        label = "Reopened" if state.reopened else "Check"
        badges.append({"label": label, "color": "#b45309", "tooltip": tooltip})
    return badges


#: Filter options beyond the statuses.
REOPENED = "Reopened"
TO_CHECK = "Needs re-review"


def filter_categories(state: FamilyState) -> List[str]:
    """The Families filter's options a family answers to."""
    cats = [status_label(state.status)]
    if state.reopened:
        cats.append(REOPENED)
    if state.reasons:
        cats.append(TO_CHECK)
    return cats


def filter_options() -> List[str]:
    return [v["label"] for v in STATUSES.values()] + [REOPENED, TO_CHECK]


def status_rank(status: str) -> int:
    """Sort key: the config's review order."""
    keys = list(STATUSES)
    return keys.index(status) if status in keys else len(keys)


def newest(paths: Iterable[Path]) -> Optional[float]:
    """The newest ``st_mtime`` among ``paths`` that exist."""
    best: Optional[float] = None
    for path in paths:
        try:
            m = path.stat().st_mtime
        except OSError:
            continue
        best = m if best is None else max(best, m)
    return best


def step_outputs(data_dir: Path, family_id: str) -> List[Tuple[str, float]]:
    """``(step, st_mtime)`` of each step snapshot (``*.inputs.tsv``) of the family.

    Blocking: walks the family directory. Each pedigree-reading step rewrites
    its snapshot when it reruns for the family.
    """
    from .pedigree_edit import EDITS_DIRNAME
    from .sharding import get_family_path
    from .tsv import inputs_step

    family_dir = get_family_path(data_dir, family_id)
    out: List[Tuple[str, float]] = []
    if not family_dir.is_dir():
        return out
    for path in family_dir.rglob("*.inputs.tsv"):
        if EDITS_DIRNAME in path.parts:
            continue
        try:
            step = inputs_step(path) or path.parent.name
            out.append((step, path.stat().st_mtime))
        except OSError:
            continue
    return out


def state_map(
    family_ids: Iterable[str],
    records: Mapping[str, StatusRecord],
    edits: Mapping[str, Any],
    diagnostics: Mapping[str, FamilyDiagnostics],
    *,
    affected: Mapping[str, bool],
    qc_high_written: Mapping[str, Any] = {},
    disagreements: Mapping[str, Sequence[str]] = {},
) -> Dict[str, FamilyState]:
    """Every family's state, from data read once for all of them."""
    return {
        fid: family_state(
            fid,
            records.get(fid),
            edit=edits.get(fid),
            diagnostics=diagnostics.get(fid),
            qc_high_written=qc_high_written.get(fid),
            affected=affected.get(fid, True),
            disagreements=disagreements.get(fid, ()),
        )
        for fid in family_ids
    }


# ---------------------------------------------------------------------------
# Reading what a state needs
# ---------------------------------------------------------------------------


def cohort_edits(cohort: Any, data_dir: Path) -> Dict[str, Any]:
    """``FID -> EditStatus``: each family's last pedigree edit, against the
    moment every pedigree-reading step of the cohort had run again
    (``None``: no cohort, so every edit waits)."""
    from .cohort_state import pedigree_run_started
    from .pedigree_edit import edit_statuses, read_edit_index

    started = None
    if cohort is not None:
        try:
            started = pedigree_run_started(cohort)
        except Exception:
            logger.warning("Failed to read %s's run record", cohort.name, exc_info=True)
    return edit_statuses(read_edit_index(data_dir), started)


def cohort_qc_high_written(cohort: Any, data_dir: Path) -> Dict[str, float]:
    """``FID -> newest write of the cohort's QC flag tables``, for the families
    with an unacknowledged high-severity flag in them (Sample QC, genome QC)."""
    from .genome_qc import cohort_flags as genome_flags
    from .genome_qc import find_cohort_gqc_labels, get_cohort_gqc_dir, gqc_path
    from .sample_qc import load_cohort_qc, qc_path
    from .sample_qc_acks import active_flags, load_acks

    flags: List[Any] = []
    paths: List[Path] = []
    try:
        qc = load_cohort_qc(cohort)
        if qc is not None:
            flags.extend(qc.flags)
            paths.append(qc_path(qc.directory, qc.prefix, qc.label, "flags"))
        flags.extend(genome_flags(cohort))
        directory = get_cohort_gqc_dir(cohort)
        labels = find_cohort_gqc_labels(directory, cohort.name)
        if labels:
            paths.append(gqc_path(directory, cohort.name, labels[0], "flags"))
    except Exception:
        logger.warning("Failed to read %s's QC flags", cohort.name, exc_info=True)
        return {}
    written = newest(paths)
    if written is None:
        return {}
    out: Dict[str, float] = {}
    for f in active_flags(flags, load_acks(data_dir)):
        if f.severity == "high":
            for fid in f.families:
                out[fid] = written
    return out


def has_affected(family: Any) -> bool:
    from .pedigree_labels import normalize_pheno

    return any(
        normalize_pheno(s.phenotype) == "2" for s in getattr(family, "samples", [])
    )


def load_family_state(
    data_dir: Path,
    family_id: str,
    family: Any,
    cohorts: Sequence[Any] = (),
    results: Optional[Sequence[Tuple[str, Any]]] = None,
    store: Any = None,
) -> FamilyState:
    """Blocking: one family's state, with every reason (the family pages).

    ``cohorts`` hold the family; the first one's run record says whether an
    edit still waits for ghfc-ngs. ``results`` are :func:`step_outputs`,
    walked here when not given (a page walks the family directory once).
    """
    from .genome_qc import family_flags as genome_family_flags
    from .genome_qc import find_cohort_gqc_labels, get_cohort_gqc_dir, gqc_path
    from .pedigree_edit import last_edits
    from .sample_qc import load_family_qc, qc_path
    from .sample_qc_acks import active_flags, load_acks

    record = family_statuses(data_dir).get(family_id)
    edit = None
    if family_id in last_edits(data_dir):
        edit = cohort_edits(cohorts[0] if cohorts else None, data_dir).get(family_id)
    qc_written = None
    try:
        flags: List[Any] = []
        paths: List[Path] = []
        qc = load_family_qc(data_dir, family_id)
        if qc is not None:
            flags.extend(qc.flags)
            paths.append(qc_path(qc.directory, qc.prefix, qc.label, "flags"))
        for cohort in cohorts:
            flags.extend(genome_family_flags(cohort, family_id))
            directory = get_cohort_gqc_dir(cohort)
            labels = find_cohort_gqc_labels(directory, cohort.name)
            if labels:
                paths.append(gqc_path(directory, cohort.name, labels[0], "flags"))
        if any(f.severity == "high" for f in active_flags(flags, load_acks(data_dir))):
            qc_written = newest(paths)
    except Exception:
        logger.warning("Failed to read the QC flags of %s", family_id, exc_info=True)
    if results is None:
        try:
            results = step_outputs(data_dir, family_id) if record is not None else []
        except Exception:
            logger.warning(
                "Failed to list the step outputs of %s", family_id, exc_info=True
            )
            results = []
    return family_state(
        family_id,
        record,
        edit=edit,
        diagnostics=diagnostics_by_family(data_dir).get(family_id),
        qc_high_written=qc_written,
        results=results,
        affected=has_affected(family) if family is not None else True,
        disagreements=(
            disagreements_by_family(store).get(family_id, ())
            if store is not None
            else ()
        ),
    )


def cohort_status_counts(
    cohort: Any, data_dir: Path, store: Any = None
) -> StatusCounts:
    """Blocking: the cohort's families by status (the home card).

    Every reason but the step outputs, which would walk every family
    directory.
    """
    states = state_map(
        cohort.families,
        family_statuses(data_dir),
        cohort_edits(cohort, data_dir),
        diagnostics_by_family(data_dir),
        affected={fid: has_affected(f) for fid, f in cohort.families.items()},
        qc_high_written=cohort_qc_high_written(cohort, data_dir),
        disagreements=disagreements_by_family(store) if store is not None else {},
    )
    return count_statuses(states.values())


def disagreements_by_family(store: Any) -> Dict[str, List[str]]:
    """``FID -> open high-severity pedigree findings``, worded. Never raises.

    Aliases (medium) are left out: a family known under another ID is no
    reason to review it again.
    """
    from . import pedigree_consistency as pc

    try:
        findings = pc.open_findings(
            pc.scan(store).findings, pc.load_acks(store.data_dir)
        )
    except Exception:
        logger.warning("Failed to compare the cohort pedigrees", exc_info=True)
        return {}
    out: Dict[str, List[str]] = {}
    for f in findings:
        if f.severity != "high":
            continue
        for fid in f.family_ids:
            others = ", ".join(
                f"{c}" if g == fid else f"{g} ({c})" for c, g in f.others(fid)
            )
            out.setdefault(fid, []).append(f"{f.info.get('label', f.code)}: {others}")
    return out
