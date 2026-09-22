"""The ghfc-ngs run record of a cohort: ``.ghfc-ngs.state.json``.

The workflow writes one of these into every cohort directory, recording when the
cohort last ran, with which pedigree and parameters, and **how far along each
step is**. The schema is documented in ``COHORT_STATE.md`` of
bourgeron-lab/ghfc-ngs and written by ``lib/CohortState.groovy``.

genetics-viz reads that record and never measures completion itself. Two
reasons, both decisive: re-deriving it would duplicate the pipeline's planning
logic and drift from it, and it is unaffordable on a network-mounted data
directory -- a cold existence check measures ~5.6 ms over CIFS, so a
2700-sample cohort costs about a minute of stat calls, against one stat and one
small read for the record.

Things the schema documentation warns about, honoured here:

- **A null step is unmeasured, not 0%.** ``ancestry`` is null unless the step
  was requested and ``extractor`` is always null, because the pipeline has no
  on-disk check for it. Conflating the two is the easiest way to mislead.
- **``schema_version`` must be checked before any field is trusted.**
- **A missing file means "not run since this feature shipped"**, not "never
  run". Absence is not evidence.
- **``status: running`` is not a liveness signal.** There is no heartbeat; a
  killed run leaves the record behind until the next run closes it out.
- **The recorded ``path`` fields are compute-cluster paths** and do not resolve
  here, so drift detection hashes the local file instead.
"""

from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import yaml

logger = logging.getLogger(__name__)

#: Written into every cohort directory by the pipeline. Nothing else may edit it.
STATE_FILENAME = ".ghfc-ngs.state.json"

#: The only schema this module knows how to read. New keys are additive by
#: design, so an unknown key is ignored; a different version is not guessed at.
SUPPORTED_SCHEMA_VERSION = 1

#: How long a ``running`` record may sit untouched before it is called suspect.
#: There is no heartbeat, so this is a judgement, not a fact.
STALE_RUNNING_AFTER = timedelta(hours=24)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "pipeline_steps.yaml"


def _load_config() -> Dict[str, Any]:
    """Load the step display metadata from YAML."""
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


_CONFIG: Dict[str, Any] = _load_config()

STEP_META: Dict[str, Dict[str, str]] = _CONFIG.get("steps", {})
#: Display order. Insertion order from the YAML.
STEP_ORDER: List[str] = list(STEP_META)
RUN_STATUS_META: Dict[str, Dict[str, str]] = _CONFIG.get("run_status", {})
DEFAULT_STATUS_META: Dict[str, str] = _CONFIG.get("default_status", {})
NO_RECORD_STATUS_META: Dict[str, str] = _CONFIG.get("no_record_status", {})


def reload_pipeline_steps_config() -> None:
    """Reload the step display metadata from YAML."""
    global _CONFIG, STEP_META, STEP_ORDER, RUN_STATUS_META
    global DEFAULT_STATUS_META, NO_RECORD_STATUS_META
    _CONFIG = _load_config()
    STEP_META = _CONFIG.get("steps", {})
    STEP_ORDER = list(STEP_META)
    RUN_STATUS_META = _CONFIG.get("run_status", {})
    DEFAULT_STATUS_META = _CONFIG.get("default_status", {})
    NO_RECORD_STATUS_META = _CONFIG.get("no_record_status", {})


# --------------------------------------------------------------------------
# Parsed record
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class StepCompletion:
    """One step's progress, as the pipeline measured it against the pedigree."""

    step: str
    done: Optional[int] = None
    total: Optional[int] = None
    pct: Optional[float] = None

    @property
    def measured(self) -> bool:
        """False when the pipeline reported null -- unmeasured, NOT zero."""
        return self.done is not None and self.total is not None

    @property
    def complete(self) -> bool:
        """True only when a measured step has reached its pedigree denominator."""
        return (
            self.measured
            and self.total is not None
            and self.total > 0
            and self.done == self.total
        )


@dataclass(frozen=True)
class RunRecord:
    """One pipeline run against one cohort. Self-contained, as the schema is."""

    status: str
    run_id: Optional[str] = None
    started_at: Optional[datetime] = None
    finished_at: Optional[datetime] = None
    cohort_name: Optional[str] = None
    pipeline_version: Optional[str] = None
    pipeline_revision: Optional[str] = None
    pipeline_commit: Optional[str] = None
    pedigree_path: Optional[str] = None
    pedigree_sha256: Optional[str] = None
    pedigree_families: Optional[int] = None
    pedigree_individuals: Optional[int] = None
    params_path: Optional[str] = None
    params_sha256: Optional[str] = None
    steps_requested: List[str] = field(default_factory=list)
    completion_measured: str = ""
    outputs_may_be_incomplete: bool = False
    reason: Optional[str] = None
    #: Empty when the record carried ``completion: null`` -- a run that stopped
    #: before a plan was built, which is every validation failure.
    completion: Dict[str, StepCompletion] = field(default_factory=dict)

    @property
    def duration(self) -> Optional[timedelta]:
        """Wall-clock duration; the schema stores the endpoints, not this."""
        if self.started_at is None or self.finished_at is None:
            return None
        return self.finished_at - self.started_at

    @property
    def is_running(self) -> bool:
        return self.status == "running"

    @property
    def is_stale_running(self) -> bool:
        """A ``running`` record old enough that the run has probably died.

        Nextflow exits without a shutdown hook, so a walltime kill or a lost
        launch node leaves this record with nothing to close it out.
        """
        if not self.is_running or self.started_at is None:
            return False
        return _now(self.started_at) - self.started_at > STALE_RUNNING_AFTER

    @property
    def measured_before_run(self) -> bool:
        """The counts predate the run, so they describe the tree at launch."""
        return self.completion_measured == "before"

    @property
    def short_commit(self) -> Optional[str]:
        """The first 8 characters of the commit that ran, when known."""
        return self.pipeline_commit[:8] if self.pipeline_commit else None


@dataclass(frozen=True)
class CohortState:
    """The whole state file: the last run, the last success, and the history."""

    path: Path
    schema_version: Optional[int] = None
    #: True when the file announces a schema this module does not know. Nothing
    #: below is trustworthy in that case.
    unsupported: bool = False
    cohort_name: Optional[str] = None
    last_run: Optional[RunRecord] = None
    last_successful_run: Optional[RunRecord] = None
    #: Up to 10 terminal records, newest first. A ``running`` record is never here.
    history: List[RunRecord] = field(default_factory=list)

    @property
    def reference_run(self) -> Optional[RunRecord]:
        """The run to compare the current files against.

        The last success when there is one, so a later failure never makes the
        cohort look undocumented; otherwise the last run.
        """
        return self.last_successful_run or self.last_run


# --------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------


def _now(reference: datetime) -> datetime:
    """Current time, made comparable with ``reference``."""
    if reference.tzinfo is None:
        return datetime.now()
    return datetime.now(timezone.utc)


def _parse_timestamp(value: Any) -> Optional[datetime]:
    """Parse an ISO-8601 timestamp with offset, tolerating anything else."""
    if not isinstance(value, str) or not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        logger.warning("Unparseable timestamp in %s: %r", STATE_FILENAME, value)
        return None


def _parse_int(value: Any) -> Optional[int]:
    return value if isinstance(value, int) and not isinstance(value, bool) else None


def _parse_completion(raw: Any) -> Dict[str, StepCompletion]:
    """Parse the ``completion`` block.

    A null block means no plan was built. A null *entry* means that step was not
    measured, which this represents as an entry with no numbers rather than as
    an absent key, so the UI can say "not measured" instead of "0%".
    """
    if not isinstance(raw, dict):
        return {}

    parsed: Dict[str, StepCompletion] = {}
    for step, value in raw.items():
        name = str(step)
        if not isinstance(value, dict):
            parsed[name] = StepCompletion(step=name)
            continue
        pct = value.get("pct")
        parsed[name] = StepCompletion(
            step=name,
            done=_parse_int(value.get("done")),
            total=_parse_int(value.get("total")),
            pct=float(pct) if isinstance(pct, (int, float)) else None,
        )
    return parsed


def _parse_record(raw: Any) -> Optional[RunRecord]:
    """Parse one run record. Unknown keys are ignored: the schema is additive."""
    if not isinstance(raw, dict):
        return None

    pipeline = raw.get("pipeline") if isinstance(raw.get("pipeline"), dict) else {}
    pedigree = raw.get("pedigree") if isinstance(raw.get("pedigree"), dict) else {}
    params = raw.get("params_file") if isinstance(raw.get("params_file"), dict) else {}

    steps_raw = raw.get("steps_requested")
    steps = (
        [str(s) for s in steps_raw if s is not None]
        if isinstance(steps_raw, (list, tuple))
        else []
    )

    return RunRecord(
        status=str(raw.get("status") or ""),
        run_id=_as_str(raw.get("run_id")),
        started_at=_parse_timestamp(raw.get("started_at")),
        finished_at=_parse_timestamp(raw.get("finished_at")),
        cohort_name=_as_str(raw.get("cohort_name")),
        pipeline_version=_as_str(pipeline.get("version")),
        pipeline_revision=_as_str(pipeline.get("revision")),
        pipeline_commit=_as_str(pipeline.get("commit_id")),
        pedigree_path=_as_str(pedigree.get("path")),
        pedigree_sha256=_as_str(pedigree.get("sha256")),
        pedigree_families=_parse_int(pedigree.get("families")),
        pedigree_individuals=_parse_int(pedigree.get("individuals")),
        params_path=_as_str(params.get("path")),
        params_sha256=_as_str(params.get("sha256")),
        steps_requested=steps,
        completion_measured=str(raw.get("completion_measured") or ""),
        outputs_may_be_incomplete=bool(raw.get("outputs_may_be_incomplete")),
        # Present only on validation failures; a task failure carries none.
        reason=_as_str(raw.get("reason")),
        completion=_parse_completion(raw.get("completion")),
    )


def _as_str(value: Any) -> Optional[str]:
    return str(value) if isinstance(value, str) and value else None


#: ``(path, mtime, size)`` -> parsed state. Keyed on the stat, so a record
#: rewritten mid-run is picked up on the next read without any invalidation
#: plumbing -- which is why the change monitor deliberately does not watch this
#: file: it is rewritten on every run and would force a full pedigree reparse.
_state_cache: Dict[Tuple[str, float, int], CohortState] = {}


def get_state_path(cohort_dir: Path) -> Path:
    """Return the state file path for a cohort directory (may not exist)."""
    return cohort_dir / STATE_FILENAME


def read_cohort_state(cohort_dir: Path) -> Optional[CohortState]:
    """Read a cohort's run record.

    Returns None when the file is absent -- which, per the schema
    documentation, means "not run since run recording shipped", not "never
    run". A malformed file is also None, with a warning: a broken record is a
    degraded state, not a reason to fail the page.
    """
    path = get_state_path(cohort_dir)
    try:
        stat = path.stat()
    except OSError:
        return None

    key = (str(path), stat.st_mtime, stat.st_size)
    cached = _state_cache.get(key)
    if cached is not None:
        return cached

    try:
        with open(path, "r") as f:
            raw = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        logger.warning("Could not read %s: %s", path, e)
        return None

    if not isinstance(raw, dict):
        logger.warning("%s is not a JSON object", path)
        return None

    version = _parse_int(raw.get("schema_version"))
    if version != SUPPORTED_SCHEMA_VERSION:
        # Deliberately parse nothing else: a field may have changed meaning.
        logger.warning(
            "%s announces schema_version %r, this build reads %d",
            path,
            raw.get("schema_version"),
            SUPPORTED_SCHEMA_VERSION,
        )
        state = CohortState(path=path, schema_version=version, unsupported=True)
        _state_cache[key] = state
        return state

    history_raw = raw.get("history")
    history = (
        [rec for rec in (_parse_record(r) for r in history_raw) if rec is not None]
        if isinstance(history_raw, (list, tuple))
        else []
    )

    state = CohortState(
        path=path,
        schema_version=version,
        unsupported=False,
        cohort_name=_as_str(raw.get("cohort_name")),
        last_run=_parse_record(raw.get("last_run")),
        last_successful_run=_parse_record(raw.get("last_successful_run")),
        history=history,
    )
    _state_cache[key] = state
    return state


def clear_state_cache() -> None:
    """Drop every cached run record."""
    _state_cache.clear()


# --------------------------------------------------------------------------
# Drift: has the input changed since the run that produced these outputs?
# --------------------------------------------------------------------------

#: ``unchanged`` | ``changed`` | ``unknown``. ``unknown`` covers both a record
#: without a checksum and a local file that cannot be hashed.
Drift = str

_sha_cache: Dict[Tuple[str, float, int], str] = {}


def sha256_file(path: Path) -> Optional[str]:
    """SHA-256 of a file's bytes, cached on its stat."""
    try:
        stat = path.stat()
    except OSError:
        return None

    key = (str(path), stat.st_mtime, stat.st_size)
    cached = _sha_cache.get(key)
    if cached is not None:
        return cached

    digest = hashlib.sha256()
    try:
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                digest.update(chunk)
    except OSError as e:
        logger.warning("Could not hash %s: %s", path, e)
        return None

    hexdigest = digest.hexdigest()
    _sha_cache[key] = hexdigest
    return hexdigest


def file_drift(recorded_sha: Optional[str], local_path: Optional[Path]) -> Drift:
    """Compare a recorded checksum against the local file.

    The ``path`` recorded in the state file is a compute-cluster path and does
    not resolve here, so the caller passes the local file and this only ever
    compares bytes.
    """
    if not recorded_sha or local_path is None:
        return "unknown"
    local_sha = sha256_file(local_path)
    if local_sha is None:
        return "unknown"
    return "unchanged" if local_sha == recorded_sha else "changed"


def clear_sha_cache() -> None:
    """Drop every cached file checksum."""
    _sha_cache.clear()


# --------------------------------------------------------------------------
# Display helpers
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class StepRow:
    """One step as the UI shows it: the measurement plus its display metadata."""

    step: str
    label: str
    unit: str
    note: str
    completion: StepCompletion
    #: In ``last_run.steps_requested``. False means the number, if any, does not
    #: describe work this run was asked to do.
    requested: bool
    #: In the params file's ``steps`` array, i.e. what the next run would do.
    in_params: bool

    @property
    def measured(self) -> bool:
        return self.completion.measured

    @property
    def complete(self) -> bool:
        return self.completion.complete

    def detail(self, *, compact: bool = False) -> str:
        """``done/total unit``, or why there is no number.

        The compact form drops the unit noun. It is what the home card uses:
        "100% / 522/522 samples" has a 121px min-content width, and a row that
        cannot shrink below that drags the whole 320px card wider with it.
        """
        c = self.completion
        if not c.measured:
            return "not measured"
        if compact or not self.unit or self.unit == "cohort":
            return f"{c.done}/{c.total}"
        return f"{c.done}/{c.total} {self.unit}"


@dataclass(frozen=True)
class StatusSummary:
    """Headline counts over a list of step rows."""

    complete: int
    measured: int
    unmeasured: int
    total: int

    @property
    def text(self) -> str:
        if self.total == 0:
            return "no steps declared"
        if self.measured == 0:
            return f"{self.total} steps, none measured"
        suffix = f", {self.unmeasured} not measured" if self.unmeasured else ""
        return f"{self.complete}/{self.measured} steps complete{suffix}"

    @property
    def all_complete(self) -> bool:
        return self.measured > 0 and self.complete == self.measured


def step_label(step: str) -> str:
    """Display label for a step, falling back to its raw pipeline name."""
    return STEP_META.get(step, {}).get("label", step)


def _step_sort_key(step: str) -> Tuple[int, str]:
    """Configured steps in configured order; anything new appended by name."""
    try:
        return (STEP_ORDER.index(step), "")
    except ValueError:
        return (len(STEP_ORDER), step)


def build_step_rows(
    state: Optional[CohortState], params_steps: List[str]
) -> List[StepRow]:
    """Build the step table for a cohort.

    With a run record, the rows are what that run was asked to do plus anything
    else it reported. Without one, they are what the params file asks for, all
    unmeasured -- so the intended pipeline shape is visible before the first
    recorded run.
    """
    record = state.last_run if state is not None else None
    requested = record.steps_requested if record is not None else []
    completion = record.completion if record is not None else {}

    steps = set(params_steps) | set(requested) | set(completion)
    rows: List[StepRow] = []
    for step in sorted(steps, key=_step_sort_key):
        meta = STEP_META.get(step, {})
        rows.append(
            StepRow(
                step=step,
                label=meta.get("label", step),
                unit=meta.get("unit", ""),
                note=meta.get("note", ""),
                completion=completion.get(step, StepCompletion(step=step)),
                requested=step in requested,
                in_params=step in params_steps,
            )
        )
    return rows


def summarize_steps(rows: List[StepRow]) -> StatusSummary:
    """Count complete / measured / unmeasured over the step rows."""
    measured = [r for r in rows if r.measured]
    return StatusSummary(
        complete=sum(1 for r in measured if r.complete),
        measured=len(measured),
        unmeasured=len(rows) - len(measured),
        total=len(rows),
    )


def record_status_meta(record: Optional[RunRecord]) -> Dict[str, str]:
    """Chip label, colour and icon for one run record."""
    if record is None:
        return dict(NO_RECORD_STATUS_META)
    meta = RUN_STATUS_META.get(record.status)
    if meta is None:
        # A status the pipeline added after this build shipped.
        return {**DEFAULT_STATUS_META, "label": record.status or "unknown"}
    if record.is_stale_running:
        # Still "running" on paper, but nothing has touched it in a day.
        return {**meta, "label": "running (stale)", "color": "orange"}
    return dict(meta)


def status_meta(state: Optional[CohortState]) -> Dict[str, str]:
    """Chip label, colour and icon for a cohort's last run."""
    if state is None:
        return dict(NO_RECORD_STATUS_META)
    if state.unsupported:
        return {**DEFAULT_STATUS_META, "label": "unreadable record"}
    return record_status_meta(state.last_run)


def format_timestamp(value: Optional[datetime]) -> str:
    """Render a timestamp for display, or an em dash when absent."""
    if value is None:
        return "—"
    return value.strftime("%Y-%m-%d %H:%M")


def format_duration(value: Optional[timedelta]) -> str:
    """Render a run duration as ``3h 12m``, or an em dash when unknown."""
    if value is None:
        return "—"
    seconds = int(value.total_seconds())
    if seconds < 0:
        return "—"
    hours, rem = divmod(seconds, 3600)
    minutes, secs = divmod(rem, 60)
    if hours:
        return f"{hours}h {minutes:02d}m"
    if minutes:
        return f"{minutes}m {secs:02d}s"
    return f"{secs}s"
