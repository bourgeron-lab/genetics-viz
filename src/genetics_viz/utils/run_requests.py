"""Requests for a ghfc-ngs run of a cohort, picked up by ghfc-ngs-swiftbar.

genetics-viz never starts a run. The VM has no cluster account, and collecting
one (a password, or a key per user) would put a credential in the app. A
curator writes a request into the cohort directory instead; the swiftbar
dashboard, polling the cluster as the person who runs the pipeline, shows it
with Preview, Start and Dismiss, and records what it did beside it.

The files (vault: ``Workflow run requests``), in
``cohorts/<C>/.ghfc-ngs.requests/``:

- ``<ID>.json``, the request (this module);
- ``<ID>.withdrawn.json``, taken back before it was handled (this module);
- ``<ID>.handled.json``, ``started`` (job, steps sent) or ``dismissed`` (note),
  written by swiftbar.

Each writer creates only its own files and never edits another's, so nothing
races. ``<ID>`` is ``YYYYMMDDTHHMMSS-<6 of [a-z0-9]>``; both sides check it
before putting it in a path. Every write goes to a temporary file in the same
directory, then a rename.

**One open request per user and cohort**: to change it, withdraw it and make a
new one. Only curators and administrators write or see requests.
"""

from __future__ import annotations

import asyncio
import fcntl
import json
import logging
import os
import re
import secrets
import string
import tempfile
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from . import cohort_state
from .auth import can_write, is_admin

logger = logging.getLogger(__name__)

REQUESTS_DIRNAME = ".ghfc-ngs.requests"
SCHEMA_VERSION = 1
ID_PATTERN = re.compile(r"^\d{8}T\d{6}-[a-z0-9]{6}$")
#: Requests settled longer ago than this leave the Workflow tab's list.
RECENT = timedelta(days=30)
#: Free text is capped here, and again by swiftbar, which only ever shows it.
REASON_MAX = 500
ORIGIN_MAX = 120
FAMILIES_MAX = 50

#: ``RunRequest.state`` values.
OPEN, WITHDRAWN, DISMISSED, STARTED = "open", "withdrawn", "dismissed", "started"


@dataclass(frozen=True)
class Handled:
    """What swiftbar did with a request (``<ID>.handled.json``)."""

    status: str  # "started" | "dismissed"
    user: str = ""
    at: Optional[datetime] = None
    note: str = ""
    job_id: Optional[str] = None
    #: The steps actually sent (``None``: a normal run), when started.
    steps_sent: Optional[Tuple[str, ...]] = None
    clean_sent: bool = False
    #: Started by another request's Start, whose run covered this one.
    covered_by: Optional[str] = None


@dataclass(frozen=True)
class RunRequest:
    """A curator's request for a run of one cohort, with what became of it."""

    id: str
    cohort: str
    #: ``None``: a normal run, the steps of the cohort's params file.
    steps: Optional[Tuple[str, ...]]
    clean_stale_families: bool
    reason: str
    requested_by: str
    requested_at: Optional[datetime]
    families: Tuple[str, ...] = ()
    origin: str = ""
    withdrawn_by: str = ""
    withdrawn_at: Optional[datetime] = None
    handled: Optional[Handled] = None
    #: The run started for it, from the run record: ``running``, ``success``,
    #: ``failed``…; ``None`` until the record names its job.
    outcome: Optional[str] = None

    @property
    def state(self) -> str:
        if self.withdrawn_at is not None or self.withdrawn_by:
            return WITHDRAWN
        if self.handled is not None:
            return DISMISSED if self.handled.status == "dismissed" else STARTED
        return OPEN

    @property
    def is_open(self) -> bool:
        return self.state == OPEN

    @property
    def settled_at(self) -> Optional[datetime]:
        if self.state == WITHDRAWN:
            return self.withdrawn_at
        if self.handled is not None:
            return self.handled.at
        return None

    @property
    def edited(self) -> bool:
        """Started with other steps, or another clean, than asked for."""
        h = self.handled
        return (
            h is not None
            and h.status == "started"
            and (
                h.steps_sent != self.steps or h.clean_sent != self.clean_stale_families
            )
        )


def steps_text(steps: Optional[Sequence[str]]) -> str:
    """``sampleqc, genomeqc``, or ``normal run`` for ``None``."""
    return "normal run" if steps is None else ", ".join(steps)


def command_line(cohort: str, steps: Optional[Sequence[str]], clean: bool) -> str:
    """The ghfc-ngs command swiftbar will send, as shown to the curator."""
    parts = ["ghfc-ngs", cohort, "--submit"]
    if steps is not None:
        parts += ["--steps", ",".join(steps)]
    if clean:
        parts.append("--clean-stale-families")
    return " ".join(parts)


# ---------------------------------------------------------------------------
# Where requests are allowed, and which steps
# ---------------------------------------------------------------------------


def requests_enabled(data_dir: Path) -> bool:
    """Whether Admin → Directories turned run requests on for ``data_dir``.

    Off without a loaded config (unit tests) and for a directory not in it.
    """
    from ..config_model import get_config, is_config_loaded

    if not is_config_loaded():
        return False
    target = _resolved(data_dir)
    return any(
        d.run_requests and _resolved(Path(d.path)) == target
        for d in get_config().data_directories
    )


def _resolved(path: Path) -> str:
    try:
        return str(Path(path).resolve())
    except OSError:
        return str(path)


def allowed_steps(cohort: Any) -> List[str]:
    """The steps a request may name: those of the cohort's params file."""
    from .pipeline_params import get_requested_steps, load_params

    params_file = getattr(cohort, "params_file", None)
    if params_file is None:
        return []
    steps: List[str] = []
    for item in get_requested_steps(load_params(params_file)):
        # ghfc-ngs accepts "a,b" as well as a list
        steps += [s.strip() for s in item.split(",") if s.strip()]
    return list(dict.fromkeys(steps))


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def requests_dir(cohort_dir: Path) -> Path:
    return Path(cohort_dir) / REQUESTS_DIRNAME


def _read_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None
    except (OSError, ValueError):
        logger.warning("Unreadable run request file %s", path, exc_info=True)
        return None
    return raw if isinstance(raw, dict) else None


def _time(value: Any) -> Optional[datetime]:
    if not isinstance(value, str) or not value:
        return None
    try:
        moment = datetime.fromisoformat(value)
    except ValueError:
        return None
    return moment if moment.tzinfo else moment.astimezone()


def _text(value: Any, limit: int) -> str:
    return str(value)[:limit] if isinstance(value, (str, int, float)) else ""


def _steps(value: Any) -> Optional[Tuple[str, ...]]:
    if not isinstance(value, list):
        return None
    return tuple(str(s) for s in value if isinstance(s, str) and s)


def _handled(raw: Optional[Dict[str, Any]]) -> Optional[Handled]:
    if raw is None or raw.get("status") not in ("started", "dismissed"):
        return None
    job = raw.get("job_id")
    return Handled(
        status=raw["status"],
        user=_text(raw.get("user"), 64),
        at=_time(raw.get("at")),
        note=_text(raw.get("note"), REASON_MAX),
        job_id=str(job) if isinstance(job, (str, int)) and str(job) else None,
        steps_sent=_steps(raw.get("steps_sent")),
        clean_sent=raw.get("clean_sent") is True,
        covered_by=_text(raw.get("covered_by"), 64) or None,
    )


def _outcome(
    job_id: Optional[str], state: Optional[cohort_state.CohortState]
) -> Optional[str]:
    """The status of the run whose Slurm job is ``job_id``, from the run record."""
    if not job_id or state is None or state.unsupported:
        return None
    for record in [state.last_run, *state.history]:
        if record is not None and record.slurm_job_id == job_id:
            if record.is_running and record.is_stale_running:
                return "died"
            return record.status
    return None


def read_requests(
    cohort_dir: Path, state: Optional[cohort_state.CohortState] = None
) -> List[RunRequest]:
    """Every request of the cohort, newest first. Blocking (a few small files).

    ``state`` is the cohort's run record, read when not given, for the outcome
    of a started request.
    """
    folder = requests_dir(cohort_dir)
    try:
        names = sorted(p.name for p in folder.iterdir())
    except OSError:
        return []
    if state is None:
        state = cohort_state.read_cohort_state(Path(cohort_dir))
    out: List[RunRequest] = []
    for name in names:
        if not name.endswith(".json") or name.count(".") != 1:
            continue
        rid = name[: -len(".json")]
        if not ID_PATTERN.match(rid):
            continue
        raw = _read_json(folder / name)
        if raw is None or raw.get("schema_version") != SCHEMA_VERSION:
            continue
        withdrawn = _read_json(folder / f"{rid}.withdrawn.json") or {}
        handled = _handled(_read_json(folder / f"{rid}.handled.json"))
        families = raw.get("families") if isinstance(raw.get("families"), list) else []
        out.append(
            RunRequest(
                id=rid,
                cohort=Path(cohort_dir).name,
                steps=_steps(raw.get("steps")),
                clean_stale_families=raw.get("clean_stale_families") is True,
                reason=_text(raw.get("reason"), REASON_MAX),
                requested_by=_text(raw.get("requested_by"), 64),
                requested_at=_time(raw.get("requested_at")),
                families=tuple(_text(f, 64) for f in families[:FAMILIES_MAX]),
                origin=_text(raw.get("origin"), ORIGIN_MAX),
                withdrawn_by=_text(withdrawn.get("user"), 64),
                withdrawn_at=_time(withdrawn.get("at")),
                handled=handled,
                outcome=_outcome(handled.job_id if handled else None, state),
            )
        )
    out.sort(key=lambda r: r.id, reverse=True)
    return out


def recent(
    requests: Sequence[RunRequest], now: Optional[datetime] = None
) -> List[RunRequest]:
    """Open requests, and those settled in the last :data:`RECENT`."""
    now = now or datetime.now(timezone.utc)
    return [
        r
        for r in requests
        if r.is_open or r.settled_at is None or now - r.settled_at <= RECENT
    ]


def count_open_requests(cohort_dir: Path) -> int:
    """How many requests of the cohort wait for the person who runs ghfc-ngs."""
    return sum(r.is_open for r in read_requests(cohort_dir, state=None))


def open_request_of(requests: Sequence[RunRequest], user: str) -> Optional[RunRequest]:
    """The user's open request, if any (there is at most one per cohort)."""
    return next((r for r in requests if r.is_open and r.requested_by == user), None)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def _check_role() -> None:
    # Checked here as well as on the control, which a crafted event bypasses.
    if not can_write():
        raise PermissionError("Only curators and administrators can request a run")


def new_id(now: Optional[datetime] = None) -> str:
    now = now or datetime.now().astimezone()
    tail = "".join(
        secrets.choice(string.ascii_lowercase + string.digits) for _ in range(6)
    )
    return f"{now.strftime('%Y%m%dT%H%M%S')}-{tail}"


def _write_json(folder: Path, name: str, content: Dict[str, Any]) -> None:
    """Write ``folder/name`` whole or not at all: a temporary file, then a rename."""
    fd, tmp = tempfile.mkstemp(dir=folder, prefix=".tmp.", suffix=".json")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            json.dump(content, f, indent=2)
            f.write("\n")
        os.chmod(tmp, 0o664)
        os.replace(tmp, folder / name)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


class _Lock:
    """One request written at a time per cohort, so two clicks cannot both
    pass the one-open-request check."""

    def __init__(self, folder: Path) -> None:
        self.path = folder / ".lock"

    def __enter__(self) -> "_Lock":
        self.f = open(self.path, "a")
        fcntl.flock(self.f.fileno(), fcntl.LOCK_EX)
        return self

    def __exit__(self, *exc: Any) -> None:
        fcntl.flock(self.f.fileno(), fcntl.LOCK_UN)
        self.f.close()


def _write_request(
    cohort: Any,
    steps: Optional[Sequence[str]],
    clean: bool,
    reason: str,
    families: Sequence[str],
    origin: str,
    user: str,
    now: Optional[datetime] = None,
) -> RunRequest:
    reason = reason.strip()
    if not reason:
        raise ValueError("Say why the run is needed")
    if not user:
        raise PermissionError("No user")
    if steps is not None:
        allowed = allowed_steps(cohort)
        steps = [s for s in dict.fromkeys(steps)]
        if not steps:
            raise ValueError("Choose at least one step, or a normal run")
        unknown = [s for s in steps if s not in allowed]
        if unknown:
            raise ValueError(
                f"Not a step of {cohort.name}'s params file: {', '.join(unknown)}"
            )
    cohort_dir = Path(cohort.pedigree_file).parent
    folder = requests_dir(cohort_dir)
    folder.mkdir(exist_ok=True)
    now = now or datetime.now().astimezone()
    with _Lock(folder):
        mine = open_request_of(read_requests(cohort_dir), user)
        if mine is not None:
            raise ValueError(
                f"You already have an open request for {cohort.name}"
                f" ({steps_text(mine.steps)}). Withdraw it first."
            )
        rid = new_id(now)
        content = {
            "schema_version": SCHEMA_VERSION,
            "id": rid,
            "cohort": cohort.name,
            "steps": list(steps) if steps is not None else None,
            "clean_stale_families": bool(clean),
            "reason": reason[:REASON_MAX],
            "families": [str(f)[:64] for f in families][:FAMILIES_MAX],
            "origin": origin[:ORIGIN_MAX],
            "requested_by": user,
            "requested_at": now.isoformat(timespec="seconds"),
        }
        _write_json(folder, f"{rid}.json", content)
    written = next(r for r in read_requests(cohort_dir) if r.id == rid)
    return written


def write_request(
    cohort: Any,
    steps: Optional[Sequence[str]],
    *,
    clean: bool = False,
    reason: str,
    user: str,
    families: Sequence[str] = (),
    origin: str = "",
    now: Optional[datetime] = None,
) -> RunRequest:
    """Ask for a run of ``cohort``: ``steps`` among its params steps, ``None``
    for a normal run. Refused to viewers, without a reason, and while the
    user has an open request for the cohort."""
    _check_role()
    return _write_request(cohort, steps, clean, reason, families, origin, user, now)


async def write_request_async(
    cohort: Any,
    steps: Optional[Sequence[str]],
    *,
    clean: bool = False,
    reason: str,
    user: str,
    families: Sequence[str] = (),
    origin: str = "",
) -> RunRequest:
    """:func:`write_request` for a page: the role is checked on the event loop
    (it reads the session), the files are written in a thread."""
    _check_role()
    return await asyncio.to_thread(
        _write_request, cohort, steps, clean, reason, families, origin, user, None
    )


def _withdraw(cohort_dir: Path, request_id: str, user: str, admin: bool) -> None:
    if not ID_PATTERN.match(request_id):
        raise ValueError("Not a request")
    folder = requests_dir(cohort_dir)
    with _Lock(folder):
        request = next(
            (r for r in read_requests(cohort_dir) if r.id == request_id), None
        )
        if request is None:
            raise ValueError("No such request")
        if not request.is_open:
            raise ValueError(f"This request is already {request.state}")
        if request.requested_by != user and not admin:
            raise PermissionError("Only its author or an administrator can withdraw it")
        _write_json(
            folder,
            f"{request_id}.withdrawn.json",
            {
                "user": user,
                "at": datetime.now().astimezone().isoformat(timespec="seconds"),
            },
        )


def withdraw_request(cohort_dir: Path, request_id: str, *, user: str) -> None:
    """Take back an open request: its author's, or any for an administrator."""
    _check_role()
    _withdraw(Path(cohort_dir), request_id, user, is_admin())


async def withdraw_request_async(
    cohort_dir: Path, request_id: str, *, user: str
) -> None:
    _check_role()
    admin = is_admin()
    await asyncio.to_thread(_withdraw, Path(cohort_dir), request_id, user, admin)
