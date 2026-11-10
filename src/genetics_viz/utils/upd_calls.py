"""UPD calls from the ghfc-ngs ``upd`` step: discovery, reading, wording.

The step analyses every sequenced member in the best mode its family allows
(trio, duo, or solo from the ``roh`` step's runs) and writes, under ``upd/`` of
each level, files named by its label (``upd_label``: ``upd_v1`` … ``upd_v4``):

- family: ``<FID>.<label>.{segments,events,members,inputs}.tsv`` and, from
  ``upd_v4``, ``<FID>.<label>.windows.tsv``;
- cohort: ``<NAME>.<label>.{segments,events,members,inputs}.tsv``, the family
  tables concatenated, with ``recurrent_families`` and ``recurrent`` on the
  events (no windows).

A new label is a new version of the analysis, written beside the older ones:
the newest-written label is shown by default (:func:`find_upd_labels`).

``windows.tsv`` holds, per analysed member and 500 kb window, the counts the
chromosome view draws: the child's het and hom-alt calls, the alleles only one
parent has (``n_mat``, ``n_pat``) and the Mendelian errors against each parent
(``n_lacks_father``, ``n_lacks_mother``), plus the depth ratio. Raw Mendelian
classes, not the HMM's states: the view checks the caller, it does not redraw it.

Every table writes ``.`` for a missing value. The pipeline grades everything;
the app only words and draws it (``config/upd.yaml``). The curator side
(validations, evidence, candidates) is :mod:`genetics_viz.utils.upd`.
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import polars as pl
import yaml

from .ancestry import NULL_VALUES
from .sharding import get_family_path
from .tsv import find_labels, read_tsv_or_none

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "upd.yaml"

UPD_SUBDIR = "upd"

#: ``upd_segments``'s column lists. Older labels lack some (``upd_v1`` has no
#: ``support_fraction`` nor ``shared_with``): missing ones are added as nulls.
SEGMENT_COLUMNS = [
    "family_id", "sample", "mode", "father", "mother", "chrom", "start", "end",
    "length_mb", "state", "origin", "source", "n_sites", "n_informative",
    "n_supporting", "n_contradicting", "n_strong", "support_fraction", "lod",
    "depth_ratio", "roh_fraction", "interpretation", "confidence", "event_id",
]  # fmt: skip
EVENT_COLUMNS = [
    "family_id", "sample", "mode", "event_id", "chrom", "start", "end",
    "length_mb", "origin", "composition", "extent", "terminal", "n_segments",
    "n_supporting", "n_strong", "support_fraction", "lod", "depth_ratio",
    "interpretation", "imprinting", "confidence", "shared_with", "notes",
]  # fmt: skip
#: Only the cohort's events table has them.
COHORT_EVENT_COLUMNS = ["recurrent_families", "recurrent"]
MEMBER_COLUMNS = [
    "family_id", "sample", "father", "mother", "mode", "father_used",
    "mother_used", "chrx", "n_sites", "error_rate", "status", "n_events",
    "n_events_high",
]  # fmt: skip
WINDOW_COLUMNS = [
    "family_id", "sample", "mode", "chrom", "start", "end", "n_sites", "n_het",
    "n_hom_alt", "n_mat", "n_pat", "n_lacks_father", "n_lacks_mother",
    "depth_ratio",
]  # fmt: skip

_COLUMNS = {
    "segments": SEGMENT_COLUMNS,
    "events": EVENT_COLUMNS + COHORT_EVENT_COLUMNS,
    "members": MEMBER_COLUMNS,
    "windows": WINDOW_COLUMNS,
}

#: Identifiers and codes, read as strings. ``family_id`` is numeric-looking
#: (``88171200009``), and so is ``chrom``.
_STRING_COLUMNS = {
    "family_id", "sample", "father", "mother", "father_used", "mother_used",
    "mode", "chrom", "state", "origin", "source", "interpretation",
    "confidence", "event_id", "composition", "extent", "terminal",
    "imprinting", "shared_with", "notes", "recurrent", "chrx", "status",
}  # fmt: skip

_INT_COLUMNS = {
    "start", "end", "n_sites", "n_informative", "n_supporting",
    "n_contradicting", "n_strong", "n_segments", "n_events", "n_events_high",
    "recurrent_families", "n_het", "n_hom_alt", "n_mat", "n_pat",
    "n_lacks_father", "n_lacks_mother",
}  # fmt: skip


# ---------------------------------------------------------------------------
# Display configuration
# ---------------------------------------------------------------------------


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


CONFIG: Dict[str, Any] = _load_config()
CONFIDENCES: Dict[str, Dict[str, Any]] = CONFIG.get("confidences", {})
DEFAULT_CONFIDENCES: List[str] = list(CONFIG.get("default_confidences", []))
ORIGINS: Dict[str, Dict[str, str]] = CONFIG.get("origins", {})
COMPOSITIONS: Dict[str, Dict[str, str]] = CONFIG.get("compositions", {})
EXTENTS: Dict[str, Dict[str, str]] = CONFIG.get("extents", {})
MODES: Dict[str, Dict[str, str]] = CONFIG.get("modes", {})
STATUSES: Dict[str, Dict[str, Any]] = CONFIG.get("statuses", {})
WINDOW_STYLE: Dict[str, Any] = CONFIG.get("windows", {})
VALIDATION_FORM: Dict[str, List[str]] = CONFIG.get("validation", {})


def confidence_info(confidence: Optional[str]) -> Dict[str, Any]:
    """Label, colour and rank of a confidence; an unknown one is raw, grey, rank 0."""
    cfg = CONFIDENCES.get(confidence or "", {})
    return {
        "key": confidence or "",
        "label": cfg.get("label", confidence or "—"),
        "color": cfg.get("color", "#9e9e9e"),
        "rank": int(cfg.get("rank", 0)),
    }


def origin_info(origin: Optional[str]) -> Dict[str, str]:
    cfg = ORIGINS.get(origin or "unknown", {})
    return {
        "label": cfg.get("label", origin or "unknown"),
        "color": cfg.get("color", "#757575"),
    }


def composition_label(code: Optional[str]) -> str:
    return COMPOSITIONS.get(code or "", {}).get("label", code or "")


def extent_label(code: Optional[str]) -> str:
    return EXTENTS.get(code or "", {}).get("label", (code or "").replace("_", " "))


def mode_info(mode: Optional[str]) -> Dict[str, str]:
    cfg = MODES.get(mode or "", {})
    return {"label": cfg.get("label", mode or "—"), "help": cfg.get("help", "")}


def member_status_label(status: Optional[str], show_detail: bool) -> str:
    """A member's ``status``, worded for the role.

    ``status`` holds one or more codes joined by ``"; "``, each optionally
    followed by ``": <detail>"`` (``parentage: S1 excluded; no_roh_tables``).
    Readers get each code's ``reader_label``; a code whose reader label is
    null (a declared parent excluded, the parents' kinship) is left out, and
    nothing left reads "Analysed".
    """
    parts = [p.strip() for p in (status or "").split(";") if p.strip()]
    words: List[str] = []
    for part in parts:
        code, _, detail = part.partition(":")
        code, detail = code.strip(), detail.strip()
        cfg = STATUSES.get(code)
        if show_detail:
            label = cfg.get("label", code) if cfg else code
            words.append(f"{label}: {detail}" if detail else label)
        elif cfg is None:
            words.append(code)
        elif cfg.get("reader_label"):
            words.append(cfg["reader_label"])
    if not words:
        analysed = STATUSES.get("analysed", {})
        return analysed.get("label" if show_detail else "reader_label") or "Analysed"
    return "; ".join(dict.fromkeys(words))


def normalise_chrom(chrom: Any) -> str:
    """``"chr8"``, ``"8"`` and ``8`` are all ``"8"``."""
    text = str(chrom if chrom is not None else "").strip()
    return text[3:] if text.lower().startswith("chr") else text


def span_text(start: Any, end: Any) -> str:
    """``80.0–101.9 Mb``; ``""`` when either end is missing."""
    try:
        return f"{int(start) / 1e6:.1f}–{int(end) / 1e6:.1f} Mb"
    except (TypeError, ValueError):
        return ""


def event_type(event: Dict[str, Any]) -> str:
    """``isodisomy · segmental · terminal`` of an event row."""
    parts = [composition_label(event.get("composition"))]
    parts.append(extent_label(event.get("extent")))
    if event.get("terminal") == "yes":
        parts.append("terminal")
    return " · ".join(p for p in parts if p)


def confidence_rank(event: Dict[str, Any]) -> int:
    return confidence_info(event.get("confidence"))["rank"]


def best_event(events: Iterable[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """The most confident event, then the highest LOD."""
    ranked = sorted(
        events, key=lambda e: (-confidence_rank(e), -(_float(e.get("lod")) or 0.0))
    )
    return ranked[0] if ranked else None


def _float(value: Any) -> Optional[float]:
    try:
        return None if value is None else float(value)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def get_family_upd_dir(data_dir: Path, family_id: str) -> Path:
    return get_family_path(data_dir, family_id) / UPD_SUBDIR


def get_cohort_upd_dir(cohort: Any) -> Path:
    return cohort.path / UPD_SUBDIR


def upd_path(directory: Path, prefix: str, label: str, kind: str) -> Path:
    """``<directory>/<prefix>.<label>.<kind>.tsv``."""
    return directory / f"{prefix}.{label}.{kind}.tsv"


def find_upd_labels(directory: Path, prefix: str) -> List[str]:
    """Labels with an events table in ``directory``, the newest written first."""
    return find_labels(directory, prefix, "events")


def probe_upd_data(data_dir: Path, family_id: str) -> bool:
    """Whether the family has upd step results."""
    return bool(find_upd_labels(get_family_upd_dir(data_dir, family_id), family_id))


def probe_cohort_upd_data(cohort: Any) -> bool:
    """Whether the cohort has upd step results."""
    return bool(find_upd_labels(get_cohort_upd_dir(cohort), cohort.name))


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_upd_table(path: Optional[Path], kind: str) -> Optional[pl.DataFrame]:
    """Read one upd TSV, ``.`` as null and identifiers as strings.

    Columns the file lacks (an older label, a stub run) are added as nulls.
    ``None`` when the path is ``None``, missing, or empty.
    """
    if path is None or not path.exists():
        return None
    df = read_tsv_or_none(
        path,
        null_values=NULL_VALUES,
        infer_schema_length=10000,
        string_columns=_STRING_COLUMNS,
    )
    if df is None:
        return None
    missing = [c for c in _COLUMNS.get(kind, []) if c not in df.columns]
    if missing:
        df = df.with_columns([pl.lit(None).alias(c) for c in missing])
    if "chrom" in df.columns:
        df = df.with_columns(pl.col("chrom").cast(pl.Utf8).str.replace(r"(?i)^chr", ""))
    return df


def _rows(df: Optional[pl.DataFrame]) -> List[Dict[str, Any]]:
    if df is None or len(df) == 0:
        return []
    rows = df.to_dicts()
    for row in rows:
        for column in _INT_COLUMNS:
            value = row.get(column)
            if isinstance(value, float) and value == value:
                row[column] = int(value)
    return rows


#: ``(path, mtime, size)`` -> rows; a rewrite by the next run is picked up.
_table_cache: Dict[Tuple[str, float, int], List[Dict[str, Any]]] = {}
_window_cache: Dict[Tuple[str, float, int], Optional[pl.DataFrame]] = {}


def _stat_key(path: Path) -> Optional[Tuple[str, float, int]]:
    try:
        stat = path.stat()
    except OSError:
        return None
    return (str(path), stat.st_mtime, stat.st_size)


def _read_cached(path: Path, kind: str) -> List[Dict[str, Any]]:
    key = _stat_key(path)
    if key is None:
        return []
    if key not in _table_cache:
        for old in [k for k in _table_cache if k[0] == key[0]]:
            del _table_cache[old]
        try:
            _table_cache[key] = _rows(read_upd_table(path, kind))
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            _table_cache[key] = []
    return _table_cache[key]


def clear_upd_cache() -> None:
    """Drop every cached table (tests, and after a data directory switch)."""
    _table_cache.clear()
    _window_cache.clear()


@dataclass
class UpdData:
    """The upd tables of one family or cohort, for one label."""

    prefix: str
    directory: Path
    label: str
    labels: List[str]
    events: List[Dict[str, Any]] = field(default_factory=list)
    segments: List[Dict[str, Any]] = field(default_factory=list)
    members: List[Dict[str, Any]] = field(default_factory=list)

    def path(self, kind: str) -> Path:
        return upd_path(self.directory, self.prefix, self.label, kind)

    @property
    def has_windows(self) -> bool:
        """``upd_v4`` and later write a windows table (the family level only)."""
        return self.path("windows").is_file()

    def member(self, sample: str) -> Optional[Dict[str, Any]]:
        return next((m for m in self.members if m.get("sample") == sample), None)

    def events_of(self, sample: str, chrom: Any = None) -> List[Dict[str, Any]]:
        wanted = normalise_chrom(chrom) if chrom is not None else None
        return [
            e
            for e in self.events
            if e.get("sample") == sample
            and (wanted is None or e.get("chrom") == wanted)
        ]

    def segments_of(self, sample: str, chrom: Any = None) -> List[Dict[str, Any]]:
        wanted = normalise_chrom(chrom) if chrom is not None else None
        return [
            s
            for s in self.segments
            if s.get("sample") == sample
            and (wanted is None or s.get("chrom") == wanted)
        ]


def load_upd(
    directory: Path, prefix: str, label: Optional[str] = None
) -> Optional[UpdData]:
    """The upd tables in ``directory`` for ``label`` (default: the newest).

    Blocking I/O over the data mount; the tables are cached on their stat.
    ``None`` when there is no events table at all.
    """
    labels = find_upd_labels(directory, prefix)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    return UpdData(
        prefix=prefix,
        directory=directory,
        label=chosen,
        labels=labels,
        events=_read_cached(upd_path(directory, prefix, chosen, "events"), "events"),
        segments=_read_cached(
            upd_path(directory, prefix, chosen, "segments"), "segments"
        ),
        members=_read_cached(upd_path(directory, prefix, chosen, "members"), "members"),
    )


def load_family_upd(
    data_dir: Path, family_id: str, label: Optional[str] = None
) -> Optional[UpdData]:
    return load_upd(get_family_upd_dir(data_dir, family_id), family_id, label)


def load_cohort_upd(cohort: Any, label: Optional[str] = None) -> Optional[UpdData]:
    return load_upd(get_cohort_upd_dir(cohort), cohort.name, label)


def load_windows(data: UpdData) -> Optional[pl.DataFrame]:
    """The label's windows table; ``None`` before ``upd_v4``. Blocking, cached."""
    path = data.path("windows")
    key = _stat_key(path)
    if key is None:
        return None
    if key not in _window_cache:
        for old in [k for k in _window_cache if k[0] == key[0]]:
            del _window_cache[old]
        try:
            _window_cache[key] = read_upd_table(path, "windows")
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            _window_cache[key] = None
    return _window_cache[key]


def windows_of(
    windows: Optional[pl.DataFrame], sample: str, chrom: Any = None
) -> List[Dict[str, Any]]:
    """One member's windows, optionally of one chromosome, in position order."""
    if windows is None or len(windows) == 0:
        return []
    sub = windows.filter(pl.col("sample") == sample)
    if chrom is not None:
        sub = sub.filter(pl.col("chrom") == normalise_chrom(chrom))
    return _rows(sub.sort(["chrom", "start"]))


# ---------------------------------------------------------------------------
# Exports, shared by the family tab and the cohort panel
# ---------------------------------------------------------------------------

#: Member columns readers do not get: which declared parent was used says
#: when one was excluded as not the parent.
_CURATOR_MEMBER_COLUMNS = {"father_used", "mother_used"}


def _identity(row: Dict[str, Any]) -> Dict[str, Any]:
    return {"FID": row.get("family_id") or "", "IID": row.get("sample") or ""}


def table_export(
    rows: Sequence[Dict[str, Any]], kind: str, show_detail: bool
) -> Tuple[List[Dict[str, Any]], List[str]]:
    """Rows and columns of one upd table, ``FID`` and ``IID`` first.

    ``kind`` is ``events``, ``segments``, ``members`` or ``windows``. Readers
    get the member status in their wording and not the parents used.
    """
    base = [c for c in _COLUMNS[kind] if c not in ("family_id", "sample")]
    # Only the cohort's events table has them; a family's reads them as nulls.
    if kind == "events" and not any(r.get("recurrent") is not None for r in rows):
        base = [c for c in base if c not in COHORT_EVENT_COLUMNS]
    if kind == "members" and not show_detail:
        base = [c for c in base if c not in _CURATOR_MEMBER_COLUMNS]
    out = []
    for row in rows:
        item = {**_identity(row), **{c: row.get(c) for c in base}}
        if kind == "members":
            item["status"] = member_status_label(row.get("status"), show_detail)
        out.append(item)
    return out, ["FID", "IID", *base]
