"""GeneTrek: a snapshot of its gene table, fetched from its API and cached.

GeneTrek (``genetrek_url`` in the application config) is the lab's gene
database. The app uses it for geneset membership (badge colours, tooltips, the
Search geneset filter) and for modes of inheritance (the MOI flag).

It is read **as a snapshot**, not per gene. Tables render thousands of gene
badges whose tooltips are computed in Python before the table is sent, so a
call per gene is out of the question. The whole table is one ~4 MB parquet
(``/api/download/database``), downloaded once per GeneTrek *release* and cached
on disk; between releases a refresh is one ``/api/health`` call. GeneTrek never
learns which genes a user looks at, only which one they click.

Failure is a degraded state, never an exception: with no snapshot at all the
badges are neutral and the profile page says why. A cached release survives
restarts with the network down.

Two traps worth knowing:

- Symbols in VEP output lag HGNC. ``KIAA1109`` is ``BLTP1`` now, so resolution
  falls back to previous and alias symbols, kept only when they name exactly one
  gene. An exact current symbol always wins.
- Geneset versions coexist (``hc_ndd_v1`` … ``hc_ndd_v7``). Their alias
  (``hc_ndd``) names the current one, so saved preferences use the alias where
  there is one: :attr:`GenesetInfo.key`.
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import tempfile
import threading
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, FrozenSet, Iterable, List, Optional, Tuple

import polars as pl
import yaml

from genetics_viz.config_model import get_genetrek_cache_dir, get_genetrek_url

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "genetrek.yaml"
with open(_CONFIG_PATH) as _f:
    GENETREK_CONFIG: Dict[str, Any] = yaml.safe_load(_f) or {}

_TIMEOUT_SECONDS = 10
#: The parquet is ~4 MB; three seconds on the lab network, so a generous cap.
_DOWNLOAD_TIMEOUT_SECONDS = 60
_USER_AGENT = "genetics-viz"

#: Identity columns GeneTrek has no id role for, read by name.
_PREV_COLUMN = "hgnc_prev_symbol"
_ALIAS_COLUMN = "hgnc_alias_symbol"
_NAME_COLUMN = "hgnc_name"
#: Per-gene values shown in the badge tooltip, when the table has them.
_INFO_COLUMNS = {"loeuf": "LOEUF", "shet_mean": "s-het", "eagle_score": "EAGLE"}
#: GeneTrek's category for mode-of-inheritance fields.
INHERITANCE_CATEGORY = "Inheritance"

_META_FILE = "fields.json"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class GenesetInfo:
    """One boolean GeneTrek column: a gene list."""

    name: str
    fullname: str = ""
    description: str = ""
    category: str = ""
    subcategory: str = ""
    visibility: str = "secondary"
    aliases: Tuple[str, ...] = ()
    derived_from: Tuple[str, ...] = ()
    n_genes: int = 0

    @property
    def key(self) -> str:
        """The name preferences are saved under: the alias if there is one.

        An alias follows the current version of a list; the versioned name
        would stay pinned to an old one.
        """
        return self.aliases[0] if self.aliases else self.name

    @property
    def label(self) -> str:
        return self.fullname or self.name

    @property
    def is_primary(self) -> bool:
        return self.visibility == "primary"

    @property
    def is_rollup(self) -> bool:
        return bool(self.derived_from)


@dataclass(frozen=True)
class FieldInfo:
    """A non-boolean GeneTrek column the app reads (inheritance sources)."""

    name: str
    fullname: str = ""
    category: str = ""


@dataclass(frozen=True)
class GeneInfo:
    """What the tooltip says about a gene beyond its genesets."""

    symbol: str
    name: str = ""
    #: ``(label, value)`` pairs for the tooltip, in display order.
    values: Tuple[Tuple[str, str], ...] = ()
    #: ``(source fullname, raw mode string)`` for each inheritance field.
    inheritance: Tuple[Tuple[str, str], ...] = ()


@dataclass
class GenetrekSnapshot:
    """An indexed copy of GeneTrek's gene table, or the reason there is none."""

    release: Optional[str] = None
    fetched_at: Optional[datetime] = None
    checked_at: Optional[datetime] = None
    error: Optional[str] = None
    #: True until the first load attempt has finished.
    loading: bool = False
    catalog: Tuple[GenesetInfo, ...] = ()
    inheritance_fields: Tuple[FieldInfo, ...] = ()
    #: symbol -> names of the genesets it belongs to.
    members: Dict[str, FrozenSet[str]] = field(default_factory=dict)
    #: geneset name -> its member symbols.
    set_members: Dict[str, FrozenSet[str]] = field(default_factory=dict)
    info: Dict[str, GeneInfo] = field(default_factory=dict)
    #: Upper-cased identifier -> current symbol, one map per layer.
    by_symbol: Dict[str, str] = field(default_factory=dict)
    by_ensg: Dict[str, str] = field(default_factory=dict)
    by_prev: Dict[str, str] = field(default_factory=dict)
    by_alias: Dict[str, str] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.release is not None and bool(self.by_symbol)

    def status_text(self) -> str:
        """One line for the profile page."""
        if self.ok:
            when = self.fetched_at.strftime("%Y-%m-%d %H:%M") if self.fetched_at else ""
            text = f"Release {self.release}, {len(self.by_symbol):,} genes, {len(self.catalog)} genesets"
            if when:
                text += f" (downloaded {when})"
            if self.error:
                text += f". Last check failed: {self.error}"
            return text
        if self.loading:
            return "Loading GeneTrek data…"
        return f"GeneTrek data unavailable: {self.error or 'not loaded'}"

    def resolve(self, identifier: Optional[str]) -> Optional[str]:
        """Map a symbol, ENSG, previous or alias symbol to the current symbol."""
        if not identifier:
            return None
        ident = identifier.strip().upper()
        if not ident:
            return None
        if ident in self.by_symbol:
            return self.by_symbol[ident]
        if ident.startswith("ENSG"):
            return self.by_ensg.get(ident.split(".", 1)[0])
        return self.by_prev.get(ident) or self.by_alias.get(ident)

    def geneset(self, key: str) -> Optional[GenesetInfo]:
        """Look a geneset up by name or alias."""
        return self._sets_by_key().get(key)

    def _sets_by_key(self) -> Dict[str, GenesetInfo]:
        cached = self.__dict__.get("_by_key")
        if cached is None:
            cached = {}
            for gs in self.catalog:
                cached[gs.name] = gs
                for alias in gs.aliases:
                    cached.setdefault(alias, gs)
            self.__dict__["_by_key"] = cached
        return cached

    def matching_identifiers(self, symbols: Iterable[str]) -> FrozenSet[str]:
        """Every upper-cased identifier that resolves into ``symbols``.

        Used to filter tables natively: a cell holding the old ``KIAA1109``
        should match a set containing ``BLTP1``, without resolving every cell.
        """
        wanted = {s.upper() for s in symbols}
        out = set(wanted)
        for layer in (self.by_prev, self.by_alias):
            out.update(k for k, v in layer.items() if v.upper() in wanted)
        return frozenset(out)


# ---------------------------------------------------------------------------
# Index building
# ---------------------------------------------------------------------------


def _split_pipe(value: Optional[str]) -> List[str]:
    if not value:
        return []
    return [v.strip() for v in str(value).split("|") if v.strip()]


def _unambiguous(
    pairs: Iterable[Tuple[str, str]], current: Dict[str, str]
) -> Dict[str, str]:
    """Keep the identifiers that name exactly one gene and no current symbol."""
    seen: Dict[str, set] = {}
    for ident, symbol in pairs:
        seen.setdefault(ident.upper(), set()).add(symbol)
    return {
        ident: next(iter(syms))
        for ident, syms in seen.items()
        if len(syms) == 1 and ident not in current
    }


def _column_values(series: pl.Series) -> List[Any]:
    """Python values, categoricals as plain strings, numbers left as numbers."""
    if series.dtype in (pl.Categorical, pl.Enum) or isinstance(series.dtype, pl.Enum):
        series = series.cast(pl.String)
    return series.to_list()


def _format_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.3g}"
    return str(value)


def build_snapshot(
    df: pl.DataFrame,
    fields: List[Dict[str, Any]],
    *,
    release: str,
    fetched_at: Optional[datetime] = None,
) -> GenetrekSnapshot:
    """Index GeneTrek's table and field list into a snapshot."""
    symbol_col = next(
        (f["name"] for f in fields if f.get("id_type") == "symbol"), "hgnc_symbol"
    )
    ensg_col = next((f["name"] for f in fields if f.get("id_type") == "ensg"), "ensg")

    catalog = tuple(
        GenesetInfo(
            name=f["name"],
            fullname=f.get("fullname") or "",
            description=f.get("description") or "",
            category=f.get("category") or "",
            subcategory=f.get("subcategory") or "",
            visibility=f.get("visibility") or "secondary",
            aliases=tuple(f.get("aliases") or ()),
            derived_from=tuple(f.get("derived_from") or ()),
            n_genes=int(f.get("n_genes") or 0),
        )
        for f in fields
        if f.get("value_type") == "boolean" and f["name"] in df.columns
    )
    inheritance_fields = tuple(
        FieldInfo(f["name"], f.get("fullname") or f["name"], f.get("category") or "")
        for f in fields
        if f.get("category") == INHERITANCE_CATEGORY
        and f.get("value_type") == "categorical"
        and f["name"] in df.columns
    )

    df = df.filter(pl.col(symbol_col).is_not_null())
    symbols: List[str] = df[symbol_col].cast(pl.String).to_list()
    by_symbol = {s.upper(): s for s in symbols}

    by_ensg: Dict[str, str] = {}
    if ensg_col in df.columns:
        for sym, ensg in zip(symbols, df[ensg_col].cast(pl.String).to_list()):
            if ensg:
                by_ensg.setdefault(ensg.split(".", 1)[0].upper(), sym)

    def layer(column: str) -> Dict[str, str]:
        if column not in df.columns:
            return {}
        values = df[column].cast(pl.String).to_list()
        return _unambiguous(
            (
                (ident, sym)
                for sym, v in zip(symbols, values)
                for ident in _split_pipe(v)
            ),
            by_symbol,
        )

    by_prev = layer(_PREV_COLUMN)
    # A previous symbol is the stronger claim: drop aliases it already covers.
    by_alias = {k: v for k, v in layer(_ALIAS_COLUMN).items() if k not in by_prev}

    set_members: Dict[str, FrozenSet[str]] = {}
    gene_sets: Dict[str, set] = {}
    for gs in catalog:
        flags = df[gs.name].fill_null(False).to_list()
        hits = frozenset(s for s, flag in zip(symbols, flags) if flag)
        set_members[gs.name] = hits
        for sym in hits:
            gene_sets.setdefault(sym, set()).add(gs.name)
    members = {sym: frozenset(names) for sym, names in gene_sets.items()}

    names = (
        df[_NAME_COLUMN].cast(pl.String).to_list()
        if _NAME_COLUMN in df.columns
        else [None] * len(symbols)
    )
    value_cols = [(c, label) for c, label in _INFO_COLUMNS.items() if c in df.columns]
    value_lists = {c: _column_values(df[c]) for c, _ in value_cols}
    inh_lists = {f.name: _column_values(df[f.name]) for f in inheritance_fields}
    info: Dict[str, GeneInfo] = {}
    for i, sym in enumerate(symbols):
        values = tuple(
            (label, _format_value(value_lists[c][i]))
            for c, label in value_cols
            if value_lists[c][i] is not None
        )
        inheritance = tuple(
            (f.fullname, inh_lists[f.name][i])
            for f in inheritance_fields
            if inh_lists[f.name][i]
        )
        if names[i] or values or inheritance:
            info[sym] = GeneInfo(sym, names[i] or "", values, inheritance)

    return GenetrekSnapshot(
        release=release,
        fetched_at=fetched_at,
        checked_at=fetched_at,
        catalog=catalog,
        inheritance_fields=inheritance_fields,
        members=members,
        set_members=set_members,
        info=info,
        by_symbol=by_symbol,
        by_ensg=by_ensg,
        by_prev=by_prev,
        by_alias=by_alias,
    )


# ---------------------------------------------------------------------------
# Global snapshot, disk cache and refresh
# ---------------------------------------------------------------------------

_snapshot: GenetrekSnapshot = GenetrekSnapshot(loading=True)
_refresh_lock = threading.Lock()
_refresh_task: Optional[asyncio.Task[None]] = None


def get_snapshot() -> GenetrekSnapshot:
    """The current snapshot. Never blocks; may be empty (see ``status_text``)."""
    return _snapshot


def set_snapshot(snapshot: GenetrekSnapshot) -> None:
    """Install a snapshot (startup, refresh, tests)."""
    global _snapshot
    _snapshot = snapshot


def clear_genetrek_cache() -> None:
    """Back to the initial "loading" state. For tests."""
    set_snapshot(GenetrekSnapshot(loading=True))


def _fetch(url: str, timeout: float) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": _USER_AGENT})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read()


def _api(path: str) -> str:
    return get_genetrek_url() + path.lstrip("/")


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def load_cached_snapshot() -> Optional[GenetrekSnapshot]:
    """Install the snapshot cached on disk, if any. Blocking (parquet read)."""
    cache_dir = get_genetrek_cache_dir()
    if cache_dir is None:
        return None
    meta_path = cache_dir / _META_FILE
    try:
        meta = json.loads(meta_path.read_text())
        release = meta["release"]
        df = pl.read_parquet(cache_dir / f"{release}.parquet")
        fetched_at = datetime.fromisoformat(meta["fetched_at"])
        snapshot = build_snapshot(
            df, meta["fields"], release=release, fetched_at=fetched_at
        )
    except FileNotFoundError:
        return None
    except (OSError, ValueError, KeyError, pl.exceptions.PolarsError) as e:
        logger.warning("Ignoring unreadable GeneTrek cache in %s: %s", cache_dir, e)
        return None
    set_snapshot(snapshot)
    logger.info("GeneTrek release %s loaded from cache", release)
    return snapshot


def refresh_snapshot(*, force: bool = False) -> GenetrekSnapshot:
    """Bring the snapshot up to GeneTrek's current release.

    Blocking (network and parquet read): call it through ``asyncio.to_thread``.
    One ``/api/health`` call when the release is unchanged. On any failure the
    previous snapshot stays in place, carrying the error.
    """
    if not _refresh_lock.acquire(blocking=False):
        return get_snapshot()  # a refresh is already running
    try:
        return _refresh(force)
    finally:
        _refresh_lock.release()


def _refresh(force: bool) -> GenetrekSnapshot:
    current = get_snapshot()
    now = datetime.now()
    try:
        health = json.loads(_fetch(_api("api/health"), _TIMEOUT_SECONDS))
        release = str(health["release"])
        if current.ok and current.release == release and not force:
            snapshot = replace(current, checked_at=now, error=None, loading=False)
            set_snapshot(snapshot)
            return snapshot

        fields_raw = json.loads(
            _fetch(
                _api("api/fields?include_hidden=true&include_counts=true"),
                _TIMEOUT_SECONDS,
            )
        )
        fields = fields_raw["fields"]
        data = _fetch(
            _api("api/download/database?file_format=parquet"), _DOWNLOAD_TIMEOUT_SECONDS
        )
        df = pl.read_parquet(data)
        snapshot = build_snapshot(df, fields, release=release, fetched_at=now)
    except Exception as e:  # noqa: BLE001 - a refresh never raises, see module doc
        logger.warning("GeneTrek refresh from %s failed: %s", get_genetrek_url(), e)
        snapshot = replace(
            current, error=f"{type(e).__name__}: {e}", loading=False, checked_at=now
        )
        set_snapshot(snapshot)
        return snapshot

    set_snapshot(snapshot)
    logger.info("GeneTrek release %s downloaded", release)
    _write_cache(release, data, fields, now)
    return snapshot


def _write_cache(
    release: str, data: bytes, fields: List[Dict[str, Any]], when: datetime
) -> None:
    cache_dir = get_genetrek_cache_dir()
    if cache_dir is None:
        return
    try:
        _atomic_write(cache_dir / f"{release}.parquet", data)
        meta = {"release": release, "fetched_at": when.isoformat(), "fields": fields}
        _atomic_write(cache_dir / _META_FILE, json.dumps(meta).encode())
        for old in cache_dir.glob("*.parquet"):
            if old.stem != release:
                old.unlink(missing_ok=True)
    except OSError as e:
        # The snapshot is already in memory; only the next cold start loses.
        logger.warning("Could not cache GeneTrek release in %s: %s", cache_dir, e)


async def _refresh_loop(interval_hours: float) -> None:
    await asyncio.to_thread(load_cached_snapshot)
    while True:
        await asyncio.to_thread(refresh_snapshot)
        await asyncio.sleep(interval_hours * 3600)


def start_refresh() -> None:
    """Load the cached snapshot, then keep it current. Safe to call twice."""
    global _refresh_task
    if _refresh_task is not None and not _refresh_task.done():
        return
    hours = float(GENETREK_CONFIG.get("refresh_hours") or 6)
    _refresh_task = asyncio.get_event_loop().create_task(_refresh_loop(hours))


def stop_refresh() -> None:
    global _refresh_task
    if _refresh_task is not None:
        _refresh_task.cancel()
        _refresh_task = None


# ---------------------------------------------------------------------------
# Links
# ---------------------------------------------------------------------------


def gene_url(symbol: str) -> str:
    """GeneTrek's page for one gene (its info dialog)."""
    return f"{get_genetrek_url()}?gene={urllib.parse.quote(symbol)}"


def search_url(tokens: Iterable[str], genome: str = "GRCh38") -> str:
    """GeneTrek's table filtered to genes and/or ``chr:start-end`` loci.

    The query travels in the fragment, which never reaches a server log: a
    locus or a list of genes from one patient's variant is not a public query.
    """
    query = ",".join(t for t in tokens if t)
    fragment = urllib.parse.urlencode({"search": query, "genome": genome}, safe=",:")
    return f"{get_genetrek_url()}#{fragment}"
