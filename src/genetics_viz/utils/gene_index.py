"""Cross-cohort gene index: every cohort's results by gene, as Parquet.

The gene page (``/gene/<symbol>``) and the variant page answer "have we seen
this gene, or this variant, before?" across every cohort of every data
directory. Reading the cohort tables per request is out of the question: they
sit on the data mount and run to hundreds of MB. So each is read once, off the
event loop, into one Parquet file sorted by gene under ``gene_index_dir``
(default ``<config dir>/.gene_index/``). A lookup is a filtered scan of those
files, small and sorted, and the row-group statistics skip most of each one.

Sources, all cohort-level (a cohort table is the union of its families'):

- ``snv``: every Wombat results table, one row per variant, sample and gene;
- ``sv``: the WisecondorX aberrations, one row per call and gene it overlaps;
- ``str``: the newest catalog's ExpansionHunter genotypes, by the locus's gene;
- ``roh``: the ``roh_genes`` table, one row per segment and gene.

A family with results but no cohort table is not indexed.

Validations and diagnostics are not indexed: they are small, already one file
per data directory, and change whenever a curator saves, so the pages read
them live.

A table is read again only when its ``(mtime, size)`` changes, and a schema
change here (:data:`INDEX_VERSION`) rebuilds everything. A refresh stats every
source, a few per cohort, on a timer (``refresh_minutes``). Failure is a
degraded state: a table that cannot be read is recorded with its error and the
rest of the index stays usable.
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import logging
import os
import re
import tempfile
import threading
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass, field, replace
from multiprocessing import get_context
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import polars as pl
import yaml

from genetics_viz.config_model import get_gene_index_dir
from genetics_viz.models import Cohort, DataStore
from genetics_viz.utils.ancestry import NULL_VALUES
from genetics_viz.utils.data import get_all_data_stores
from genetics_viz.utils.roh import find_cohort_roh_files
from genetics_viz.utils.sample_counts import (
    cohort_wcx_file,
    cohort_wombat_sources,
    family_header_sources,
    scan_cohort_wombat,
)
from genetics_viz.utils.str_eh import (
    eh_path,
    find_eh_labels,
    get_cohort_eh_dir,
    read_eh_table,
)
from genetics_viz.utils.tsv import read_tsv_or_none
from genetics_viz.utils.wisecondorx import infer_sv_type, parse_wisecondorx_bed

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "gene_index.yaml"
with open(_CONFIG_PATH) as _f:
    GENE_INDEX_CONFIG: Dict[str, Any] = yaml.safe_load(_f) or {}

#: Bumped whenever a schema below changes: every table is then read again.
INDEX_VERSION = 2

KINDS = ("snv", "sv", "str", "roh")

#: Columns of every kind's Parquet files, in order. Every file of a kind has
#: exactly these, so the files of all cohorts scan as one table.
SCHEMAS: Dict[str, Dict[str, Any]] = {
    "snv": {
        "gene_key": pl.Utf8,
        "gene": pl.Utf8,
        "cohort": pl.Utf8,
        "config": pl.Utf8,
        "family_id": pl.Utf8,
        "sample": pl.Utf8,
        "variant": pl.Utf8,
        "chrom": pl.Utf8,
        "pos": pl.Int64,
        "ref": pl.Utf8,
        "alt": pl.Utf8,
        "consequence": pl.Utf8,
        "impact": pl.Utf8,
        "clin_sig": pl.Utf8,
        "gnomad": pl.Float64,
        "sample_gt": pl.Utf8,
        "father_gt": pl.Utf8,
        "mother_gt": pl.Utf8,
    },
    "sv": {
        "gene_key": pl.Utf8,
        "gene": pl.Utf8,
        "cohort": pl.Utf8,
        "family_id": pl.Utf8,
        "sample": pl.Utf8,
        "locus": pl.Utf8,
        "variant": pl.Utf8,
        "chrom": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "sv_type": pl.Utf8,
        "call": pl.Utf8,
        "ratio": pl.Float64,
        "zscore": pl.Float64,
        "exonic": pl.Boolean,
    },
    "str": {
        "gene_key": pl.Utf8,
        "gene": pl.Utf8,
        "cohort": pl.Utf8,
        "label": pl.Utf8,
        "family_id": pl.Utf8,
        "sample": pl.Utf8,
        "locus": pl.Utf8,
        "variant": pl.Utf8,
        "chrom": pl.Utf8,
        "pos": pl.Int64,
        "end": pl.Int64,
        "motif": pl.Utf8,
        "genotype": pl.Utf8,
        "allele1": pl.Int64,
        "allele2": pl.Int64,
        "status": pl.Utf8,
        "severity": pl.Utf8,
        "disease": pl.Utf8,
        "inheritance": pl.Utf8,
    },
    "roh": {
        "gene_key": pl.Utf8,
        "gene": pl.Utf8,
        "cohort": pl.Utf8,
        "label": pl.Utf8,
        "family_id": pl.Utf8,
        "sample": pl.Utf8,
        "chrom": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "length_mb": pl.Float64,
        "roh_class": pl.Utf8,
        "gene_fraction": pl.Float64,
    },
}

#: The Wombat columns the SNV index reads; the hundreds of others are skipped.
_WOMBAT_COLUMNS = (
    "#CHROM", "POS", "REF", "ALT", "sample", "sample_gt", "father_gt",
    "mother_gt", "VEP_Consequence", "VEP_IMPACT", "VEP_SYMBOL", "VEP_CLIN_SIG",
    "fafmax_faf95_max_genomes",
)  # fmt: skip

#: The ``roh_genes`` columns the ROH index reads (the cohort table is ~93 MB).
_ROH_GENES_COLUMNS = (
    "family_id", "IID", "chrom", "roh_start", "roh_end", "roh_length_mb",
    "roh_class", "gene", "gene_fraction_in_roh",
)  # fmt: skip

_MANIFEST = "manifest.json"
#: Set to build in the app's own process: the tests, where spawning a process
#: per table would cost seconds each.
IN_PROCESS_ENV = "GENETICS_VIZ_GENE_INDEX_IN_PROCESS"
#: The error of a table whose build process died (the kernel's OOM killer
#: leaves no exception to catch). Kept until the table changes.
KILLED = (
    "The process reading this table was stopped (out of memory?). "
    "It is skipped until the table changes."
)

#: ``classify_cnv``'s label for a segment that is not a call.
BELOW_THRESHOLD = "Below threshold"

_VARIANT_RE = re.compile(
    r"^(?:chr)?([0-9]{1,2}|X|Y|M|MT)[:\-_ ](\d+)[:\-_ ]([ACGTN]+)[:\-_ >/]([ACGTN*]+)$",
    re.IGNORECASE,
)


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Source:
    """One cohort table the index reads."""

    kind: str
    cohort: str
    #: The Wombat config, the STR catalog label or the ROH bundle and tag.
    label: str
    path: Path


def discover_sources(store: DataStore) -> List[Source]:
    """Every indexable cohort table of a data directory (blocking: it lists)."""
    sources: List[Source] = []
    for name, cohort in sorted(store.cohorts.items()):
        try:
            sources.extend(_cohort_sources(store.data_dir, cohort))
        except OSError as e:
            logger.warning("Cannot list the results of cohort %s: %s", name, e)
    return sources


def _cohort_sources(data_dir: Path, cohort: Cohort) -> List[Source]:
    name = cohort.name
    out = [
        Source("snv", name, src.wombat_config, src.path)
        for src in cohort_wombat_sources(data_dir, name)
    ]
    wcx = cohort_wcx_file(data_dir, name)
    if wcx is not None:
        out.append(Source("sv", name, "wisecondorx", wcx))
    eh_dir = get_cohort_eh_dir(cohort)
    labels = find_eh_labels(eh_dir, name)
    if labels:  # newest catalog only, as the STR tabs open on
        out.append(
            Source(
                "str", name, labels[0], eh_path(eh_dir, name, labels[0], "genotypes")
            )
        )
    roh = find_cohort_roh_files(cohort)
    if roh is not None and roh.paths.get("roh_genes") is not None:
        out.append(
            Source(
                "roh",
                name,
                f"{roh.bundle_version}_{roh.filter_tag}",
                roh.paths["roh_genes"],
            )
        )
    return out


# ---------------------------------------------------------------------------
# Building one table's frame
# ---------------------------------------------------------------------------


def _conform(df: pl.DataFrame, kind: str) -> pl.DataFrame:
    """Exactly the kind's columns, in order and typed; missing ones are null."""
    exprs = []
    for col, dtype in SCHEMAS[kind].items():
        if col in df.columns:
            exprs.append(pl.col(col).cast(dtype, strict=False))
        else:
            exprs.append(pl.lit(None, dtype=dtype).alias(col))
    return df.select(exprs)


def _with_chr(col: str) -> pl.Expr:
    """``22`` and ``chr22`` both as ``chr22``, so lookups need one spelling."""
    value = pl.col(col).cast(pl.Utf8)
    return (
        pl.when(value.str.starts_with("chr"))
        .then(value)
        .otherwise(pl.lit("chr") + value)
    )


def _explode_genes(df: pl.DataFrame, column: str) -> pl.DataFrame:
    """One row per gene of a ``,``/``&`` joined cell; rows with none dropped."""
    return (
        df.with_columns(
            pl.col(column)
            .cast(pl.Utf8)
            .fill_null("")
            .str.replace_all("&", ",")
            .str.split(",")
            .alias("gene")
        )
        .explode("gene")
        .with_columns(pl.col("gene").str.strip_chars())
        .filter(pl.col("gene").is_not_null() & (pl.col("gene") != ""))
        .unique(maintain_order=True)
        .with_columns(pl.col("gene").str.to_uppercase().alias("gene_key"))
    )


def _sample_families(cohort: Cohort) -> Dict[str, str]:
    """Sample -> family, the first family a sample sits in."""
    out: Dict[str, str] = {}
    for family in cohort.families.values():
        for sample in family.samples:
            out.setdefault(sample.sample_id, family.family_id)
    return out


def _family_from_pedigree(df: pl.DataFrame, cohort: Cohort) -> pl.DataFrame:
    families = _sample_families(cohort)
    return df.with_columns(
        pl.col("sample")
        .cast(pl.Utf8)
        .replace_strict(families, default=None, return_dtype=pl.Utf8)
        .alias("family_id")
    )


def _collect(lf: pl.LazyFrame) -> pl.DataFrame:
    """Collect with the streaming engine (polars 2.x name, then 1.x)."""
    try:
        return lf.collect(engine="streaming")
    except (TypeError, ValueError):
        return lf.collect(streaming=True)


def _snv_frame(source: Source, data_dir: Path, cohort: Cohort) -> pl.DataFrame:
    wombat_config = source.label
    vep_config = re.match(
        rf"{re.escape(cohort.name)}\.rare\.([^.]+)\.", source.path.name
    )
    header_sources = (
        family_header_sources(data_dir, cohort, vep_config.group(1), wombat_config)
        if vep_config
        else ()
    )
    lf = scan_cohort_wombat(
        source.path, header_sources, _WOMBAT_COLUMNS, null_values=[".", ""]
    )
    if lf is None:
        return _conform(pl.DataFrame(), "snv")
    columns = lf.collect_schema().names()
    if "VEP_SYMBOL" not in columns:
        return _conform(pl.DataFrame(), "snv")
    keys = ["#CHROM", "POS", "REF", "ALT", "sample"]
    first = [c for c in columns if c not in (*keys, "VEP_SYMBOL", "VEP_Consequence")]
    # Wombat writes a row per transcript: one per variant, sample and gene is
    # kept, the gene's consequences joined. Grouped before anything else, so
    # the streaming engine never holds a whole table.
    df = _collect(
        lf.with_columns(
            pl.col("VEP_SYMBOL")
            .fill_null("")
            .str.replace_all("&", ",")
            .str.split(",")
            .alias("gene")
        )
        .explode("gene")
        .with_columns(pl.col("gene").str.strip_chars())
        .filter(pl.col("gene").is_not_null() & (pl.col("gene") != ""))
        .group_by([*keys, "gene"])
        .agg(
            pl.col("VEP_Consequence")
            .drop_nulls()
            .unique()
            .sort()
            .str.join(",")
            .alias("consequence")
            if "VEP_Consequence" in columns
            else pl.lit(None, dtype=pl.Utf8).alias("consequence"),
            *[pl.col(c).drop_nulls().first() for c in first],
        )
    )
    if df.is_empty():
        return _conform(pl.DataFrame(), "snv")
    df = _family_from_pedigree(df, cohort)
    return _conform(
        df.with_columns(
            pl.col("gene").str.to_uppercase().alias("gene_key"),
            pl.lit(cohort.name).alias("cohort"),
            pl.lit(wombat_config).alias("config"),
            (
                pl.col("#CHROM")
                + ":"
                + pl.col("POS")
                + ":"
                + pl.col("REF")
                + ":"
                + pl.col("ALT")
            ).alias("variant"),
            _with_chr("#CHROM").alias("chrom"),
            pl.col("POS").alias("pos"),
            pl.col("REF").alias("ref"),
            pl.col("ALT").alias("alt"),
            *[
                pl.col(src).alias(dst)
                for src, dst in (
                    ("VEP_IMPACT", "impact"),
                    ("VEP_CLIN_SIG", "clin_sig"),
                    ("fafmax_faf95_max_genomes", "gnomad"),
                )
                if src in df.columns
            ],
        ),
        "snv",
    )


def _sv_frame(source: Source, cohort: Cohort) -> pl.DataFrame:
    df = parse_wisecondorx_bed(source.path)
    if df is None or len(df) == 0 or "genic_symbol" not in df.columns:
        return _conform(pl.DataFrame(), "sv")
    # Not a call: the SVs tab and Search hide these by default, and they are
    # most of a WisecondorX table (871k of CANDY_mpx's 1M gene rows).
    if "wisecondorX" in df.columns:
        df = df.filter(pl.col("wisecondorX") != BELOW_THRESHOLD)
    if df.is_empty():
        return _conform(pl.DataFrame(), "sv")
    # The key validations and diagnostics use, built as the Search page does
    sv_types = [
        infer_sv_type(r)
        for r in df.select(
            [c for c in ("wisecondorX", "type", "ratio") if c in df.columns]
        ).to_dicts()
    ]
    df = df.with_columns(pl.Series("sv_type", sv_types, dtype=pl.Utf8))
    if "family_id" not in df.columns:
        df = _family_from_pedigree(df, cohort)
    df = df.with_row_index("_row")
    # One row per call and gene, exploded apart from the call's other columns:
    # exploded with them, a call over 80 genes copied its exonic gene list 80
    # times, and one cohort table took 13 GB.
    genes = _explode_genes(df.select("_row", "genic_symbol"), "genic_symbol")
    if "exonic_symbol" in df.columns:
        exonic = (
            _explode_genes(df.select("_row", "exonic_symbol"), "exonic_symbol")
            .select("_row", "gene_key")
            .unique()
            .with_columns(pl.lit(True).alias("exonic"))
        )
        genes = genes.join(exonic, on=["_row", "gene_key"], how="left")
    else:
        genes = genes.with_columns(pl.lit(None, dtype=pl.Boolean).alias("exonic"))
    calls = df.select(
        "_row",
        pl.lit(cohort.name).alias("cohort"),
        "family_id",
        "sample",
        pl.col("chr:start-end").alias("locus"),
        (pl.col("chr:start-end") + ":" + pl.col("sv_type")).alias("variant"),
        pl.col("chr").alias("chrom"),
        "start",
        "end",
        "sv_type",
        pl.col("wisecondorX").alias("call")
        if "wisecondorX" in df.columns
        else pl.lit(None, dtype=pl.Utf8).alias("call"),
        *[c for c in ("ratio", "zscore") if c in df.columns],
    )
    return _conform(
        genes.select("_row", "gene", "gene_key", pl.col("exonic").fill_null(False))
        .join(calls, on="_row")
        .sort("_row"),
        "sv",
    )


def _str_frame(source: Source, cohort: Cohort) -> pl.DataFrame:
    df = read_eh_table(source.path)
    if df is None or len(df) == 0 or "locus" not in df.columns:
        return _conform(pl.DataFrame(), "str")
    # The genotype table has no gene column; jmsa50 loci are named by their gene
    df = df.with_columns(
        pl.col("locus").str.split("_").list.first().alias("gene"),
    ).with_columns(pl.col("gene").str.to_uppercase().alias("gene_key"))
    return _conform(
        df.with_columns(
            pl.lit(cohort.name).alias("cohort"),
            pl.lit(source.label).alias("label"),
            # str_validation.str_variant_key
            (
                pl.col("chrom").cast(pl.Utf8)
                + ":"
                + pl.col("pos").cast(pl.Utf8)
                + "-"
                + pl.col("end").cast(pl.Utf8)
                + ":"
                + pl.col("locus")
            ).alias("variant"),
        ),
        "str",
    )


def _roh_frame(source: Source, cohort: Cohort) -> pl.DataFrame:
    with open(source.path) as f:
        header = f.readline().rstrip("\n").split("\t")
    wanted = [c for c in _ROH_GENES_COLUMNS if c in header]
    if "gene" not in wanted:
        return _conform(pl.DataFrame(), "roh")
    df = read_tsv_or_none(
        source.path, columns=wanted, infer_schema_length=0, null_values=NULL_VALUES
    )
    if df is None or len(df) == 0:
        return _conform(pl.DataFrame(), "roh")
    df = df.rename({"IID": "sample"}) if "IID" in df.columns else df
    if "family_id" not in df.columns:
        df = _family_from_pedigree(df, cohort)
    df = _explode_genes(df, "gene")
    return _conform(
        df.with_columns(
            pl.lit(cohort.name).alias("cohort"),
            pl.lit(source.label).alias("label"),
            _with_chr("chrom").alias("chrom"),
            pl.col("roh_start").alias("start"),
            pl.col("roh_end").alias("end"),
            pl.col("roh_length_mb").alias("length_mb"),
            pl.col("gene_fraction_in_roh").alias("gene_fraction")
            if "gene_fraction_in_roh" in df.columns
            else pl.lit(None, dtype=pl.Float64).alias("gene_fraction"),
        ),
        "roh",
    )


def build_frame(source: Source, data_dir: Path, cohort: Cohort) -> pl.DataFrame:
    """One source table, conformed to its kind's schema (blocking)."""
    if source.kind == "snv":
        return _snv_frame(source, data_dir, cohort)
    if source.kind == "sv":
        return _sv_frame(source, cohort)
    if source.kind == "str":
        return _str_frame(source, cohort)
    return _roh_frame(source, cohort)


# ---------------------------------------------------------------------------
# On-disk index and its state
# ---------------------------------------------------------------------------


@dataclass
class IndexStatus:
    """What the gene pages say about the index."""

    building: bool = False
    #: Sources checked so far in the running refresh, and how many there are.
    done: int = 0
    total: int = 0
    checked_at: Optional[datetime] = None
    #: Source path -> why it could not be indexed.
    errors: Dict[str, str] = field(default_factory=dict)


#: Data directory path -> {source path -> manifest entry}. Replaced whole, so a
#: lookup on the event loop never sees a half-updated mapping.
_manifests: Dict[str, Dict[str, Dict[str, Any]]] = {}
_status: Dict[str, IndexStatus] = {}
_build_lock = threading.Lock()
_refresh_task: Optional[asyncio.Task[None]] = None


def _dir_key(data_dir: str) -> str:
    """The data directory's subdirectory: its path would not make a file name."""
    return hashlib.sha1(data_dir.encode()).hexdigest()[:16]


def _store_root(data_dir: str) -> Optional[Path]:
    root = get_gene_index_dir()
    return None if root is None else root / _dir_key(data_dir)


def _ensure_dir(path: Path) -> None:
    # The index holds genotypes: readable by the app's user only.
    path.mkdir(parents=True, exist_ok=True, mode=0o700)


def _atomic_write_bytes(path: Path, data: bytes) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "wb") as fh:
            fh.write(data)
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _write_parquet(df: pl.DataFrame, path: Path) -> None:
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    os.close(fd)
    try:
        df.sort("gene_key").write_parquet(
            tmp,
            statistics=True,
            row_group_size=int(GENE_INDEX_CONFIG.get("row_group_size") or 20000),
        )
        os.replace(tmp, path)
    except BaseException:
        Path(tmp).unlink(missing_ok=True)
        raise


def _file_name(source: Source) -> str:
    """``<kind>/<cohort>.<label>.<hash>.parquet``: readable, and unique per path."""
    safe = re.sub(r"[^A-Za-z0-9_.-]", "_", f"{source.cohort}.{source.label}")
    digest = hashlib.sha1(str(source.path).encode()).hexdigest()[:8]
    return f"{source.kind}/{safe}.{digest}.parquet"


def _load_manifest(root: Path) -> Dict[str, Dict[str, Any]]:
    try:
        raw = json.loads((root / _MANIFEST).read_text())
    except FileNotFoundError:
        return {}
    except (OSError, ValueError) as e:
        logger.warning("Ignoring unreadable gene index manifest in %s: %s", root, e)
        return {}
    if raw.get("version") != INDEX_VERSION:
        return {}
    return dict(raw.get("entries") or {})


def _save_manifest(root: Path, data_dir: str, entries: Dict[str, Dict]) -> None:
    payload = {"version": INDEX_VERSION, "data_dir": data_dir, "entries": entries}
    _atomic_write_bytes(root / _MANIFEST, json.dumps(payload, indent=1).encode())


def _stat_key(path: Path) -> Optional[Tuple[float, int]]:
    try:
        st = path.stat()
    except OSError:
        return None
    return st.st_mtime, st.st_size


def refresh_store(store: DataStore) -> None:
    """Bring one data directory's index up to date (blocking).

    The manifest left on disk by a previous run is published first, so a
    restart serves lookups from it while the sources are checked again.
    """
    data_dir = str(store.data_dir)
    root = _store_root(data_dir)
    if root is None:
        return
    status = _status.setdefault(data_dir, IndexStatus())
    entries = dict(_manifests.get(data_dir) or _load_manifest(root))
    _manifests[data_dir] = dict(entries)

    status.building = True
    status.done = 0
    try:
        sources = discover_sources(store)
        status.total = len(sources)
        wanted = {str(s.path) for s in sources}
        changed = False
        for stale in [p for p in entries if p not in wanted]:
            (root / entries.pop(stale).get("file", "")).unlink(missing_ok=True)
            changed = True
        if changed:
            _ensure_dir(root)
            _save_manifest(root, data_dir, entries)
        for source in sources:
            # Saved after every table, so a restart keeps what was read.
            if _refresh_source(store, root, source, entries, status):
                _save_manifest(root, data_dir, entries)
                _manifests[data_dir] = dict(entries)
            status.done += 1
        _manifests[data_dir] = dict(entries)
    finally:
        status.building = False
        status.checked_at = datetime.now()


def _refresh_source(
    store: DataStore,
    root: Path,
    source: Source,
    entries: Dict[str, Dict[str, Any]],
    status: IndexStatus,
) -> bool:
    """Re-read one source if it changed. Whether its entry changed."""
    key = str(source.path)
    stat = _stat_key(source.path)
    old = entries.get(key)
    if (
        old is not None
        and stat is not None
        and (old.get("mtime"), old.get("size")) == stat
        and (old.get("error") or (root / old.get("file", "")).exists())
    ):
        if old.get("error"):
            status.errors[key] = old["error"]
        return False
    entry: Dict[str, Any] = {
        "kind": source.kind,
        "cohort": source.cohort,
        "label": source.label,
        "mtime": stat[0] if stat else None,
        "size": stat[1] if stat else None,
        "indexed_at": datetime.now().isoformat(timespec="seconds"),
    }
    cohort = store.cohorts.get(source.cohort)
    try:
        if cohort is None or stat is None:
            raise FileNotFoundError(key)
        name = _file_name(source)
        _ensure_dir((root / name).parent)
        rows = _build(source, store.data_dir, cohort, root / name)
        entry.update(file=name, rows=rows)
        status.errors.pop(key, None)
    except Exception as e:  # noqa: BLE001 - one bad table must not stop the rest
        logger.warning("Gene index: cannot read %s: %s", key, e)
        entry["error"] = (
            str(e) if isinstance(e, ChildKilled) else f"{type(e).__name__}: {e}"
        )
        status.errors[key] = entry["error"]
        if old is not None and old.get("file"):
            (root / old["file"]).unlink(missing_ok=True)
    entries[key] = entry
    return True


def _build_to_file(source: Source, data_dir: Path, cohort: Cohort, path: Path) -> int:
    """Read one source and write its Parquet file. The number of rows."""
    frame = build_frame(source, data_dir, cohort)
    _write_parquet(frame, path)
    return len(frame)


def _child_init() -> None:
    """A build process: lower priority than the app, the OOM killer's first pick."""
    try:
        os.nice(10)
    except OSError:
        pass
    try:  # Linux only; raising one's own score needs no privilege
        with open("/proc/self/oom_score_adj", "w") as f:
            f.write("1000")
    except OSError:
        pass


def _build(source: Source, data_dir: Path, cohort: Cohort, path: Path) -> int:
    """:func:`_build_to_file` in a process of its own, unless :data:`IN_PROCESS_ENV`.

    The largest cohort tables take gigabytes to read. In a separate process an
    out-of-memory kill costs that one table (:data:`KILLED`), not the app: in
    v0.26.0 it restarted the app every 45 s, each restart beginning the index
    again. The process gets half the cores, so the pages stay responsive.
    """
    if os.environ.get(IN_PROCESS_ENV):
        return _build_to_file(source, data_dir, cohort, path)
    light = replace(cohort, _dataframe=None)  # the pedigree frame is not needed
    threads = str(max(1, (os.cpu_count() or 2) // 2))
    previous = os.environ.get("POLARS_MAX_THREADS")
    os.environ["POLARS_MAX_THREADS"] = threads  # read by the child's polars only
    try:
        with ProcessPoolExecutor(
            max_workers=1, mp_context=get_context("spawn"), initializer=_child_init
        ) as pool:
            return pool.submit(_build_to_file, source, data_dir, light, path).result()
    except BrokenProcessPool as e:
        raise ChildKilled(KILLED) from e
    finally:
        if previous is None:
            os.environ.pop("POLARS_MAX_THREADS", None)
        else:
            os.environ["POLARS_MAX_THREADS"] = previous


class ChildKilled(RuntimeError):
    """The build process of a table died without an exception."""


def refresh_all(*, wait: bool = False) -> None:
    """Refresh every registered data directory (blocking).

    Returns at once when a refresh is already running, unless ``wait``.
    """
    if not _build_lock.acquire(blocking=wait):
        return
    try:
        for store in list(get_all_data_stores().values()):
            try:
                refresh_store(store)
            except Exception:  # noqa: BLE001 - see the module doc
                logger.exception("Gene index refresh failed for %s", store.data_dir)
    finally:
        _build_lock.release()


async def _refresh_loop(minutes: float) -> None:
    while True:
        await asyncio.to_thread(refresh_all)
        await asyncio.sleep(minutes * 60)


def start_gene_index() -> None:
    """Build or update the index off the event loop, then keep it current."""
    global _refresh_task
    if _refresh_task is not None and not _refresh_task.done():
        return
    minutes = float(GENE_INDEX_CONFIG.get("refresh_minutes") or 30)
    _refresh_task = asyncio.get_event_loop().create_task(_refresh_loop(minutes))


def stop_gene_index() -> None:
    global _refresh_task
    if _refresh_task is not None:
        _refresh_task.cancel()
        _refresh_task = None


def clear_gene_index() -> None:
    """Forget the in-memory manifests and status (tests). The files stay."""
    _manifests.clear()
    _status.clear()


def index_status() -> IndexStatus:
    """The status of every registered data directory, merged."""
    merged = IndexStatus()
    stamps = []
    for path_str in get_all_data_stores():
        st = _status.get(path_str)
        if st is None:
            # Not checked yet in this run: the first refresh is pending.
            merged.building = True
            continue
        merged.building |= st.building
        merged.done += st.done if st.building else st.total
        merged.total += st.total
        merged.errors.update(st.errors)
        if st.checked_at is not None:
            stamps.append(st.checked_at)
    merged.checked_at = min(stamps) if stamps else None
    return merged


def indexed_tables() -> int:
    """How many source tables the published manifests hold."""
    return sum(
        1
        for path_str in get_all_data_stores()
        for e in (_manifests.get(path_str) or {}).values()
        if e.get("file")
    )


# ---------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------


def _empty(kind: str) -> pl.DataFrame:
    return pl.DataFrame(schema={**SCHEMAS[kind], "data_dir": pl.Utf8})


def _scan(kind: str, predicate: pl.Expr) -> pl.DataFrame:
    """Rows of every data directory's ``kind`` files matching ``predicate``."""
    frames: List[pl.DataFrame] = []
    for path_str in get_all_data_stores():
        root = _store_root(path_str)
        entries = _manifests.get(path_str) or {}
        if root is None:
            continue
        files = [
            root / e["file"]
            for e in entries.values()
            if e.get("kind") == kind and e.get("file")
        ]
        files = [f for f in files if f.exists()]
        if not files:
            continue
        try:
            found = pl.scan_parquet(files).filter(predicate).collect()
        except (OSError, pl.exceptions.PolarsError) as e:
            # A file replaced by a refresh between the listing and the scan
            logger.warning("Gene index lookup in %s failed: %s", root, e)
            continue
        frames.append(found.with_columns(pl.lit(path_str).alias("data_dir")))
    if not frames:
        return _empty(kind)
    return pl.concat(frames, how="vertical_relaxed")


def lookup_gene(kind: str, gene_keys: Iterable[str]) -> pl.DataFrame:
    """Every indexed row of ``kind`` on one of ``gene_keys`` (upper-cased symbols)."""
    keys = sorted({k.upper() for k in gene_keys if k})
    if not keys:
        return _empty(kind)
    return _scan(kind, pl.col("gene_key").is_in(keys))


@dataclass(frozen=True)
class VariantQuery:
    """A small variant, as typed: ``chr7-32000000-C-T``, ``7:32000000:C:T``."""

    chrom: str
    pos: int
    ref: str
    alt: str

    @property
    def label(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


def parse_variant(text: Optional[str]) -> Optional[VariantQuery]:
    """A :class:`VariantQuery`, or ``None`` when ``text`` is not a small variant."""
    m = _VARIANT_RE.match((text or "").strip())
    if not m:
        return None
    chrom = m.group(1).upper()
    chrom = "M" if chrom == "MT" else chrom
    return VariantQuery(
        f"chr{chrom}", int(m.group(2)), m.group(3).upper(), m.group(4).upper()
    )


def lookup_variant(query: VariantQuery) -> pl.DataFrame:
    """Every indexed SNV row of the variant: one per carrier, table and gene."""
    return _scan(
        "snv",
        (pl.col("chrom") == query.chrom)
        & (pl.col("pos") == query.pos)
        & (pl.col("ref") == query.ref)
        & (pl.col("alt") == query.alt),
    )
