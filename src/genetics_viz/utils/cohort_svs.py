"""The cohort SVs tab without UI: large copy-number calls across families.

Reads the cohort's WisecondorX aberrations (``cohorts/<C>/svs/wisecondorx/
<C>_aberrations[.annotated].bed``) **lazily**: a cohort table can run to
gigabytes (CANDY_mpx's took 13.4 GB read whole), and 87 % of its rows are
below threshold. Only the columns used are scanned, and only calls past the
permissive thresholds and the smallest size offered are collected. From those:

- the calls of the individuals shown, filtered by strength and size;
- whole-chromosome gains and losses, by the family report's rule
  (``family_report.whole_chromosome_events``), so the two always agree;
- regions where same-direction calls of two or more families overlap.

Vedette has no cohort table yet, so it is not here.
"""

import logging
from pathlib import Path
from typing import Any, Collection, Dict, List, Optional, Tuple

import polars as pl

from .wisecondorx import WISECONDORX_CONFIG

logger = logging.getLogger(__name__)

Key = Tuple[str, str]

#: Minimum call sizes the tab offers, in kb, and the default.
SIZE_OPTIONS_KB = (100, 500, 1000, 5000, 10000)
DEFAULT_MIN_KB = 1000
#: Robust first: the same labels and colours as the family SVs tab.
STRENGTHS = ("robust", "permissive")
_LEVELS = {
    ("robust", "loss"): WISECONDORX_CONFIG["robust_loss"],
    ("robust", "gain"): WISECONDORX_CONFIG["robust_gain"],
    ("permissive", "loss"): WISECONDORX_CONFIG["permissive_loss"],
    ("permissive", "gain"): WISECONDORX_CONFIG["permissive_gain"],
}
_NUMERIC = ("start", "end", "ratio", "zscore")
_TEXT = ("chr", "barcode", "sample", "family_id", "exonic_symbol", "genic_symbol")

#: ``(path, mtime, size)`` -> the calls past the permissive thresholds.
_cache: Dict[Tuple[str, float, int], pl.DataFrame] = {}


def clear_cache() -> None:
    _cache.clear()


def _past(level: Dict[str, Any], gain: bool) -> pl.Expr:
    r, z = pl.col("ratio"), pl.col("zscore")
    if gain:
        return (r >= level["ratio_threshold"]) & (z >= level["zscore_threshold"])
    return (r <= level["ratio_threshold"]) & (z <= level["zscore_threshold"])


def scan_calls(path: Path, min_kb: int = SIZE_OPTIONS_KB[0]) -> pl.DataFrame:
    """Blocking: the calls past the permissive thresholds and ``min_kb``.

    Columns: ``chr`` (with its ``chr`` prefix), ``start``, ``end``, ``ratio``,
    ``zscore``, ``sample``, ``family_id`` (``""`` when the table has none),
    ``exonic_symbol``, ``direction`` (gain/loss) and ``strength``
    (robust/permissive). Cached on the file's stat.
    """
    stat = path.stat()
    key = (str(path), stat.st_mtime, stat.st_size)
    if key in _cache:
        return _cache[key]
    with open(path, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
    wanted = [c for c in (*_NUMERIC, *_TEXT) if c in header]
    if not {"chr", "start", "end", "ratio", "zscore"} <= set(wanted) or not (
        {"barcode", "sample"} & set(wanted)
    ):
        raise ValueError(f"{path.name}: not a WisecondorX aberrations table")
    lazy = (
        pl.scan_csv(
            path,
            separator="\t",
            infer_schema=False,
            quote_char=None,
            truncate_ragged_lines=True,
        )
        .select(wanted)
        .with_columns(pl.col(c).cast(pl.Float64, strict=False) for c in _NUMERIC)
        .filter((pl.col("end") - pl.col("start")) >= min_kb * 1000)
        .filter(
            _past(_LEVELS[("permissive", "gain")], True)
            | _past(_LEVELS[("permissive", "loss")], False)
        )
    )
    if "sample" not in wanted:
        lazy = lazy.rename({"barcode": "sample"})
    for column in ("family_id", "exonic_symbol"):
        if column not in wanted:
            lazy = lazy.with_columns(pl.lit("").alias(column))
    gain = pl.col("ratio") > 0
    lazy = lazy.with_columns(
        pl.when(pl.col("chr").str.starts_with("chr"))
        .then(pl.col("chr"))
        .otherwise(pl.lit("chr") + pl.col("chr"))
        .alias("chr"),
        pl.when(gain).then(pl.lit("gain")).otherwise(pl.lit("loss")).alias("direction"),
        pl.when(
            _past(_LEVELS[("robust", "gain")], True)
            | _past(_LEVELS[("robust", "loss")], False)
        )
        .then(pl.lit("robust"))
        .otherwise(pl.lit("permissive"))
        .alias("strength"),
    ).select(
        "chr", "start", "end", "ratio", "zscore", "sample", "family_id",
        "exonic_symbol", "direction", "strength",
    )  # fmt: skip
    calls = lazy.collect(engine="streaming")
    _cache.clear()  # one cohort table at a time: they can be large
    _cache[key] = calls
    return calls


def call_rows(
    calls: pl.DataFrame,
    families_of: Dict[str, List[str]],
    validations: Optional[Dict[Tuple[str, str], List[Tuple]]] = None,
) -> List[Dict[str, Any]]:
    """One table row per call and family of its sample.

    ``families_of`` maps a sample to its families in the cohort (a sample can
    sit in two); the table's own ``family_id`` wins when set. ``validations``
    is ``load_validation_map(svs.tsv)``: the call's standing status, else
    ``TODO``.
    """
    from .validation_status import STATUS_COLORS, UNVALIDATED, consensus

    out: List[Dict[str, Any]] = []
    for r in calls.to_dicts():
        start, end = int(r["start"]), int(r["end"])
        sv_type = "dup" if r["direction"] == "gain" else "del"
        locus = f"{r['chr']}:{start}-{end}"
        level = _LEVELS[(r["strength"], r["direction"])]
        sample = str(r["sample"])
        fids = [str(r["family_id"])] if r["family_id"] else families_of.get(sample, [])
        for fid in fids or [""]:
            row = {
                "_key": f"{fid}:{sample}:{locus}",
                "FID": fid,
                "sample": sample,
                "locus": locus,
                "chr": r["chr"],
                "start": start,
                "end": end,
                "size_mb": round((end - start) / 1e6, 2),
                "call": level["label"],
                "call_badges": [{"label": level["label"], "color": level["color"]}],
                "direction": r["direction"],
                "strength": r["strength"],
                "ratio": round(r["ratio"], 3),
                "zscore": round(r["zscore"], 2),
                "exonic_symbol": r["exonic_symbol"] or "",
                "_genes": [g for g in (r["exonic_symbol"] or "").split(",") if g],
                "Variant": f"{locus}:{sv_type}",
            }
            if validations is not None:
                entries = [
                    e
                    for e in validations.get((row["Variant"], sample), [])
                    if e[3] != "1"
                ]
                status = (
                    consensus((e[0], e[1]) for e in entries)[0] if entries else ""
                ) or UNVALIDATED
                row["Validation"] = status
                row["Validation_badges"] = [
                    {"label": status, "color": STATUS_COLORS.get(status, "#9e9e9e")}
                ]
            out.append(row)
    return out


def filter_calls(
    rows: List[Dict[str, Any]],
    keys: Optional[Collection[Key]],
    *,
    min_kb: int = DEFAULT_MIN_KB,
    strengths: Collection[str] = ("robust",),
) -> List[Dict[str, Any]]:
    """The calls of the individuals in ``keys`` (``None``: all) the filters leave.

    Sorted by chromosome order, then start.
    """
    from .cytobands import CHROM_ORDER, norm_chrom

    order = {c: i for i, c in enumerate(CHROM_ORDER)}
    wanted = set(strengths)
    out = [
        r
        for r in rows
        if (keys is None or (r["FID"], r["sample"]) in keys)
        and r["end"] - r["start"] >= min_kb * 1000
        and (not wanted or r["strength"] in wanted)
    ]
    out.sort(
        key=lambda r: (order.get(norm_chrom(r["chr"]), 99), r["start"], r["sample"])
    )
    return out


def aneuploidy_rows(rows: List[Dict[str, Any]], min_fraction: float) -> List[Dict]:
    """Whole-chromosome gains and losses of the calls in ``rows``.

    The family report's rule, on the calls shown: same-direction calls summed
    per individual and chromosome against its assessable span.
    """
    from .family_report import whole_chromosome_events

    if not rows:
        return []
    frame = pl.DataFrame(
        {
            "chr": [r["chr"] for r in rows],
            "start": [float(r["start"]) for r in rows],
            "end": [float(r["end"]) for r in rows],
            "ratio": [float(r["ratio"]) for r in rows],
            "sample": [f"{r['FID']}\t{r['sample']}" for r in rows],
        }
    )
    out = []
    for key, events in whole_chromosome_events(frame, min_fraction).items():
        fid, sample = key.split("\t", 1)
        out.append(
            {"_key": key, "FID": fid, "sample": sample, "events": ", ".join(events)}
        )
    out.sort(key=lambda r: (r["FID"], r["sample"]))
    return out


def recurrent_regions(
    rows: List[Dict[str, Any]], minimum_families: int = 2
) -> List[Dict[str, Any]]:
    """Regions where same-direction calls of ``minimum_families`` families overlap.

    Calls are merged per chromosome and direction while they overlap (a chain
    of overlaps is one region). Most families first.
    """
    from .cytobands import CHROM_ORDER, norm_chrom

    order = {c: i for i, c in enumerate(CHROM_ORDER)}
    groups: Dict[Tuple[str, str], List[Dict[str, Any]]] = {}
    for r in rows:
        groups.setdefault((r["chr"], r["direction"]), []).append(r)
    regions = []
    for (chrom, direction), calls in groups.items():
        calls = sorted(calls, key=lambda r: r["start"])
        current: List[Dict[str, Any]] = []
        end = -1
        for call in calls + [None]:
            if call is not None and call["start"] <= end:
                current.append(call)
                end = max(end, call["end"])
                continue
            if current:
                fams = {c["FID"] for c in current}
                if len(fams) >= minimum_families:
                    start = min(c["start"] for c in current)
                    genes = dict.fromkeys(g for c in current for g in c["_genes"])
                    regions.append(
                        {
                            "_key": f"{chrom}:{start}-{end}:{direction}",
                            "region": f"{chrom}:{start}-{end}",
                            "chr": chrom,
                            "start": start,
                            "direction": direction,
                            "size_mb": round((end - start) / 1e6, 2),
                            "n_families": len(fams),
                            "n_samples": len(
                                {(c["FID"], c["sample"]) for c in current}
                            ),
                            "families": ", ".join(sorted(fams)),
                            "exonic_symbol": ",".join(genes),
                            "_genes": list(genes),
                        }
                    )
            if call is not None:
                current, end = [call], call["end"]
    regions.sort(
        key=lambda g: (
            -g["n_families"],
            order.get(norm_chrom(g["chr"]), 99),
            g["start"],
        )
    )
    return regions
