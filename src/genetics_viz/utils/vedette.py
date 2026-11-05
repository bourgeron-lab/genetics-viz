"""Vedette SV calls: discovery, reading, one row per carrier, filters.

Vedette (Alex's SV caller) writes, per family, under ``svs/vedette/`` (format:
``docs/EXPORT.md`` of the vedette repository, from v1.4.0):

- ``<FID>.vedette.tsv``: one row per row of Vedette's events page: cohort sites
  carried by the same members collapsed into one event (``row_id``,
  ``site_ids``), plus sub-threshold, mosaic and chromosomal rows (``kind``).
  Coded columns (``impact_level``, ``support``, ``own_call``), complete gene
  lists (``exonic_genes``, ``genic_genes``, ``;``-separated), ``carriers``,
  then ``GT_ CN_ GQ_ INH_ DNQ_<id>`` per member. A missing value is an empty
  field. ``start``/``end`` are the first and last affected bases, 1-based,
  both included; ``vcf_pos = start - 1``.
- ``<FID>.vedette.qc.tsv``: one row per member, the per-sample QC (wave score,
  carrier load, de novo candidates, karyotype, ``flags``).
- ``<FID>.vedette.json``: ``schema_version``, ``vedette_version``,
  ``model_version``, cohort, member order.

and per sample, under ``samples/.../<ID>/svs/vedette/``, the sample's own calls
before merging (``<ID>.filtered.vcf.gz``, every FILTER kept) with their
evidence (SOURCE, IMPRECISE, CIPOS, SQ, PE, SR, CG). VCFs from v1.4.0 carry
``##vedette_coords=padded`` (``POS`` the base before the event, ``END`` its
last base); older ones have ``POS`` on the first affected base and ``END`` one
past the last.

Exports before v1.4.0 (one row per PASS site, ``site_id``, French ``impact``
label, truncated ``genes``, an ``inheritance`` column, no JSON) still open:
columns are read by name, a missing one leaves its field empty, unknown ones
are kept.

Row ids are renumbered when the cohort is merged again, so validations are
keyed by locus and type (``chr:start-end:del``), as WisecondorX's, and keep the
row id only as a note. Wording and colours are in ``config/vedette.yaml``.
"""

from __future__ import annotations

import gzip
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import polars as pl
import yaml

from genetics_viz.utils.genesets import cell_matches
from genetics_viz.utils.sharding import get_family_path, get_sample_path

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "vedette.yaml"


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


VEDETTE_CONFIG: Dict[str, Any] = _load_config()

#: Row columns the table knows (v1.4 and earlier); any other non-member
#: column is kept raw.
ROW_COLUMNS = (
    "row_id", "kind", "site_ids", "site_id", "contig", "start", "end", "vcf_pos",
    "svtype", "size", "filter", "rescue", "n_sites", "sqmax", "cohort_af",
    "cohort_carrier_frac", "cohort_ncarrier", "ref_carrier_frac", "lcr_ovl",
    "segdup_ovl", "support", "own_call", "mosaic_cn", "impact_score",
    "impact_level", "impact_reason", "impact", "exonic_genes", "exonic_gene_ids",
    "genic_genes", "genic_gene_ids", "carriers", "qc_alerts", "genes", "inheritance",
)  # fmt: skip
_MEMBER_RE = re.compile(r"^(GT|CN|GQ|INH|DNQ)_(.+)$")
_IMPACT_RE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d+)?)\s*(?:[—–-]\s*)?([^·]*)(?:·\s*(.*))?$"
)
_MORE_GENES_RE = re.compile(r"^\(\+(\d+)\)$")
_NULLS = ("", ".", "NA", "nan")
_ALT_GT = ("0/1", "1/0", "1/1", "0|1", "1|0", "1|1", "1", "./1", "1/.")
_KARYOTYPE_PREFIX = "karyotype_"
#: VCF header key of Vedette files in the VCF-spec convention (from v1.4.0).
COORDS_KEY = "vedette_coords"

UNKNOWN_COLOR = "#9e9e9e"


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------


def family_vedette_dir(data_dir: Path, family_id: str) -> Path:
    return get_family_path(data_dir, family_id) / VEDETTE_CONFIG.get(
        "family_subdir", "svs/vedette"
    )


def family_vedette_file(data_dir: Path, family_id: str) -> Path:
    suffix = VEDETTE_CONFIG.get("family_suffix", ".vedette.tsv")
    return family_vedette_dir(data_dir, family_id) / f"{family_id}{suffix}"


def sample_vedette_vcf(data_dir: Path, sample_id: str) -> Path:
    suffix = VEDETTE_CONFIG.get("sample_vcf_suffix", ".filtered.vcf.gz")
    subdir = VEDETTE_CONFIG.get("family_subdir", "svs/vedette")
    return get_sample_path(data_dir, sample_id) / subdir / f"{sample_id}{suffix}"


def probe_vedette(data_dir: Path, family_id: str) -> bool:
    return family_vedette_file(data_dir, family_id).exists()


def _family_file(data_dir: Path, family_id: str, key: str, default: str) -> Path:
    suffix = VEDETTE_CONFIG.get(key, default)
    return family_vedette_dir(data_dir, family_id) / f"{family_id}{suffix}"


def _vcf_source(path: Path) -> str:
    """The last ``##source=`` of a VCF header (``vedette-cohort-1.1.0``), or ``""``."""
    version = ""
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                if not line.startswith("##"):
                    break
                if line.startswith("##source="):
                    version = line.split("=", 1)[1].strip()
    except (OSError, EOFError, UnicodeDecodeError):
        logger.warning("Cannot read the header of %s", path, exc_info=True)
    return version


@dataclass
class ExportInfo:
    """``<FID>.vedette.json``: what wrote the family export."""

    schema_version: Optional[int] = None
    vedette_version: str = ""
    model_version: str = ""
    cohort: str = ""
    #: Member order of the ``GT_/CN_/GQ_/INH_/DNQ_<id>`` columns.
    members: List[str] = field(default_factory=list)
    generated: str = ""
    #: ``##source=`` of the family VCF, for exports without a JSON (< v1.4.0).
    legacy_source: str = ""

    @property
    def caller(self) -> str:
        """What validations record as ``CallerVersion``."""
        if self.vedette_version:
            return f"vedette-{self.vedette_version}"
        return self.legacy_source

    @property
    def newer_schema(self) -> bool:
        supported = int(VEDETTE_CONFIG.get("supported_schema_version", 1))
        return self.schema_version is not None and self.schema_version > supported


def read_export_info(data_dir: Path, family_id: str) -> ExportInfo:
    """The family's JSON; without one, the caller version from the family VCF."""
    path = _family_file(data_dir, family_id, "family_json_suffix", ".vedette.json")
    if path.exists():
        try:
            raw = json.loads(path.read_text())
        except (OSError, ValueError):
            logger.warning("Cannot read %s", path, exc_info=True)
            raw = {}
        return ExportInfo(
            schema_version=_int(raw.get("schema_version")),
            vedette_version=str(raw.get("vedette_version") or ""),
            model_version=str(raw.get("model_version") or ""),
            cohort=str(raw.get("cohort") or ""),
            members=[str(m) for m in raw.get("members") or []],
            generated=str(raw.get("generated") or ""),
        )
    vcf = _family_file(data_dir, family_id, "family_vcf_suffix", ".vedette.vcf.gz")
    return ExportInfo(legacy_source=_vcf_source(vcf) if vcf.exists() else "")


def caller_version(data_dir: Path, family_id: str) -> str:
    """``vedette-<version>`` from the JSON, else the family VCF's ``##source=``."""
    return read_export_info(data_dir, family_id).caller


@dataclass
class MemberQC:
    """One row of ``<FID>.vedette.qc.tsv``."""

    sample: str
    wave_score: Optional[float] = None
    pass_sites_carried: Optional[int] = None
    carrier_excess: Optional[float] = None
    denovo_candidates: Optional[int] = None
    karyotype: str = ""
    flags: List[str] = field(default_factory=list)
    flags_text: str = ""


def read_member_qc(data_dir: Path, family_id: str) -> List[MemberQC]:
    """The family's per-member QC, ``[]`` when the file is missing or empty."""
    path = _family_file(data_dir, family_id, "family_qc_suffix", ".vedette.qc.tsv")
    if not path.exists():
        return []
    try:
        df = pl.read_csv(
            path, separator="\t", infer_schema=False, null_values=list(_NULLS),
            quote_char=None,
        )  # fmt: skip
    except (pl.exceptions.NoDataError, OSError):
        return []
    out = []
    for r in df.iter_rows(named=True):
        if not r.get("sample"):
            continue
        out.append(
            MemberQC(
                sample=r["sample"],
                wave_score=_float(r.get("wave_score")),
                pass_sites_carried=_int(r.get("pass_sites_carried")),
                carrier_excess=_float(r.get("carrier_excess")),
                denovo_candidates=_int(r.get("denovo_candidates")),
                karyotype=r.get("karyotype") or "",
                flags=_split_list(r.get("flags")),
                flags_text=r.get("flags_text") or "",
            )
        )
    return out


def qc_flag_label(code: str) -> str:
    return VEDETTE_CONFIG.get("qc_flags", {}).get(code, code)


_TABLE_CACHE: Dict[Tuple[str, float, int], pl.DataFrame] = {}


def clear_vedette_cache() -> None:
    _TABLE_CACHE.clear()


def read_vedette_table(path: Path) -> pl.DataFrame:
    """The family table, every column a string, cached on ``(path, mtime, size)``."""
    stat = path.stat()
    key = (str(path), stat.st_mtime, stat.st_size)
    cached = _TABLE_CACHE.get(key)
    if cached is not None:
        return cached
    try:
        df = pl.read_csv(
            path,
            separator="\t",
            infer_schema=False,
            null_values=list(_NULLS),
            quote_char=None,
        )
    except pl.exceptions.NoDataError:
        df = pl.DataFrame()
    _TABLE_CACHE.clear()  # one family at a time is enough
    _TABLE_CACHE[key] = df
    return df


def table_members(columns: Iterable[str], order: Sequence[str] = ()) -> List[str]:
    """Member IDs: the JSON's order, else as the ``GT_<id>`` (etc.) columns come."""
    found: List[str] = []
    for col in columns:
        m = _MEMBER_RE.match(col)
        if m and m.group(2) not in found:
            found.append(m.group(2))
    ordered = [m for m in order if m in found]
    return ordered + [m for m in found if m not in ordered]


# ---------------------------------------------------------------------------
# Cell parsing
# ---------------------------------------------------------------------------


def _float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _int(value: Any) -> Optional[int]:
    f = _float(value)
    return None if f is None else int(f)


def _split_list(value: Any) -> List[str]:
    """A ``;``-separated list cell (empty fields dropped)."""
    return [v.strip() for v in str(value or "").split(";") if v.strip()]


def impact_level_code(text: Optional[str]) -> str:
    """The configured level whose ``match`` names ``text`` (else ``text`` itself)."""
    t = (text or "").strip().lower()
    if not t:
        return ""
    for code, cfg in VEDETTE_CONFIG.get("impact_levels", {}).items():
        if t == code or t in (m.lower() for m in cfg.get("match", [])):
            return code
    return t


def parse_impact(text: Optional[str]) -> Tuple[Optional[float], str, str]:
    """``(score, level code, reason)`` from ``+0.90 — élevé · DEL exonique RERE``."""
    if not text:
        return None, "", ""
    m = _IMPACT_RE.match(text)
    if not m:
        return None, "", text.strip()
    return float(m.group(1)), impact_level_code(m.group(2)), (m.group(3) or "").strip()


def parse_inheritance(text: Optional[str]) -> Dict[str, str]:
    """``{sample: code}`` from ``S3:denovo; S4:maternal``."""
    out: Dict[str, str] = {}
    for part in (text or "").replace(",", ";").split(";"):
        sample, sep, code = part.strip().rpartition(":")
        if sep and sample:
            out[sample.strip()] = code.strip()
    return out


def split_site_genes(text: Optional[str]) -> Tuple[List[str], int]:
    """Gene symbols of a ``genes`` cell, and the ``(+N)`` more it left out."""
    genes: List[str] = []
    more = 0
    for part in (text or "").split(","):
        part = part.strip()
        if not part:
            continue
        # The last symbol may carry the count: "PANK4 (+5)".
        tokens = part.split()
        if tokens and _MORE_GENES_RE.match(tokens[-1]):
            more += int(_MORE_GENES_RE.match(tokens[-1]).group(1))
            tokens = tokens[:-1]
        symbol = " ".join(tokens).strip()
        if symbol and symbol not in genes:
            genes.append(symbol)
    return genes, more


def is_carrier(gt: Optional[str]) -> bool:
    return (gt or "") in _ALT_GT


def variant_key(contig: str, start: Any, end: Any, svtype: str) -> str:
    """The validation and diagnostic key, as WisecondorX's: ``chr:start-end:del``."""
    return f"{contig}:{start}-{end}:{(svtype or '').lower()}"


def level_info(code: str) -> Dict[str, Any]:
    cfg = VEDETTE_CONFIG.get("impact_levels", {}).get(code)
    if cfg is None:
        return {"label": code, "color": UNKNOWN_COLOR, "rank": -9}
    return {"label": cfg.get("label", code), "color": cfg["color"], "rank": cfg["rank"]}


def inheritance_info(code: str) -> Dict[str, str]:
    """Label and colour of an inheritance code, ``karyotype_<K>`` included."""
    cfg = VEDETTE_CONFIG.get("inheritance", {}).get(code)
    if cfg is not None:
        return {"label": cfg.get("label", code), "color": cfg["color"]}
    if code.startswith(_KARYOTYPE_PREFIX):
        k = VEDETTE_CONFIG.get("karyotype_inheritance", {})
        karyotype = code[len(_KARYOTYPE_PREFIX) :]
        return {
            "label": k.get("label", "karyotype {K}").replace("{K}", karyotype),
            "color": k.get("color", UNKNOWN_COLOR),
        }
    return {"label": code, "color": UNKNOWN_COLOR}


def _coded(section: str, code: str) -> Dict[str, str]:
    """Label and colour of a ``kinds`` / ``own_calls`` / ``support`` code."""
    cfg = VEDETTE_CONFIG.get(section, {}).get(code)
    if cfg is None:
        return {"label": code, "color": UNKNOWN_COLOR}
    return {"label": cfg.get("label", code), "color": cfg["color"]}


def kind_info(code: str) -> Dict[str, str]:
    return _coded("kinds", code)


def own_call_info(code: str) -> Dict[str, str]:
    return _coded("own_calls", code)


def support_info(code: str) -> Dict[str, str]:
    return _coded("support", code)


def type_color(svtype: str) -> str:
    return VEDETTE_CONFIG.get("types", {}).get(svtype, {}).get("color", UNKNOWN_COLOR)


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


def _impact(row: Dict[str, Any]) -> Tuple[Optional[float], str, str]:
    """``(score, level code, reason)``: the v1.4 columns, else the label parsed."""
    if row.get("impact_score") is not None or row.get("impact_level"):
        return (
            _float(row.get("impact_score")),
            impact_level_code(row.get("impact_level")),
            row.get("impact_reason") or "",
        )
    return parse_impact(row.get("impact"))


def _genes(row: Dict[str, Any]) -> Tuple[List[str], List[str], str]:
    """``(exonic, every gene, display)``: the v1.4 lists, else the ``genes`` cell.

    Before v1.4 the list was cut after six names and did not tell exonic genes
    apart; the display keeps its ``(+N)``.
    """
    if "exonic_genes" in row or "genic_genes" in row:
        exonic = _split_list(row.get("exonic_genes"))
        genic = [g for g in _split_list(row.get("genic_genes")) if g not in exonic]
        genes = exonic + genic
        return exonic, genes, ", ".join(genes)
    genes, more = split_site_genes(row.get("genes"))
    return [], genes, ", ".join(genes) + (f" (+{more})" if more else "")


def _inheritance(row: Dict[str, Any], members: Sequence[str]) -> Dict[str, str]:
    """``{member: code}``: the ``INH_<id>`` columns, else ``inheritance``."""
    if any(f"INH_{m}" in row for m in members):
        return {m: row[f"INH_{m}"] for m in members if row.get(f"INH_{m}")}
    return parse_inheritance(row.get("inheritance"))


def _carriers(row: Dict[str, Any], members: Sequence[str]) -> List[str]:
    """The ``carriers`` column (Vedette's: GT ≥ 0/1 and GQ ≥ its minimum);
    before v1.4, every member with an alternate genotype."""
    if "carriers" in row:
        return [m for m in _split_list(row.get("carriers"))]
    return [m for m in members if is_carrier(row.get(f"GT_{m}"))]


def carrier_rows(
    df: pl.DataFrame,
    *,
    labels: Optional[Dict[str, str]] = None,
    members: Sequence[str] = (),
) -> List[Dict[str, Any]]:
    """One row per family row and carrier.

    ``labels`` names members in the Family badges' tooltips (``father``…);
    ``members`` is the JSON's member order. Badge and validation fields are
    added later, by the tab.
    """
    if df.is_empty():
        return []
    members = table_members(df.columns, members)
    extra = [c for c in df.columns if c not in ROW_COLUMNS and not _MEMBER_RE.match(c)]
    gt_colors = VEDETTE_CONFIG.get("genotype_colors", {})
    labels = labels or {}
    rows: List[Dict[str, Any]] = []
    for site in df.iter_rows(named=True):
        contig = site.get("contig") or ""
        start, end = _int(site.get("start")), _int(site.get("end"))
        svtype = (site.get("svtype") or "").upper()
        if start is None or end is None:
            continue
        row_id = site.get("row_id") or site.get("site_id") or ""
        kind = site.get("kind") or "site"
        score, level, reason = _impact(site)
        exonic, genes, genes_display = _genes(site)
        inheritance = _inheritance(site, members)
        size = _int(site.get("size"))
        support, own = site.get("support") or "", site.get("own_call") or ""
        family = []
        for m in members:
            gt = site.get(f"GT_{m}") or ""
            cn, gq = site.get(f"CN_{m}"), site.get(f"GQ_{m}")
            who = f"{m} ({labels[m]})" if labels.get(m) else m
            tip = f"{who}: GT {gt or '-'}, CN {cn or '-'}, GQ {gq or '-'}"
            if inheritance.get(m):
                tip += f", {inheritance_info(inheritance[m])['label']}"
                if site.get(f"DNQ_{m}"):
                    tip += f" (DNQ {site[f'DNQ_{m}']})"
            family.append(
                {
                    "label": f"{gt or '-'}·{cn or '-'}",
                    "color": gt_colors.get(gt, gt_colors.get("missing", "#f3f4f6")),
                    "textColor": "white" if is_carrier(gt) else "black",
                    "tooltip": tip,
                }
            )
        base = {
            "site_id": row_id,
            "site_ids": _split_list(site.get("site_ids"))
            or ([row_id] if row_id else []),
            "kind": kind_info(kind)["label"],
            "_kind": kind,
            "_kind_color": kind_info(kind)["color"],
            "chr:start-end": f"{contig}:{start}-{end}",
            "contig": contig,
            "start": start,
            "end": end,
            "svtype": svtype,
            "_type_color": type_color(svtype),
            "size": size if size is not None else end - start + 1,
            "filter": site.get("filter") or "",
            "rescue": site.get("rescue") or "",
            "n_sites": _int(site.get("n_sites")),
            "sqmax": _float(site.get("sqmax")),
            "cohort_af": _float(site.get("cohort_af")),
            "cohort_carrier_frac": _float(site.get("cohort_carrier_frac")),
            "cohort_ncarrier": _int(site.get("cohort_ncarrier")),
            "ref_carrier_frac": _float(site.get("ref_carrier_frac")),
            "lcr_ovl": _float(site.get("lcr_ovl")),
            "segdup_ovl": _float(site.get("segdup_ovl")),
            "support": support_info(support)["label"] if support else "",
            "_support": support,
            "_support_color": support_info(support)["color"],
            "own_call": own_call_info(own)["label"] if own else "",
            "_own_call": own,
            "_own_call_color": own_call_info(own)["color"],
            "mosaic_cn": _float(site.get("mosaic_cn")),
            "qc_alerts": site.get("qc_alerts") or "",
            "genes": genes_display,
            "_genes": genes,
            "_exonic": exonic,
            "impact_score": score,
            "impact_level": level_info(level)["label"] if level else "",
            "_impact_level": level,
            "_impact_color": level_info(level)["color"],
            "impact_reason": reason,
            "family": " ".join(b["label"] for b in family),
            "family_badges": family,
            "_variant_key": variant_key(contig, start, end, svtype),
            **{c: site.get(c) for c in extra},
        }
        for m in _carriers(site, members):
            code = inheritance.get(m, "")
            gt = site.get(f"GT_{m}") or ""
            rows.append(
                {
                    **base,
                    "id": f"{row_id}|{m}",
                    "sample": m,
                    "GT": gt,
                    "CN": _int(site.get(f"CN_{m}")),
                    "GQ": _int(site.get(f"GQ_{m}")),
                    "DNQ": _int(site.get(f"DNQ_{m}")),
                    "inheritance": inheritance_info(code)["label"] if code else "",
                    "_inheritance": code,
                    "_inheritance_color": inheritance_info(code)["color"],
                }
            )
    return rows


# ---------------------------------------------------------------------------
# Filters
# ---------------------------------------------------------------------------


def impact_level_codes() -> List[str]:
    """Configured levels, most severe first."""
    levels = VEDETTE_CONFIG.get("impact_levels", {})
    return sorted(levels, key=lambda c: -levels[c].get("rank", 0))


def inheritance_codes() -> List[str]:
    return list(VEDETTE_CONFIG.get("inheritance", {}))


def kind_codes() -> List[str]:
    return list(VEDETTE_CONFIG.get("kinds", {}))


def own_call_codes() -> List[str]:
    return list(VEDETTE_CONFIG.get("own_calls", {}))


@dataclass
class VedetteFilters:
    """The tab's filters. ``None`` / empty lists filter nothing, and a row
    missing the value a filter reads (no GQ on a mosaic row) passes it."""

    types: List[str] = field(default_factory=lambda: ["DEL", "DUP"])
    kinds: List[str] = field(default_factory=kind_codes)
    #: ``own_call`` codes; empty = any.
    own_calls: List[str] = field(default_factory=list)
    min_size: Optional[int] = None
    max_size: Optional[int] = None
    #: Cohort carrier fraction.
    max_freq: Optional[float] = None
    min_gq: Optional[int] = None
    impact_levels: List[str] = field(default_factory=impact_level_codes)
    #: Inheritance codes; ``""`` stands for "none given" (founders).
    inheritance: List[str] = field(default_factory=list)
    genesets: List[str] = field(default_factory=list)
    validations: List[str] = field(default_factory=list)

    _PRESET_FIELDS = (
        "types", "kinds", "own_calls", "min_size", "max_size", "max_freq",
        "min_gq", "impact_levels", "inheritance",
    )  # fmt: skip

    @classmethod
    def preset(cls, name: str) -> "VedetteFilters":
        """The filters a ``presets`` entry of the config sets."""
        cfg = {
            **VEDETTE_CONFIG.get("defaults", {}),
            **VEDETTE_CONFIG.get("presets", {}).get(name, {}),
        }
        f = cls()
        f.types = list(cfg.get("types") or f.types)
        f.kinds = list(cfg.get("kinds") or f.kinds)
        f.own_calls = list(cfg.get("own_calls") or [])
        f.min_size = cfg.get("min_size")
        f.max_size = cfg.get("max_size")
        f.max_freq = cfg.get("max_freq")
        f.min_gq = cfg.get("min_gq")
        f.impact_levels = list(cfg.get("impact_levels") or f.impact_levels)
        f.inheritance = list(cfg.get("inheritance") or [])
        return f

    def matches_preset(self, name: str) -> bool:
        p = VedetteFilters.preset(name)
        for attr in self._PRESET_FIELDS:
            mine, theirs = getattr(self, attr), getattr(p, attr)
            if isinstance(mine, list):
                if set(mine) != set(theirs):
                    return False
            elif mine != theirs:
                return False
        return True


def preset_names() -> List[str]:
    return list(VEDETTE_CONFIG.get("presets", {}))


def default_preset() -> str:
    names = preset_names()
    return names[0] if names else ""


def apply_filters(
    rows: Sequence[Dict[str, Any]],
    f: VedetteFilters,
    *,
    members: Optional[Iterable[str]] = None,
    gene_symbols: Optional[frozenset] = None,
) -> List[Dict[str, Any]]:
    """Rows of the selected members passing every filter but validation.

    The validation filter needs the row's status, added by the tab after
    loading the curation files; it applies ``filter_by_status`` itself.
    """
    wanted_members = set(members) if members is not None else None
    types, kinds, owns = set(f.types), set(f.kinds), set(f.own_calls)
    all_kinds = kinds >= set(kind_codes())
    levels = set(f.impact_levels)
    all_levels = levels >= set(impact_level_codes())
    inherit = set(f.inheritance)
    out = []
    for r in rows:
        if wanted_members is not None and r["sample"] not in wanted_members:
            continue
        if types and r["svtype"] not in types:
            continue
        if not all_kinds and r.get("_kind", "site") not in kinds:
            continue
        if owns and r.get("_own_call", "") not in owns:
            continue
        size = r.get("size") or 0
        if f.min_size is not None and size < f.min_size:
            continue
        if f.max_size is not None and size > f.max_size:
            continue
        freq = r.get("cohort_carrier_frac")
        if f.max_freq is not None and freq is not None and freq > f.max_freq:
            continue
        gq = r.get("GQ")
        if f.min_gq is not None and gq is not None and gq < f.min_gq:
            continue
        if not all_levels and r.get("_impact_level") not in levels:
            continue
        if inherit and r.get("_inheritance", "") not in inherit:
            continue
        if gene_symbols is not None and not cell_matches(
            ",".join(r.get("_genes", [])), gene_symbols
        ):
            continue
        out.append(r)
    return out


#: Link keys of the tab, prefixed ``v``.
LINK_KEYS = (
    "vpre", "vtype", "vkind", "vown", "vsize", "vfreq", "vgq", "vimp", "vinh",
    "vgs", "vval",
)  # fmt: skip


def _range(lo: Optional[int], hi: Optional[int]) -> Optional[str]:
    if lo is None and hi is None:
        return None
    return f"{'' if lo is None else lo}-{'' if hi is None else hi}"


def filters_to_params(f: VedetteFilters) -> Dict[str, Any]:
    """What a link records: the preset the filters match, else every value.

    Genesets and validation statuses are written either way: no preset sets them.
    """
    out: Dict[str, Any] = {}
    preset = next((p for p in preset_names() if f.matches_preset(p)), None)
    if preset is not None:
        out["vpre"] = preset
    else:
        out["vpre"] = "custom"
        out["vtype"] = f.types
        out["vkind"] = f.kinds
        out["vown"] = f.own_calls
        out["vsize"] = _range(f.min_size, f.max_size)
        out["vfreq"] = f.max_freq
        out["vgq"] = f.min_gq
        out["vimp"] = f.impact_levels
        out["vinh"] = [c or "-none" for c in f.inheritance]
    out["vgs"] = f.genesets
    out["vval"] = f.validations
    return out


def filters_from_params(get: Any, get_list: Any) -> VedetteFilters:
    """The filters a link sets (``get`` / ``get_list`` as :class:`PageLink`'s)."""
    preset = get("vpre") or default_preset()
    if preset != "custom":
        f = VedetteFilters.preset(
            preset if preset in preset_names() else default_preset()
        )
    else:
        f = VedetteFilters()
        if (types := get_list("vtype")) is not None:
            f.types = [t for t in types if t]
        if (kinds := get_list("vkind")) is not None:
            f.kinds = [k for k in kinds if k]
        if (owns := get_list("vown")) is not None:
            f.own_calls = [o for o in owns if o]
        if raw := get("vsize"):
            lo, _, hi = raw.partition("-")
            f.min_size, f.max_size = _int(lo), _int(hi)
        f.max_freq = _float(get("vfreq"))
        f.min_gq = _int(get("vgq"))
        if (levels := get_list("vimp")) is not None:
            f.impact_levels = levels
        if (inh := get_list("vinh")) is not None:
            f.inheritance = ["" if c == "-none" else c for c in inh]
    f.genesets = get_list("vgs") or []
    f.validations = get_list("vval") or []
    return f


# ---------------------------------------------------------------------------
# The carrier's own call
# ---------------------------------------------------------------------------


@dataclass
class OwnCall:
    """The carrier's own call matching a row, from its ``filtered.vcf.gz``.

    ``start``/``end`` are the first and last affected bases, 1-based, as the
    family export's, whichever convention the VCF was written in.
    """

    call_id: str
    start: int
    end: int
    filters: List[str]
    source: str
    precise: bool
    cipos: Optional[Tuple[int, int]]
    ciend: Optional[Tuple[int, int]]
    model: str
    overlap: float
    fields: Dict[str, Any]

    @property
    def passed(self) -> bool:
        return self.filters in ([], ["PASS"])


def _info(rec: Any, key: str, default: Any = None) -> Any:
    """An INFO value; ``default`` when the header does not declare the key
    (pysam raises on those)."""
    if key not in rec.header.info:
        return default
    return rec.info.get(key, default)


def _reciprocal(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    inter = min(a[1], b[1]) - max(a[0], b[0])
    if inter <= 0:
        return 0.0
    return inter / max(a[1] - a[0], b[1] - b[0], 1)


def find_own_call(
    vcf: Path, chrom: str, start: int, end: int, svtype: str
) -> Optional[OwnCall]:
    """The sample's call of this type overlapping the row best, if any.

    ``start``/``end`` are 1-based inclusive (the family export's). The VCF's
    records are read in their own convention: ``##vedette_coords=padded``
    (from v1.4.0) has ``POS`` on the base before the event and ``END`` on its
    last base; older files have ``POS`` on the first affected base and ``END``
    one past the last. Blocking, but one tabix query; ``None`` without a VCF,
    an index or a match.
    """
    if not vcf.exists() or not (
        Path(f"{vcf}.tbi").exists() or Path(f"{vcf}.csi").exists()
    ):
        return None
    import pysam

    min_overlap = float(VEDETTE_CONFIG.get("own_call_min_overlap", 0.5))
    best: Optional[OwnCall] = None
    try:
        with pysam.VariantFile(str(vcf)) as f:
            padded = any(r.key == COORDS_KEY for r in f.header.records)
            try:
                records = f.fetch(chrom, max(0, start - 2), end + 1)
            except ValueError:
                return None
            for rec in records:
                rtype = str(_info(rec, "SVTYPE", "")).upper()
                if rtype != svtype.upper():
                    continue
                # pysam's stop is the END field for a symbolic allele.
                r_start, r_end = (
                    (rec.pos + 1, rec.stop) if padded else (rec.pos, rec.stop - 1)
                )
                ov = _reciprocal((r_start - 1, r_end), (start - 1, end))
                if ov < min_overlap or (best is not None and ov <= best.overlap):
                    continue
                sample = rec.samples[0] if len(rec.samples) else {}
                fields = {}
                for key in ("GT", "GQ", "SQ", "CN", "CNQ", "RD", "PE", "SR", "CG"):
                    try:
                        value = sample.get(key)
                    except (KeyError, ValueError):
                        continue
                    if key == "GT" and value is not None:
                        value = "/".join("." if a is None else str(a) for a in value)
                    if value is not None:
                        fields[key] = value
                cipos, ciend = _info(rec, "CIPOS"), _info(rec, "CIEND")
                best = OwnCall(
                    call_id=rec.id or "",
                    start=r_start,
                    end=r_end,
                    filters=list(rec.filter.keys()),
                    source=str(_info(rec, "SOURCE", "")),
                    precise=not _info(rec, "IMPRECISE", False),
                    cipos=tuple(cipos) if cipos else None,
                    ciend=tuple(ciend) if ciend else None,
                    model=str(_info(rec, "MODEL_VER", "")),
                    overlap=ov,
                    fields=fields,
                )
    except (OSError, ValueError):
        logger.warning("Cannot read %s", vcf, exc_info=True)
        return None
    return best


def assume_precise(size: Optional[int]) -> bool:
    """Without an own call: short sites have base-pair breakpoints."""
    return (size or 0) < int(VEDETTE_CONFIG.get("assume_precise_below_bp", 1000))
