"""Vedette SV calls: discovery, reading, one row per carrier, filters.

Vedette (Alex's SV caller) writes, per family, under ``svs/vedette/``:

- ``<FID>.vedette.tsv``: one row per PASS site. ``site_id contig start end
  svtype size filter cohort_af cohort_carrier_frac cohort_ncarrier genes
  impact``, then ``GT_<id> CN_<id> GQ_<id>`` per member, then
  ``inheritance`` (``S3:denovo; S4:maternal``, children only).
- ``<FID>.vedette.vcf.gz``: every site, PASS or not. Read here only for the
  caller version (``##source=``).

and per sample, under ``samples/.../<ID>/svs/vedette/``, the sample's own calls
before merging (``<ID>.filtered.vcf.gz``, every FILTER kept) with their
evidence (SOURCE, IMPRECISE, CIPOS, SQ, PE, SR, CG). A carrier may have no own
call: the cohort step re-genotypes every sample at every site.

``site_id`` is positional (``vedette_site_N``) and renumbered whenever the
cohort is merged again, so validations are keyed by locus and type
(``chr:start-end:del``), as WisecondorX's, and keep the site id only as a note.

Vedette changes quickly: columns are read by name, a missing one leaves its
field empty, unknown ones are kept. ``impact`` packs a score, a level and a
reason into one string (``+0.90 — élevé · DEL exonique RERE``); separate
``impact_score`` / ``impact_level`` / ``impact_reason`` columns win when present.
Wording and colours are in ``config/vedette.yaml``.
"""

from __future__ import annotations

import gzip
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

#: The site columns the table knows; any other non-member column is kept raw.
SITE_COLUMNS = (
    "site_id", "contig", "start", "end", "svtype", "size", "filter", "cohort_af",
    "cohort_carrier_frac", "cohort_ncarrier", "genes", "impact", "impact_score",
    "impact_level", "impact_reason", "inheritance",
)  # fmt: skip
_MEMBER_FIELDS = ("GT", "CN", "GQ")
_MEMBER_RE = re.compile(r"^(GT|CN|GQ)_(.+)$")
_IMPACT_RE = re.compile(
    r"^\s*([+-]?\d+(?:\.\d+)?)\s*(?:[—–-]\s*)?([^·]*)(?:·\s*(.*))?$"
)
_MORE_GENES_RE = re.compile(r"^\(\+(\d+)\)$")
_NULLS = ("", ".", "NA", "nan")
_ALT_GT = ("0/1", "1/0", "1/1", "0|1", "1|0", "1|1", "1", "./1", "1/.")

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


def caller_version(data_dir: Path, family_id: str) -> str:
    """The last ``##source=`` of the family VCF (``vedette-cohort-1.1.0``), or ``""``."""
    suffix = VEDETTE_CONFIG.get("family_vcf_suffix", ".vedette.vcf.gz")
    path = family_vedette_dir(data_dir, family_id) / f"{family_id}{suffix}"
    if not path.exists():
        return ""
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


def table_members(columns: Iterable[str]) -> List[str]:
    """Member IDs in column order, from the ``GT_<id>`` columns."""
    out: List[str] = []
    for col in columns:
        m = _MEMBER_RE.match(col)
        if m and m.group(2) not in out:
            out.append(m.group(2))
    return out


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
        return {"label": code, "color": UNKNOWN_COLOR, "rank": -1}
    return {"label": cfg.get("label", code), "color": cfg["color"], "rank": cfg["rank"]}


def inheritance_info(code: str) -> Dict[str, str]:
    cfg = VEDETTE_CONFIG.get("inheritance", {}).get(code)
    if cfg is None:
        return {"label": code, "color": UNKNOWN_COLOR}
    return {"label": cfg.get("label", code), "color": cfg["color"]}


def type_color(svtype: str) -> str:
    return VEDETTE_CONFIG.get("types", {}).get(svtype, {}).get("color", UNKNOWN_COLOR)


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


def carrier_rows(
    df: pl.DataFrame,
    *,
    labels: Optional[Dict[str, str]] = None,
) -> List[Dict[str, Any]]:
    """One row per site and carrier (a member with an alternate genotype).

    ``labels`` names members in the Family badges' tooltips (``father``…).
    Badge and validation fields are added later, by the tab.
    """
    if df.is_empty():
        return []
    members = table_members(df.columns)
    extra = [c for c in df.columns if c not in SITE_COLUMNS and not _MEMBER_RE.match(c)]
    gt_colors = VEDETTE_CONFIG.get("genotype_colors", {})
    labels = labels or {}
    rows: List[Dict[str, Any]] = []
    for site in df.iter_rows(named=True):
        contig = site.get("contig") or ""
        start, end = _int(site.get("start")), _int(site.get("end"))
        svtype = (site.get("svtype") or "").upper()
        if start is None or end is None:
            continue
        if site.get("impact_score") is not None or site.get("impact_level"):
            score = _float(site.get("impact_score"))
            level = impact_level_code(site.get("impact_level"))
            reason = site.get("impact_reason") or ""
        else:
            score, level, reason = parse_impact(site.get("impact"))
        inheritance = parse_inheritance(site.get("inheritance"))
        genes, more_genes = split_site_genes(site.get("genes"))
        size = _int(site.get("size"))
        family = []
        for m in members:
            gt = site.get(f"GT_{m}") or "./."
            cn, gq = site.get(f"CN_{m}"), site.get(f"GQ_{m}")
            who = f"{m} ({labels[m]})" if labels.get(m) else m
            tip = f"{who}: GT {gt}, CN {cn or '-'}, GQ {gq or '-'}"
            if inheritance.get(m):
                tip += f", {inheritance_info(inheritance[m])['label']}"
            family.append(
                {
                    "label": f"{gt}·{cn or '-'}",
                    "color": gt_colors.get(gt, gt_colors.get("missing", "#f3f4f6")),
                    "textColor": "black" if not is_carrier(gt) else "white",
                    "tooltip": tip,
                }
            )
        base = {
            "site_id": site.get("site_id") or "",
            "chr:start-end": f"{contig}:{start}-{end}",
            "contig": contig,
            "start": start,
            "end": end,
            "svtype": svtype,
            "_type_color": type_color(svtype),
            "size": size if size is not None else end - start,
            "filter": site.get("filter") or "",
            "cohort_af": _float(site.get("cohort_af")),
            "cohort_carrier_frac": _float(site.get("cohort_carrier_frac")),
            "cohort_ncarrier": _int(site.get("cohort_ncarrier")),
            "genes": ", ".join(genes) + (f" (+{more_genes})" if more_genes else ""),
            "_genes": genes,
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
        for m in members:
            gt = site.get(f"GT_{m}")
            if not is_carrier(gt):
                continue
            code = inheritance.get(m, "")
            rows.append(
                {
                    **base,
                    "id": f"{base['site_id']}|{m}",
                    "sample": m,
                    "GT": gt,
                    "CN": _int(site.get(f"CN_{m}")),
                    "GQ": _int(site.get(f"GQ_{m}")),
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


@dataclass
class VedetteFilters:
    """The tab's filters. ``None`` / empty lists filter nothing."""

    types: List[str] = field(default_factory=lambda: ["DEL", "DUP"])
    min_size: Optional[int] = None
    max_size: Optional[int] = None
    max_af: Optional[float] = None
    min_gq: Optional[int] = None
    impact_levels: List[str] = field(default_factory=impact_level_codes)
    #: Inheritance codes; ``""`` stands for "none given" (founders).
    inheritance: List[str] = field(default_factory=list)
    genesets: List[str] = field(default_factory=list)
    validations: List[str] = field(default_factory=list)

    @classmethod
    def preset(cls, name: str) -> "VedetteFilters":
        """The filters a ``presets`` entry of the config sets."""
        cfg = {
            **VEDETTE_CONFIG.get("defaults", {}),
            **VEDETTE_CONFIG.get("presets", {}).get(name, {}),
        }
        f = cls()
        f.types = list(cfg.get("types") or f.types)
        f.min_size = cfg.get("min_size")
        f.max_size = cfg.get("max_size")
        f.max_af = cfg.get("max_af")
        f.min_gq = cfg.get("min_gq")
        f.impact_levels = list(cfg.get("impact_levels") or f.impact_levels)
        f.inheritance = list(cfg.get("inheritance") or [])
        return f

    def matches_preset(self, name: str) -> bool:
        p = VedetteFilters.preset(name)
        return (
            set(self.types) == set(p.types)
            and (self.min_size, self.max_size, self.max_af, self.min_gq)
            == (p.min_size, p.max_size, p.max_af, p.min_gq)
            and set(self.impact_levels) == set(p.impact_levels)
            and set(self.inheritance) == set(p.inheritance)
        )


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
    types = set(f.types)
    levels = set(f.impact_levels)
    all_levels = levels >= set(impact_level_codes())
    inherit = set(f.inheritance)
    out = []
    for r in rows:
        if wanted_members is not None and r["sample"] not in wanted_members:
            continue
        if types and r["svtype"] not in types:
            continue
        size = r.get("size") or 0
        if f.min_size is not None and size < f.min_size:
            continue
        if f.max_size is not None and size > f.max_size:
            continue
        if f.max_af is not None and (r.get("cohort_af") or 0.0) > f.max_af:
            continue
        if f.min_gq is not None and (r.get("GQ") or 0) < f.min_gq:
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
LINK_KEYS = ("vpre", "vtype", "vsize", "vaf", "vgq", "vimp", "vinh", "vgs", "vval")


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
        out["vsize"] = _range(f.min_size, f.max_size)
        out["vaf"] = f.max_af
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
        if raw := get("vsize"):
            lo, _, hi = raw.partition("-")
            f.min_size, f.max_size = _int(lo), _int(hi)
        f.max_af = _float(get("vaf"))
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
    """The carrier's own call matching a site, from its ``filtered.vcf.gz``."""

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


def _reciprocal(a: Tuple[int, int], b: Tuple[int, int]) -> float:
    inter = min(a[1], b[1]) - max(a[0], b[0])
    if inter <= 0:
        return 0.0
    return inter / max(a[1] - a[0], b[1] - b[0], 1)


def find_own_call(
    vcf: Path, chrom: str, start: int, end: int, svtype: str
) -> Optional[OwnCall]:
    """The sample's call of this type overlapping the site best, if any.

    Blocking, but one tabix query; ``None`` without a VCF, an index or a match.
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
            try:
                records = f.fetch(chrom, max(0, start - 1), end + 1)
            except ValueError:
                return None
            for rec in records:
                rtype = str(rec.info.get("SVTYPE", "")).upper()
                if rtype != svtype.upper():
                    continue
                r_end = int(rec.info.get("END", rec.stop))
                ov = _reciprocal((rec.pos, r_end), (start, end))
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
                cipos = rec.info.get("CIPOS") if "CIPOS" in rec.header.info else None
                ciend = rec.info.get("CIEND") if "CIEND" in rec.header.info else None
                best = OwnCall(
                    call_id=rec.id or "",
                    start=rec.pos,
                    end=r_end,
                    filters=list(rec.filter.keys()),
                    source=str(rec.info.get("SOURCE", "")),
                    precise="IMPRECISE" not in rec.info,
                    cipos=tuple(cipos) if cipos else None,
                    ciend=tuple(ciend) if ciend else None,
                    model=str(rec.info.get("MODEL_VER", "")),
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
