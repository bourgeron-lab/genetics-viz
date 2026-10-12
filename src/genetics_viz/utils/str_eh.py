"""Short tandem repeats (STR) from ExpansionHunter: discovery, reading, shaping.

The ghfc-ngs ``expansionhunter`` step writes, under ``str/expansionhunter/`` of
each level, files named by the catalog label (``expansionhunter_catalog_name``):

- sample: ``<id>.<label>.{vcf,json}`` and ``<id>.<label>_realigned.bam``;
- family: ``<FID>.<label>.genotypes.tsv``, one row per locus and member;
- cohort: ``<NAME>.<label>.genotypes.tsv`` (the family tables plus
  ``cohort_percentile``), ``.loci.tsv`` (per-locus founders' distribution and
  status counts) and ``.flags.tsv`` (the flagged genotypes).

Every table writes ``.`` for a missing value. ``allele1`` is the shorter
allele; a male's chrX locus has one allele and ``allele2`` is ``.``.
``father`` and ``mother`` are filled in only for parents genotyped in the
family, so a member with neither is a *founder* of the family's data - the
cohort statistics count only founders, so a family does not weigh an allele
several times over.

Classes, statuses and wording come from ``config/str.yaml``; the pipeline
assigns them, the app only draws them.
"""

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import polars as pl
import yaml

from .sharding import get_family_path, get_sample_path
from .tsv import read_tsv_or_none
from .validation_badges import STATUS_CFG

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "str.yaml"

EH_TOOL = "expansionhunter"

#: Statuses the pipeline flags (``flags.tsv``), most severe first.
FLAGGED = ("pathogenic", "carrier", "intermediate")

#: Columns that hold identifiers or labels and must stay strings. Left to
#: inference, ``family_id`` comes out as an integer (``88181200007``) and the
#: read then fails on ``AU-ANG-FER-MPX3``; ``genotype`` reads ``0`` and ``1``
#: for a male's chrX locus.
_STRING_COLUMNS = {
    "family_id", "sample", "father", "mother", "sex", "phenotype",
    "locus", "chrom", "motif", "filter", "genotype",
    "allele1_ci", "allele2_ci", "allele1_support", "allele2_support",
    "disease", "inheritance", "allele1_class", "allele2_class", "status",
}  # fmt: skip


# ---------------------------------------------------------------------------
# Display configuration
# ---------------------------------------------------------------------------


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


_CONFIG: Dict[str, Any] = _load_config()
STATUSES: Dict[str, Dict[str, Any]] = _CONFIG.get("statuses", {})
SUPPORT: Dict[str, Dict[str, str]] = _CONFIG.get("support", {})
COMMON_FRACTION: float = float(_CONFIG.get("common_in_cohort_fraction", 0.05))
REVIEWER_FILENAME: str = _CONFIG.get(
    "reviewer_filename", "{sample}.{label}.{locus}.svg"
)


def reload_str_config() -> None:
    """Reload the STR display configuration from YAML."""
    global _CONFIG, STATUSES, SUPPORT, COMMON_FRACTION, REVIEWER_FILENAME
    _CONFIG = _load_config()
    STATUSES = _CONFIG.get("statuses", {})
    SUPPORT = _CONFIG.get("support", {})
    COMMON_FRACTION = float(_CONFIG.get("common_in_cohort_fraction", 0.05))
    REVIEWER_FILENAME = _CONFIG.get("reviewer_filename", REVIEWER_FILENAME)


def status_info(status: Optional[str]) -> Dict[str, Any]:
    """Label, colour and rank of a class or status; ``None`` is unclassified."""
    key = status if status in STATUSES else "unclassified"
    info = STATUSES.get(key, {})
    return {
        "key": key,
        "label": info.get("label", key),
        "color": info.get("color", "#9e9e9e"),
        "rank": int(info.get("rank", -1)),
    }


def support_code(support: Optional[str]) -> str:
    """Short code for an allele's read support (``F``, ``IR``; spanning is blank)."""
    return SUPPORT.get(support or "", {}).get("code", "")


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def get_str_dir(base: Path, tool: str = EH_TOOL) -> Path:
    """Return ``<base>/str/<tool>/`` (may not exist).

    ``base`` is a family, sample or cohort directory, from
    :func:`get_family_path`, :func:`get_sample_path` or ``cohort.path``.
    """
    return base / "str" / tool


def get_family_eh_dir(data_dir: Path, family_id: str) -> Path:
    return get_str_dir(get_family_path(data_dir, family_id))


def get_sample_eh_dir(data_dir: Path, sample_id: str) -> Path:
    return get_str_dir(get_sample_path(data_dir, sample_id))


def get_cohort_eh_dir(cohort: Any) -> Path:
    return get_str_dir(cohort.path)


def eh_path(directory: Path, prefix: str, label: str, kind: str) -> Path:
    """``<directory>/<prefix>.<label>.<kind>.tsv``."""
    return directory / f"{prefix}.{label}.{kind}.tsv"


def find_eh_labels(directory: Path, prefix: str) -> List[str]:
    """Catalog labels with a genotype table in ``directory``, newest first.

    The label is the part between the prefix and ``.genotypes.tsv``; a lazy
    capture would cut a dotted label short, so the whole middle is taken.
    """
    if not directory.is_dir():
        return []
    pattern = re.compile(rf"^{re.escape(prefix)}\.(.+)\.genotypes\.tsv$")
    found: List[Tuple[float, str]] = []
    try:
        for entry in directory.iterdir():
            m = pattern.match(entry.name)
            if m:
                found.append((entry.stat().st_mtime, m.group(1)))
    except OSError:
        logger.warning("Cannot list %s", directory, exc_info=True)
        return []
    return [label for _, label in sorted(found, key=lambda t: (-t[0], t[1]))]


def probe_eh_data(data_dir: Path, family_id: str) -> bool:
    """Whether the family has an ExpansionHunter genotype table."""
    return bool(find_eh_labels(get_family_eh_dir(data_dir, family_id), family_id))


def probe_cohort_eh_data(cohort: Any) -> bool:
    """Whether the cohort has an ExpansionHunter genotype table."""
    return bool(find_eh_labels(get_cohort_eh_dir(cohort), cohort.name))


def reviewer_path(data_dir: Path, sample_id: str, label: str, locus: str) -> Path:
    """Where the REViewer image of one sample at one locus would be."""
    name = REVIEWER_FILENAME.format(sample=sample_id, label=label, locus=locus)
    return get_sample_eh_dir(data_dir, sample_id) / name


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_eh_table(path: Optional[Path]) -> Optional[pl.DataFrame]:
    """Read one ExpansionHunter TSV, ``.`` as null and identifiers as strings.

    Returns ``None`` when the path is ``None``, missing, or holds nothing.
    """
    if path is None or not path.exists():
        return None
    return read_tsv_or_none(
        path,
        null_values=["."],
        infer_schema_length=10000,
        string_columns=_STRING_COLUMNS,
    )


#: ``(path, mtime, size)`` -> frame. The cohort genotype table is several MB
#: and every family page of the cohort reads it; keyed on the stat, a rewrite
#: by the next pipeline run is picked up by itself.
_table_cache: Dict[Tuple[str, float, int], Optional[pl.DataFrame]] = {}


def read_eh_table_cached(path: Path) -> Optional[pl.DataFrame]:
    """:func:`read_eh_table`, cached on the file's stat."""
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (str(path), stat.st_mtime, stat.st_size)
    if key not in _table_cache:
        for old in [k for k in _table_cache if k[0] == key[0]]:
            del _table_cache[old]
        _table_cache[key] = read_eh_table(path)
    return _table_cache[key]


def clear_eh_cache() -> None:
    """Drop every cached table (tests, and after a data directory switch)."""
    _table_cache.clear()


@dataclass
class CohortContext:
    """A cohort's ExpansionHunter tables for one catalog label, either may be None."""

    name: str
    label: str
    loci: Optional[pl.DataFrame]
    genotypes: Optional[pl.DataFrame]

    @property
    def available(self) -> bool:
        return self.loci is not None or self.genotypes is not None

    @property
    def n_founders(self) -> Optional[int]:
        if self.loci is None or "n_founders" not in self.loci.columns:
            return None
        values = self.loci["n_founders"].drop_nulls()
        return int(values.max()) if len(values) else None


def load_cohort_context(cohort: Any, label: str) -> CohortContext:
    """Read the cohort's loci and genotype tables for ``label``.

    Blocking I/O over the data mount - call it through ``asyncio.to_thread``.
    A missing or unreadable table is ``None``: cohort context is optional.
    """
    directory = get_cohort_eh_dir(cohort)
    tables: Dict[str, Optional[pl.DataFrame]] = {}
    for kind in ("loci", "genotypes"):
        path = eh_path(directory, cohort.name, label, kind)
        try:
            tables[kind] = read_eh_table_cached(path)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            tables[kind] = None
    return CohortContext(cohort.name, label, tables["loci"], tables["genotypes"])


# ---------------------------------------------------------------------------
# Shaping
# ---------------------------------------------------------------------------


def _int(value: Any) -> Optional[int]:
    if value is None:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _fmt(value: Any) -> str:
    """A threshold or statistic for display: ``.`` for missing, no ``.0``."""
    if value is None:
        return "–"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _coverage(value: Any) -> str:
    try:
        return f"{float(value):.1f}x"
    except (TypeError, ValueError):
        return "no coverage"


def thresholds_text(row: Dict[str, Any]) -> str:
    """``≤20 · I 21 · P ≥24``, leaving out a threshold the locus lacks."""
    parts = []
    if row.get("normal_max") is not None:
        parts.append(f"N ≤{_fmt(row['normal_max'])}")
    if row.get("intermediate_min") is not None:
        parts.append(f"I ≥{_fmt(row['intermediate_min'])}")
    if row.get("pathogenic_min") is not None:
        parts.append(f"P ≥{_fmt(row['pathogenic_min'])}")
    return " · ".join(parts) or "no thresholds"


def longest_allele(row: Dict[str, Any]) -> Optional[int]:
    values = [
        v for v in (_int(row.get("allele1")), _int(row.get("allele2"))) if v is not None
    ]
    return max(values) if values else None


def allele_badges(row: Dict[str, Any], percentile: Any = None) -> List[Dict[str, str]]:
    """One badge per allele: repeat count plus support code, coloured by class."""
    badges = []
    for i in (1, 2):
        count = _int(row.get(f"allele{i}"))
        if count is None:
            continue
        support = row.get(f"allele{i}_support")
        code = support_code(support)
        info = status_info(row.get(f"allele{i}_class"))
        reads = (
            f"{row.get(f'allele{i}_spanning') or 0} spanning,"
            f" {row.get(f'allele{i}_flanking') or 0} flanking,"
            f" {row.get(f'allele{i}_inrepeat') or 0} in-repeat"
        )
        tip = [
            f"{count} units (CI {row.get(f'allele{i}_ci') or '–'})",
            info["label"],
            (support or "").lower() or "no support",
            reads,
            _coverage(row.get("locus_coverage")),
        ]
        if row.get("filter") and row.get("filter") != "PASS":
            tip.append(row["filter"])
        if percentile is not None:
            tip.append(f"cohort pct {_fmt(percentile)}")
        badges.append(
            {
                "label": f"{count}·{code}" if code else str(count),
                "color": info["color"],
                "tooltip": " · ".join(tip),
            }
        )
    if row.get("filter") and row.get("filter") != "PASS":
        badges.append(
            {
                "label": row["filter"],
                "color": "#ffffff",
                "tooltip": "ExpansionHunter filter: depth under 10x or few"
                " reads spanning a breakend",
            }
        )
    return badges


def validation_badge(status: str) -> Dict[str, str]:
    """The curators' aggregated status as a cell badge: symbol, status colour."""
    cfg = STATUS_CFG.get(status, {})
    return {
        "label": cfg.get("symbol") or status,
        "color": cfg.get("color", "#9e9e9e"),
        "tooltip": f"Curator validation: {status}",
    }


def is_transmission(row: Dict[str, Any]) -> bool:
    """A child's longest allele above both parents', at a non-normal status.

    An expansion on transmission - or an artefact; either deserves a look.
    """
    parents_max = _int(row.get("parents_max_allele"))
    longest = longest_allele(row)
    if parents_max is None or longest is None:
        return False
    return longest > parents_max and status_info(row.get("status"))["rank"] >= 1


def common_in_cohort(loci_row: Optional[Dict[str, Any]]) -> Optional[float]:
    """Share of the cohort flagged at this locus, when it passes the configured fraction."""
    if not loci_row:
        return None
    n = _int(loci_row.get("n_samples"))
    if not n:
        return None
    flagged = sum(_int(loci_row.get(f"n_{s}")) or 0 for s in FLAGGED)
    share = flagged / n
    return share if share >= COMMON_FRACTION else None


@dataclass(frozen=True)
class Member:
    """A family member as the genotype table describes them."""

    sample: str
    sex: Optional[str]
    phenotype: Optional[str]
    father: Optional[str]
    mother: Optional[str]

    @property
    def column(self) -> str:
        """Header text: ID, then sex and affection glyphs."""
        sex = {"male": "♂", "female": "♀"}.get(self.sex or "", "?")
        pheno = {"affected": "●", "unaffected": "○"}.get(self.phenotype or "", "")
        return f"{self.sample} {sex}{pheno}".strip()


def family_members(fam: pl.DataFrame) -> List[Member]:
    """Members in the table's order (pedigree order), each once."""
    seen: Dict[str, Member] = {}
    for r in fam.select(["sample", "sex", "phenotype", "father", "mother"]).to_dicts():
        if r["sample"] not in seen:
            seen[r["sample"]] = Member(
                r["sample"], r["sex"], r["phenotype"], r["father"], r["mother"]
            )
    return list(seen.values())


def member_key(index: int) -> str:
    """Column id of a member. Sample IDs may hold dots, which TanStack reads as paths."""
    return f"m{index}"


def family_locus_rows(
    fam: pl.DataFrame,
    members: List[Member],
    selected: List[str],
    context: Optional[CohortContext] = None,
    validations: Optional[Dict[Tuple[str, str], str]] = None,
) -> List[Dict[str, Any]]:
    """Pivot a family genotype table to one row per locus, most severe first.

    Each selected member gets ``m<i>_badges`` (one badge per allele) and
    ``m<i>_max`` (the longest allele, for sorting), ``i`` being their index in
    ``members`` so the column ids stay stable when the selection changes. The
    worst status, flags and cohort columns are computed over the selected
    members only. ``validations`` maps ``(locus, sample)`` to the curators'
    aggregated status, drawn as one more badge in the member's cell.
    """
    loci_by_name: Dict[str, Dict[str, Any]] = {}
    if context is not None and context.loci is not None:
        loci_by_name = {r["locus"]: r for r in context.loci.to_dicts()}
    percentiles: Dict[Tuple[str, str], Any] = {}
    if (
        context is not None
        and context.genotypes is not None
        and "cohort_percentile" in context.genotypes.columns
    ):
        family_ids = fam["family_id"].unique().to_list()
        sub = context.genotypes.filter(pl.col("family_id").is_in(family_ids))
        percentiles = {
            (r["sample"], r["locus"]): r["cohort_percentile"]
            for r in sub.select(["sample", "locus", "cohort_percentile"]).to_dicts()
        }

    index = {m.sample: i for i, m in enumerate(members)}
    wanted = [s for s in selected if s in index]

    by_locus: Dict[str, Dict[str, Dict[str, Any]]] = {}
    order: List[str] = []
    first: Dict[str, Dict[str, Any]] = {}
    for r in fam.to_dicts():
        locus = r["locus"]
        if locus not in by_locus:
            by_locus[locus] = {}
            order.append(locus)
            first[locus] = r
        by_locus[locus][r["sample"]] = r

    rows: List[Dict[str, Any]] = []
    for position, locus in enumerate(order):
        info = first[locus]
        row: Dict[str, Any] = {
            "id": locus,
            "order": position,
            "locus": locus,
            "chrom": info.get("chrom"),
            "pos": info.get("pos"),
            "end": info.get("end"),
            "position": f"{info.get('chrom')}:{info.get('pos')}-{info.get('end')}",
            "disease": info.get("disease") or "",
            "inheritance": info.get("inheritance") or "",
            "motif": info.get("motif") or "",
            "ref_units": info.get("ref_units"),
            "thresholds": thresholds_text(info),
        }
        worst = status_info(None)
        flags: List[Dict[str, str]] = []
        for sample in wanted:
            key = member_key(index[sample])
            g = by_locus[locus].get(sample)
            if g is None:
                row[f"{key}_badges"] = []
                row[f"{key}_max"] = None
                row[key] = "–"
                continue
            row[f"{key}_badges"] = allele_badges(g, percentiles.get((sample, locus)))
            validated = (validations or {}).get((locus, sample))
            if validated:
                row[f"{key}_badges"].append(validation_badge(validated))
            row[f"{key}_max"] = longest_allele(g)
            row[key] = ""
            s = status_info(g.get("status"))
            if s["rank"] > worst["rank"]:
                worst = s
            if is_transmission(g):
                flags.append(
                    {
                        "label": f"↑ {sample}",
                        "color": "#5e35b1",
                        "tooltip": f"{sample}'s longest allele ({longest_allele(g)})"
                        f" is above both parents' ({g.get('parents_max_allele')}):"
                        " an expansion on transmission, or an artefact",
                    }
                )
        row["worst_rank"] = worst["rank"]
        row["worst_key"] = worst["key"]
        row["worst"] = worst["label"]
        row["worst_badges"] = [
            {"label": worst["label"], "color": worst["color"], "tooltip": ""}
        ]

        loci_row = loci_by_name.get(locus)
        if loci_row:
            row["founders_median"] = loci_row.get("founders_median")
            row["founders_p999"] = loci_row.get("founders_p99.9")
            row["founders_max"] = loci_row.get("founders_max")
            row["cohort_flagged"] = sum(
                _int(loci_row.get(f"n_{s}")) or 0 for s in FLAGGED
            )
            share = common_in_cohort(loci_row)
            if share is not None:
                flags.insert(
                    0,
                    {
                        "label": "common in cohort",
                        "color": "#455a64",
                        "tooltip": f"{share:.0%} of the cohort is pathogenic, carrier"
                        " or intermediate here: suspect the locus, catalog or"
                        " thresholds before the individual",
                    },
                )
        else:
            row["founders_median"] = row["founders_p999"] = row["founders_max"] = None
            row["cohort_flagged"] = None
        row["flags_badges"] = flags
        row["flags"] = ""
        rows.append(row)

    rows.sort(key=lambda r: (-r["worst_rank"], r["order"]))
    return rows


def member_details(
    fam: pl.DataFrame, locus: str, context: Optional[CohortContext] = None
) -> List[Dict[str, Any]]:
    """Every member's full row at one locus, with the cohort percentile when known."""
    rows = fam.filter(pl.col("locus") == locus).to_dicts()
    if context is not None and context.genotypes is not None:
        fids = {r["family_id"] for r in rows}
        sub = context.genotypes.filter(
            (pl.col("locus") == locus) & pl.col("family_id").is_in(list(fids))
        )
        pct = {r["sample"]: r.get("cohort_percentile") for r in sub.to_dicts()}
        for r in rows:
            r["cohort_percentile"] = pct.get(r["sample"])
    return rows


def locus_distribution(genotypes: pl.DataFrame, locus: str) -> Dict[int, int]:
    """Founders' longest allele at one locus -> number of founders.

    The same quantity ``loci.tsv`` summarises and ``cohort_percentile`` ranks,
    so the plot agrees with both.
    """
    sub = genotypes.filter(
        (pl.col("locus") == locus)
        & pl.col("father").is_null()
        & pl.col("mother").is_null()
    ).select(["allele1", "allele2"])
    counts: Dict[int, int] = {}
    for r in sub.to_dicts():
        longest = longest_allele(r)
        if longest is not None:
            counts[longest] = counts.get(longest, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Cohort panel
# ---------------------------------------------------------------------------


def restrict_to(
    genotypes: pl.DataFrame, keys: List[Tuple[Optional[str], str]]
) -> pl.DataFrame:
    """The genotype rows of the given ``(family_id, sample)`` pairs.

    Keyed by both: a sample listed under two families has a row in each.
    """
    if not keys:
        return genotypes.clear()
    wanted = pl.DataFrame(
        {"family_id": [k[0] for k in keys], "sample": [k[1] for k in keys]},
        schema={"family_id": pl.String, "sample": pl.String},
    ).unique()
    return genotypes.join(wanted, on=["family_id", "sample"], how="semi")


def cohort_locus_rows(
    context: CohortContext, keys: List[Tuple[Optional[str], str]]
) -> List[Dict[str, Any]]:
    """One row per locus: how many of the selected individuals reach each status.

    The founders' statistics come from ``loci.tsv`` and describe the whole
    cohort; the counts follow the selection. Most flagged first.
    """
    if context.genotypes is None:
        return []
    sub = restrict_to(context.genotypes, keys)
    counts: Dict[str, Dict[str, int]] = {}
    for r in sub.group_by(["locus", "status"]).len().to_dicts():
        counts.setdefault(r["locus"], {})[r["status"]] = r["len"]
    loci = (
        context.loci.to_dicts()
        if context.loci is not None
        else [
            r
            for r in context.genotypes.unique(
                "locus", keep="first", maintain_order=True
            )
            .select(["locus", "disease", "inheritance"])
            .to_dicts()
        ]
    )
    rows = []
    for position, loc in enumerate(loci):
        c = counts.get(loc["locus"], {})
        flagged = sum(c.get(s, 0) for s in FLAGGED)
        share = common_in_cohort(loc)
        rows.append(
            {
                "id": loc["locus"],
                "order": position,
                "locus": loc["locus"],
                "disease": loc.get("disease") or "",
                "inheritance": loc.get("inheritance") or "",
                "n_pathogenic": c.get("pathogenic", 0),
                "n_carrier": c.get("carrier", 0),
                "n_intermediate": c.get("intermediate", 0),
                "n_flagged": flagged,
                "founders_p999": loc.get("founders_p99.9"),
                "founders_max": loc.get("founders_max"),
                "thresholds": thresholds_text(loc),
                "flags": "",
                "flags_badges": (
                    [
                        {
                            "label": "common",
                            "color": "#455a64",
                            "tooltip": f"{share:.0%} of the whole cohort is flagged"
                            " here: suspect the locus, catalog or thresholds",
                        }
                    ]
                    if share is not None
                    else []
                ),
            }
        )
    rows.sort(key=lambda r: (-r["n_flagged"], r["order"]))
    return rows


def cohort_flagged_rows(
    context: CohortContext, keys: List[Tuple[Optional[str], str]]
) -> List[Dict[str, Any]]:
    """The selected individuals' pathogenic, carrier and intermediate genotypes.

    Taken from the genotype table rather than ``flags.tsv``: it is already in
    memory and carries the same rows, so a third file is not read.
    """
    if context.genotypes is None:
        return []
    sub = restrict_to(context.genotypes, keys).filter(pl.col("status").is_in(FLAGGED))
    rows = []
    for r in sub.to_dicts():
        info = status_info(r.get("status"))
        rows.append(
            {
                "id": f"{r['family_id']}|{r['sample']}|{r['locus']}",
                "sample": r["sample"],
                "FID": r["family_id"],
                "locus": r["locus"],
                "status": info["label"],
                "status_rank": info["rank"],
                "status_badges": [
                    {"label": info["label"], "color": info["color"], "tooltip": ""}
                ],
                "alleles": "",
                "alleles_badges": allele_badges(r, r.get("cohort_percentile")),
                "allele_values": [
                    a
                    for a in (_int(r.get("allele1")), _int(r.get("allele2")))
                    if a is not None
                ],
                "longest": longest_allele(r),
                "pct": r.get("cohort_percentile"),
            }
        )
    rows.sort(key=lambda r: (-r["status_rank"], r["locus"], -(r["pct"] or 0)))
    return rows
