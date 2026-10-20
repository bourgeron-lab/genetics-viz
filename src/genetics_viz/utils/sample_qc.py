"""Sample QC from the ghfc-ngs ``sampleqc`` step: discovery, reading, shaping.

The step runs somalier on every CRAM and checks the result against the
pedigree. It writes, under ``qc/`` of each level, files named by the sites
label (``sampleqc_label``, e.g. ``somalier_hg38``):

- family: ``<FID>.<label>.{samples,pairs,flags,inputs}.tsv``;
- cohort: ``<NAME>.<label>.{samples,pairs,flags,inputs}.tsv`` and ``.html``
  (somalier's own relate report).

``samples`` has one row per member with sex, karyotype and quality metrics;
``pairs`` one row per pair of members (KING kinship, IBS0/IBS2, expected and
observed relationship, status); ``flags`` one row per finding. The cohort
tables are the family ones concatenated plus what only the cohort sees: the
cross-family pairs observed related, and their flags (``level`` ``cohort``,
``family_id`` ``"FA,FB"``). A family's own files therefore never hold its
members' relatives in other families.

Every table writes ``.`` for a missing value. The pipeline assigns every
flag, severity and relationship class; the app only words and draws them
(``config/sample_qc.yaml``).
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import (
    Any,
    Collection,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Set,
    Tuple,
)

import polars as pl
import yaml

from .ancestry import NULL_VALUES
from .roh import to_float
from .sample_counts import group_stats
from .sharding import get_family_path
from .tsv import find_labels, read_tsv_or_none

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "sample_qc.yaml"

QC_SUBDIR = "qc"

#: Published column lists (``sampleqc_classify``). A stub run writes cut-down
#: headers, so missing columns are added as nulls rather than trusted.
SAMPLE_COLUMNS = [
    "family_id", "sample", "father", "mother", "pedigree_sex", "inferred_sex",
    "karyotype", "sex_evidence", "x_copies", "y_copies", "somalier_x_copies",
    "somalier_y_copies", "x_het_rate", "x_sites", "depth_mean", "n_sites_called",
    "het_hom_alt_ratio", "contamination_charr", "charr_sites", "somalier_charr",
    "p_middling_ab", "trio_informative_sites", "trio_mendelian_error_rate", "flags",
]  # fmt: skip
PAIR_COLUMNS = [
    "family_id_a", "sample_a", "family_id_b", "sample_b", "n_sites", "ibs0",
    "ibs2", "ibs0_ibs2_ratio", "king_kinship", "somalier_relatedness",
    "hom_concordance", "x_ibs0", "x_ibs2", "expected_kinship",
    "expected_relationship", "observed_relationship", "po_informative_sites",
    "po_opposite_hom_rate", "status",
]  # fmt: skip
FLAG_COLUMNS = [
    "level", "family_id", "sample_a", "sample_b", "flag", "severity", "evidence",
]  # fmt: skip

#: Identifiers and labels, read as strings. ``family_id`` is numeric-looking
#: (``88171200005``) and a cohort flag's is ``"FA,FB"``.
_STRING_COLUMNS = {
    "family_id", "family_id_a", "family_id_b", "sample", "sample_a", "sample_b",
    "father", "mother", "pedigree_sex", "inferred_sex", "karyotype",
    "sex_evidence", "flags", "expected_relationship", "observed_relationship",
    "status", "level", "flag", "severity", "evidence",
}  # fmt: skip

_MISSING_PARENT = {"", "0", "-9", "."}


# ---------------------------------------------------------------------------
# Display configuration
# ---------------------------------------------------------------------------


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


_CONFIG: Dict[str, Any] = _load_config()


def _apply_config(config: Dict[str, Any]) -> None:
    global CATEGORIES, SEVERITIES, FLAGS, KARYOTYPE_FLAG, UNKNOWN_FLAG
    global DENOVO_RISK, RELATIONSHIPS, PAIR_STATUS, KARYOTYPE_COLORS
    global THRESHOLDS, KARYOTYPE_POINTS, ACKNOWLEDGEMENT
    CATEGORIES = config.get("categories", {})
    SEVERITIES = config.get("severities", {})
    FLAGS = config.get("flags", {})
    KARYOTYPE_FLAG = config.get("karyotype_flag", {})
    UNKNOWN_FLAG = config.get("unknown_flag", {})
    DENOVO_RISK = frozenset(config.get("denovo_risk", []))
    RELATIONSHIPS = config.get("relationships", {})
    PAIR_STATUS = config.get("pair_status", {})
    KARYOTYPE_COLORS = config.get("karyotypes", {})
    THRESHOLDS = {k: float(v) for k, v in config.get("thresholds", {}).items()}
    KARYOTYPE_POINTS = {
        k: (float(v[0]), float(v[1]))
        for k, v in config.get("karyotype_points", {}).items()
    }
    ACKNOWLEDGEMENT = config.get("acknowledgement", {})


CATEGORIES: Dict[str, Dict[str, Any]] = {}
SEVERITIES: Dict[str, Dict[str, Any]] = {}
FLAGS: Dict[str, Dict[str, Any]] = {}
KARYOTYPE_FLAG: Dict[str, Any] = {}
UNKNOWN_FLAG: Dict[str, Any] = {}
DENOVO_RISK: frozenset = frozenset()
RELATIONSHIPS: Dict[str, Dict[str, str]] = {}
PAIR_STATUS: Dict[str, Dict[str, str]] = {}
KARYOTYPE_COLORS: Dict[str, str] = {}
THRESHOLDS: Dict[str, float] = {}
KARYOTYPE_POINTS: Dict[str, Tuple[float, float]] = {}
ACKNOWLEDGEMENT: Dict[str, Any] = {}
_apply_config(_CONFIG)


def reload_sample_qc_config() -> None:
    """Reload the sample QC display configuration from YAML."""
    _apply_config(_load_config())


@dataclass(frozen=True)
class FlagInfo:
    """How one flag code is worded and coloured."""

    code: str
    category: str
    label: str
    reader_label: str
    explain: str

    @property
    def category_label(self) -> str:
        return CATEGORIES.get(self.category, {}).get("label", self.category)

    @property
    def category_color(self) -> str:
        return CATEGORIES.get(self.category, {}).get("color", "#9e9e9e")


def flag_info(code: str) -> FlagInfo:
    """Wording of a flag code; ``KARYOTYPE_<k>`` by prefix, unknown codes raw."""
    if code in FLAGS:
        cfg = FLAGS[code]
        return FlagInfo(
            code=code,
            category=cfg.get("category", "sample_quality"),
            label=cfg.get("label", code),
            reader_label=cfg.get("reader_label", cfg.get("label", code)),
            explain=cfg.get("explain", ""),
        )
    prefix = KARYOTYPE_FLAG.get("prefix", "KARYOTYPE_")
    if code.startswith(prefix):
        karyotype = code[len(prefix) :]
        return FlagInfo(
            code=code,
            category=KARYOTYPE_FLAG.get("category", "karyotype"),
            label=KARYOTYPE_FLAG.get("label", "{karyotype}").format(
                karyotype=karyotype
            ),
            reader_label=KARYOTYPE_FLAG.get("reader_label", code),
            explain=KARYOTYPE_FLAG.get("explain", "").format(karyotype=karyotype),
        )
    return FlagInfo(
        code=code,
        category=UNKNOWN_FLAG.get("category", "sample_quality"),
        label=code,
        reader_label=UNKNOWN_FLAG.get("reader_label", code),
        explain="",
    )


def severity_info(severity: Optional[str]) -> Dict[str, Any]:
    """Label, colour and rank of a severity; anything unknown ranks as info."""
    key = severity if severity in SEVERITIES else "info"
    cfg = SEVERITIES.get(key, {})
    return {
        "key": key,
        "label": cfg.get("label", key),
        "color": cfg.get("color", "#757575"),
        "rank": int(cfg.get("rank", 0)),
    }


def relationship_info(name: Optional[str]) -> Dict[str, str]:
    """Label and colour of a relationship class (``None`` is not measured)."""
    if not name:
        return {"label": "—", "color": "#bdbdbd"}
    cfg = RELATIONSHIPS.get(name, {})
    return {"label": cfg.get("label", name), "color": cfg.get("color", "#9e9e9e")}


def status_info(status: Optional[str]) -> Dict[str, str]:
    cfg = PAIR_STATUS.get(status or "", {})
    return {
        "label": cfg.get("label", status or "—"),
        "color": cfg.get("color", "#9e9e9e"),
    }


def karyotype_color(karyotype: Optional[str]) -> str:
    if karyotype in ("XX", "XY"):
        return KARYOTYPE_COLORS.get(karyotype, "#9e9e9e")
    if karyotype in (None, "", "unknown"):
        return "#9e9e9e"
    return KARYOTYPE_COLORS.get("other", "#6a1b9a")


def category_order() -> List[str]:
    return list(CATEGORIES)


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def get_family_qc_dir(data_dir: Path, family_id: str) -> Path:
    return get_family_path(data_dir, family_id) / QC_SUBDIR


def get_cohort_qc_dir(cohort: Any) -> Path:
    return cohort.path / QC_SUBDIR


def qc_path(directory: Path, prefix: str, label: str, kind: str) -> Path:
    """``<directory>/<prefix>.<label>.<kind>.tsv``."""
    return directory / f"{prefix}.{label}.{kind}.tsv"


def find_qc_labels(directory: Path, prefix: str) -> List[str]:
    """Sites labels with a samples table in ``directory``, newest first."""
    return find_labels(directory, prefix, "samples")


def probe_qc_data(data_dir: Path, family_id: str) -> bool:
    """Whether the family has sample QC results."""
    return bool(find_qc_labels(get_family_qc_dir(data_dir, family_id), family_id))


def probe_cohort_qc_data(cohort: Any) -> bool:
    """Whether the cohort has sample QC results."""
    return bool(find_qc_labels(get_cohort_qc_dir(cohort), cohort.name))


def cohort_report_path(cohort: Any, label: str) -> Optional[Path]:
    """somalier's HTML relate report of the cohort, when it was written."""
    path = get_cohort_qc_dir(cohort) / f"{cohort.name}.{label}.html"
    return path if path.is_file() else None


def cohort_report_url(cohort: Any, label: str, static_prefix: str) -> Optional[str]:
    """URL of :func:`cohort_report_path` under the data directory's static prefix."""
    if cohort_report_path(cohort, label) is None:
        return None
    return (
        f"{static_prefix}/cohorts/{cohort.name}/{QC_SUBDIR}/{cohort.name}.{label}.html"
    )


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_qc_table(
    path: Optional[Path], columns: Optional[List[str]] = None
) -> Optional[pl.DataFrame]:
    """Read one sample QC TSV, ``.`` as null and identifiers as strings.

    Columns of ``columns`` the file lacks (a stub run's cut-down header) are
    added as nulls. ``None`` when the path is ``None``, missing, or empty.
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
    missing = [c for c in columns or [] if c not in df.columns]
    if missing:
        df = df.with_columns([pl.lit(None).alias(c) for c in missing])
    return df


#: ``(path, mtime, size)`` -> frame; a rewrite by the next run is picked up.
_table_cache: Dict[Tuple[str, float, int], Optional[pl.DataFrame]] = {}


def _read_cached(path: Path, columns: List[str]) -> Optional[pl.DataFrame]:
    try:
        stat = path.stat()
    except OSError:
        return None
    key = (str(path), stat.st_mtime, stat.st_size)
    if key not in _table_cache:
        for old in [k for k in _table_cache if k[0] == key[0]]:
            del _table_cache[old]
        try:
            _table_cache[key] = read_qc_table(path, columns)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            _table_cache[key] = None
    return _table_cache[key]


def clear_qc_cache() -> None:
    """Drop every cached table (tests, and after a data directory switch)."""
    _table_cache.clear()


@dataclass(frozen=True)
class FlagRow:
    """One row of a ``flags.tsv``."""

    level: str
    #: The ``family_id`` cell: one FID, or ``"FA,FB"`` for a cross-family pair.
    family_id: str
    sample_a: str
    sample_b: Optional[str]
    code: str
    severity: str
    evidence: str

    @property
    def families(self) -> Tuple[str, ...]:
        return tuple(f for f in self.family_id.split(",") if f)

    @property
    def keys(self) -> List[Tuple[str, str]]:
        """``(family_id, sample)`` of every sample the flag is about.

        A cross-family pair's ``sample_a`` is in the first family and
        ``sample_b`` in the second; otherwise both are in the one family.
        """
        fams = self.families or ("",)
        out = [(fams[0], self.sample_a)]
        if self.sample_b:
            out.append((fams[1] if len(fams) > 1 else fams[0], self.sample_b))
        return out

    @property
    def samples(self) -> List[str]:
        return [s for _, s in self.keys]

    @property
    def is_cross_family(self) -> bool:
        return len(set(self.families)) > 1

    @property
    def ack_key(self) -> Tuple[str, str, str, str]:
        """``(code, family_id, sample_a, sample_b)``: what a curator acknowledges.

        The family files and the cohort file write a family flag alike, so one
        acknowledgement covers it on both pages.
        """
        return (self.code, self.family_id, self.sample_a, self.sample_b or "")

    @property
    def info(self) -> FlagInfo:
        return flag_info(self.code)

    @property
    def rank(self) -> int:
        return severity_info(self.severity)["rank"]

    def other_family(self, family_id: str) -> Optional[str]:
        """The family on the other side of a cross-family flag."""
        others = [f for f in self.families if f != family_id]
        return others[0] if others and self.is_cross_family else None


def flag_rows(df: Optional[pl.DataFrame]) -> List[FlagRow]:
    """Every flag of a flags table, most severe first, in file order otherwise."""
    if df is None or len(df) == 0:
        return []
    rows = []
    for r in df.to_dicts():
        if not r.get("flag") or not r.get("sample_a"):
            continue
        rows.append(
            FlagRow(
                level=r.get("level") or "family",
                family_id=r.get("family_id") or "",
                sample_a=r["sample_a"],
                sample_b=r.get("sample_b") or None,
                code=r["flag"],
                severity=r.get("severity") or "info",
                evidence=r.get("evidence") or "",
            )
        )
    return sorted(rows, key=lambda f: -f.rank)


@dataclass
class QcData:
    """The sample QC tables of one family or cohort, for one sites label."""

    prefix: str
    directory: Path
    label: str
    labels: List[str]
    samples: Optional[pl.DataFrame]
    pairs: Optional[pl.DataFrame]
    flags: List[FlagRow] = field(default_factory=list)

    def sample_rows(self) -> List[Dict[str, Any]]:
        return [] if self.samples is None else self.samples.to_dicts()

    def pair_rows(self) -> List[Dict[str, Any]]:
        return [] if self.pairs is None else self.pairs.to_dicts()


def load_qc(
    directory: Path, prefix: str, label: Optional[str] = None
) -> Optional[QcData]:
    """The QC tables in ``directory`` for ``label`` (default: the newest).

    Blocking I/O over the data mount; the tables are cached on their stat.
    ``None`` when there is no samples table at all.
    """
    labels = find_qc_labels(directory, prefix)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    samples = _read_cached(
        qc_path(directory, prefix, chosen, "samples"), SAMPLE_COLUMNS
    )
    pairs = _read_cached(qc_path(directory, prefix, chosen, "pairs"), PAIR_COLUMNS)
    flags = _read_cached(qc_path(directory, prefix, chosen, "flags"), FLAG_COLUMNS)
    return QcData(
        prefix=prefix,
        directory=directory,
        label=chosen,
        labels=labels,
        samples=samples,
        pairs=pairs,
        flags=flag_rows(flags),
    )


def load_family_qc(
    data_dir: Path, family_id: str, label: Optional[str] = None
) -> Optional[QcData]:
    return load_qc(get_family_qc_dir(data_dir, family_id), family_id, label)


def load_cohort_qc(cohort: Any, label: Optional[str] = None) -> Optional[QcData]:
    return load_qc(get_cohort_qc_dir(cohort), cohort.name, label)


# ---------------------------------------------------------------------------
# Shaping
# ---------------------------------------------------------------------------


def flags_by_key(flags: Iterable[FlagRow]) -> Dict[Tuple[str, str], List[FlagRow]]:
    """``(family_id, sample) -> flags`` about that sample, a pair counting twice."""
    out: Dict[Tuple[str, str], List[FlagRow]] = {}
    for f in flags:
        for key in f.keys:
            out.setdefault(key, []).append(f)
    return out


def flags_by_sample(flags: Iterable[FlagRow]) -> Dict[str, List[FlagRow]]:
    """``sample -> flags``, for one family's tables."""
    out: Dict[str, List[FlagRow]] = {}
    for f in flags:
        for sample in dict.fromkeys(f.samples):
            out.setdefault(sample, []).append(f)
    return out


def family_flags_with_cohort(
    family_qc: Optional[QcData],
    cohort_qc: Optional[QcData],
    family_id: str,
) -> List[FlagRow]:
    """The family's own flags, then the cohort-level flags about its members.

    The cohort file repeats the family flags; only its ``cohort`` level ones
    are added, so nothing is shown twice.
    """
    own = list(family_qc.flags) if family_qc else []
    if cohort_qc is None:
        return own
    extra = [
        f
        for f in cohort_qc.flags
        if f.level == "cohort" and any(fid == family_id for fid, _ in f.keys)
    ]
    return sorted(own + extra, key=lambda f: -f.rank)


def count_by_severity(
    flags: Iterable[FlagRow], restrict_to: Optional[Set[Tuple[str, str]]] = None
) -> Dict[str, int]:
    """How many flags of each severity, optionally only those about ``restrict_to``."""
    counts = {key: 0 for key in SEVERITIES}
    for f in flags:
        if restrict_to is not None and not any(k in restrict_to for k in f.keys):
            continue
        key = severity_info(f.severity)["key"]
        counts[key] = counts.get(key, 0) + 1
    return counts


def severity_summary(counts: Dict[str, int]) -> str:
    """``"2 high · 1 medium"``; empty when nothing above info."""
    parts = [f"{counts[k]} {k}" for k in ("high", "medium") if counts.get(k)]
    return " · ".join(parts)


def flag_title(f: FlagRow, show_detail: bool) -> str:
    """The flag's label for the user's role."""
    return f.info.label if show_detail else f.info.reader_label


def flag_tooltip(f: FlagRow, show_detail: bool) -> str:
    """Multi-line hover text: what it means, then the pipeline's evidence."""
    if not show_detail:
        return (
            f"{f.info.reader_label} ({severity_info(f.severity)['label'].lower()})."
            "\nA curator can see the details."
        )
    lines = [f"{f.info.label} — {severity_info(f.severity)['label'].lower()}"]
    if f.info.explain:
        lines.append(f.info.explain)
    if f.evidence:
        lines.append(f"Evidence: {f.evidence}")
    return "\n".join(lines)


def flag_badge(f: FlagRow, show_detail: bool) -> Dict[str, str]:
    """``badge_list`` entry of one flag, coloured by severity."""
    return {
        "label": flag_title(f, show_detail),
        "color": severity_info(f.severity)["color"],
        "tooltip": flag_tooltip(f, show_detail),
    }


def flag_badges(flags: Iterable[FlagRow], show_detail: bool) -> List[Dict[str, str]]:
    """One badge per distinct flag label (two pair flags read as one)."""
    seen: Dict[str, Dict[str, str]] = {}
    for f in sorted(flags, key=lambda f: -f.rank):
        title = flag_title(f, show_detail)
        if title in seen:
            if show_detail and f.evidence:
                seen[title]["tooltip"] += f"\nEvidence: {f.evidence}"
            continue
        seen[title] = flag_badge(f, show_detail)
    return list(seen.values())


def category_badges(
    flags: Iterable[FlagRow], show_detail: bool
) -> List[Dict[str, str]]:
    """One badge per category, coloured by its worst severity, in config order.

    The category names are neutral, so the badges read the same for every
    role; only the tooltip carries the flags, worded for the role.
    """
    by_cat: Dict[str, List[FlagRow]] = {}
    for f in flags:
        by_cat.setdefault(f.info.category, []).append(f)
    order = category_order()
    badges = []
    for cat in sorted(by_cat, key=lambda c: order.index(c) if c in order else 99):
        group = sorted(by_cat[cat], key=lambda f: -f.rank)
        label = CATEGORIES.get(cat, {}).get("label", cat)
        lines = dict.fromkeys(
            f"• {flag_title(f, show_detail)}"
            + (f": {f.evidence}" if show_detail and f.evidence else "")
            for f in group
        )
        badges.append(
            {
                "label": label,
                "color": severity_info(group[0].severity)["color"],
                # Always several lines, so the floating tooltip shows it whole.
                "tooltip": "\n".join([label, *lines]),
            }
        )
    return badges


def denovo_risk(flags: Iterable[FlagRow], samples: Iterable[str]) -> List[FlagRow]:
    """The flags about ``samples`` that make their de novo calls unreliable."""
    wanted = set(samples)
    return [
        f
        for f in flags
        if f.code in DENOVO_RISK and any(s in wanted for s in f.samples)
    ]


def is_founder(row: Dict[str, Any]) -> bool:
    return (row.get("father") or "0") in _MISSING_PARENT and (
        row.get("mother") or "0"
    ) in _MISSING_PARENT


def mendelian_fence(rows: Iterable[Dict[str, Any]]) -> Optional[float]:
    """Upper outlier fence (Q3 + 1.5 IQR) of the trio Mendelian error rates.

    Drawn only: the pipeline raises no flag on the rate itself.
    """
    rates = [
        r
        for r in (to_float(row.get("trio_mendelian_error_rate")) for row in rows)
        if r is not None
    ]
    if len(rates) < 4:
        return None
    return group_stats("mendelian", rates).high_fence


@dataclass(frozen=True)
class PairRow:
    """One ``pairs.tsv`` row, typed."""

    family_a: str
    sample_a: str
    family_b: str
    sample_b: str
    kinship: Optional[float]
    ratio: Optional[float]
    n_sites: Optional[float]
    expected_kinship: Optional[float]
    expected: Optional[str]
    observed: Optional[str]
    status: Optional[str]
    po_rate: Optional[float]

    @property
    def cross_family(self) -> bool:
        return self.family_a != self.family_b

    @property
    def key(self) -> str:
        return f"{self.family_a}|{self.sample_a}|{self.family_b}|{self.sample_b}"

    @property
    def declared(self) -> bool:
        """The pedigree declares the two related (or they are the same family)."""
        return (self.expected or "unrelated") != "unrelated"

    @property
    def flag_key(self) -> FrozenSet[Tuple[str, str]]:
        """The pair as :func:`pair_severities` keys it."""
        return frozenset(
            {(self.family_a, self.sample_a), (self.family_b, self.sample_b)}
        )


def pair_severities(flags: Iterable[FlagRow]) -> Dict[FrozenSet[Tuple[str, str]], str]:
    """The worst severity of the flags naming each pair, keyed as
    :attr:`PairRow.flag_key`.

    The pipeline marks a pair ``mismatch`` whatever its flag's severity, so a
    pair of MZ twins, of consanguineous parents or of third-degree relatives
    across families is one; the severity tells them from a contradiction.
    """
    out: Dict[FrozenSet[Tuple[str, str]], FlagRow] = {}
    for f in flags:
        if not f.sample_b:
            continue
        key = frozenset(f.keys)
        if key not in out or f.rank > out[key].rank:
            out[key] = f
    return {k: severity_info(f.severity)["key"] for k, f in out.items()}


def pair_rows(df: Optional[pl.DataFrame]) -> List[PairRow]:
    if df is None or len(df) == 0:
        return []
    out = []
    for r in df.to_dicts():
        if not r.get("sample_a") or not r.get("sample_b"):
            continue
        out.append(
            PairRow(
                family_a=r.get("family_id_a") or "",
                sample_a=r["sample_a"],
                family_b=r.get("family_id_b") or r.get("family_id_a") or "",
                sample_b=r["sample_b"],
                kinship=to_float(r.get("king_kinship")),
                ratio=to_float(r.get("ibs0_ibs2_ratio")),
                n_sites=to_float(r.get("n_sites")),
                expected_kinship=to_float(r.get("expected_kinship")),
                expected=r.get("expected_relationship"),
                observed=r.get("observed_relationship"),
                status=r.get("status"),
                po_rate=to_float(r.get("po_opposite_hom_rate")),
            )
        )
    return out


#: A sample with this many relatives in other families is called a hub. Fewer
#: is often one real family split over several family IDs.
HUB_MINIMUM = 5


@dataclass
class GraphNode:
    family_id: str
    sample: str
    #: Cross-family relatives of this sample.
    degree: int = 0


@dataclass
class CrossFamilyGraph:
    nodes: Dict[Tuple[str, str], GraphNode]
    edges: List[PairRow]

    def hubs(self, minimum: int = HUB_MINIMUM) -> List[GraphNode]:
        """Samples related to at least ``minimum`` others outside their family."""
        return sorted(
            (n for n in self.nodes.values() if n.degree >= minimum),
            key=lambda n: -n.degree,
        )


def cross_family_graph(
    pairs: Iterable[PairRow],
    restrict_to: Optional[Set[Tuple[str, str]]] = None,
) -> CrossFamilyGraph:
    """The cross-family related pairs as a graph of samples.

    With ``restrict_to``, an edge is kept when either end is in it: a shown
    sample's relative in a filtered-out family is still its relative.
    """
    nodes: Dict[Tuple[str, str], GraphNode] = {}
    edges: List[PairRow] = []
    for p in pairs:
        if not p.cross_family:
            continue
        a, b = (p.family_a, p.sample_a), (p.family_b, p.sample_b)
        if restrict_to is not None and a not in restrict_to and b not in restrict_to:
            continue
        edges.append(p)
        for fid, sample in (a, b):
            node = nodes.setdefault((fid, sample), GraphNode(fid, sample))
            node.degree += 1
    return CrossFamilyGraph(nodes=nodes, edges=edges)


def worst_by_individual(flags: Iterable[FlagRow]) -> Dict[str, int]:
    """How many individuals have each severity as their worst flag."""
    counts = {key: 0 for key in SEVERITIES}
    for key_flags in flags_by_key(flags).values():
        worst = max(key_flags, key=lambda f: f.rank)
        sev = severity_info(worst.severity)["key"]
        counts[sev] = counts.get(sev, 0) + 1
    return counts


def cohort_headline(
    cohort: Any, acknowledged: Collection[Tuple[str, str, str, str]] = ()
) -> Optional[Dict[str, int]]:
    """Individuals by worst severity, for the home card; ``None`` when QC has not run.

    Flags whose :attr:`FlagRow.ack_key` is in ``acknowledged`` do not count.
    Blocking I/O - call it through ``asyncio.to_thread``.
    """
    try:
        qc = load_cohort_qc(cohort)
    except Exception:
        logger.warning("Failed to read sample QC of %s", cohort.name, exc_info=True)
        return None
    if qc is None:
        return None
    return worst_by_individual(f for f in qc.flags if f.ack_key not in acknowledged)


#: Categories a sex / karyotype plot marks as flagged.
SEX_CATEGORIES = frozenset({"sex", "karyotype"})
QUALITY_CATEGORIES = frozenset({"sample_quality", "mendelian"})


def sample_point(row: Dict[str, Any]) -> Dict[str, Any]:
    """One samples-table row as a point for the sex and quality plots."""
    return {
        "name": row.get("sample") or "",
        "family": row.get("family_id") or "",
        "x": to_float(row.get("x_copies")),
        "y": to_float(row.get("y_copies")),
        "het": to_float(row.get("x_het_rate")),
        "pedigree_sex": row.get("pedigree_sex"),
        "inferred_sex": row.get("inferred_sex"),
        "karyotype": row.get("karyotype"),
        "evidence": row.get("sex_evidence"),
        "charr": to_float(row.get("contamination_charr")),
        "mendel": to_float(row.get("trio_mendelian_error_rate")),
    }


def has_category(flags: Iterable[FlagRow], categories: Iterable[str]) -> bool:
    wanted = set(categories)
    return any(f.info.category in wanted for f in flags)


def is_typical_karyotype(karyotype: Optional[str]) -> bool:
    return karyotype in ("XX", "XY")
