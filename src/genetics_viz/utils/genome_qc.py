"""Genome QC from the ghfc-ngs ``genomeqc`` step: discovery, reading, wording.

The step says what kind of DNA and library each sample is, which ``sampleqc``
cannot: a genome from an EBV-transformed cell line has the person's genotypes,
so every relatedness, sex and Mendelian check passes, while its clonal somatic
mutations come out as extra de novo calls. It writes, under the label
``genomeqc_label`` (``gqc_v1``):

- sample: ``samples/…/<ID>/qc/genomeqc/<ID>.<label>.tsv`` (one row: EBV and
  mtDNA copies per cell, read length, chr20 depth and alignment metrics) and
  ``.stats.txt`` (samtools' own output, for MultiQC);
- family: ``families/…/<FID>/qc/<FID>.<label>.{variants,inputs}.tsv`` (one row
  per member: SNV and indel profile from the family BCF). **No flags**;
- cohort: ``cohorts/<NAME>/qc/<NAME>.<label>.{samples,flags,inputs}.tsv``:
  every sample's metrics and its member row, its ``family_id``, ``library``,
  ``dna_source``, ``seq_class`` and ``flags``; one row per flag.

Flags exist at cohort level only: outliers are judged against the samples of
the same ``seq_class`` (read length class × PCR), 20 of them at least. The
step writes one row per **sample**, with one ``family_id`` even for a sample
in two families.

The family and cohort files sit in the same ``qc/`` directories as
``sampleqc``'s and share its file names, so a label is the step's only when
its ``inputs.tsv`` says ``#step genomeqc``.

Every table writes ``.`` for a missing value. The pipeline measures, classes
and flags; the app words and draws (``config/genome_qc.yaml`` for metrics and
classes, ``config/sample_qc.yaml`` for the flags, which are shown, counted and
acknowledged as Sample QC flags).
"""

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set

import polars as pl
import yaml

from .ancestry import NULL_VALUES
from .sample_qc import FlagRow
from .sharding import get_family_path, get_sample_path
from .tsv import find_labels, inputs_step, read_tsv_or_none

logger = logging.getLogger(__name__)

_CONFIG_PATH = Path(__file__).parent.parent / "config" / "genome_qc.yaml"

GQC_STEP = "genomeqc"
QC_SUBDIR = "qc"
SAMPLE_SUBDIR = ("qc", "genomeqc")

#: ``FlagRow.level`` of a genome QC flag (Sample QC's are family or cohort).
GENOME_LEVEL = "genome"

#: The pipeline's ``low`` severity (``EBV_TRACE``) is Sample QC's ``info``.
SEVERITY_ALIASES = {"low": "info"}

#: Published column lists (``genomeqc_tables``). Missing columns (an older
#: label, a stub run) are added as nulls.
SAMPLE_COLUMNS = [
    "sample", "read_length", "mean_read_length", "chr20_reads", "chr20_depth",
    "ebv_reads", "ebv_copies_per_cell", "chrm_reads", "mt_copies_per_cell",
    "duplicate_share", "mapq0_share", "properly_paired_share", "error_rate",
    "mean_base_quality", "insert_median", "insert_mad", "read_indels_per_kb",
]  # fmt: skip
FAMILY_COLUMNS = [
    "family_id", "sample", "n_snv", "n_snv_het", "n_snv_hom", "het_hom_ratio",
    "ti_tv", "private_snvs", "low_vaf_het_share", "n_indel", "ins_del_ratio",
    "indel_1bp_share", "ins_over_10bp_share", "indel_tr_share",
]  # fmt: skip
COHORT_COLUMNS = (
    ["sample", "family_id", "seq_class", "library", "dna_source"]
    + SAMPLE_COLUMNS[1:]
    + FAMILY_COLUMNS[2:]
    + ["flags"]
)
FLAG_COLUMNS = ["sample", "family_id", "seq_class", "flag", "severity", "evidence"]

#: Identifiers and labels, read as strings: ``family_id`` is numeric-looking
#: (``88171200005``).
_STRING_COLUMNS = {
    "sample", "family_id", "seq_class", "library", "dna_source", "flags", "flag",
    "severity", "evidence",
}  # fmt: skip

CELL_LINE = "CELL_LINE"


# ---------------------------------------------------------------------------
# Display configuration
# ---------------------------------------------------------------------------


def _load_config() -> Dict[str, Any]:
    with open(_CONFIG_PATH, "r") as f:
        return yaml.safe_load(f) or {}


THRESHOLDS: Dict[str, float] = {}
SEQ_CLASSES: Dict[str, Dict[str, str]] = {}
GROUPS: Dict[str, str] = {}
METRICS: Dict[str, Dict[str, Any]] = {}
CAVEATS: Dict[str, str] = {}


def _apply_config(config: Dict[str, Any]) -> None:
    global THRESHOLDS, SEQ_CLASSES, GROUPS, METRICS, CAVEATS
    THRESHOLDS = {k: float(v) for k, v in config.get("thresholds", {}).items()}
    SEQ_CLASSES = config.get("seq_classes", {})
    GROUPS = config.get("groups", {})
    METRICS = config.get("metrics", {})
    CAVEATS = {
        k: " ".join(str(v).split()) for k, v in config.get("caveats", {}).items()
    }


_apply_config(_load_config())


def reload_genome_qc_config() -> None:
    """Reload the genome QC display configuration from YAML."""
    _apply_config(_load_config())


@dataclass(frozen=True)
class MetricInfo:
    """How one column of the samples table is named, formatted and drawn."""

    key: str
    group: str
    label: str
    format: str
    digits: int
    #: Direction of the pipeline's class test: ``high``, ``low``, ``both``.
    side: Optional[str]
    log: bool
    #: Fixed thresholds (values of ``THRESHOLDS``) drawn as lines.
    lines: List[float]
    plot: bool
    explain: str
    #: The flags the pipeline raises on this column.
    flags: List[str] = field(default_factory=list)
    #: Hidden from readers (the EBV columns).
    curator_only: bool = False

    def fmt(self, value: Any) -> str:
        """The value as shown: ``—`` when missing."""
        if value is None or value == "":
            return "—"
        if self.format == "text":
            return str(value)
        try:
            number = float(value)
        except (TypeError, ValueError):
            return str(value)
        if self.format == "int":
            return f"{number:,.0f}"
        if self.format == "percent":
            return f"{number * 100:.{self.digits}f}%"
        return f"{number:,.{self.digits}f}"


def metric_info(key: str) -> MetricInfo:
    """Wording of a column; an unknown one is shown raw, as a number."""
    cfg = METRICS.get(key, {})
    return MetricInfo(
        key=key,
        group=cfg.get("group", ""),
        label=cfg.get("label", key),
        format=cfg.get("format", "float"),
        digits=int(cfg.get("digits", 2)),
        side=cfg.get("side"),
        log=bool(cfg.get("log", False)),
        lines=[THRESHOLDS[k] for k in cfg.get("lines", []) if k in THRESHOLDS],
        plot=bool(cfg.get("plot", cfg.get("format", "float") != "text")),
        explain=cfg.get("explain", ""),
        flags=list(cfg.get("flags", [])),
        curator_only=bool(cfg.get("curator_only", False)),
    )


def metrics(group: Optional[str] = None, curator: bool = True) -> List[MetricInfo]:
    """Every configured metric in file order, or one group's.

    ``curator`` false leaves out the curator-only ones (the EBV columns).
    """
    return [
        info
        for info in (
            metric_info(k)
            for k, cfg in METRICS.items()
            if group is None or cfg.get("group") == group
        )
        if curator or not info.curator_only
    ]


def seq_class_info(name: Optional[str]) -> Dict[str, str]:
    """Label, colour and explanation of a ``seq_class``; unknown ones raw."""
    if not name:
        name = "unknown"
    cfg = SEQ_CLASSES.get(name)
    if cfg is None:
        return {"key": name, "label": name, "color": "#9e9e9e", "explain": ""}
    return {
        "key": name,
        "label": cfg.get("label", name),
        "color": cfg.get("color", "#9e9e9e"),
        "explain": cfg.get("explain", ""),
    }


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def get_cohort_gqc_dir(cohort: Any) -> Path:
    return cohort.path / QC_SUBDIR


def get_family_gqc_dir(data_dir: Path, family_id: str) -> Path:
    return get_family_path(data_dir, family_id) / QC_SUBDIR


def get_sample_gqc_dir(data_dir: Path, sample_id: str) -> Path:
    return get_sample_path(data_dir, sample_id).joinpath(*SAMPLE_SUBDIR)


def gqc_path(directory: Path, prefix: str, label: str, kind: str) -> Path:
    """``<directory>/<prefix>.<label>.<kind>.tsv``."""
    return directory / f"{prefix}.{label}.{kind}.tsv"


def _step_labels(directory: Path, prefix: str, kind: str) -> List[str]:
    """Labels of ``kind`` tables written by ``genomeqc``, newest first.

    Strict: the step always writes ``inputs.tsv``, and a label without one
    may be ``sampleqc``'s (whose cohort files share the ``samples`` kind).
    """
    return [
        label
        for label in find_labels(directory, prefix, kind)
        if inputs_step(gqc_path(directory, prefix, label, "inputs")) == GQC_STEP
    ]


def find_cohort_gqc_labels(directory: Path, prefix: str) -> List[str]:
    """Labels of the cohort samples tables in ``directory``, newest first."""
    return _step_labels(directory, prefix, "samples")


def find_family_gqc_labels(directory: Path, prefix: str) -> List[str]:
    """Labels of the family variants tables in ``directory``, newest first."""
    return _step_labels(directory, prefix, "variants")


def find_sample_gqc_labels(directory: Path, sample_id: str) -> List[str]:
    """Labels of ``<ID>.<label>.tsv`` in a sample's ``qc/genomeqc/``.

    The directory is the step's own, so every ``.tsv`` there is one of its
    sample tables. Newest first.
    """
    if not directory.is_dir():
        return []
    head, tail = f"{sample_id}.", ".tsv"
    found = []
    try:
        for entry in directory.iterdir():
            name = entry.name
            if name.startswith(head) and name.endswith(tail):
                label = name[len(head) : -len(tail)]
                if label:
                    found.append((entry.stat().st_mtime, label))
    except OSError:
        logger.warning("Cannot list %s", directory, exc_info=True)
        return []
    return [label for _, label in sorted(found, key=lambda t: (-t[0], t[1]))]


def probe_cohort_gqc(cohort: Any) -> bool:
    """Whether the cohort has genome QC results."""
    return bool(find_cohort_gqc_labels(get_cohort_gqc_dir(cohort), cohort.name))


def probe_family_gqc(data_dir: Path, family_id: str, cohort: Any = None) -> bool:
    """Whether the family has genome QC results: its own table, or the cohort's."""
    if find_family_gqc_labels(get_family_gqc_dir(data_dir, family_id), family_id):
        return True
    return cohort is not None and probe_cohort_gqc(cohort)


def sample_stats_path(
    data_dir: Path, sample_id: str, label: Optional[str] = None
) -> Optional[Path]:
    """The sample's ``samtools stats`` output (curators only), if any."""
    directory = get_sample_gqc_dir(data_dir, sample_id)
    labels = find_sample_gqc_labels(directory, sample_id)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    path = directory / f"{sample_id}.{chosen}.stats.txt"
    return path if path.exists() else None


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_gqc_table(
    path: Optional[Path], columns: Optional[List[str]] = None
) -> Optional[pl.DataFrame]:
    """Read one genome QC TSV, ``.`` as null and identifiers as strings.

    Columns of ``columns`` the file lacks are added as nulls. ``None`` when
    the path is ``None``, missing, or empty.
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


_table_cache: Dict[tuple, Optional[pl.DataFrame]] = {}


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
            _table_cache[key] = read_gqc_table(path, columns)
        except Exception:
            logger.warning("Failed to read %s", path, exc_info=True)
            _table_cache[key] = None
    return _table_cache[key]


def clear_gqc_cache() -> None:
    """Drop every cached table (tests, and after a data directory switch)."""
    _table_cache.clear()


def flag_rows(df: Optional[pl.DataFrame]) -> List[FlagRow]:
    """A flags table as Sample QC ``FlagRow``s, most severe first.

    One sample each (``sample_b`` empty), level ``genome``, the pipeline's
    ``low`` severity as ``info``. Their ``ack_key`` is the Sample QC one, so
    ``validations/sampleqc.tsv`` acknowledges them unchanged.
    """
    if df is None or len(df) == 0:
        return []
    rows = []
    for r in df.to_dicts():
        if not r.get("flag") or not r.get("sample"):
            continue
        severity = r.get("severity") or "info"
        rows.append(
            FlagRow(
                level=GENOME_LEVEL,
                family_id=r.get("family_id") or "",
                sample_a=r["sample"],
                sample_b=None,
                code=r["flag"],
                severity=SEVERITY_ALIASES.get(severity, severity),
                evidence=r.get("evidence") or "",
            )
        )
    return sorted(rows, key=lambda f: -f.rank)


@dataclass
class GenomeQc:
    """The genome QC of a cohort, or of one family, for one label.

    ``source`` is ``cohort`` (the cohort tables: metrics, ``seq_class`` and
    flags) or ``family`` (no cohort table: the family's variants table and its
    members' sample tables, without class nor flags).
    """

    prefix: str
    label: str
    labels: List[str]
    source: str
    samples: Optional[pl.DataFrame]
    flags: List[FlagRow] = field(default_factory=list)

    @property
    def flagged(self) -> bool:
        """Whether flags were computed (cohort tables): no flag means clean."""
        return self.source == "cohort"

    def sample_rows(self) -> List[Dict[str, Any]]:
        return [] if self.samples is None else self.samples.to_dicts()

    def row_for(self, sample: str) -> Optional[Dict[str, Any]]:
        if self.samples is None or "sample" not in self.samples.columns:
            return None
        rows = self.samples.filter(pl.col("sample") == sample).to_dicts()
        return rows[0] if rows else None

    def flags_for(self, sample: str) -> List[FlagRow]:
        return [f for f in self.flags if f.sample_a == sample]

    def cell_lines(self) -> Set[str]:
        """Samples flagged ``CELL_LINE``."""
        return {f.sample_a for f in self.flags if f.code == CELL_LINE}


def load_cohort_gqc(cohort: Any, label: Optional[str] = None) -> Optional[GenomeQc]:
    """The cohort's genome QC for ``label`` (default: the newest).

    Blocking I/O over the data mount; tables are cached on their stat.
    ``None`` when the cohort has no genome QC samples table.
    """
    directory = get_cohort_gqc_dir(cohort)
    labels = find_cohort_gqc_labels(directory, cohort.name)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    samples = _read_cached(
        gqc_path(directory, cohort.name, chosen, "samples"), COHORT_COLUMNS
    )
    flags = _read_cached(
        gqc_path(directory, cohort.name, chosen, "flags"), FLAG_COLUMNS
    )
    return GenomeQc(
        prefix=cohort.name,
        label=chosen,
        labels=labels,
        source="cohort",
        samples=samples,
        flags=flag_rows(flags),
    )


def cohort_flags(cohort: Any, label: Optional[str] = None) -> List[FlagRow]:
    """The cohort's genome QC flags (the flags table only, a few KB).

    For the views that add them to the Sample QC flags: family header and
    member icons, review badges, Issues, the cohort column and pills, the
    home figure, the report. Blocking I/O, cached on the file's stat; ``[]``
    without a cohort or its tables, or when they cannot be read.
    """
    if cohort is None:
        return []
    try:
        directory = get_cohort_gqc_dir(cohort)
        labels = find_cohort_gqc_labels(directory, cohort.name)
        if not labels:
            return []
        chosen = label if label in labels else labels[0]
        return flag_rows(
            _read_cached(
                gqc_path(directory, cohort.name, chosen, "flags"), FLAG_COLUMNS
            )
        )
    except Exception:
        logger.warning(
            "Failed to read the genome QC flags of %s", cohort.name, exc_info=True
        )
        return []


def family_flags(cohort: Any, family_id: str) -> List[FlagRow]:
    """:func:`cohort_flags` of one family."""
    return [f for f in cohort_flags(cohort) if f.family_id == family_id]


def cell_lines(cohorts: Iterable[Any], family_id: str) -> Dict[str, FlagRow]:
    """Sample -> ``CELL_LINE`` flag of the family's members, over ``cohorts``.

    A family may sit in several cohorts; any of their genome QC tables
    counts. Acknowledged or not: the caveat stays. Blocking I/O (the flags
    tables, a few KB each, cached).
    """
    out: Dict[str, FlagRow] = {}
    for cohort in cohorts:
        for f in family_flags(cohort, family_id):
            if f.code == CELL_LINE:
                out.setdefault(f.sample_a, f)
    return out


def caveat_text(context: str, samples: Iterable[str]) -> str:
    """The caveat on cell-line ``samples``' results in ``context``; ``""`` without."""
    names = sorted(set(samples))
    template = CAVEATS.get(context, "")
    if not names or not template:
        return ""
    return template.format(samples=", ".join(names))


def _family_from_cohort(cohort_gqc: GenomeQc, family_id: str) -> GenomeQc:
    samples = cohort_gqc.samples
    if samples is not None:
        samples = samples.filter(pl.col("family_id") == family_id)
    return GenomeQc(
        prefix=family_id,
        label=cohort_gqc.label,
        labels=cohort_gqc.labels,
        source="cohort",
        samples=samples,
        flags=[f for f in cohort_gqc.flags if f.family_id == family_id],
    )


def _family_own(
    data_dir: Path, family_id: str, label: Optional[str]
) -> Optional[GenomeQc]:
    directory = get_family_gqc_dir(data_dir, family_id)
    labels = find_family_gqc_labels(directory, family_id)
    if not labels:
        return None
    chosen = label if label in labels else labels[0]
    members = _read_cached(
        gqc_path(directory, family_id, chosen, "variants"), FAMILY_COLUMNS
    )
    if members is None or "sample" not in members.columns:
        return GenomeQc(family_id, chosen, labels, "family", members)
    rows = []
    for sid in members.get_column("sample").drop_nulls().to_list():
        path = get_sample_gqc_dir(data_dir, sid) / f"{sid}.{chosen}.tsv"
        row = _read_cached(path, SAMPLE_COLUMNS)
        if row is not None and len(row):
            rows.append(row.select(SAMPLE_COLUMNS))
    if rows:
        per_sample = pl.concat(rows, how="vertical_relaxed")
        members = members.join(per_sample, on="sample", how="left")
    # The cohort table's shape: class, library, source and flags are null.
    missing = [c for c in COHORT_COLUMNS if c not in members.columns]
    members = members.with_columns(
        [
            pl.lit(None, dtype=pl.Utf8 if c in _STRING_COLUMNS else pl.Float64).alias(c)
            for c in missing
        ]
    )
    return GenomeQc(family_id, chosen, labels, "family", members)


def load_family_gqc(
    data_dir: Path,
    family_id: str,
    cohort: Any = None,
    label: Optional[str] = None,
) -> Optional[GenomeQc]:
    """One family's genome QC: from the cohort tables when there are any.

    Only the cohort tables carry ``seq_class`` and the flags. Without them
    (a standalone family, a cohort not gathered yet) the family's variants
    table is joined to its members' sample tables, with no flag. Rows are
    the cohort's whose ``family_id`` is the family: the step keeps one row,
    and one family, per sample. Blocking I/O; ``None`` without any table.
    """
    if cohort is not None:
        cohort_gqc = load_cohort_gqc(cohort, label)
        if cohort_gqc is not None:
            return _family_from_cohort(cohort_gqc, family_id)
    return _family_own(data_dir, family_id, label)
