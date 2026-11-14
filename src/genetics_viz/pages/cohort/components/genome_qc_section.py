"""The "DNA and library" section of the Sample QC tabs (ghfc-ngs ``genomeqc``).

What kind of DNA and library each sample is: EBV genomes per cell (a cell
line), read length and PCR (``seq_class``), chr20 alignment metrics and each
member's variant profile. The flags themselves are listed with the Sample QC
flags (Issues, flagged individuals); this section shows the metrics behind
them, one row per ``seq_class``, since the pipeline judges a sample against
its own class.

Readers see no EBV column or chart (EBV in the person is health information);
the ``CELL_LINE`` flag says what they need.
"""

from typing import Any, Callable, Dict, Iterable, List, Optional, Set, Tuple

from nicegui import ui

from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.genome_qc import (
    CELL_LINE,
    GenomeQc,
    MetricInfo,
    metric_info,
    metrics,
    seq_class_info,
)
from genetics_viz.utils.genome_qc_plots import chart_height, class_strip_option
from genetics_viz.utils.roh import to_float
from genetics_viz.utils.sample_qc import FlagRow, flag_badges, flag_info
from genetics_viz.utils.sample_qc_acks import (
    Acks,
    acknowledged_badges,
    split_acknowledged,
)

Key = Tuple[str, str]

#: The table's metric columns, in order (the others are in the download).
TABLE_METRICS = [
    "ebv_copies_per_cell", "mt_copies_per_cell", "chr20_depth", "duplicate_share",
    "error_rate", "insert_median", "read_indels_per_kb", "het_hom_ratio", "ti_tv",
    "low_vaf_het_share", "n_indel",
]  # fmt: skip

_NUM = {"cellType": "number", "sorting": "numerical", "sortable": True}

CELL_LINE_TEXT = (
    "Its inherited genotypes are sound. Its de novo calls, UPD isodisomies and"
    " CNVs may come from the culture and need confirmation on primary DNA"
    " (blood or saliva)."
)


def _cell(metric: MetricInfo, value: Any) -> Optional[float]:
    """A number for the table: percentages in %, counts whole, others as read.

    Not rounded to the metric's digits: an EBV count of 0.0005 per cell
    would read 0.
    """
    number = to_float(value)
    if number is None:
        return None
    if metric.format == "percent":
        return round(number * 100, 4)
    if metric.format == "int":
        return round(number)
    return number


def _header(metric: MetricInfo) -> str:
    return f"{metric.label} (%)" if metric.format == "percent" else metric.label


def _table_metrics(curator: bool) -> List[MetricInfo]:
    return [m for m in map(metric_info, TABLE_METRICS) if curator or not m.curator_only]


def _class_badge(cls: Optional[str]) -> Dict[str, str]:
    info = seq_class_info(cls)
    return {"label": info["label"], "color": info["color"], "tooltip": info["explain"]}


def _flag_cell(flags: List[FlagRow], acks: Acks, curator: bool) -> List[Dict]:
    active, acked = split_acknowledged(flags, acks)
    return flag_badges(active, curator) + acknowledged_badges(acked, acks, curator)


def _table_rows(
    rows: Iterable[Dict[str, Any]],
    by_sample: Dict[str, List[FlagRow]],
    acks: Acks,
    curator: bool,
) -> List[Dict[str, Any]]:
    out = []
    for r in rows:
        sample = r.get("sample") or ""
        row: Dict[str, Any] = {
            "family_id": r.get("family_id") or "",
            "sample": sample,
            "_key": f"{r.get('family_id') or ''}:{sample}",
            "seq_class": "",
            "seq_class_badges": [_class_badge(r.get("seq_class"))]
            if r.get("seq_class")
            else [],
            "flags": "",
            "flags_badges": _flag_cell(by_sample.get(sample, []), acks, curator),
        }
        for m in _table_metrics(curator):
            row[m.key] = _cell(m, r.get(m.key))
        out.append(row)
    return out


def _columns(curator: bool, with_family: bool) -> List[Dict[str, Any]]:
    columns: List[Dict[str, Any]] = []
    if with_family:
        columns.append({"id": "family_id", "header": "Family", "sortable": True})
    columns += [
        {"id": "sample", "header": "Sample", "sortable": True},
        {"id": "seq_class", "header": "Class", "cellType": "badge_list"},
    ]
    columns += [
        {"id": m.key, "header": _header(m), **_NUM} for m in _table_metrics(curator)
    ]
    columns.append({"id": "flags", "header": "Flags", "cellType": "badge_list"})
    return columns


def _flagged_on(metric: MetricInfo, flags: Iterable[FlagRow], acks: Acks) -> Set[str]:
    """Samples with an unacknowledged flag raised on this metric."""
    return {f.sample_a for f in flags if f.code in metric.flags and f not in acks}


def _plot_options(curator: bool) -> Dict[str, str]:
    return {m.key: m.label for m in metrics(curator=curator) if m.plot}


def _render_chart(
    rows: List[Dict[str, Any]],
    reference: List[Dict[str, Any]],
    flags: List[FlagRow],
    acks: Acks,
    curator: bool,
    highlighted: Set[str],
    mark: str,
) -> None:
    """A metric selector and its chart."""
    options = _plot_options(curator)
    if not options:
        return
    default = (
        "ebv_copies_per_cell" if "ebv_copies_per_cell" in options else "chr20_depth"
    )
    state = {"metric": default if default in options else next(iter(options))}

    @ui.refreshable
    def chart() -> None:
        metric = metric_info(state["metric"])
        ui.echart(
            class_strip_option(
                metric,
                rows,
                highlighted=highlighted,
                flagged=_flagged_on(metric, flags, acks),
                reference=reference,
            )
        ).classes("w-full").style(f"height: {chart_height(rows)}px").mark(
            f"{mark}-chart"
        )
        if metric.explain:
            ui.label(metric.explain).classes("text-xs text-gray-500")

    def on_metric(e: Any) -> None:
        state["metric"] = e.value
        chart.refresh()

    ui.select(
        options, value=state["metric"], label="Metric", on_change=on_metric
    ).props("outlined dense").classes("w-64").mark(f"{mark}-metric")
    chart()


def _footnote(curator: bool) -> None:
    ui.label(
        "Measured on each CRAM's chr20 (about 2% of the genome) and the family's"
        " variant calls (GQ ≥ 20, DP ≥ 10). Each row of a chart is one class of"
        " read length and library: the pipeline compares a sample with its own"
        " class only (solid line: median; dotted: 5 scaled MADs, with 20 samples"
        " or more)."
        + (
            " EBV genomes per cell: blood carries almost none; a cell line 1 or more"
            " (dashed lines at 0.01 and 1)."
            if curator
            else ""
        )
    ).classes("text-xs text-gray-500")


def render_cell_line_badge(
    flags: Iterable[FlagRow], on_click: Optional[Callable[[], Any]] = None
) -> None:
    """The family header's "Cell line" badge, acknowledged or not.

    The QC badge counts the flag until a curator acknowledges it; this one
    stays, because the caveat on the de novo, UPD and CNV results does.
    """
    samples = sorted({f.sample_a for f in flags if f.code == CELL_LINE})
    if not samples:
        return
    badge = (
        ui.badge(f"Cell line: {', '.join(samples)}", color="#880e4f")
        .classes("cursor-pointer")
        .mark("cell-line-badge")
    )
    if on_click is not None:
        badge.on("click", on_click)
    badge.tooltip(f"{flag_info(CELL_LINE).label}. {CELL_LINE_TEXT}")


def render_cell_line_callout(
    rows: List[Dict[str, Any]], flags: List[FlagRow], curator: bool
) -> None:
    """A standing warning for every cell-line member, acknowledged or not."""
    by_sample = {f.sample_a: f for f in flags if f.code == CELL_LINE}
    if not by_sample:
        return
    ebv = {r.get("sample"): r.get("ebv_copies_per_cell") for r in rows}
    with (
        ui.card()
        .classes("w-full bg-red-50 border border-red-200")
        .mark("genome-qc-cell-line")
    ):
        for sample, flag in by_sample.items():
            label = flag_info(CELL_LINE).label
            detail = ""
            if curator and ebv.get(sample) is not None:
                detail = f" ({metric_info('ebv_copies_per_cell').fmt(ebv[sample])} EBV genomes per cell)"
            with ui.row().classes("items-start gap-2 no-wrap"):
                ui.icon("biotech", color="red-8").classes("mt-0.5")
                ui.label(f"{sample}: {label}{detail}. {CELL_LINE_TEXT}").classes(
                    "text-sm text-red-900"
                )


def render_family_genome_qc(
    gqc: Optional[GenomeQc],
    reference: List[Dict[str, Any]],
    members: List[str],
    acks: Acks,
    curator: bool,
    stats_links: Dict[str, str],
    loading: bool,
) -> None:
    """The family tab's section: its members over the cohort, by class.

    Args:
        gqc: The family's genome QC (``load_family_gqc``), ``None`` without.
        reference: The cohort's rows, drawn grey behind the members.
        members: The selected members.
        stats_links: Sample -> URL of its samtools stats (curators only).
        loading: The tables are still being read.
    """
    if loading:
        with ui.row().classes("items-center gap-2"):
            ui.spinner(size="sm")
            ui.label("Reading the genome QC tables…").classes("text-xs text-gray-500")
        return
    if gqc is None:
        ui.label("No genome QC (ghfc-ngs genomeqc step) for this family.").classes(
            "text-gray-500 italic"
        )
        return
    member_set = set(members)
    rows = [r for r in gqc.sample_rows() if r.get("sample") in member_set]
    flags = [f for f in gqc.flags if f.sample_a in member_set]
    by_sample: Dict[str, List[FlagRow]] = {}
    for f in flags:
        by_sample.setdefault(f.sample_a, []).append(f)

    render_cell_line_callout(rows, flags, curator)
    if not gqc.flagged:
        ui.label(
            "Read-length classes and flags are set over the whole cohort, and this"
            " family's cohort has no genome QC table: metrics only."
        ).classes("text-xs text-gray-500")
    missing = [m for m in members if m not in {r.get("sample") for r in rows}]
    if missing:
        ui.label(f"Not measured: {', '.join(missing)}.").classes(
            "text-xs text-gray-500"
        )
    with ui.element("div").classes("w-full min-w-0"):
        DataTable(
            columns=_columns(curator, with_family=False),
            rows=_table_rows(rows, by_sample, acks, curator),
            row_key="sample",
            pagination={"rowsPerPage": 10},
        )
    links = {s: u for s, u in stats_links.items() if s in member_set}
    if curator and links:
        with ui.row().classes("items-center gap-3"):
            ui.label("samtools stats (chr20):").classes("text-xs text-gray-600")
            for sample, url in sorted(links.items()):
                ui.link(f"{sample} ↗", url, new_tab=True).classes("text-xs")
    if rows:
        background = [r for r in reference if r.get("sample") not in member_set]
        _render_chart(
            background + rows,
            reference or rows,
            flags + [f for f in gqc.flags if f.sample_a not in member_set],
            acks,
            curator,
            highlighted=member_set,
            mark="genome-qc",
        )
    _footnote(curator)


def render_cohort_genome_qc(
    gqc: Optional[GenomeQc],
    keys: Set[Key],
    acks: Acks,
    curator: bool,
    loading: bool,
) -> None:
    """The cohort tab's section: the shown individuals, by class."""
    if loading:
        with ui.row().classes("items-center gap-2"):
            ui.spinner(size="sm")
            ui.label("Reading the genome QC tables…").classes("text-xs text-gray-500")
        return
    if gqc is None:
        ui.label("No genome QC (ghfc-ngs genomeqc step) for this cohort.").classes(
            "text-gray-500 italic"
        )
        return
    reference = gqc.sample_rows()
    rows = [
        r
        for r in reference
        if ((r.get("family_id") or ""), (r.get("sample") or "")) in keys
    ]
    shown = {r.get("sample") for r in rows}
    flags = [f for f in gqc.flags if f.sample_a in shown]
    by_sample: Dict[str, List[FlagRow]] = {}
    for f in flags:
        by_sample.setdefault(f.sample_a, []).append(f)

    counts: Dict[str, int] = {}
    for r in rows:
        cls = r.get("seq_class") or "unknown"
        counts[cls] = counts.get(cls, 0) + 1
    with ui.row().classes("items-center gap-2"):
        ui.label(
            f"{len(rows)} of {len(keys)} shown individuals measured ({gqc.label})."
        ).classes("text-xs text-gray-500")
        for cls, n in counts.items():
            info = seq_class_info(cls)
            ui.badge(f"{info['label']}: {n}", color=info["color"]).tooltip(
                info["explain"]
            )
        cell_lines = sum(1 for f in flags if f.code == CELL_LINE)
        if cell_lines:
            ui.badge(
                f"{flag_info(CELL_LINE).label}: {cell_lines}", color="#c62828"
            ).tooltip(CELL_LINE_TEXT).mark("cohort-genome-qc-cell-lines")
    if rows:
        _render_chart(
            rows,
            reference,
            flags,
            acks,
            curator,
            highlighted=set(),
            mark="cohort-genome-qc",
        )
    with ui.element("div").classes("w-full min-w-0"):
        DataTable(
            columns=_columns(curator, with_family=True),
            rows=_table_rows(rows, by_sample, acks, curator),
            row_key="_key",
            pagination={"rowsPerPage": 25},
        )
    _footnote(curator)
