"""Cohort-wide ROH sub-tab for the cohort page's right-hand panel.

Follows the filters of the individuals table. Reads only the cohort ``froh``
table - the segment and gene tables run to tens of megabytes at cohort level
and the family page already shows them per family. Blocks:

1. violins of a metric per reference group, the cohort beside the unrelated
   reference genomes of that region, with the fixed thresholds as lines; or
   the cohort alone by phenotype or sex;
2. the population-vs-recent plane of one reference group;
3. counts of grades and patterns per reference group;
4. the individuals worth a look, linked to their family;
5. a button opening the genome-wide ROH islands of the cohort's founders,
   which do not follow the filters (the pipeline counts founders only).

The violins are Plotly, as for the cohort PGS: ECharts has no violin series.
"""

import asyncio
import logging
import math
from collections import Counter
from typing import Any, Callable, Dict, List, Optional, Tuple

import polars as pl
from nicegui import ui

from genetics_viz.components.download_button import Download, download_button
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.utils.ancestry import Variant, get_region_color
from genetics_viz.utils.auth import can_write
from genetics_viz.utils.export import export_filename, person, rows_to_tsv
from genetics_viz.utils.pedigree_labels import (
    PHENO_COLORS,
    PHENO_ORDER,
    SEX_COLORS,
    SEX_ORDER,
    normalize_pheno,
    normalize_sex,
    pheno_label,
    sex_label,
)
from genetics_viz.utils.roh import (
    froh_export,
    EXCESS_PERCENTILE,
    F_BACKGROUND,
    GRADE_HELP,
    GRADES,
    PATTERNS,
    POOLED,
    find_cohort_roh_files,
    flag_badges,
    get_cohort_roh_dir,
    get_reference_roh_dir,
    grade_badge,
    grade_for,
    is_flagged,
    load_reference_froh,
    load_reference_hotspots,
    metric_value,
    ordinal,
    pattern_badge,
    probe_cohort_roh_data,
    quantile,
    read_roh_table,
    reference_values,
    threshold_grades,
    to_float,
)
from genetics_viz.utils.roh_plots import islands_option, pattern_scatter_option

logger = logging.getLogger(__name__)

__all__ = ["probe_cohort_roh_data", "render_cohort_roh_tab"]

_METRICS = {
    "froh_1.5mb": "F_ROH (ROH ≥ 1.5 Mb)",
    "f_long": "F long (≥ 10 cM, recent)",
    F_BACKGROUND: "F short + intermediate (population)",
}

_GROUP_REFERENCE = "Reference group"
_GROUP_PHENOTYPE = "Phenotype"
_GROUP_SEX = "Sex"

_REFERENCE_COLOR = "#9e9e9e"

#: The tool's floor on f_long for a "recent" reading (MIN_LONG_F_FOR_RECENT).
_MIN_LONG_F = 0.005


def _rgba(hex_color: str, alpha: float) -> str:
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def _table_box() -> ui.element:
    """A block box the table can scroll inside.

    ``nicegui-card`` and ``nicegui-column`` size their children by content, so
    a bare table as wide as its nine columns pushes the right-hand panel past
    the page instead of scrolling.
    """
    return ui.element("div").classes("w-full min-w-0")


def _group_name(group: str) -> str:
    return "pooled" if group == POOLED else group


# ---------------------------------------------------------------------------
# 1. Violins
# ---------------------------------------------------------------------------


def _violin_figure(
    paired: List[Tuple[Dict, Dict]],
    reference: Optional[pl.DataFrame],
    metric: str,
    group_by: str,
    log_x: bool,
    show_relationship: bool,
) -> Any:
    """Horizontal violins: one row per group, values along x.

    Horizontal because the panel is narrow and there are up to nine reference
    groups: rotated tick labels under vertical violins were unreadable. Split
    by reference group, the cohort is the upper half of each violin and its
    unrelated reference genomes the lower half.
    """
    import plotly.graph_objects as go

    fig = go.Figure()
    common = dict(
        orientation="h",
        box_visible=True,
        meanline_visible=True,
        width=0.9,
        spanmode="hard",
    )

    if group_by == _GROUP_REFERENCE:
        by_group: Dict[str, List[float]] = {}
        for row, _ in paired:
            value = metric_value(row, metric)
            if value is not None:
                by_group.setdefault(row.get("ref_group") or POOLED, []).append(value)
        # Reverse: plotly draws the first category at the bottom.
        order = sorted(by_group, reverse=True)
        labels = {g: f"{_group_name(g)} (n={len(by_group[g])})" for g in order}
        xs, ys = [], []
        for g in order:
            xs += by_group[g]
            ys += [labels[g]] * len(by_group[g])
        fig.add_trace(
            go.Violin(
                x=xs,
                y=ys,
                name="Cohort",
                side="positive",
                line_color="#1e88e5",
                fillcolor=_rgba("#1e88e5", 0.35),
                points="all",
                jitter=0.3,
                pointpos=0.5,
                marker=dict(size=3),
                **common,
            )
        )
        if reference is not None:
            rx, ry = [], []
            for g in order:
                values = reference_values(reference, g, metric)
                rx += values
                ry += [labels[g]] * len(values)
            fig.add_trace(
                go.Violin(
                    x=rx,
                    y=ry,
                    name="Reference",
                    side="negative",
                    line_color=_REFERENCE_COLOR,
                    fillcolor=_rgba(_REFERENCE_COLOR, 0.35),
                    points=False,
                    **common,
                )
            )
        fig.update_layout(violinmode="overlay", violingap=0.1, showlegend=True)
        n_rows = len(order)
    else:
        groups: Dict[str, List[float]] = {}
        for row, ind in paired:
            value = metric_value(row, metric)
            if value is None:
                continue
            if group_by == _GROUP_SEX:
                key = normalize_sex(ind.get("Sex"))
            else:
                key = normalize_pheno(ind.get("Phenotype"))
            groups.setdefault(key, []).append(value)
        order_all = SEX_ORDER if group_by == _GROUP_SEX else PHENO_ORDER
        present = [g for g in order_all if g in groups]
        for g in reversed(present):
            label = sex_label(g) if group_by == _GROUP_SEX else pheno_label(g)
            palette = SEX_COLORS if group_by == _GROUP_SEX else PHENO_COLORS
            color = palette.get(g, "#7f7f7f")
            name = f"{label} (n={len(groups[g])})"
            fig.add_trace(
                go.Violin(
                    x=groups[g],
                    y=[name] * len(groups[g]),
                    name=name,
                    line_color=color,
                    fillcolor=_rgba(color, 0.35),
                    points="all",
                    jitter=0.3,
                    pointpos=0,
                    marker=dict(size=3),
                    **common,
                )
            )
        fig.update_layout(showlegend=False)
        n_rows = len(present)

    if metric == "froh_1.5mb":
        for grade in threshold_grades():
            text = grade.grade
            if show_relationship:
                text += f" ({grade.category})"
            fig.add_vline(
                # Shapes on a log axis take log10 coordinates, not data values.
                x=math.log10(grade.threshold) if log_x else grade.threshold,
                line_dash="dash",
                line_color=grade.color,
                line_width=1,
                annotation_text=f"{text} > {grade.threshold:g}",
                annotation_font_size=10,
                annotation_font_color=grade.color,
                annotation_position="top",
                annotation_textangle=-90,
            )

    fig.update_layout(
        title=dict(text=_METRICS[metric], font=dict(size=13)),
        xaxis_title=_METRICS[metric],
        xaxis_type="log" if log_x else "linear",
        margin=dict(l=10, r=10, t=60, b=40),
        height=140 + 70 * max(n_rows, 1),
        legend=dict(orientation="h", y=1.12, x=1, xanchor="right"),
    )
    fig.update_yaxes(automargin=True)
    return fig


# ---------------------------------------------------------------------------
# 3. Counts and 4. flagged individuals
# ---------------------------------------------------------------------------


def _render_counts(rows: List[Dict]) -> None:
    by_group: Dict[str, List[Dict]] = {}
    for r in rows:
        by_group.setdefault(r.get("ref_group") or POOLED, []).append(r)
    grades = [g for g in GRADES if g.category != "none"]
    table = []
    for group in sorted(by_group):
        members = by_group[group]
        cats = Counter(r.get("froh_category") or "none" for r in members)
        pats = Counter(r.get("homozygosity_pattern") or "typical" for r in members)
        entry: Dict[str, Any] = {"group": _group_name(group), "n": len(members)}
        for g in grades:
            entry[f"g_{g.category}"] = cats.get(g.category, 0)
        for p in PATTERNS:
            if p != "typical":
                entry[f"p_{p}"] = pats.get(p, 0)
        table.append(entry)

    num = {"cellType": "number", "sorting": "numerical", "sortable": True}
    columns: List[Dict[str, Any]] = [
        {"id": "group", "header": "Reference group", "sortable": True},
        {"id": "n", "header": "n", **num},
    ]
    columns += [{"id": f"g_{g.category}", "header": g.grade, **num} for g in grades]
    columns += [
        {"id": f"p_{p}", "header": info["label"], **num}
        for p, info in PATTERNS.items()
        if p != "typical"
    ]
    with _table_box():
        DataTable(
            columns=columns,
            rows=table,
            row_key="group",
            pagination={"rowsPerPage": 10},
        )
    ui.label("Grade columns, then pattern columns. " + GRADE_HELP).classes(
        "text-xs text-gray-500"
    )


def _grade_rank(category: Optional[str]) -> int:
    """Sort key of a category: the strongest grade highest."""
    order = [g.category for g in GRADES]
    return len(order) - order.index(category) if category in order else 0


def _render_flagged(
    rows: List[Dict], cohort_name: str, show_relationship: bool
) -> None:
    flagged = [r for r in rows if is_flagged(r)]
    # Highest F_ROH first, matching the sort indicator: the table does not
    # reorder its initial rows by itself.
    flagged.sort(key=lambda r: -(to_float(r.get("froh_1.5mb")) or 0))
    ui.label(f"{len(flagged)} of {len(rows)} individuals flagged").classes(
        "text-xs text-gray-500"
    )
    if not flagged:
        return
    table = []
    for r in flagged:
        pct = to_float(r.get("pct_froh_1.5mb"))
        table.append(
            {
                "key": f"{r.get('family_id')}|{r.get('IID')}",
                "IID": r.get("IID"),
                "FID": r.get("family_id"),
                "region": r.get("region") or "—",
                "region_color": get_region_color(r.get("region")),
                "froh": to_float(r.get("froh_1.5mb")),
                "pct": ordinal(pct),
                "pct_value": pct,
                "f_long": to_float(r.get("f_long")),
                "pattern": r.get("homozygosity_pattern"),
                "pattern_badges": [pattern_badge(r.get("homozygosity_pattern"))],
                "grade": r.get("froh_category"),
                "grade_badges": [
                    grade_badge(r.get("froh_category"), show_relationship)
                ],
                "grade_rank": _grade_rank(r.get("froh_category")),
                "flags": "",
                "flags_badges": flag_badges(r),
            }
        )
    num = {"cellType": "number", "sorting": "numerical", "sortable": True}
    with _table_box():
        DataTable(
            columns=[
                {"id": "IID", "header": "Sample", "sortable": True},
                {
                    "id": "FID",
                    "header": "Family",
                    "cellType": "link",
                    "href": f"/cohort/{cohort_name}/family/{{FID}}",
                    "sortable": True,
                },
                {
                    "id": "region",
                    "header": "Region",
                    "cellType": "badge",
                    "colorField": "region_color",
                    "sortable": True,
                },
                {"id": "froh", "header": "F_ROH", **num},
                {
                    "id": "pct",
                    "header": "Pct",
                    "sortable": True,
                    "sortField": "pct_value",
                    "sorting": "numerical",
                },
                {"id": "f_long", "header": "F long", **num},
                {"id": "pattern", "header": "Pattern", "cellType": "badge_list"},
                {
                    "id": "grade",
                    "header": "Grade",
                    "cellType": "badge_list",
                    "sortable": True,
                    "sortField": "grade_rank",
                    "sorting": "numerical",
                },
                {"id": "flags", "header": "Flags", "cellType": "badge_list"},
            ],
            rows=table,
            row_key="key",
            pagination={"rowsPerPage": 15},
            initial_sorting=[{"id": "froh", "desc": True}],
        )


# ---------------------------------------------------------------------------
# 5. Islands dialog
# ---------------------------------------------------------------------------


async def _open_islands(path: Any, bundle_version: str) -> None:
    with ui.dialog().props("full-width") as dialog, ui.card().classes("w-full"):
        with ui.row().classes("w-full items-center"):
            ui.label("ROH islands among the cohort's founders").classes(
                "text-lg font-semibold text-blue-700"
            )
            ui.space()
            ui.button(icon="close", on_click=dialog.close).props("flat round dense")
        ui.label(
            "Per 100 kb window, the fraction of the cohort's founders (both parents"
            " absent from their family) with a ROH there. Founders only, so a family"
            " does not count its shared segments several times. Shaded: reference ROH"
            " hotspots of the same group, windows homozygous in many outbred people."
            " This view does not follow the individuals-table filters."
        ).classes("text-xs text-gray-500")
        body = ui.column().classes("w-full")
        with body:
            ui.spinner(size="lg")
    dialog.open()

    try:
        windows = await asyncio.to_thread(read_roh_table, path)
    except Exception:
        logger.warning("Failed to read %s", path, exc_info=True)
        windows = None
    hotspots = load_reference_hotspots(bundle_version)

    body.clear()
    with body:
        if windows is None or len(windows) == 0:
            ui.label("No ROH windows in the cohort table.").classes(
                "text-gray-500 italic"
            )
            return
        groups = windows["group"].unique().to_list()
        groups = [POOLED] + sorted(g for g in groups if g != POOLED)
        founders = {
            r["group"]: r["n_founders"]
            for r in windows.group_by("group")
            .agg(pl.col("n_founders").max())
            .to_dicts()
        }
        if hotspots is None:
            ui.label(
                "Reference hotspots unavailable for this bundle; showing the cohort only."
            ).classes("text-orange-600 text-sm")
        selected = {"group": POOLED}

        @ui.refreshable
        def plot() -> None:
            group = selected["group"]
            sub = windows.filter(pl.col("group") == group)
            rows = [
                (r["chrom"], int(r["start"]), int(r["end"]), float(r["freq"]))
                for r in sub.select("chrom", "start", "end", "freq").iter_rows(
                    named=True
                )
            ]
            hot = []
            if hotspots is not None:
                hot = [
                    (h["chrom"], int(h["start"]), int(h["end"]))
                    for h in hotspots.filter(pl.col("group") == group).to_dicts()
                ]
            ui.echart(
                islands_option(
                    rows,
                    hot,
                    title=f"{_group_name(group)} — {founders.get(group, 0)} founders",
                )
            ).classes("w-full h-[28rem]")

        ui.select(
            options=groups,
            value=POOLED,
            label="Group",
            on_change=lambda e: (selected.update({"group": e.value}), plot.refresh()),
        ).props("outlined dense").classes("w-48")
        plot()


# ---------------------------------------------------------------------------
# Tab
# ---------------------------------------------------------------------------


def render_cohort_roh_tab(
    store: Any,
    cohort: Any,
    filtered_state: Dict[str, Any],
    refreshers: List[Callable[[], None]],
) -> None:
    """Render the cohort-wide ROH sub-tab.

    Args:
        store: DataStore instance
        cohort: Cohort instance
        filtered_state: Dict with 'individuals' key holding the filtered rows
        refreshers: List to append the refresh callback to
    """
    roh_dir = get_cohort_roh_dir(cohort)
    files = find_cohort_roh_files(cohort)
    if files is None or files.path("froh") is None:
        ui.label(f"No cohort-level ROH results in: {roh_dir}").classes(
            "text-gray-500 italic"
        )
        return

    show_relationship = can_write()
    state: Dict[str, Any] = {"files": files, "by_sample": None, "reference": None}

    def load_variant(variant: Optional[Variant] = None) -> None:
        selected = find_cohort_roh_files(cohort, variant) or files
        state["files"] = selected
        try:
            froh = read_roh_table(selected.path("froh"))
        except Exception:
            logger.warning(
                "Failed to read cohort ROH for %s", cohort.name, exc_info=True
            )
            froh = None
        # Keyed by (family, sample): a sample listed under two families has a
        # row in each, scored separately.
        state["by_sample"] = (
            {(r.get("family_id"), r["IID"]): r for r in froh.to_dicts() if r.get("IID")}
            if froh is not None
            else None
        )
        state["reference"] = load_reference_froh(selected.bundle_version)

    load_variant()
    if not state["by_sample"]:
        ui.label("No individuals in the cohort ROH table").classes(
            "text-gray-500 italic"
        )
        return

    controls: Dict[str, Any] = {
        "metric": "froh_1.5mb",
        "group_by": _GROUP_REFERENCE,
        "log_x": False,
        "groups": [],
        "scatter_group": None,
    }

    def shown() -> List[Tuple[Dict, Dict]]:
        """``(froh row, individual)`` of the shown individuals with ROH results."""
        by_sample = state["by_sample"] or {}
        allowed = set(controls["groups"])
        paired: List[Tuple[Dict, Dict]] = []
        seen = set()
        for ind in filtered_state.get("individuals", []):
            key = (ind.get("FID"), ind["Sample ID"])
            # Some pedigrees repeat a row; count each person once.
            if key in seen:
                continue
            seen.add(key)
            row = by_sample.get(key)
            if row is None:
                continue
            if allowed and (row.get("ref_group") or POOLED) not in allowed:
                continue
            paired.append((row, ind))
        return paired

    def export() -> str:
        paired = []
        for row, ind in shown():
            who = person(ind["FID"], ind["Sample ID"], ind["Sex"], ind["Phenotype"])
            paired.append((row, who))
        rows, columns = froh_export(paired, show_relationship)
        return rows_to_tsv(rows, columns) if rows else ""

    with ui.card().classes("w-full min-w-0"):
        with ui.row().classes("items-center gap-2 w-full"):
            ui.label("Runs of homozygosity").classes(
                "text-lg font-semibold text-blue-700"
            )
            provenance = ui.label(state["files"].label).classes("text-xs text-gray-500")
            ui.space()
            download_button(
                [
                    Download(
                        "froh",
                        "F_ROH, pattern and grade of the individuals shown",
                        lambda: export_filename(
                            cohort.name, "roh", tag=state["files"].file_tag
                        ),
                        export,
                    )
                ],
                mark="download-roh",
            )

        variants = state["files"].all_variants
        if len(variants) > 1:

            def on_variant(event: Any) -> None:
                bundle, tag = event.value.split(" / ", 1)
                load_variant((bundle, tag))
                provenance.text = state["files"].label
                render_body.refresh()

            ui.select(
                options=[f"{b} / {t}" for b, t in variants],
                value=f"{state['files'].bundle_version} / {state['files'].filter_tag}",
                label="Bundle / filters",
                on_change=on_variant,
            ).props("outlined dense").classes("w-full mb-2")

        windows_path = state["files"].path("roh_windows")
        islands_btn = ui.button(
            "ROH islands (genome-wide)",
            icon="travel_explore",
            on_click=lambda: _open_islands(
                state["files"].path("roh_windows"), state["files"].bundle_version
            ),
        ).props("outline dense no-caps")
        if windows_path is None:
            islands_btn.disable()
            islands_btn.tooltip("No roh_windows table for this cohort")

        def update(**kw: Any) -> None:
            controls.update(kw)
            render_body.refresh()

        with ui.row().classes("items-center gap-2 w-full mt-2"):
            ui.select(
                options=_METRICS,
                value=controls["metric"],
                label="Metric",
                on_change=lambda e: update(metric=e.value),
            ).props("outlined dense").classes("flex-1 min-w-0")
            ui.select(
                options=[_GROUP_REFERENCE, _GROUP_PHENOTYPE, _GROUP_SEX],
                value=controls["group_by"],
                label="Group by",
                on_change=lambda e: update(group_by=e.value),
            ).props("outlined dense").classes("flex-1 min-w-0")
            ui.checkbox(
                "Log scale", value=False, on_change=lambda e: update(log_x=e.value)
            ).props("dense")

        all_groups = sorted(
            {r.get("ref_group") or POOLED for r in state["by_sample"].values()}
        )
        ui.select(
            options=all_groups,
            value=[],
            label="Reference groups (all when empty)",
            multiple=True,
            on_change=lambda e: update(groups=e.value or []),
        ).props("use-chips outlined dense").classes("w-full mt-1")

        @ui.refreshable
        def render_body() -> None:
            reference = state["reference"]
            paired = shown()

            if reference is None:
                ref_dir = get_reference_roh_dir(state["files"].bundle_version)
                detail = (
                    "set 'ancestry_reference_dir' in the application config"
                    if ref_dir is None
                    else f"expected at {ref_dir}"
                )
                ui.label(
                    f"ROH reference unavailable for bundle {state['files'].bundle_version}"
                    f" — {detail}. The cohort is shown without the reference."
                ).classes("text-orange-600 text-sm")

            if not paired:
                ui.label(
                    "No individuals with ROH results match the current filters"
                ).classes("text-gray-500 italic")
                return
            rows = [r for r, _ in paired]
            ui.label(f"{len(paired)} individuals with ROH results").classes(
                "text-xs text-gray-500"
            )

            fig = _violin_figure(
                paired,
                reference,
                controls["metric"],
                controls["group_by"],
                controls["log_x"],
                show_relationship,
            )
            ui.plotly(fig).classes("w-full")
            caption = (
                "Each region beside its unrelated reference genomes (grey)."
                if controls["group_by"] == _GROUP_REFERENCE
                else "Cohort only: groups mix ancestry regions, so restrict to one"
                " reference group before comparing them."
            )
            if controls["log_x"]:
                caption += " Log scale hides zero values."
            ui.label(caption).classes("text-xs text-gray-500")

            # Population vs recent, for one reference group at a time: the
            # percentile lines belong to a group.
            group_counts = Counter(r.get("ref_group") or POOLED for r in rows)
            groups = [g for g, _ in group_counts.most_common()]
            if controls["scatter_group"] not in group_counts:
                controls["scatter_group"] = groups[0]
            group = controls["scatter_group"]

            ui.label("Population vs recent").classes(
                "text-base font-semibold text-blue-700 mt-4"
            )
            ui.select(
                options=groups,
                value=group,
                label="Reference group",
                on_change=lambda e: update(scatter_group=e.value),
            ).props("outlined dense").classes("w-48")
            ref_bg = (
                reference_values(reference, group, F_BACKGROUND)
                if reference is not None
                else []
            )
            ref_long = (
                reference_values(reference, group, "f_long")
                if reference is not None
                else []
            )
            q = EXCESS_PERCENTILE / 100
            y_line = quantile(ref_long, q)
            points = []
            for r in rows:
                if (r.get("ref_group") or POOLED) != group:
                    continue
                bg, long = metric_value(r, F_BACKGROUND), to_float(r.get("f_long"))
                if bg is not None and long is not None:
                    points.append(
                        {
                            "name": r["IID"],
                            "x": bg,
                            "y": long,
                            "color": grade_for(r.get("froh_category")).color,
                        }
                    )
            ui.echart(
                pattern_scatter_option(
                    list(zip(ref_bg, ref_long)),
                    points,
                    title=f"{_group_name(group)}: {len(points)} cohort vs {len(ref_bg)} reference",
                    x_line=quantile(ref_bg, q),
                    y_line=max(y_line, _MIN_LONG_F) if y_line is not None else None,
                    show_labels=False,
                    symbol_size=8,
                )
            ).classes("w-full h-80")

            ui.label("Counts").classes("text-base font-semibold text-blue-700 mt-4")
            _render_counts(rows)

            ui.label("Flagged individuals").classes(
                "text-base font-semibold text-blue-700 mt-4"
            )
            _render_flagged(rows, cohort.name, show_relationship)

        render_body()
        refreshers.append(render_body.refresh)
