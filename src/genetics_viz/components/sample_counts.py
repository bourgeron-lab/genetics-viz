"""Variants-per-sample section of the family Stats dialogs (Wombat and SVs).

Two parts: a bar chart of each member's count in the current table view, and
each member placed on the cohort's per-sample distribution. The second loads
off-thread, so the dialog opens at once. See ``utils/sample_counts.py``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

from nicegui import run, ui

from genetics_viz.models import Cohort
from genetics_viz.utils.sample_counts import (
    CountKind,
    GroupStats,
    MemberComparison,
    MissingHeaderError,
    cohort_sample_counts,
    compare_to_cohort,
    histogram_bins,
)

logger = logging.getLogger(__name__)

_FLAG_COLORS = {"high": "#dc2626", "low": "#d97706"}
_FLAG_LABELS = {"high": "High", "low": "Low"}
_MEMBER_COLOR = "#1e40af"
_HIST_COLOR = "#94a3b8"


def sample_count_bar_options(
    per_member: Dict[str, Dict[str, int]],
    categories: List[str],
    colors: Dict[str, str],
) -> Dict[str, Any]:
    """ECharts options: one stacked bar per member, one stack per category."""
    members = list(per_member)
    return {
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "shadow"}},
        "legend": {"data": categories, "top": 0},
        "grid": {"top": 30, "containLabel": True},
        "xAxis": {"type": "category", "data": members, "name": "Sample"},
        "yAxis": {"type": "value", "name": "Count", "minInterval": 1},
        "series": [
            {
                "name": cat,
                "type": "bar",
                "stack": "total",
                "data": [per_member[m].get(cat, 0) for m in members],
                "itemStyle": {"color": colors.get(cat, _HIST_COLOR)},
            }
            for cat in categories
        ],
    }


def render_sample_count_bars(
    per_member: Dict[str, Dict[str, int]],
    categories: List[str],
    colors: Dict[str, str],
    title: str,
) -> None:
    """The current-view bar chart, with its heading."""
    ui.label(title).classes("text-lg font-semibold text-gray-800 mt-4")
    if not per_member:
        ui.label("No family member selected").classes("text-gray-500 italic")
        return
    ui.echart(sample_count_bar_options(per_member, categories, colors)).classes(
        "w-full h-56"
    ).mark("sample-count-bars")


def _histogram_options(
    group: GroupStats, members: List[MemberComparison], unit: str
) -> Dict[str, Any]:
    """A group's per-sample count histogram, one marker line per member's bin.

    Members in the same bin share a line, so their labels do not overlap; the
    line is red if any of them is flagged.
    """
    starts, width = histogram_bins(group.counts)
    freq = [0] * len(starts)
    for c in group.counts:
        freq[min(c // width, len(starts) - 1)] += 1
    labels = [str(s) if width == 1 else f"{s}-{s + width - 1}" for s in starts]

    by_bin: Dict[int, List[MemberComparison]] = {}
    for m in members:
        by_bin.setdefault(min(m.count // width, len(labels) - 1), []).append(m)
    mark_lines = []
    for idx, in_bin in sorted(by_bin.items()):
        flags = [m.flag for m in in_bin if m.flag]
        mark_lines.append(
            {
                "xAxis": labels[idx],
                "label": {
                    "formatter": ", ".join(f"{m.sample_id}: {m.count}" for m in in_bin),
                    "position": "end",
                },
                "lineStyle": {
                    "color": _FLAG_COLORS[flags[0]] if flags else _MEMBER_COLOR,
                    "width": 2,
                    "type": "solid",
                },
            }
        )
    return {
        "title": {
            "text": f"Cohort {group.name} (n={group.n})",
            "left": "center",
            "textStyle": {"fontSize": 13},
        },
        "tooltip": {"trigger": "axis", "axisPointer": {"type": "shadow"}},
        "grid": {"top": 55, "left": 30, "right": 30, "containLabel": True},
        "xAxis": {
            "type": "category",
            "data": labels,
            "name": unit,
            "nameLocation": "middle",
            "nameGap": 28,
        },
        "yAxis": {
            "type": "value",
            "name": "Samples",
            "nameLocation": "middle",
            "nameGap": 30,
            "minInterval": 1,
        },
        "series": [
            {
                "type": "bar",
                "name": "Samples",
                "data": freq,
                "itemStyle": {"color": _HIST_COLOR},
                "barCategoryGap": "10%",
                "markLine": {"symbol": "none", "data": mark_lines},
            }
        ],
    }


def _fmt(v: float) -> str:
    return f"{v:g}"


def _render_comparison(
    comparisons: List[MemberComparison], missing: List[str], unit: str
) -> None:
    """The per-member table and the group histograms."""
    if comparisons:
        with (
            ui.element("div")
            .classes("grid gap-x-6 gap-y-1 text-sm items-center mt-2")
            .style("grid-template-columns: repeat(6, max-content)")
        ):
            for head in ("Sample", "Group", unit, "Group median [Q1-Q3]", "Pct", ""):
                ui.label(head).classes("font-semibold text-gray-600")
            for m in comparisons:
                ui.label(m.sample_id).classes("font-mono")
                ui.label(m.group.name)
                ui.label(str(m.count)).classes("font-semibold")
                ui.label(
                    f"{_fmt(m.group.median)} [{_fmt(m.group.q1)}-{_fmt(m.group.q3)}]"
                )
                ui.label(f"{m.percentile:.0f}%")
                if m.flag:
                    ui.badge(_FLAG_LABELS[m.flag]).style(
                        f"background-color: {_FLAG_COLORS[m.flag]}"
                    ).tooltip(
                        f"Outside Q1-Q3 by more than 1.5 × IQR "
                        f"(fence {_fmt(m.group.low_fence)} to "
                        f"{_fmt(m.group.high_fence)})"
                    )
                else:
                    ui.label("")

    if missing:
        ui.label(f"Not in the cohort pedigree: {', '.join(missing)}").classes(
            "text-gray-500 italic text-sm mt-1"
        )

    groups: Dict[str, List[MemberComparison]] = {}
    for m in comparisons:
        groups.setdefault(m.group.name, []).append(m)
    if groups:
        with ui.row().classes("w-full gap-4 flex-wrap mt-2"):
            for members in groups.values():
                with ui.column().classes("flex-1 min-w-[400px]"):
                    ui.echart(
                        _histogram_options(members[0].group, members, unit)
                    ).classes("w-full h-64").mark("sample-count-histogram")


def render_cohort_comparison(
    path: Optional[Path],
    kind: CountKind,
    cohort: Optional[Cohort],
    member_ids: List[str],
    unit: str = "Variants",
    header_source: Optional[Path] = None,
) -> None:
    """Compare members with their cohort group, loading the counts off-thread.

    ``header_source`` is the family's own file of the same config, lending
    its header to a cohort file written without one.

    Degraded states (no cohort file, unknown cohort, unreadable or headerless
    file) are a grey line, never an exception.
    """
    ui.label("Compared with the cohort (unfiltered)").classes(
        "text-lg font-semibold text-gray-800 mt-4"
    )
    if cohort is None:
        ui.label("The cohort is not loaded; no comparison available.").classes(
            "text-gray-500 italic"
        )
        return
    if path is None:
        ui.label(
            "No cohort-level file for this configuration; no comparison available."
        ).classes("text-gray-500 italic")
        return
    if not member_ids:
        ui.label("No family member selected").classes("text-gray-500 italic")
        return

    ui.label(
        f"Raw counts from {path.name}, the table's filters not applied. "
        "Founders and non-founders are compared within their own group; "
        "samples without a call count as 0."
    ).classes("text-xs text-gray-500")

    body = ui.column().classes("w-full")
    with body:
        with ui.row().classes("items-center gap-2"):
            ui.spinner(size="md", color="blue")
            ui.label("Counting the cohort's variants...").classes(
                "text-sm text-gray-500"
            )

    async def load() -> None:
        try:
            counts = await run.io_bound(cohort_sample_counts, path, kind, header_source)
        except MissingHeaderError as exc:
            body.clear()
            with body:
                ui.label(f"{exc}; no comparison available.").classes(
                    "text-gray-500 italic"
                )
            return
        except Exception:
            logger.warning("Failed to count %s", path, exc_info=True)
            body.clear()
            with body:
                ui.label(f"Could not read {path.name}").classes("text-gray-500 italic")
            return
        comparisons, missing = compare_to_cohort(counts, cohort, member_ids)
        body.clear()
        with body:
            _render_comparison(comparisons, missing, unit)

    ui.timer(0.1, load, once=True)
