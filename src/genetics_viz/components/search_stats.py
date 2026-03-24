"""Search page statistics dialog — chromosome distribution, consequence/validation charts, ideogram."""

import re
import statistics
from collections import Counter
from typing import Any, Dict, List

from nicegui import ui

from genetics_viz.utils.cytobands import (
    CHROM_ORDER,
    CHROM_SIZES_MB,
    CYTOBANDS,
    GIESTAIN_COLORS,
    VALIDATION_COLORS,
    norm_chrom,
)
from genetics_viz.utils.vep import (
    VEP_CONSEQUENCES,
    format_consequence_display,
    get_highest_consequence_term,
)
from genetics_viz.utils.wisecondorx import WISECONDORX_CONFIG, infer_sv_type

_SV_VARIANT_PATTERN = re.compile(r"^chr[^:]+:\d+-\d+$")

# Colors for SV gain/loss in charts
_SV_CALL_COLORS = {
    "GAIN": WISECONDORX_CONFIG["robust_gain"]["color"],
    "LOSS": WISECONDORX_CONFIG["robust_loss"]["color"],
    "Unknown": "#94a3b8",
}


def show_stats_dialog(rows: List[Dict[str, Any]]) -> None:
    """Create and open a variant statistics dialog for the given result rows.

    Processes both SNV/Indel and SV (WisecondorX) rows, deduplicates,
    and displays chromosome distribution (stacked bar), consequence
    and validation pie charts, outlier analysis, and an ideogram view.
    """
    # Parse Variant column into #CHROM/POS/REF/ALT if not already present
    # (validation file data uses combined Variant column: chr:pos:ref:alt)
    for r in rows:
        if "#CHROM" not in r and "Variant" in r:
            parts = str(r["Variant"]).split(":")
            if len(parts) == 4:
                r["#CHROM"], r["POS"], r["REF"], r["ALT"] = parts

    # Separate SVs from SNV/Indels
    sv_rows_raw: List[Dict[str, Any]] = []
    snv_rows_raw: List[Dict[str, Any]] = []
    for r in rows:
        variant_str = str(r.get("Variant", ""))
        if r.get("_source_type") == "wisecondorx" or _SV_VARIANT_PATTERN.match(
            variant_str
        ):
            sv_rows_raw.append(r)
        else:
            snv_rows_raw.append(r)

    # Deduplicate SNVs by (#CHROM, POS, REF, ALT)
    seen: set = set()
    unique_snvs: list = []
    for r in snv_rows_raw:
        key = (
            r.get("#CHROM", ""),
            r.get("POS", ""),
            r.get("REF", ""),
            r.get("ALT", ""),
        )
        if key not in seen:
            seen.add(key)
            unique_snvs.append(r)

    # Classify SNV vs Indel
    for r in unique_snvs:
        ref = str(r.get("REF", ""))
        alt = str(r.get("ALT", ""))
        r["_is_snv"] = len(ref) == 1 and len(alt) == 1

    # Deduplicate SVs by (Variant, sample)
    sv_seen: set = set()
    unique_svs: list = []
    for r in sv_rows_raw:
        key = (r.get("Variant", ""), r.get("sample", ""))
        if key not in sv_seen:
            sv_seen.add(key)
            unique_svs.append(r)

    # Parse SV coordinates and classify gain/loss
    for r in unique_svs:
        variant = str(r.get("Variant", ""))
        m = _SV_VARIANT_PATTERN.match(variant)
        if m:
            parts = variant.split(":")
            r["_sv_chrom"] = norm_chrom(parts[0])
            range_parts = parts[1].split("-")
            try:
                r["_sv_start_mb"] = float(range_parts[0]) / 1_000_000
                r["_sv_end_mb"] = float(range_parts[1]) / 1_000_000
            except (ValueError, TypeError, IndexError):
                r["_sv_start_mb"] = 0
                r["_sv_end_mb"] = 0
        else:
            r["_sv_chrom"] = ""
            r["_sv_start_mb"] = 0
            r["_sv_end_mb"] = 0
        # Classify gain/loss using shared inference logic
        sv_type = infer_sv_type(r)
        r["_sv_call"] = "GAIN" if sv_type == "dup" else "LOSS"

    chrom_order = CHROM_ORDER
    chrom_sizes_mb = CHROM_SIZES_MB
    validation_colors = VALIDATION_COLORS

    # Filter state
    type_filter = {"snv": True, "indel": True, "sv": True}
    show_ideogram: Dict[str, bool] = {"value": False}
    _containers: Dict[str, Any] = {"charts": None, "ideo": None}

    with ui.dialog().props("full-width") as stats_dialog, ui.card().classes("w-full"):
        with ui.column().classes("w-full p-4"):
            # Header
            with ui.row().classes("items-center justify-between w-full mb-2"):
                with ui.row().classes("items-center gap-3"):
                    ui.label("Variant Statistics").classes(
                        "text-xl font-bold text-blue-900"
                    )
                    subtitle_label = ui.label("").classes("text-sm text-gray-500")
                    ideogram_btn = ui.button(
                        "Ideogram",
                    ).props("outline color=blue size=sm dense no-caps")
                    snv_cb = (
                        ui.checkbox("SNVs", value=True)
                        .props("dense")
                        .classes("text-sm")
                    )
                    indel_cb = (
                        ui.checkbox("Indels", value=True)
                        .props("dense")
                        .classes("text-sm")
                    )
                    sv_cb = (
                        ui.checkbox("SVs", value=True).props("dense").classes("text-sm")
                    )
                ui.button(
                    icon="close",
                    on_click=lambda: stats_dialog.close(),
                ).props("flat round")

            @ui.refreshable
            def render_stats_content():
                # Filter SNV/Indel variants by type
                filtered_snvs = [
                    r
                    for r in unique_snvs
                    if (type_filter["snv"] and r["_is_snv"])
                    or (type_filter["indel"] and not r["_is_snv"])
                ]
                # Filter SVs
                filtered_svs = list(unique_svs) if type_filter["sv"] else []
                # Combined for validation/outlier stats
                filtered_all = filtered_snvs + filtered_svs

                snv_n = sum(1 for r in filtered_snvs if r["_is_snv"])
                indel_n = sum(1 for r in filtered_snvs if not r["_is_snv"])
                sv_n = len(filtered_svs)
                count_parts = []
                if snv_n:
                    count_parts.append(f"{snv_n} SNVs")
                if indel_n:
                    count_parts.append(f"{indel_n} Indels")
                if sv_n:
                    count_parts.append(f"{sv_n} SVs")
                subtitle_label.text = (
                    f"{len(filtered_all)} unique variants ({', '.join(count_parts)})"
                    if count_parts
                    else "0 variants"
                )

                # Chromosome distribution stacked by validation
                chrom_validation: Dict[str, Dict[str, int]] = {
                    c: {} for c in chrom_order
                }
                # Count SNVs/Indels per chromosome
                for r in filtered_snvs:
                    chrom = norm_chrom(r.get("#CHROM", ""))
                    status = r.get("Validation", "") or "TODO"
                    if chrom in chrom_validation:
                        chrom_validation[chrom][status] = (
                            chrom_validation[chrom].get(status, 0) + 1
                        )
                # Count SVs per chromosome (each SV = 1 count)
                for r in filtered_svs:
                    chrom = r.get("_sv_chrom", "")
                    status = r.get("Validation", "") or "TODO"
                    if chrom in chrom_validation:
                        chrom_validation[chrom][status] = (
                            chrom_validation[chrom].get(status, 0) + 1
                        )

                all_statuses: List[str] = []
                for c in chrom_order:
                    for s in chrom_validation[c]:
                        if s not in all_statuses:
                            all_statuses.append(s)

                # Consequence distribution (SNVs use VEP, SVs use gain/loss)
                consequence_counts: Counter = Counter()
                for r in filtered_snvs:
                    consequence_counts[
                        get_highest_consequence_term(str(r.get("VEP_Consequence", "")))
                    ] += 1
                for r in filtered_svs:
                    consequence_counts[r.get("_sv_call", "Unknown")] += 1

                # Validation distribution (all variant types combined)
                validation_counts = Counter(
                    r.get("Validation", "") or "TODO" for r in filtered_all
                )

                # Scatter data for ideogram (SNVs)
                scatter_data: List[List[Any]] = []
                for r in filtered_snvs:
                    chrom = norm_chrom(r.get("#CHROM", ""))
                    pos = r.get("POS", 0)
                    try:
                        pos_mb = round(float(pos) / 1_000_000, 2)
                    except (ValueError, TypeError):
                        continue
                    if chrom in chrom_sizes_mb:
                        status = r.get("Validation", "") or "TODO"
                        scatter_data.append([pos_mb, chrom, status])

                # SV range data for ideogram
                sv_scatter_data: List[List[Any]] = []
                for r in filtered_svs:
                    chrom = r.get("_sv_chrom", "")
                    start_mb = r.get("_sv_start_mb", 0)
                    end_mb = r.get("_sv_end_mb", 0)
                    call = r.get("_sv_call", "Unknown")
                    if chrom in chrom_sizes_mb:
                        sv_scatter_data.append([start_mb, end_mb, chrom, call])

                # Outlier analysis data
                all_samples = {r.get("sample", "") for r in filtered_all}
                all_samples.discard("")
                has_multiple_samples = len(all_samples) > 1

                # --- Charts container ---
                _containers["charts"] = ui.column().classes("w-full")
                _containers["charts"].set_visibility(not show_ideogram["value"])
                with _containers["charts"]:
                    ui.label("Variants per Chromosome").classes(
                        "text-lg font-semibold text-gray-800 mt-2"
                    )
                    stacked_series = [
                        {
                            "name": status,
                            "type": "bar",
                            "stack": "total",
                            "data": [
                                chrom_validation[c].get(status, 0) for c in chrom_order
                            ],
                            "itemStyle": {
                                "color": validation_colors.get(status, "#94a3b8")
                            },
                        }
                        for status in all_statuses
                    ]
                    ui.echart(
                        {
                            "tooltip": {
                                "trigger": "axis",
                                "axisPointer": {"type": "shadow"},
                            },
                            "legend": {"data": all_statuses, "top": 0},
                            "grid": {"top": 30},
                            "xAxis": {
                                "type": "category",
                                "data": chrom_order,
                                "name": "Chromosome",
                            },
                            "yAxis": {"type": "value", "name": "Count"},
                            "series": stacked_series,
                        }
                    ).classes("w-full h-64")

                    with ui.row().classes("w-full gap-4 flex-wrap mt-4"):
                        # Consequence pie chart
                        with ui.column().classes("flex-1 min-w-[400px]"):
                            ui.label(
                                "Consequence Distribution (highest per variant)"
                            ).classes("text-lg font-semibold text-gray-800")
                            cons_data = []
                            for cons, count in consequence_counts.most_common():
                                if cons in _SV_CALL_COLORS:
                                    color = _SV_CALL_COLORS[cons]
                                    display_name = cons
                                else:
                                    color = VEP_CONSEQUENCES.get(cons, ("", "#6b7280"))[
                                        1
                                    ]
                                    display_name = format_consequence_display(cons)
                                cons_data.append(
                                    {
                                        "name": display_name,
                                        "value": count,
                                        "itemStyle": {"color": color},
                                    }
                                )
                            ui.echart(
                                {
                                    "tooltip": {"trigger": "item"},
                                    "series": [
                                        {
                                            "type": "pie",
                                            "radius": "70%",
                                            "data": cons_data,
                                            "label": {"formatter": "{b}: {c} ({d}%)"},
                                        }
                                    ],
                                }
                            ).classes("w-full h-80")

                        # Validation pie chart
                        with ui.column().classes("flex-1 min-w-[400px]"):
                            ui.label("Validation Status Distribution").classes(
                                "text-lg font-semibold text-gray-800"
                            )
                            val_data = [
                                {
                                    "name": st,
                                    "value": cnt,
                                    "itemStyle": {
                                        "color": validation_colors.get(st, "#6b7280")
                                    },
                                }
                                for st, cnt in validation_counts.most_common()
                            ]
                            ui.echart(
                                {
                                    "tooltip": {"trigger": "item"},
                                    "series": [
                                        {
                                            "type": "pie",
                                            "radius": "70%",
                                            "data": val_data,
                                            "label": {"formatter": "{b}: {c} ({d}%)"},
                                        }
                                    ],
                                }
                            ).classes("w-full h-80")

                        # Outlier analysis (only for multi-sample data)
                        if has_multiple_samples:
                            _render_outlier_section(filtered_all, validation_colors)

                # --- Ideogram container ---
                _containers["ideo"] = ui.column().classes("w-full")
                _containers["ideo"].set_visibility(show_ideogram["value"])
                with _containers["ideo"]:
                    svg_w = 1800
                    lbl_w = 50
                    plot_w = svg_w - lbl_w - 20
                    row_h = 16
                    tri_h = 6
                    row_gap = tri_h + 4
                    svg_h = len(chrom_order) * (row_h + row_gap) + 60
                    max_mb = max(chrom_sizes_mb.values())

                    svg_parts = [
                        f'<svg viewBox="0 0 {svg_w} {svg_h}" '
                        f'xmlns="http://www.w3.org/2000/svg" '
                        f'preserveAspectRatio="xMinYMin meet" '
                        f'style="font-family: sans-serif; width: 100%; height: auto;">'
                    ]

                    axis_y = len(chrom_order) * (row_h + row_gap)
                    for mb_val in range(0, 260, 50):
                        gx = lbl_w + (mb_val / max_mb) * plot_w
                        svg_parts.append(
                            f'<line x1="{gx:.1f}" y1="0" '
                            f'x2="{gx:.1f}" y2="{axis_y}" '
                            f'stroke="#e5e7eb" stroke-width="0.5" '
                            f'stroke-dasharray="3,3"/>'
                        )
                        svg_parts.append(
                            f'<text x="{gx:.1f}" y="{axis_y + 14}" '
                            f'text-anchor="middle" font-size="12" '
                            f'fill="#6b7280">{mb_val}</text>'
                        )
                    svg_parts.append(
                        f'<text x="{svg_w / 2}" y="{axis_y + 32}" '
                        f'text-anchor="middle" font-size="13" '
                        f'fill="#6b7280">Position (Mb)</text>'
                    )

                    for ci, chrom in enumerate(chrom_order):
                        bar_y = ci * (row_h + row_gap) + tri_h
                        bands = CYTOBANDS.get(chrom, [])
                        cs = chrom_sizes_mb.get(chrom, 0)
                        total_w = (cs / max_mb) * plot_w

                        svg_parts.append(
                            f'<text x="{lbl_w - 6}" y="{bar_y + row_h * 0.75}" '
                            f'text-anchor="end" font-size="12" '
                            f'fill="#374151">{chrom}</text>'
                        )
                        for band in bands:
                            bx = lbl_w + (band["start"] / max_mb) * plot_w
                            bw = max(
                                ((band["end"] - band["start"]) / max_mb) * plot_w,
                                0.5,
                            )
                            color = GIESTAIN_COLORS.get(band["stain"], "#e5e7eb")
                            svg_parts.append(
                                f'<rect x="{bx:.1f}" y="{bar_y}" '
                                f'width="{bw:.1f}" height="{row_h}" '
                                f'fill="{color}"/>'
                            )
                        svg_parts.append(
                            f'<rect x="{lbl_w}" y="{bar_y}" '
                            f'width="{total_w:.1f}" height="{row_h}" '
                            f'fill="none" stroke="#9ca3af" '
                            f'stroke-width="0.5" rx="3"/>'
                        )

                    # SNV/Indel scatter triangles
                    for sd in scatter_data:
                        v_mb, v_chrom, v_status = sd
                        if v_chrom not in chrom_order:
                            continue
                        v_idx = chrom_order.index(v_chrom)
                        bar_y = v_idx * (row_h + row_gap) + tri_h
                        vx = lbl_w + (v_mb / max_mb) * plot_w
                        v_color = validation_colors.get(v_status, "#94a3b8")
                        tw = 5
                        svg_parts.append(
                            f'<polygon points="{vx - tw:.1f},{bar_y - tri_h} '
                            f"{vx + tw:.1f},{bar_y - tri_h} "
                            f'{vx:.1f},{bar_y}" '
                            f'fill="{v_color}" opacity="0.9"/>'
                        )
                        svg_parts.append(
                            f'<line x1="{vx:.1f}" y1="{bar_y}" '
                            f'x2="{vx:.1f}" y2="{bar_y + row_h}" '
                            f'stroke="{v_color}" stroke-width="1.5" '
                            f'opacity="0.85"/>'
                        )

                    # SV range lines (colored by gain/loss)
                    for sd in sv_scatter_data:
                        s_mb, e_mb, s_chrom, s_call = sd
                        if s_chrom not in chrom_order:
                            continue
                        s_idx = chrom_order.index(s_chrom)
                        bar_y = s_idx * (row_h + row_gap) + tri_h
                        sx = lbl_w + (s_mb / max_mb) * plot_w
                        ex = lbl_w + (e_mb / max_mb) * plot_w
                        line_w = max(ex - sx, 2)
                        sv_color = _SV_CALL_COLORS.get(s_call, "#94a3b8")
                        line_y = bar_y - tri_h / 2
                        svg_parts.append(
                            f'<rect x="{sx:.1f}" y="{line_y - 2:.1f}" '
                            f'width="{line_w:.1f}" height="4" '
                            f'fill="{sv_color}" opacity="0.7" rx="1"/>'
                        )

                    # Legend
                    legend_y = axis_y + 40
                    legend_x = lbl_w
                    # Validation status legend
                    for v_status in all_statuses:
                        v_color = validation_colors.get(v_status, "#94a3b8")
                        svg_parts.append(
                            f'<rect x="{legend_x}" y="{legend_y}" '
                            f'width="12" height="12" rx="2" fill="{v_color}"/>'
                        )
                        svg_parts.append(
                            f'<text x="{legend_x + 16}" y="{legend_y + 10}" '
                            f'font-size="12" fill="#374151">{v_status}</text>'
                        )
                        legend_x += len(v_status) * 8 + 32
                    # SV gain/loss legend (if SVs present)
                    if sv_scatter_data:
                        legend_x += 20
                        svg_parts.append(
                            f'<text x="{legend_x}" y="{legend_y + 10}" '
                            f'font-size="12" fill="#374151" '
                            f'font-weight="bold">SVs:</text>'
                        )
                        legend_x += 40
                        for sv_label, sv_color in [
                            ("GAIN", _SV_CALL_COLORS["GAIN"]),
                            ("LOSS", _SV_CALL_COLORS["LOSS"]),
                        ]:
                            svg_parts.append(
                                f'<rect x="{legend_x}" y="{legend_y}" '
                                f'width="20" height="12" rx="2" '
                                f'fill="{sv_color}" opacity="0.7"/>'
                            )
                            svg_parts.append(
                                f'<text x="{legend_x + 24}" '
                                f'y="{legend_y + 10}" '
                                f'font-size="12" fill="#374151">'
                                f"{sv_label}</text>"
                            )
                            legend_x += len(sv_label) * 8 + 40

                    svg_parts.append("</svg>")
                    ui.html(
                        "\n".join(svg_parts),
                        sanitize=False,
                    ).classes("w-full")

            render_stats_content()

            # Toggle ideogram/charts view
            def toggle_ideogram(_e=None):
                show_ideogram["value"] = not show_ideogram["value"]
                if _containers["charts"]:
                    _containers["charts"].set_visibility(not show_ideogram["value"])
                if _containers["ideo"]:
                    _containers["ideo"].set_visibility(show_ideogram["value"])
                if show_ideogram["value"]:
                    ideogram_btn.props(remove="outline", add="unelevated")
                else:
                    ideogram_btn.props(remove="unelevated", add="outline")
                ideogram_btn.update()

            ideogram_btn.on_click(toggle_ideogram)

            # SNV / Indel / SV filter handler
            def on_type_filter_change(_e=None):
                type_filter["snv"] = snv_cb.value
                type_filter["indel"] = indel_cb.value
                type_filter["sv"] = sv_cb.value
                render_stats_content.refresh()

            snv_cb.on_value_change(on_type_filter_change)
            indel_cb.on_value_change(on_type_filter_change)
            sv_cb.on_value_change(on_type_filter_change)

    stats_dialog.open()


def _render_outlier_section(
    filtered_all: List[Dict[str, Any]],
    validation_colors: Dict[str, str],
) -> None:
    """Render the sample outlier analysis section.

    Shows median variants/sample, IQR, outlier count, and a top-10
    bar chart with outlier samples highlighted in red.
    """
    sample_counts = Counter(
        r.get("sample", "") for r in filtered_all if r.get("sample")
    )
    if not sample_counts:
        return

    counts = sorted(sample_counts.values())
    n_samples = len(counts)
    median_val = statistics.median(counts)

    # Compute Q1 and Q3 via lower/upper halves
    lower_half = counts[: n_samples // 2]
    upper_half = counts[(n_samples + 1) // 2 :]
    q1 = statistics.median(lower_half) if lower_half else median_val
    q3 = statistics.median(upper_half) if upper_half else median_val
    iqr = q3 - q1
    upper_threshold = q3 + 1.5 * iqr
    n_outliers = sum(1 for c in counts if c > upper_threshold)

    with ui.column().classes("flex-1 min-w-[400px]"):
        ui.label("Sample Outliers").classes("text-lg font-semibold text-gray-800")
        with ui.column().classes("gap-1"):
            ui.label(f"Median: {median_val:.0f} variants/sample").classes("text-sm")
            ui.label(f"IQR: {iqr:.1f} (Q1={q1:.1f}, Q3={q3:.1f})").classes("text-sm")
            outlier_cls = "text-sm font-semibold"
            if n_outliers:
                outlier_cls += " text-red-600"
            ui.label(
                f"Outliers (>{upper_threshold:.1f}): {n_outliers} sample(s)"
            ).classes(outlier_cls)

        # Top 10 samples bar chart
        top_10 = sample_counts.most_common(10)
        if top_10:
            ui.echart(
                {
                    "tooltip": {"trigger": "axis"},
                    "grid": {
                        "top": 10,
                        "bottom": 30,
                        "left": 100,
                        "right": 20,
                    },
                    "xAxis": {"type": "value", "name": "Count"},
                    "yAxis": {
                        "type": "category",
                        "data": [s for s, _ in reversed(top_10)],
                        "axisLabel": {"fontSize": 10},
                    },
                    "series": [
                        {
                            "type": "bar",
                            "data": [
                                {
                                    "value": c,
                                    "itemStyle": {
                                        "color": (
                                            "#ef4444"
                                            if c > upper_threshold
                                            else "#3b82f6"
                                        )
                                    },
                                }
                                for _, c in reversed(top_10)
                            ],
                        }
                    ],
                }
            ).classes("w-full h-64")
