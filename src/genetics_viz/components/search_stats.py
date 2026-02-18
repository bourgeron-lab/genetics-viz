"""Search page statistics dialog — chromosome distribution, consequence/validation charts, ideogram."""

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


def show_stats_dialog(rows: List[Dict[str, Any]]) -> None:
    """Create and open a variant statistics dialog for the given result rows.

    Filters out WisecondorX rows, deduplicates by variant coordinates,
    and displays chromosome distribution (stacked bar), consequence
    and validation pie charts, and an ideogram view.
    """
    # Parse Variant column into #CHROM/POS/REF/ALT if not already present
    # (validation file data uses combined Variant column: chr:pos:ref:alt)
    for r in rows:
        if "#CHROM" not in r and "Variant" in r:
            parts = str(r["Variant"]).split(":")
            if len(parts) == 4:
                r["#CHROM"], r["POS"], r["REF"], r["ALT"] = parts

    # Skip SVS rows (they lack #CHROM/POS/REF/ALT)
    snv_rows = [
        r for r in rows
        if r.get("_source_type") != "wisecondorx"
    ]

    # Deduplicate by (#CHROM, POS, REF, ALT)
    seen: set = set()
    unique_variants: list = []
    for r in snv_rows:
        key = (
            r.get("#CHROM", ""),
            r.get("POS", ""),
            r.get("REF", ""),
            r.get("ALT", ""),
        )
        if key not in seen:
            seen.add(key)
            unique_variants.append(r)

    # Classify variant type
    for r in unique_variants:
        ref = str(r.get("REF", ""))
        alt = str(r.get("ALT", ""))
        r["_is_snv"] = len(ref) == 1 and len(alt) == 1

    chrom_order = CHROM_ORDER
    chrom_sizes_mb = CHROM_SIZES_MB
    validation_colors = VALIDATION_COLORS

    # Filter state
    type_filter = {"snv": True, "indel": True}
    show_ideogram: Dict[str, bool] = {"value": False}
    _containers: Dict[str, Any] = {"charts": None, "ideo": None}

    with ui.dialog().props(
        "full-width"
    ) as stats_dialog, ui.card().classes("w-full"):
        with ui.column().classes("w-full p-4"):
            # Header
            with ui.row().classes(
                "items-center justify-between w-full mb-2"
            ):
                with ui.row().classes("items-center gap-3"):
                    ui.label("Variant Statistics").classes(
                        "text-xl font-bold text-blue-900"
                    )
                    subtitle_label = ui.label("").classes(
                        "text-sm text-gray-500"
                    )
                    ideogram_btn = ui.button(
                        "Ideogram",
                    ).props(
                        "outline color=blue size=sm dense no-caps"
                    )
                    snv_cb = ui.checkbox(
                        "SNVs", value=True
                    ).props("dense").classes("text-sm")
                    indel_cb = ui.checkbox(
                        "Indels", value=True
                    ).props("dense").classes("text-sm")
                ui.button(
                    icon="close",
                    on_click=lambda: stats_dialog.close(),
                ).props("flat round")

            @ui.refreshable
            def render_stats_content():
                # Filter variants by type
                filtered = [
                    r for r in unique_variants
                    if (type_filter["snv"] and r["_is_snv"])
                    or (type_filter["indel"] and not r["_is_snv"])
                ]
                snv_n = sum(1 for r in filtered if r["_is_snv"])
                indel_n = len(filtered) - snv_n
                subtitle_label.text = (
                    f"{len(filtered)} unique variants "
                    f"({snv_n} SNVs, {indel_n} Indels)"
                )

                # Chromosome distribution stacked by validation
                chrom_validation: Dict[str, Dict[str, int]] = {
                    c: {} for c in chrom_order
                }
                for r in filtered:
                    chrom = norm_chrom(r.get("#CHROM", ""))
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

                # Consequence distribution
                consequence_counts = Counter(
                    get_highest_consequence_term(
                        str(r.get("VEP_Consequence", ""))
                    )
                    for r in filtered
                )
                # Validation distribution
                validation_counts = Counter(
                    r.get("Validation", "") or "TODO"
                    for r in filtered
                )

                # Scatter data for ideogram
                scatter_data: List[List[Any]] = []
                for r in filtered:
                    chrom = norm_chrom(r.get("#CHROM", ""))
                    pos = r.get("POS", 0)
                    try:
                        pos_mb = round(float(pos) / 1_000_000, 2)
                    except (ValueError, TypeError):
                        continue
                    if chrom in chrom_sizes_mb:
                        status = r.get("Validation", "") or "TODO"
                        scatter_data.append([pos_mb, chrom, status])

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
                                chrom_validation[c].get(status, 0)
                                for c in chrom_order
                            ],
                            "itemStyle": {
                                "color": validation_colors.get(
                                    status, "#94a3b8"
                                )
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
                            cons_data = [
                                {
                                    "name": format_consequence_display(cons),
                                    "value": count,
                                    "itemStyle": {
                                        "color": VEP_CONSEQUENCES.get(
                                            cons, ("", "#6b7280")
                                        )[1]
                                    },
                                }
                                for cons, count in consequence_counts.most_common()
                            ]
                            ui.echart(
                                {
                                    "tooltip": {"trigger": "item"},
                                    "series": [{
                                        "type": "pie",
                                        "radius": "70%",
                                        "data": cons_data,
                                        "label": {"formatter": "{b}: {c} ({d}%)"},
                                    }],
                                }
                            ).classes("w-full h-80")

                        # Validation pie chart
                        with ui.column().classes("flex-1 min-w-[400px]"):
                            ui.label(
                                "Validation Status Distribution"
                            ).classes("text-lg font-semibold text-gray-800")
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
                                    "series": [{
                                        "type": "pie",
                                        "radius": "70%",
                                        "data": val_data,
                                        "label": {"formatter": "{b}: {c} ({d}%)"},
                                    }],
                                }
                            ).classes("w-full h-80")

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
                            f'{vx + tw:.1f},{bar_y - tri_h} '
                            f'{vx:.1f},{bar_y}" '
                            f'fill="{v_color}" opacity="0.9"/>'
                        )
                        svg_parts.append(
                            f'<line x1="{vx:.1f}" y1="{bar_y}" '
                            f'x2="{vx:.1f}" y2="{bar_y + row_h}" '
                            f'stroke="{v_color}" stroke-width="1.5" '
                            f'opacity="0.85"/>'
                        )

                    legend_y = axis_y + 40
                    legend_x = lbl_w
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
                    _containers["charts"].set_visibility(
                        not show_ideogram["value"]
                    )
                if _containers["ideo"]:
                    _containers["ideo"].set_visibility(
                        show_ideogram["value"]
                    )
                if show_ideogram["value"]:
                    ideogram_btn.props(
                        remove="outline", add="unelevated"
                    )
                else:
                    ideogram_btn.props(
                        remove="unelevated", add="outline"
                    )
                ideogram_btn.update()

            ideogram_btn.on_click(toggle_ideogram)

            # SNV / Indel filter handler
            def on_type_filter_change(_e=None):
                type_filter["snv"] = snv_cb.value
                type_filter["indel"] = indel_cb.value
                render_stats_content.refresh()

            snv_cb.on_value_change(on_type_filter_change)
            indel_cb.on_value_change(on_type_filter_change)

    stats_dialog.open()
