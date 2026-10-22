"""SV dialog: the SV details card and the suggested curated coordinates."""

import re
from typing import Any, Callable, Dict, List, Optional

from nicegui import ui

from genetics_viz.utils.gene_scoring import (
    GeneScorer,
    build_gene_badges,
    gene_badges_html,
    split_genes,
    tier_legend_html,
)
from genetics_viz.utils.genetrek import search_url
from genetics_viz.utils.wisecondorx import WISECONDORX_CONFIG

#: SV columns holding genes, rendered as badges (exonic ones bordered).
_GENE_COLUMNS = ("genic_symbol", "exonic_symbol", "genic_ensg", "exonic_ensg")

#: WisecondorX call kinds, in the order the legend and the colouring use.
_CALL_KINDS = ("robust_loss", "permissive_loss", "robust_gain", "permissive_gain")

#: Fields shown on the card's first row.
_PRIMARY_FIELDS = ["sample", "ratio", "zscore", "type", "call", "gene"]

_SV_KEY_PATTERN = re.compile(r"^(chr[^:]+):(\d+)-(\d+):(\w+)$")

# Priority constants for suggestion sorting
PRIORITY_PARENT = 0
PRIORITY_FAMILY = 1
PRIORITY_COHORT = 2

_PRIORITY_COLORS = {
    PRIORITY_PARENT: "green",
    PRIORITY_FAMILY: "blue",
    PRIORITY_COHORT: "grey",
}
_PRIORITY_LABELS = {
    PRIORITY_PARENT: "parent",
    PRIORITY_FAMILY: "family",
    PRIORITY_COHORT: "cohort",
}


def build_sv_suggestions(
    all_validation_map: Dict,
    chrom: str,
    start_pos: int,
    end_pos: int,
    current_sample: str,
    sample_parents: Dict[str, Optional[str]],
    family_members: List[str],
    max_suggestions: int = 5,
) -> List[Dict[str, Any]]:
    """Build a list of suggested curated coordinates from overlapping validated SVs.

    Scans all validations for SVs on the same chromosome that overlap with
    the current SV.  Returns suggestions sorted by priority (parents first,
    then family, then cohort) and overlap percentage.
    """
    _missing = {"", "0", "-9", "-", None}
    father_id = sample_parents.get("father")
    mother_id = sample_parents.get("mother")
    family_set = set(family_members)
    current_length = end_pos - start_pos

    if current_length <= 0:
        return []

    # coord_key -> suggestion dict
    merged: Dict[tuple, Dict[str, Any]] = {}

    for (variant_key, sample_id), validations in all_validation_map.items():
        m = _SV_KEY_PATTERN.match(variant_key)
        if not m:
            continue
        other_chrom, other_start_s, other_end_s, _other_type = m.groups()
        if other_chrom != chrom:
            continue
        other_start = int(other_start_s)
        other_end = int(other_end_s)

        # Filter to "present" or "different", non-ignored validations
        good = [
            v for v in validations if v[0] in ("present", "different") and v[3] != "1"
        ]
        if not good:
            continue

        # Get suggested coordinates: prefer curated, fall back to original
        # Use the most recent validation with curated coords
        curated = [v for v in good if v[4] or v[5]]
        if curated:
            curated.sort(key=lambda v: v[6], reverse=True)
            best = curated[0]
            sug_start = int(best[4]) if best[4] else other_start
            sug_end = int(best[5]) if best[5] else other_end
        else:
            sug_start = other_start
            sug_end = other_end

        # Check overlap between suggested coords and current SV
        overlap_start = max(start_pos, sug_start)
        overlap_end = min(end_pos, sug_end)
        if overlap_start >= overlap_end:
            continue
        overlap_length = overlap_end - overlap_start
        sug_length = sug_end - sug_start
        if sug_length <= 0:
            continue
        # Overlap as % of current SV and as % of suggested SV
        overlap_pct_current = overlap_length / current_length * 100
        overlap_pct_suggested = overlap_length / sug_length * 100

        # Skip exact same coordinates as current SV (not useful)
        if sug_start == start_pos and sug_end == end_pos:
            continue

        # Classify sample priority and relationship
        if sample_id == father_id and father_id not in _missing:
            priority = PRIORITY_PARENT
            relationship = "father"
        elif sample_id == mother_id and mother_id not in _missing:
            priority = PRIORITY_PARENT
            relationship = "mother"
        elif sample_id in family_set:
            priority = PRIORITY_FAMILY
            relationship = "family"
        else:
            priority = PRIORITY_COHORT
            relationship = "cohort"

        coord_key = (sug_start, sug_end)
        if coord_key not in merged:
            merged[coord_key] = {
                "start": sug_start,
                "end": sug_end,
                "overlap_pct_current": overlap_pct_current,
                "overlap_pct_suggested": overlap_pct_suggested,
                "samples": [],
                "priority": priority,
                "inheritance": "unknown",
            }
        entry = merged[coord_key]
        # Update priority to best (lowest) seen
        entry["priority"] = min(entry["priority"], priority)
        # Update overlap to max seen
        entry["overlap_pct_current"] = max(
            entry["overlap_pct_current"], overlap_pct_current
        )
        entry["overlap_pct_suggested"] = max(
            entry["overlap_pct_suggested"], overlap_pct_suggested
        )
        # Set inheritance from parent relationship
        if relationship == "father" and entry["inheritance"] == "unknown":
            entry["inheritance"] = "paternal"
        elif relationship == "mother" and entry["inheritance"] == "unknown":
            entry["inheritance"] = "maternal"
        # Add sample if not already listed
        if sample_id not in {s["id"] for s in entry["samples"]}:
            entry["samples"].append({"id": sample_id, "priority": priority})

    # Sort: by priority (parent < family < cohort), then overlap % desc
    suggestions = sorted(
        merged.values(),
        key=lambda s: (
            s["priority"],
            -min(s["overlap_pct_current"], s["overlap_pct_suggested"]),
        ),
    )
    return suggestions[:max_suggestions]


def render_suggestions(
    suggestions: List[Dict[str, Any]],
    chrom: str,
    on_apply: Callable[[Dict[str, Any]], None],
) -> None:
    """One button per suggested position; ``on_apply`` takes the suggestion."""
    if not suggestions:
        return
    with ui.row().classes("w-full mt-4 items-center gap-2 flex-wrap"):
        ui.icon("lightbulb", color="amber").classes("text-lg")
        ui.label("Suggested Coordinates").classes("text-sm font-semibold text-gray-600")
        for sug in suggestions:
            samples_text = ", ".join(
                s["id"] for s in sorted(sug["samples"], key=lambda x: x["priority"])
            )
            pct_cur = sug["overlap_pct_current"]
            pct_sug = sug["overlap_pct_suggested"]
            label = (
                f"{chrom}:{sug['start']}-{sug['end']}"
                f"  ({pct_cur:.0f}% ours | {pct_sug:.0f}% theirs)"
            )
            tooltip_text = f"{samples_text} [{_PRIORITY_LABELS[sug['priority']]}]"
            color = _PRIORITY_COLORS[sug["priority"]]
            ui.button(
                label,
                on_click=lambda _, s=sug: on_apply(s),
                icon="content_paste",
            ).props(f"dense outline color={color} size=sm no-caps").tooltip(
                tooltip_text
            )


def _call_badge_html(color: str, text: Any) -> str:
    return (
        f'<span class="q-badge" style="background-color: {color}; color: white;'
        f' padding: 2px 6px; border-radius: 4px; font-size: 12px;">{text}</span>'
    )


def _render_call(value: Any) -> None:
    """The CNV call as its configured badge, or plain text."""
    call_value = str(value) if value else "-"
    for kind in _CALL_KINDS:
        cfg = WISECONDORX_CONFIG[kind]
        if call_value == cfg["label"]:
            ui.html(_call_badge_html(cfg["color"], call_value), sanitize=False)
            return
    ui.label(call_value).classes("text-base text-gray-600 font-medium")


def _render_threshold_value(value: Any, metric: str) -> None:
    """A ratio or z-score, coloured by the call it reaches (``metric``)."""
    key = f"{metric}_threshold"
    try:
        x = float(value) if value else 0
        robust_loss, permissive_loss, robust_gain, permissive_gain = (
            WISECONDORX_CONFIG[k] for k in _CALL_KINDS
        )
        if x <= robust_loss[key]:
            ui.label(str(value)).style(
                f"color: {robust_loss['color']}; font-weight: bold"
            )
        elif x <= permissive_loss[key]:
            ui.label(str(value)).style(
                f"color: {permissive_loss['color']}; font-weight: 600"
            )
        elif x >= robust_gain[key]:
            ui.label(str(value)).style(
                f"color: {robust_gain['color']}; font-weight: bold"
            )
        elif x >= permissive_gain[key]:
            ui.label(str(value)).style(
                f"color: {permissive_gain['color']}; font-weight: 600"
            )
        else:
            ui.label(str(value)).classes("text-base text-gray-900 font-medium")
    except (TypeError, ValueError):
        ui.label(str(value) if value is not None else "-").classes(
            "text-base text-gray-900 font-medium"
        )


def _render_sv_gene_badges(gene_value_str: Any, gene_scorer: GeneScorer) -> None:
    """Render gene badges with color coding and exonic border."""
    if not gene_value_str or str(gene_value_str) in ["-", ""]:
        ui.label("-").classes("text-base text-gray-900 font-medium")
        return

    symbols, exonic = [], set()
    for item in split_genes(gene_value_str, ","):
        gene_name, _, gene_type = item.partition(":")
        symbols.append(gene_name.strip())
        if gene_type.strip() == "exonic":
            exonic.add(gene_name.strip())
    badges = build_gene_badges(symbols, gene_scorer, cap=None, bordered=exonic)
    ui.html(gene_badges_html(badges), sanitize=False)


def _render_gene_column(value: Any, exonic: bool, gene_scorer: GeneScorer) -> None:
    """A genic/exonic symbol or ENSG column as badges."""
    idents = [g for g in split_genes(value, ",") if g != "-"]
    if not idents:
        ui.label("-").classes("text-sm text-gray-800 break-all")
        return
    badges = build_gene_badges(
        idents,
        gene_scorer,
        cap=None,
        sort=False,
        bordered=idents if exonic else (),
    )
    ui.html(gene_badges_html(badges), sanitize=False)


def render_sv_details(
    sv_data: Dict[str, Any], locus: str, gene_scorer: GeneScorer
) -> None:
    """The "SV Details" heading and card: legends, call fields, gene columns."""
    ui.label("SV Details").classes("text-xl font-semibold mb-2")

    with ui.card().classes("w-full mb-4"):
        with ui.column().classes("p-4 gap-4"):
            # Legends
            with ui.row().classes("gap-4 items-center flex-wrap"):
                # Gene badge legend: the user's tiers, exonic border
                ui.html(
                    tier_legend_html(gene_scorer, exonic_border=True),
                    sanitize=False,
                )
                ui.link(
                    "Region in GeneTrek",
                    search_url([locus]),
                    new_tab=True,
                ).classes("text-xs").tooltip(
                    "Open GeneTrek's table filtered to the genes this SV overlaps"
                )

                ui.separator().props("vertical").classes("mx-2")

                # CNV call legend (from config)
                with ui.row().classes("gap-2 items-center"):
                    ui.label("CNV calls:").classes("text-xs font-semibold")
                    for kind in _CALL_KINDS:
                        cfg = WISECONDORX_CONFIG[kind]
                        ui.html(
                            _call_badge_html(cfg["color"], cfg["label"]),
                            sanitize=False,
                        )

            ui.separator()

            with ui.row().classes("gap-6 flex-wrap items-center"):
                for key in _PRIMARY_FIELDS:
                    if key not in sv_data:
                        continue
                    value = sv_data[key]
                    with ui.column().classes("gap-0"):
                        ui.label(key).classes("text-xs font-semibold text-gray-500")
                        if key == "call":
                            _render_call(value)
                        elif key in ("ratio", "zscore"):
                            _render_threshold_value(value, key)
                        elif key == "gene":
                            _render_sv_gene_badges(value, gene_scorer)
                        else:
                            ui.label(str(value) if value is not None else "-").classes(
                                "text-base text-gray-900 font-medium"
                            )

            # See more / see less section for gene annotation details
            other_fields = {
                k: v
                for k, v in sv_data.items()
                if k not in _PRIMARY_FIELDS
                and k not in ["chr:start-end"]
                and k in _GENE_COLUMNS
            }
            if not other_fields:
                return

            ui.separator()

            show_more = {"value": False}

            def toggle_more():
                show_more["value"] = not show_more["value"]
                more_button.text = "See less ▲" if show_more["value"] else "See more ▼"
                details_container.set_visibility(show_more["value"])

            more_button = (
                ui.button("See more ▼", on_click=toggle_more)
                .props("flat dense")
                .classes("text-sm text-blue-600")
            )

            with ui.column().classes("gap-2 mt-2") as details_container:
                with ui.element("div").classes("grid grid-cols-2 gap-4"):
                    for key, value in other_fields.items():
                        with ui.column().classes("gap-0"):
                            ui.label(key).classes("text-xs font-semibold text-gray-500")
                            _render_gene_column(
                                value, key.startswith("exonic"), gene_scorer
                            )

            details_container.set_visibility(False)
