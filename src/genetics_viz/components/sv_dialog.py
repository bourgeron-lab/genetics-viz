"""Shared component for SV visualization in a dialog with IGV.js."""

import csv
import fcntl
from genetics_viz.utils.auth import can_write, get_current_user
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml
from nicegui import ui

from genetics_viz.utils.data import get_data_store, get_static_prefix
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.sharding import get_sample_path, get_sample_url
from genetics_viz.utils.wisecondorx import infer_sv_type


# Load WisecondorX thresholds and colors from YAML
def _load_wisecondorx_config():
    config_path = (
        Path(__file__).parent.parent / "config" / "wisecondorx_thresholds.yaml"
    )
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


WISECONDORX_CONFIG = _load_wisecondorx_config()


_SV_KEY_PATTERN = re.compile(r"^(chr[^:]+):(\d+)-(\d+):(\w+)$")

# Priority constants for suggestion sorting
_PRIORITY_PARENT = 0
_PRIORITY_FAMILY = 1
_PRIORITY_COHORT = 2


def _build_sv_suggestions(
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
    parent_ids = {
        pid
        for pid in (sample_parents.get("father"), sample_parents.get("mother"))
        if pid and pid not in ("", "0", "-9", "-")
    }
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
        overlap_pct = overlap_length / min(current_length, sug_length) * 100

        # Skip exact same coordinates as current SV (not useful)
        if sug_start == start_pos and sug_end == end_pos:
            continue

        # Classify sample priority
        if sample_id in parent_ids:
            priority = _PRIORITY_PARENT
        elif sample_id in family_set:
            priority = _PRIORITY_FAMILY
        else:
            priority = _PRIORITY_COHORT

        coord_key = (sug_start, sug_end)
        if coord_key not in merged:
            merged[coord_key] = {
                "start": sug_start,
                "end": sug_end,
                "overlap_pct": overlap_pct,
                "samples": [],
                "priority": priority,
            }
        entry = merged[coord_key]
        # Update priority to best (lowest) seen
        entry["priority"] = min(entry["priority"], priority)
        # Update overlap to max seen
        entry["overlap_pct"] = max(entry["overlap_pct"], overlap_pct)
        # Add sample if not already listed
        if sample_id not in {s["id"] for s in entry["samples"]}:
            entry["samples"].append({"id": sample_id, "priority": priority})

    # Sort: by priority (parent < family < cohort), then overlap % desc
    suggestions = sorted(
        merged.values(), key=lambda s: (s["priority"], -s["overlap_pct"])
    )
    return suggestions[:max_suggestions]


def show_sv_dialog(
    cohort_name: str,
    family_id: str,
    chrom: str,
    start: str,
    end: str,
    sample: str,
    sv_data: Dict[str, Any],
    on_validation_saved: Optional[callable] = None,
) -> None:
    """Show SV visualization dialog with IGV viewer.

    Args:
        cohort_name: Cohort name
        family_id: Family ID
        chrom: Chromosome
        start: Start position
        end: End position
        sample: Sample ID
        sv_data: Additional SV data to display
        on_validation_saved: Optional callback to invoke after validation is saved
    """
    store = get_data_store()

    # Calculate expanded locus (x1.6 centered on original)
    start_pos = int(start)
    end_pos = int(end)
    center = (start_pos + end_pos) // 2
    original_length = end_pos - start_pos
    expanded_length = int(original_length * 1.6)
    expanded_start = center - (expanded_length // 2)
    expanded_end = center + (expanded_length // 2)

    # Ensure non-negative positions
    expanded_start = max(0, expanded_start)

    locus = f"{chrom}:{start}-{end}"
    expanded_locus = f"{chrom}:{expanded_start}-{expanded_end}"

    # Check for curated boundaries from validations
    from genetics_viz.components.validation_loader import load_validation_map

    validation_file = store.data_dir / "validations" / "svs.tsv"
    validation_map = load_validation_map(validation_file, family_id)
    all_validation_map = load_validation_map(validation_file, family_id=None)

    # Determine SV type for variant key
    sv_type = infer_sv_type(sv_data)
    variant_key = f"{chrom}:{start}-{end}:{sv_type}"
    map_key = (variant_key, sample)

    # Look for present validations with curated boundaries
    # Keep original coordinates for title and validation lookups
    # Use curated coordinates only for ROI display and IGV locus
    curated_start = None
    curated_end = None
    curated_locus = None
    curated_expanded_locus = None
    if map_key in validation_map:
        validations = validation_map[map_key]
        present_with_curated = [
            v
            for v in validations
            if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
        ]
        if present_with_curated:
            # Sort by timestamp (most recent first)
            present_with_curated.sort(key=lambda v: v[6], reverse=True)
            most_recent = present_with_curated[0]
            if most_recent[4]:  # CuratedStart
                curated_start = most_recent[4]
            if most_recent[5]:  # CuratedEnd
                curated_end = most_recent[5]

            # Build curated locus string if any curated boundaries exist
            if curated_start or curated_end:
                display_start = curated_start if curated_start else start
                display_end = curated_end if curated_end else end
                curated_locus = f"{chrom}:{display_start}-{display_end}"

                # Calculate expanded locus based on curated boundaries
                curated_start_pos = int(display_start)
                curated_end_pos = int(display_end)
                curated_center = (curated_start_pos + curated_end_pos) // 2
                curated_length = curated_end_pos - curated_start_pos
                curated_expanded_length = int(curated_length * 1.6)
                curated_expanded_start = curated_center - (curated_expanded_length // 2)
                curated_expanded_end = curated_center + (curated_expanded_length // 2)
                curated_expanded_start = max(0, curated_expanded_start)
                curated_expanded_locus = (
                    f"{chrom}:{curated_expanded_start}-{curated_expanded_end}"
                )

    # Use curated coordinates for IGV display if available
    igv_locus = curated_expanded_locus if curated_expanded_locus else expanded_locus
    roi_start_coord = curated_start if curated_start else start
    roi_end_coord = curated_end if curated_end else end

    with (
        ui.dialog().props("maximized") as dialog,
        ui.card().classes("w-full h-full"),
    ):
        with ui.column().classes("w-full h-full p-6"):
            # Header with close button
            with ui.row().classes("items-center justify-between w-full mb-4"):
                with ui.row().classes("items-center gap-2"):
                    ui.label(f"🧬 SV: {locus} - {sample}").classes(
                        "text-2xl font-bold text-blue-900"
                    )
                    if curated_locus:
                        ui.label(f"(curated: {curated_locus})").classes(
                            "text-lg font-semibold text-green-600"
                        )
                ui.button(icon="close", on_click=lambda: dialog.close()).props(
                    "flat round"
                )

            # Get family members
            cohort = store.get_cohort(cohort_name)
            family_members: List[str] = []
            sample_parents: Dict[str, Optional[str]] = {
                "father": None,
                "mother": None,
            }

            if cohort:
                members_data = cohort.get_family_members(family_id)
                family_members = [m["Sample ID"] for m in members_data]
                # Find current sample's parents
                for member in members_data:
                    if member["Sample ID"] == sample:
                        sample_parents["father"] = member.get("Father")
                        sample_parents["mother"] = member.get("Mother")
                        break

            # Track additional samples
            additional_samples: Dict[str, List[str]] = {"value": []}

            # Automatically add parents if available
            for parent_type, parent_id in sample_parents.items():
                if (
                    parent_id
                    and parent_id != "-"
                    and parent_id != "0"
                    and parent_id != sample
                ):
                    bedgraph_file = (
                        get_sample_path(store.data_dir, parent_id)
                        / "sequences"
                        / f"{parent_id}.by1000.bedgraph.gz"
                    )
                    if bedgraph_file.exists():
                        additional_samples["value"].append(parent_id)

            def get_relationship_label(sample_id: str) -> str:
                """Get relationship label for a sample."""
                if sample_id == sample:
                    return "(carrier)"

                if sample_id == sample_parents["father"]:
                    return "(father)"
                if sample_id == sample_parents["mother"]:
                    return "(mother)"

                # Check if it's a sibling
                if cohort and sample_parents["father"] and sample_parents["mother"]:
                    members_data = cohort.get_family_members(family_id)
                    for member in members_data:
                        if member["Sample ID"] == sample_id:
                            member_father = member.get("Father")
                            member_mother = member.get("Mother")
                            if (
                                member_father == sample_parents["father"]
                                and member_mother == sample_parents["mother"]
                                and member_father
                                and member_mother
                                and member_father != "-"
                                and member_father != "0"
                                and member_mother != "-"
                                and member_mother != "0"
                            ):
                                return "(sibling)"
                            break

                return ""

            # Additional samples section with add menu
            with ui.row().classes("items-center gap-4 mb-2"):
                # Menu to add samples
                with ui.button("Add Samples", icon="add").props(
                    "outline color=blue size=sm"
                ):
                    with ui.menu():
                        ui.menu_item("Add Parents", on_click=lambda: add_parents())
                        ui.menu_item("Add Family", on_click=lambda: add_family())
                        ui.separator()
                        with ui.row().classes("items-center gap-2 px-4 py-2"):
                            barcode_input = (
                                ui.input("Barcode").classes("flex-grow").props("dense")
                            )
                            ui.button(
                                "Add",
                                icon="add",
                                on_click=lambda: add_sample(barcode_input.value),
                            ).props("flat dense size=sm")

            # Display additional samples
            additional_samples_container = ui.column().classes("gap-1 mb-4")

            def refresh_additional_samples():
                additional_samples_container.clear()
                with additional_samples_container:
                    if additional_samples["value"]:
                        ui.label("Additional Samples:").classes(
                            "text-sm font-semibold text-gray-700"
                        )
                        for add_sample_id in additional_samples["value"]:
                            with ui.row().classes("items-center gap-2"):
                                label_text = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                                ui.label(label_text).classes("text-sm text-gray-600")

                                def make_remove_handler(sid: str):
                                    return lambda: remove_sample(sid)

                                ui.button(
                                    icon="delete",
                                    on_click=make_remove_handler(add_sample_id),
                                ).props("flat dense size=xs color=red")

            def add_sample(sample_id: str):
                if (
                    sample_id
                    and sample_id not in additional_samples["value"]
                    and sample_id != sample
                ):
                    bedgraph_file = (
                        get_sample_path(store.data_dir, sample_id)
                        / "sequences"
                        / f"{sample_id}.by1000.bedgraph.gz"
                    )
                    if bedgraph_file.exists():
                        additional_samples["value"].append(sample_id)
                        refresh_additional_samples()
                        refresh_igv()
                        refresh_igv_vaf()
                        refresh_igv_cram()
                    else:
                        ui.notify(
                            f"Bedgraph file not found for sample: {sample_id}",
                            type="warning",
                        )

            def add_parents():
                added = []
                for parent_type, parent_id in sample_parents.items():
                    if parent_id and parent_id != "-" and parent_id != "0":
                        if (
                            parent_id not in additional_samples["value"]
                            and parent_id != sample
                        ):
                            bedgraph_file = (
                                get_sample_path(store.data_dir, parent_id)
                                / "sequences"
                                / f"{parent_id}.by1000.bedgraph.gz"
                            )
                            if bedgraph_file.exists():
                                additional_samples["value"].append(parent_id)
                                added.append(parent_id)
                if added:
                    refresh_additional_samples()
                    refresh_igv()
                    refresh_igv_vaf()
                    refresh_igv_cram()
                    ui.notify(f"Added parents: {', '.join(added)}", type="positive")
                else:
                    ui.notify("No parents to add or files not found", type="warning")

            def add_family():
                added = []
                for member_id in family_members:
                    if (
                        member_id not in additional_samples["value"]
                        and member_id != sample
                    ):
                        bedgraph_file = (
                            get_sample_path(store.data_dir, member_id)
                            / "sequences"
                            / f"{member_id}.by1000.bedgraph.gz"
                        )
                        if bedgraph_file.exists():
                            additional_samples["value"].append(member_id)
                            added.append(member_id)
                if added:
                    refresh_additional_samples()
                    refresh_igv()
                    refresh_igv_vaf()
                    refresh_igv_cram()
                    ui.notify(f"Added {len(added)} family members", type="positive")
                else:
                    ui.notify("No additional family members to add", type="warning")

            def remove_sample(sample_id: str):
                if sample_id in additional_samples["value"]:
                    additional_samples["value"].remove(sample_id)
                    refresh_additional_samples()
                    refresh_igv()
                    refresh_igv_vaf()
                    refresh_igv_cram()

            refresh_additional_samples()

            # Helper function to render gene badges with color coding
            gene_scorer = get_gene_scorer()

            def render_sv_gene_badges(gene_value_str: str):
                """Render gene badges with color coding and exonic border."""
                if not gene_value_str or str(gene_value_str) in ["-", ""]:
                    return ui.label("-").classes("text-base text-gray-900 font-medium")

                with ui.row().classes("gap-1 flex-wrap"):
                    for item in str(gene_value_str).split(","):
                        if ":" in item:
                            gene_name, gene_type = item.split(":", 1)
                            gene_name = gene_name.strip()
                            is_exonic = gene_type.strip() == "exonic"
                        else:
                            gene_name = item.strip()
                            is_exonic = False

                        color = gene_scorer.get_gene_color(gene_name)
                        text_color = "black" if color == "#ffffff" else "white"
                        tooltip_text = gene_scorer.get_gene_tooltip(gene_name)
                        border_style = "border: 2px solid black;" if is_exonic else ""

                        ui.html(
                            f'<span class="q-badge" style="background-color: {color}; color: {text_color}; {border_style} padding: 2px 6px; border-radius: 4px; font-size: 12px;" title="{tooltip_text}">{gene_name}</span>',
                            sanitize=False,
                        )

            # SV details card
            ui.label("SV Details").classes("text-xl font-semibold mb-2")

            with ui.card().classes("w-full mb-4"):
                with ui.column().classes("p-4 gap-4"):
                    # Legends
                    with ui.row().classes("gap-4 items-center flex-wrap"):
                        # Gene badge legend
                        with ui.row().classes("gap-2 items-center"):
                            ui.label("Gene badges:").classes("text-xs font-semibold")
                            ui.html(
                                '<span class="q-badge" style="background-color: #ffffff; color: black; border: 2px solid black; padding: 2px 6px; border-radius: 4px; font-size: 10px;">Exonic</span>',
                                sanitize=False,
                            )
                            ui.label("(black border)").classes("text-xs text-gray-500")
                            ui.html(
                                '<span class="q-badge" style="background-color: #ffffff; color: black; padding: 2px 6px; border-radius: 4px; font-size: 10px;">Genic</span>',
                                sanitize=False,
                            )
                            ui.label("(no border)").classes("text-xs text-gray-500")
                            ui.html(
                                '<span class="q-badge" style="background-color: #8b0000; color: white; padding: 2px 6px; border-radius: 4px; font-size: 10px;">Color</span>',
                                sanitize=False,
                            )
                            ui.label("(geneset score)").classes("text-xs text-gray-500")

                        ui.separator().props("vertical").classes("mx-2")

                        # CNV call legend (from config)
                        robust_loss = WISECONDORX_CONFIG["robust_loss"]
                        permissive_loss = WISECONDORX_CONFIG["permissive_loss"]
                        robust_gain = WISECONDORX_CONFIG["robust_gain"]
                        permissive_gain = WISECONDORX_CONFIG["permissive_gain"]

                        with ui.row().classes("gap-2 items-center"):
                            ui.label("CNV calls:").classes("text-xs font-semibold")
                            ui.html(
                                f'<span class="q-badge" style="background-color: {robust_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{robust_loss["label"]}</span>',
                                sanitize=False,
                            )
                            ui.html(
                                f'<span class="q-badge" style="background-color: {permissive_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{permissive_loss["label"]}</span>',
                                sanitize=False,
                            )
                            ui.html(
                                f'<span class="q-badge" style="background-color: {robust_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{robust_gain["label"]}</span>',
                                sanitize=False,
                            )
                            ui.html(
                                f'<span class="q-badge" style="background-color: {permissive_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{permissive_gain["label"]}</span>',
                                sanitize=False,
                            )

                    ui.separator()

                    # Primary fields with badge styling
                    primary_fields = [
                        "sample",
                        "ratio",
                        "zscore",
                        "type",
                        "call",
                        "gene",
                    ]

                    with ui.row().classes("gap-6 flex-wrap items-center"):
                        for key in primary_fields:
                            if key in sv_data:
                                value = sv_data[key]
                                with ui.column().classes("gap-0"):
                                    ui.label(key).classes(
                                        "text-xs font-semibold text-gray-500"
                                    )

                                    # Special rendering for specific fields
                                    if key == "call":
                                        # CNV call badge (from config)
                                        call_value = str(value) if value else "-"
                                        robust_loss = WISECONDORX_CONFIG["robust_loss"]
                                        permissive_loss = WISECONDORX_CONFIG[
                                            "permissive_loss"
                                        ]
                                        robust_gain = WISECONDORX_CONFIG["robust_gain"]
                                        permissive_gain = WISECONDORX_CONFIG[
                                            "permissive_gain"
                                        ]

                                        if call_value == robust_loss["label"]:
                                            ui.html(
                                                f'<span class="q-badge" style="background-color: {robust_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{call_value}</span>',
                                                sanitize=False,
                                            )
                                        elif call_value == permissive_loss["label"]:
                                            ui.html(
                                                f'<span class="q-badge" style="background-color: {permissive_loss["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{call_value}</span>',
                                                sanitize=False,
                                            )
                                        elif call_value == robust_gain["label"]:
                                            ui.html(
                                                f'<span class="q-badge" style="background-color: {robust_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{call_value}</span>',
                                                sanitize=False,
                                            )
                                        elif call_value == permissive_gain["label"]:
                                            ui.html(
                                                f'<span class="q-badge" style="background-color: {permissive_gain["color"]}; color: white; padding: 2px 6px; border-radius: 4px; font-size: 12px;">{call_value}</span>',
                                                sanitize=False,
                                            )
                                        else:
                                            ui.label(call_value).classes(
                                                "text-base text-gray-600 font-medium"
                                            )

                                    elif key == "ratio":
                                        # Colored ratio value (from config)
                                        try:
                                            r = float(value) if value else 0
                                            robust_loss = WISECONDORX_CONFIG[
                                                "robust_loss"
                                            ]
                                            permissive_loss = WISECONDORX_CONFIG[
                                                "permissive_loss"
                                            ]
                                            robust_gain = WISECONDORX_CONFIG[
                                                "robust_gain"
                                            ]
                                            permissive_gain = WISECONDORX_CONFIG[
                                                "permissive_gain"
                                            ]

                                            if r <= robust_loss["ratio_threshold"]:
                                                ui.label(str(value)).style(
                                                    f"color: {robust_loss['color']}; font-weight: bold"
                                                )
                                            elif (
                                                r <= permissive_loss["ratio_threshold"]
                                            ):
                                                ui.label(str(value)).style(
                                                    f"color: {permissive_loss['color']}; font-weight: 600"
                                                )
                                            elif r >= robust_gain["ratio_threshold"]:
                                                ui.label(str(value)).style(
                                                    f"color: {robust_gain['color']}; font-weight: bold"
                                                )
                                            elif (
                                                r >= permissive_gain["ratio_threshold"]
                                            ):
                                                ui.label(str(value)).style(
                                                    f"color: {permissive_gain['color']}; font-weight: 600"
                                                )
                                            else:
                                                ui.label(str(value)).classes(
                                                    "text-base text-gray-900 font-medium"
                                                )
                                        except:
                                            ui.label(
                                                str(value) if value is not None else "-"
                                            ).classes(
                                                "text-base text-gray-900 font-medium"
                                            )

                                    elif key == "zscore":
                                        # Colored zscore value (from config)
                                        try:
                                            z = float(value) if value else 0
                                            robust_loss = WISECONDORX_CONFIG[
                                                "robust_loss"
                                            ]
                                            permissive_loss = WISECONDORX_CONFIG[
                                                "permissive_loss"
                                            ]
                                            robust_gain = WISECONDORX_CONFIG[
                                                "robust_gain"
                                            ]
                                            permissive_gain = WISECONDORX_CONFIG[
                                                "permissive_gain"
                                            ]

                                            if z <= robust_loss["zscore_threshold"]:
                                                ui.label(str(value)).style(
                                                    f"color: {robust_loss['color']}; font-weight: bold"
                                                )
                                            elif (
                                                z <= permissive_loss["zscore_threshold"]
                                            ):
                                                ui.label(str(value)).style(
                                                    f"color: {permissive_loss['color']}; font-weight: 600"
                                                )
                                            elif z >= robust_gain["zscore_threshold"]:
                                                ui.label(str(value)).style(
                                                    f"color: {robust_gain['color']}; font-weight: bold"
                                                )
                                            elif (
                                                z >= permissive_gain["zscore_threshold"]
                                            ):
                                                ui.label(str(value)).style(
                                                    f"color: {permissive_gain['color']}; font-weight: 600"
                                                )
                                            else:
                                                ui.label(str(value)).classes(
                                                    "text-base text-gray-900 font-medium"
                                                )
                                        except:
                                            ui.label(
                                                str(value) if value is not None else "-"
                                            ).classes(
                                                "text-base text-gray-900 font-medium"
                                            )

                                    elif key == "gene":
                                        # Gene badges with color coding and border for exonic
                                        render_sv_gene_badges(value)

                                    else:
                                        # Default rendering
                                        ui.label(
                                            str(value) if value is not None else "-"
                                        ).classes("text-base text-gray-900 font-medium")

                    # See more / see less section for gene annotation details
                    gene_annotation_fields = [
                        "genic_symbol",
                        "exonic_symbol",
                        "genic_ensg",
                        "exonic_ensg",
                    ]
                    other_fields = {
                        k: v
                        for k, v in sv_data.items()
                        if k not in primary_fields
                        and k not in ["chr:start-end"]
                        and k in gene_annotation_fields
                    }

                    if other_fields:
                        ui.separator()

                        show_more = {"value": False}

                        def toggle_more():
                            show_more["value"] = not show_more["value"]
                            more_button.text = (
                                "See less ▲" if show_more["value"] else "See more ▼"
                            )
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
                                        ui.label(key).classes(
                                            "text-xs font-semibold text-gray-500"
                                        )
                                        # Render gene/symbol fields with color-coded badges
                                        if key in ["genic_symbol", "exonic_symbol"]:
                                            if value and str(value) not in ["-", ""]:
                                                is_exonic = "exonic" in key
                                                with ui.row().classes(
                                                    "gap-1 flex-wrap"
                                                ):
                                                    for item in str(value).split(","):
                                                        gene_name = item.strip()
                                                        color = (
                                                            gene_scorer.get_gene_color(
                                                                gene_name
                                                            )
                                                        )
                                                        text_color = (
                                                            "black"
                                                            if color == "#ffffff"
                                                            else "white"
                                                        )
                                                        tooltip_text = gene_scorer.get_gene_tooltip(
                                                            gene_name
                                                        )
                                                        border_style = (
                                                            "border: 2px solid black;"
                                                            if is_exonic
                                                            else ""
                                                        )
                                                        ui.html(
                                                            f'<span class="q-badge" style="background-color: {color}; color: {text_color}; {border_style} padding: 2px 6px; border-radius: 4px; font-size: 12px;" title="{tooltip_text}">{gene_name}</span>',
                                                            sanitize=False,
                                                        )
                                            else:
                                                ui.label("-").classes(
                                                    "text-sm text-gray-800 break-all"
                                                )
                                        elif key in ["genic_ensg", "exonic_ensg"]:
                                            # Render ENSG IDs with color-coded badges
                                            if value and str(value) not in ["-", ""]:
                                                is_exonic = "exonic" in key
                                                with ui.row().classes(
                                                    "gap-1 flex-wrap"
                                                ):
                                                    for item in str(value).split(","):
                                                        ensg_id = item.strip()
                                                        color = (
                                                            gene_scorer.get_gene_color(
                                                                ensg_id
                                                            )
                                                        )
                                                        text_color = (
                                                            "black"
                                                            if color == "#ffffff"
                                                            else "white"
                                                        )
                                                        tooltip_text = gene_scorer.get_gene_tooltip(
                                                            ensg_id
                                                        )
                                                        border_style = (
                                                            "border: 2px solid black;"
                                                            if is_exonic
                                                            else ""
                                                        )
                                                        ui.html(
                                                            f'<span class="q-badge" style="background-color: {color}; color: {text_color}; {border_style} padding: 2px 6px; border-radius: 4px; font-size: 12px;" title="{tooltip_text}">{ensg_id}</span>',
                                                            sanitize=False,
                                                        )
                                            else:
                                                ui.label("-").classes(
                                                    "text-sm text-gray-800 break-all"
                                                )
                                        else:
                                            # Default rendering for other fields
                                            ui.label(
                                                str(value) if value is not None else "-"
                                            ).classes("text-sm text-gray-800 break-all")

                        details_container.set_visibility(False)

            # Calculate dynamic height for IGV containers
            def calculate_igv_height():
                """Calculate height based on number of tracks: 200px base + (tracks * 100px) + 50px buffer."""
                num_tracks = 1 + len(additional_samples["value"])
                return 200 + (num_tracks * 100) + 50

            def calculate_vaf_height():
                """Calculate height based on number of tracks: 200px base + (tracks * 100px) + 50px buffer."""
                num_tracks = 1 + len(additional_samples["value"])
                return 200 + (num_tracks * 100) + 50

            def calculate_cram_height():
                """Calculate height for CRAM tracks: 200px base + (tracks * 250px) + 50px buffer."""
                num_tracks = 1 + len(additional_samples["value"])
                return 200 + (num_tracks * 250) + 50

            # Build IGV tracks function for bedgraph
            def build_igv_tracks():
                tracks = []

                # Main sample track
                bedgraph_file = (
                    get_sample_path(store.data_dir, sample)
                    / "sequences"
                    / f"{sample}.by1000.bedgraph.gz"
                )
                if bedgraph_file.exists():
                    sample_url_seg = get_sample_url(store.data_dir, sample)
                    main_label = f"{sample} {get_relationship_label(sample)}".strip()
                    tracks.append(
                        {
                            "name": main_label,
                            "type": "wig",
                            "format": "bedgraph",
                            "url": f"{get_static_prefix()}/{sample_url_seg}/sequences/{sample}.by1000.bedgraph.gz",
                            "indexURL": f"{get_static_prefix()}/{sample_url_seg}/sequences/{sample}.by1000.bedgraph.gz.tbi",
                            "height": 100,
                            "autoscaleGroup": "cnv",
                        }
                    )

                # Additional samples tracks
                for add_sample_id in additional_samples["value"]:
                    bedgraph_file = (
                        get_sample_path(store.data_dir, add_sample_id)
                        / "sequences"
                        / f"{add_sample_id}.by1000.bedgraph.gz"
                    )
                    if bedgraph_file.exists():
                        add_url_seg = get_sample_url(store.data_dir, add_sample_id)
                        add_label = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                        tracks.append(
                            {
                                "name": add_label,
                                "type": "wig",
                                "format": "bedgraph",
                                "url": f"{get_static_prefix()}/{add_url_seg}/sequences/{add_sample_id}.by1000.bedgraph.gz",
                                "indexURL": f"{get_static_prefix()}/{add_url_seg}/sequences/{add_sample_id}.by1000.bedgraph.gz.tbi",
                                "height": 100,
                                "autoscaleGroup": "cnv",
                            }
                        )

                return tracks

            # Build IGV tracks function for VAF bedgraph
            def build_vaf_tracks():
                tracks = []

                # Main sample track
                vaf_file = (
                    get_sample_path(store.data_dir, sample)
                    / "sequences"
                    / f"{sample}.vaf.bedgraph.gz"
                )
                if vaf_file.exists():
                    sample_url_seg = get_sample_url(store.data_dir, sample)
                    main_label = f"{sample} {get_relationship_label(sample)}".strip()
                    tracks.append(
                        {
                            "name": main_label,
                            "type": "wig",
                            "format": "bedgraph",
                            "graphType": "points",
                            "url": f"{get_static_prefix()}/{sample_url_seg}/sequences/{sample}.vaf.bedgraph.gz",
                            "indexURL": f"{get_static_prefix()}/{sample_url_seg}/sequences/{sample}.vaf.bedgraph.gz.tbi",
                            "height": 100,
                            "min": 0,
                            "max": 1,
                            "autoscale": False,
                        }
                    )

                # Additional samples tracks
                for add_sample_id in additional_samples["value"]:
                    vaf_file = (
                        get_sample_path(store.data_dir, add_sample_id)
                        / "sequences"
                        / f"{add_sample_id}.vaf.bedgraph.gz"
                    )
                    if vaf_file.exists():
                        add_url_seg = get_sample_url(store.data_dir, add_sample_id)
                        add_label = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                        tracks.append(
                            {
                                "name": add_label,
                                "type": "wig",
                                "format": "bedgraph",
                                "graphType": "points",
                                "url": f"{get_static_prefix()}/{add_url_seg}/sequences/{add_sample_id}.vaf.bedgraph.gz",
                                "indexURL": f"{get_static_prefix()}/{add_url_seg}/sequences/{add_sample_id}.vaf.bedgraph.gz.tbi",
                                "height": 100,
                                "min": 0,
                                "max": 1,
                                "autoscale": False,
                            }
                        )

                return tracks

            # Build CRAM tracks function for read-level view
            def build_cram_tracks():
                tracks = []

                # Main sample CRAM track
                cram_file = (
                    get_sample_path(store.data_dir, sample)
                    / "sequences"
                    / f"{sample}.GRCh38_GIABv3.cram"
                )
                if cram_file.exists():
                    sample_url_seg = get_sample_url(store.data_dir, sample)
                    main_label = f"{sample} {get_relationship_label(sample)}".strip()
                    tracks.append(
                        {
                            "name": main_label,
                            "type": "alignment",
                            "format": "cram",
                            "url": f"{get_static_prefix()}/{sample_url_seg}/sequences/{sample}.GRCh38_GIABv3.cram",
                            "indexURL": f"{get_static_prefix()}/{sample_url_seg}/sequences/{sample}.GRCh38_GIABv3.cram.crai",
                            "height": 250,
                            "displayMode": "SQUISHED",
                            "viewAsPairs": True,
                            "autoscale": False,
                            "autoscaleGroup": "cram",
                        }
                    )

                # Additional samples CRAM tracks
                for add_sample_id in additional_samples["value"]:
                    cram_file = (
                        get_sample_path(store.data_dir, add_sample_id)
                        / "sequences"
                        / f"{add_sample_id}.GRCh38_GIABv3.cram"
                    )
                    if cram_file.exists():
                        add_url_seg = get_sample_url(store.data_dir, add_sample_id)
                        add_label = f"{add_sample_id} {get_relationship_label(add_sample_id)}".strip()
                        tracks.append(
                            {
                                "name": add_label,
                                "type": "alignment",
                                "format": "cram",
                                "url": f"{get_static_prefix()}/{add_url_seg}/sequences/{add_sample_id}.GRCh38_GIABv3.cram",
                                "indexURL": f"{get_static_prefix()}/{add_url_seg}/sequences/{add_sample_id}.GRCh38_GIABv3.cram.crai",
                                "height": 250,
                                "displayMode": "SQUISHED",
                                "viewAsPairs": True,
                                "autoscale": False,
                                "autoscaleGroup": "cram",
                            }
                        )

                return tracks

            # Suggested coordinates from overlapping validated SVs
            suggestions = _build_sv_suggestions(
                all_validation_map=all_validation_map,
                chrom=chrom,
                start_pos=start_pos,
                end_pos=end_pos,
                current_sample=sample,
                sample_parents=sample_parents,
                family_members=family_members,
            )

            if suggestions:
                _PRIORITY_COLORS = {
                    _PRIORITY_PARENT: "green",
                    _PRIORITY_FAMILY: "blue",
                    _PRIORITY_COHORT: "grey",
                }
                _PRIORITY_LABELS = {
                    _PRIORITY_PARENT: "parent",
                    _PRIORITY_FAMILY: "family",
                    _PRIORITY_COHORT: "cohort",
                }

                with ui.row().classes("w-full mt-4 items-center gap-2 flex-wrap"):
                    ui.icon("lightbulb", color="amber").classes("text-lg")
                    ui.label("Suggested Coordinates").classes(
                        "text-sm font-semibold text-gray-600"
                    )
                    for sug in suggestions:
                        sample_labels = []
                        for s in sorted(sug["samples"], key=lambda x: x["priority"]):
                            sample_labels.append(s["id"])
                        samples_text = ", ".join(sample_labels)
                        pct = sug["overlap_pct"]
                        label = (
                            f"{chrom}:{sug['start']}-{sug['end']}  ({pct:.0f}% overlap)"
                        )
                        tooltip_text = (
                            f"{samples_text} [{_PRIORITY_LABELS[sug['priority']]}]"
                        )
                        color = _PRIORITY_COLORS[sug["priority"]]

                        def make_apply(s):
                            def on_click():
                                curated_inputs["start"].value = str(s["start"])
                                curated_inputs["end"].value = str(s["end"])
                                roi_coords["start"] = s["start"]
                                roi_coords["end"] = s["end"]
                                _update_roi()
                                ui.notify(
                                    f"Curated position set to"
                                    f" {chrom}:{s['start']}-{s['end']}",
                                    type="positive",
                                )

                            return on_click

                        ui.button(
                            label,
                            on_click=make_apply(sug),
                            icon="content_paste",
                        ).props(f"dense outline color={color} size=sm no-caps").tooltip(
                            tooltip_text
                        )

            # IGV.js viewer container for bedgraph (CNV view)
            cnv_expansion = (
                ui.expansion("CNV Coverage View", icon="analytics", value=True)
                .classes("w-full mt-4 border border-gray-400 rounded-lg")
                .props("header-class='text-lg font-bold bg-gray-100'")
            )

            with cnv_expansion:
                igv_container = (
                    ui.element("div")
                    .props('id="igv-div"')
                    .classes("w-full border border-gray-300 rounded-lg")
                    .style(f"height: {calculate_igv_height()}px")
                )

                @ui.refreshable
                def refresh_igv():
                    """Refresh IGV viewer with current tracks."""
                    igv_container.clear()

                    # Update container height
                    try:
                        new_height = calculate_igv_height()
                        igv_container.style(f"height: {new_height}px")
                    except RuntimeError:
                        pass

                    tracks = build_igv_tracks()

                    # Determine ROI color based on call type
                    call_value = sv_data.get("call", "")
                    roi_color = "rgba(255,0,0,0.2)"  # Default red for loss
                    if "GAIN" in str(call_value).upper():
                        roi_color = "rgba(0,128,0,0.2)"  # Green for gain

                    # Create ROI description
                    roi_description = f"{chrom}:{start}-{end} - {call_value}"

                    # IGV.js initialization script
                    igv_script = f"""
                    (async function() {{
                        const igvDiv = document.getElementById("igv-div");
                        if (!igvDiv) return;

                        // Clear any existing IGV instance
                        if (window.igvBrowser) {{
                            igv.removeBrowser(window.igvBrowser);
                        }}

                        const options = {{
                            genome: "hg38",
                            locus: "{igv_locus}",
                            tracks: {json.dumps(tracks)},
                            roi: [
                                {{
                                    name: "{roi_description}",
                                    color: "{roi_color}",
                                    features: [
                                        {{
                                            chr: "{chrom}",
                                            start: {roi_start_coord},
                                            end: {roi_end_coord},
                                            name: "{call_value}"
                                        }}
                                    ]
                                }}
                            ]
                        }};

                        window.igvBrowser = await igv.createBrowser(igvDiv, options);
                        
                        // Trigger resize to ensure proper width
                        setTimeout(() => {{
                            if (window.igvBrowser && window.igvBrowser.visibilityChange) {{
                                window.igvBrowser.visibilityChange();
                            }}
                        }}, 100);
                    }})();
                    """

                    ui.run_javascript(igv_script, timeout=10.0)

                refresh_igv()

            # Handle expansion panel opening to resize IGV
            def on_cnv_expansion_change(e):
                if e.value:  # Panel was opened
                    ui.run_javascript("""
                        setTimeout(() => {
                            if (window.igvBrowser && window.igvBrowser.visibilityChange) {
                                window.igvBrowser.visibilityChange();
                            }
                        }, 100);
                    """)

            cnv_expansion.on_value_change(on_cnv_expansion_change)

            # IGV.js viewer container for VAF bedgraph
            vaf_expansion = (
                ui.expansion("VAF View", icon="water_drop", value=False)
                .classes("w-full mt-4 border border-gray-400 rounded-lg")
                .props("header-class='text-lg font-bold bg-gray-100'")
            )

            with vaf_expansion:
                igv_vaf_container = (
                    ui.element("div")
                    .props('id="igv-vaf-div"')
                    .classes("w-full border border-gray-300 rounded-lg")
                    .style(f"height: {calculate_vaf_height()}px")
                )

                @ui.refreshable
                def refresh_igv_vaf():
                    """Refresh VAF IGV viewer with current tracks."""
                    igv_vaf_container.clear()

                    # Update container height
                    try:
                        new_height = calculate_vaf_height()
                        igv_vaf_container.style(f"height: {new_height}px")
                    except RuntimeError:
                        pass

                    tracks = build_vaf_tracks()

                    # Determine ROI color based on call type
                    call_value = sv_data.get("call", "")
                    roi_color = "rgba(255,0,0,0.2)"  # Default red for loss
                    if "GAIN" in str(call_value).upper():
                        roi_color = "rgba(0,128,0,0.2)"  # Green for gain

                    # Create ROI description
                    roi_description = f"{chrom}:{start}-{end} - {call_value}"

                    # IGV.js initialization script
                    igv_vaf_script = f"""
                    (async function() {{
                        const igvDiv = document.getElementById("igv-vaf-div");
                        if (!igvDiv) return;

                        // Clear any existing IGV instance
                        if (window.igvVafBrowser) {{
                            igv.removeBrowser(window.igvVafBrowser);
                        }}

                        const options = {{
                            genome: "hg38",
                            locus: "{igv_locus}",
                            tracks: {json.dumps(tracks)},
                            roi: [
                                {{
                                    name: "{roi_description}",
                                    color: "{roi_color}",
                                    features: [
                                        {{
                                            chr: "{chrom}",
                                            start: {roi_start_coord},
                                            end: {roi_end_coord},
                                            name: "{call_value}"
                                        }}
                                    ]
                                }}
                            ]
                        }};

                        window.igvVafBrowser = await igv.createBrowser(igvDiv, options);

                        // Trigger resize to ensure proper width
                        setTimeout(() => {{
                            if (window.igvVafBrowser && window.igvVafBrowser.visibilityChange) {{
                                window.igvVafBrowser.visibilityChange();
                            }}
                        }}, 100);
                    }})();
                    """

                    ui.run_javascript(igv_vaf_script, timeout=10.0)

                refresh_igv_vaf()

            # Handle expansion panel opening to resize IGV
            def on_vaf_expansion_change(e):
                if e.value:  # Panel was opened
                    ui.run_javascript("""
                        setTimeout(() => {
                            if (window.igvVafBrowser && window.igvVafBrowser.visibilityChange) {
                                window.igvVafBrowser.visibilityChange();
                            }
                        }, 100);
                    """)

            vaf_expansion.on_value_change(on_vaf_expansion_change)

            # Second IGV.js viewer container for CRAM (read-level split view)
            cram_expansion = (
                ui.expansion("Read-Level Split View", icon="dna")
                .classes("w-full mt-6 border border-gray-400 rounded-lg")
                .props("default-opened header-class='text-lg font-bold bg-gray-100'")
            )

            with cram_expansion:
                igv_cram_container = (
                    ui.element("div")
                    .props('id="igv-cram-div"')
                    .classes("w-full border border-gray-300 rounded-lg")
                    .style(f"height: {calculate_cram_height()}px")
                )

                @ui.refreshable
                def refresh_igv_cram():
                    """Refresh CRAM IGV viewer with current tracks."""
                    igv_cram_container.clear()

                    # Update container height
                    try:
                        new_height = calculate_cram_height()
                        igv_cram_container.style(f"height: {new_height}px")
                    except RuntimeError:
                        pass

                    cram_tracks = build_cram_tracks()

                    # Split view locus format: "chr:(start-1500)-(start+1500) chr:(end-1500)-(end+1500)"
                    # Use curated coordinates if available
                    split_start = int(roi_start_coord)
                    split_end = int(roi_end_coord)
                    start_window_start = max(0, split_start - 1500)
                    start_window_end = split_start + 1500
                    end_window_start = max(0, split_end - 1500)
                    end_window_end = split_end + 1500
                    split_locus = f"{chrom}:{start_window_start}-{start_window_end} {chrom}:{end_window_start}-{end_window_end}"

                    # Determine ROI color based on call type (same as CNV view)
                    call_value = sv_data.get("call", "")
                    roi_color = "rgba(255,0,0,0.2)"  # Default red for loss
                    if "GAIN" in str(call_value).upper():
                        roi_color = "rgba(0,128,0,0.2)"  # Green for gain

                    # Create ROI description
                    roi_description = f"{chrom}:{start}-{end} - {call_value}"

                    # IGV.js initialization script for CRAM viewer
                    igv_cram_script = f"""
                    (async function() {{
                        const igvDiv = document.getElementById("igv-cram-div");
                        if (!igvDiv) return;

                        // Clear any existing IGV instance
                        if (window.igvCramBrowser) {{
                            igv.removeBrowser(window.igvCramBrowser);
                        }}

                        const options = {{
                            genome: "hg38",
                            locus: "{split_locus}",
                            showCenterGuide: true,
                            tracks: {json.dumps(cram_tracks)},
                            roi: [
                                {{
                                    name: "{roi_description}",
                                    color: "{roi_color}",
                                    features: [
                                        {{
                                            chr: "{chrom}",
                                            start: {roi_start_coord},
                                            end: {roi_end_coord},
                                            name: "{call_value}"
                                        }}
                                    ]
                                }}
                            ]
                        }};

                        window.igvCramBrowser = await igv.createBrowser(igvDiv, options);
                        
                        // Trigger resize to ensure proper width
                        setTimeout(() => {{
                            if (window.igvCramBrowser && window.igvCramBrowser.visibilityChange) {{
                                window.igvCramBrowser.visibilityChange();
                            }}
                        }}, 100);
                    }})();
                    """

                    ui.run_javascript(igv_cram_script, timeout=10.0)

                refresh_igv_cram()

                # Store references to curated input fields (will be set later)
                curated_inputs = {"start": None, "end": None}

                # Store current ROI coordinates - use curated boundaries if available
                roi_start = int(curated_start) if curated_start else int(start)
                roi_end = int(curated_end) if curated_end else int(end)
                roi_coords = {"start": roi_start, "end": roi_end}

                def _update_roi():
                    """Update ROI in both IGV instances from roi_coords."""
                    call_value = sv_data.get("call", "")
                    roi_color = (
                        "rgba(255,0,0,0.2)"
                        if "GAIN" not in str(call_value).upper()
                        else "rgba(0,128,0,0.2)"
                    )
                    roi_description = f"{chrom}:{roi_coords['start']}-{roi_coords['end']} - {call_value}"
                    update_roi_script = f"""
                    const roiConfig = [
                        {{
                            name: "{roi_description}",
                            color: "{roi_color}",
                            features: [
                                {{
                                    chr: "{chrom}",
                                    start: {roi_coords["start"]},
                                    end: {roi_coords["end"]},
                                    name: "{call_value}"
                                }}
                            ]
                        }}
                    ];
                    if (window.igvBrowser) {{
                        window.igvBrowser.clearROIs();
                        window.igvBrowser.loadROI(roiConfig);
                    }}
                    if (window.igvCramBrowser) {{
                        window.igvCramBrowser.clearROIs();
                        window.igvCramBrowser.loadROI(roiConfig);
                    }}
                    """
                    ui.run_javascript(update_roi_script)

                # Buttons for new start/end positions
                async def on_new_start_click():
                    """Handle New Start button click - show current loci."""
                    loci = await ui.run_javascript(
                        """
                        return window.igvCramBrowser ? window.igvCramBrowser.currentLoci() : null;
                    """,
                        timeout=2.0,
                    )
                    if loci and isinstance(loci, list) and len(loci) > 0:
                        # Parse locus string: chr:start-end
                        locus_str = loci[0]
                        try:
                            # Split by ':' to get chr and range
                            parts = locus_str.split(":")
                            if len(parts) == 2:
                                # Split range by '-' to get start and end
                                range_parts = parts[1].split("-")
                                if len(range_parts) == 2:
                                    start_pos = float(range_parts[0])
                                    end_pos = float(range_parts[1])
                                    mean_pos = round((start_pos + end_pos) / 2)

                                    # Update curated start input
                                    if curated_inputs["start"]:
                                        curated_inputs["start"].value = str(mean_pos)

                                    # Update ROI in both IGV instances
                                    roi_coords["start"] = mean_pos
                                    _update_roi()

                                    ui.notify(f"New Start - {mean_pos}", type="info")
                                else:
                                    ui.notify(f"New Start - {locus_str}", type="info")
                            else:
                                ui.notify(f"New Start - {locus_str}", type="info")
                        except (ValueError, IndexError):
                            ui.notify(f"New Start - {locus_str}", type="info")
                    else:
                        ui.notify("IGV browser not initialized", type="warning")

                async def on_new_end_click():
                    """Handle New End button click - show current loci."""
                    loci = await ui.run_javascript(
                        """
                        return window.igvCramBrowser ? window.igvCramBrowser.currentLoci() : null;
                    """,
                        timeout=2.0,
                    )
                    if loci and isinstance(loci, list) and len(loci) > 1:
                        # Parse locus string: chr:start-end
                        locus_str = loci[1]
                        try:
                            # Split by ':' to get chr and range
                            parts = locus_str.split(":")
                            if len(parts) == 2:
                                # Split range by '-' to get start and end
                                range_parts = parts[1].split("-")
                                if len(range_parts) == 2:
                                    start_pos = float(range_parts[0])
                                    end_pos = float(range_parts[1])
                                    mean_pos = round((start_pos + end_pos) / 2)

                                    # Update curated end input
                                    if curated_inputs["end"]:
                                        curated_inputs["end"].value = str(mean_pos)

                                    # Update ROI in both IGV instances
                                    roi_coords["end"] = mean_pos
                                    _update_roi()

                                    ui.notify(f"New End - {mean_pos}", type="info")
                                else:
                                    ui.notify(f"New End - {locus_str}", type="info")
                            else:
                                ui.notify(f"New End - {locus_str}", type="info")
                        except (ValueError, IndexError):
                            ui.notify(f"New End - {locus_str}", type="info")
                    else:
                        ui.notify("IGV browser not initialized", type="warning")

                with ui.row().classes("w-full items-center justify-center gap-0 mt-4"):
                    ui.button("New Start", on_click=on_new_start_click).props(
                        "color=primary"
                    )
                    ui.element("div").classes("w-[40%]")  # 40% spacer
                    ui.button("New End", on_click=on_new_end_click).props(
                        "color=primary"
                    )

            # Handle expansion panel opening to resize IGV
            def on_cram_expansion_change(e):
                if e.value:  # Panel was opened
                    ui.run_javascript("""
                        setTimeout(() => {
                            if (window.igvCramBrowser && window.igvCramBrowser.visibilityChange) {
                                window.igvCramBrowser.visibilityChange();
                            }
                        }, 100);
                    """)

            cram_expansion.on_value_change(on_cram_expansion_change)

            # SV Validation section
            # Determine SV type for variant key (dup for gain, del for loss)
            sv_type = infer_sv_type(sv_data)
            variant_key = f"{chrom}:{start}-{end}:{sv_type}"

            validation_file = store.data_dir / "validations" / "svs.tsv"

            ui.label("SV Validation").classes("text-xl font-semibold mb-2 mt-4")

            with ui.card().classes("w-full p-4 mb-4"):
                with ui.column().classes("gap-4"):
                    with ui.row().classes("items-center gap-4 w-full flex-wrap"):
                        ui.label("User:").classes("font-semibold")
                        user_input = (
                            ui.input("Username")
                            .props("outlined dense readonly")
                            .classes("w-48")
                        )
                        user_input.value = get_current_user()

                        ui.label("Inheritance:").classes("font-semibold ml-4")
                        inheritance_select = (
                            ui.select(
                                [
                                    "unknown",
                                    "de novo",
                                    "paternal",
                                    "maternal",
                                    "not paternal",
                                    "not maternal",
                                    "either",
                                    "homozygous",
                                ],
                                value="unknown",
                            )
                            .props("outlined dense")
                            .classes("w-40")
                        )

                        ui.label("Validation:").classes("font-semibold ml-4")
                        validation_select = (
                            ui.select(
                                [
                                    "present",
                                    "absent",
                                    "uncertain",
                                    "different",
                                ],
                                value="present",
                            )
                            .props("outlined dense")
                            .classes("w-40")
                        )

                        ui.label("Curated Position:").classes("font-semibold")
                        curated_start_input = (
                            ui.input("Curated start")
                            .props("outlined dense")
                            .classes("w-48")
                        )
                        curated_inputs["start"] = curated_start_input

                        ui.label("-").classes("font-semibold ml-4")
                        curated_end_input = (
                            ui.input("Curated end")
                            .props("outlined dense")
                            .classes("w-48")
                        )
                        curated_inputs["end"] = curated_end_input

                        def _refresh_roi_from_curated():
                            """Update ROI from curated position inputs."""
                            new_start = (
                                curated_start_input.value.strip()
                                if curated_start_input.value
                                else ""
                            )
                            new_end = (
                                curated_end_input.value.strip()
                                if curated_end_input.value
                                else ""
                            )
                            if not new_start and not new_end:
                                ui.notify(
                                    "Enter curated start and/or end first",
                                    type="warning",
                                )
                                return
                            try:
                                if new_start:
                                    roi_coords["start"] = int(new_start)
                                if new_end:
                                    roi_coords["end"] = int(new_end)
                            except ValueError:
                                ui.notify(
                                    "Curated positions must be integers",
                                    type="warning",
                                )
                                return
                            _update_roi()
                            ui.notify(
                                f"ROI updated: {chrom}:{roi_coords['start']}-{roi_coords['end']}",
                                type="positive",
                            )

                        ui.button(
                            icon="refresh",
                            on_click=_refresh_roi_from_curated,
                        ).props("flat round dense size=sm color=blue").tooltip(
                            "Refresh IGV ROI to curated coordinates"
                        )

                    with ui.row().classes("items-center gap-4 w-full"):
                        ui.label("Comment:").classes("font-semibold")
                        comment_input = (
                            ui.input("Optional comment")
                            .props("outlined dense")
                            .classes("flex-grow")
                        )

                        if can_write():
                            ui.button(
                                "Save Validation",
                                icon="save",
                                on_click=lambda: save_sv_validation(),
                            ).props("color=blue")
                        else:
                            ui.label("Read-only access").classes(
                                "text-sm text-gray-400 italic"
                            )

                    def save_sv_validation():
                        """Save an SV validation."""
                        if not can_write():
                            ui.notify("Permission denied", type="negative")
                            return
                        user = user_input.value.strip()
                        inheritance = inheritance_select.value
                        validation_status = validation_select.value
                        curated_start = (
                            curated_start_input.value.strip()
                            if curated_start_input.value
                            else ""
                        )
                        curated_end = (
                            curated_end_input.value.strip()
                            if curated_end_input.value
                            else ""
                        )
                        comment = (
                            comment_input.value.strip() if comment_input.value else ""
                        )

                        if not user:
                            ui.notify("Please enter a username", type="warning")
                            return

                        timestamp = datetime.now().isoformat()

                        try:
                            # Ensure directory exists
                            validation_file.parent.mkdir(parents=True, exist_ok=True)

                            # Check if file exists
                            file_exists = validation_file.exists()

                            # Append validation
                            with open(validation_file, "a") as f:
                                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                                try:
                                    # Only write header if file is new
                                    if not file_exists:
                                        f.write(
                                            "FID\tVariant\tSample\tUser\tInheritance\tValidation\tCuratedStart\tCuratedEnd\tComment\tIgnore\tTimestamp\n"
                                        )
                                    f.write(
                                        f"{family_id}\t{variant_key}\t{sample}\t{user}\t{inheritance}\t{validation_status}\t{curated_start}\t{curated_end}\t{comment}\t0\t{timestamp}\n"
                                    )

                                    # If validation is "present" with maternal/paternal/either inheritance,
                                    # also save validations for the parent(s)
                                    parent_samples_added = []
                                    if (
                                        validation_status == "present"
                                        and inheritance
                                        in ["maternal", "paternal", "either"]
                                    ):
                                        # Load existing validations to check if parents already have validations
                                        existing_validations = set()
                                        if file_exists:
                                            # Re-read the file to get existing validations
                                            with open(validation_file, "r") as check_f:
                                                reader = csv.DictReader(
                                                    check_f, delimiter="\t"
                                                )
                                                for row in reader:
                                                    if (
                                                        row.get("FID") == family_id
                                                        and row.get("Variant")
                                                        == variant_key
                                                    ):
                                                        existing_validations.add(
                                                            row.get("Sample")
                                                        )

                                        # Determine which parent(s) to add validations for
                                        parents_to_add = []
                                        if (
                                            inheritance == "maternal"
                                            and sample_parents["mother"]
                                        ):
                                            if (
                                                sample_parents["mother"]
                                                not in existing_validations
                                            ):
                                                parents_to_add.append(
                                                    ("mother", sample_parents["mother"])
                                                )
                                        elif (
                                            inheritance == "paternal"
                                            and sample_parents["father"]
                                        ):
                                            if (
                                                sample_parents["father"]
                                                not in existing_validations
                                            ):
                                                parents_to_add.append(
                                                    ("father", sample_parents["father"])
                                                )
                                        elif inheritance == "either":
                                            if (
                                                sample_parents["mother"]
                                                and sample_parents["mother"]
                                                not in existing_validations
                                            ):
                                                parents_to_add.append(
                                                    ("mother", sample_parents["mother"])
                                                )
                                            if (
                                                sample_parents["father"]
                                                and sample_parents["father"]
                                                not in existing_validations
                                            ):
                                                parents_to_add.append(
                                                    ("father", sample_parents["father"])
                                                )

                                        # Write parent validations
                                        for parent_type, parent_id in parents_to_add:
                                            parent_comment = (
                                                f"(inherited from {sample})"
                                            )
                                            f.write(
                                                f"{family_id}\t{variant_key}\t{parent_id}\t{user}\tunknown\tpresent\t{curated_start}\t{curated_end}\t{parent_comment}\t0\t{timestamp}\n"
                                            )
                                            parent_samples_added.append(parent_id)

                                finally:
                                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)

                            # Notify user about saved validations
                            notify_message = f"SV validation saved: {validation_status}"
                            if parent_samples_added:
                                notify_message += f" (also saved for {', '.join(parent_samples_added)})"
                            ui.notify(notify_message, type="positive")

                            # Reload validation history
                            load_sv_validation_history()

                            # Call the callback if provided, then close dialog
                            if on_validation_saved:
                                on_validation_saved()
                            dialog.close()

                        except Exception as e:
                            ui.notify(f"Error saving validation: {e}", type="negative")
                            import traceback

                            print(traceback.format_exc())

            # Validation history container
            validation_history_container = ui.column().classes("w-full mb-4")

            def load_sv_validation_history():
                """Load and display SV validation history."""
                validations = []

                if validation_file.exists():
                    with open(validation_file, "r") as f:
                        reader = csv.DictReader(f, delimiter="\t")
                        for row in reader:
                            if (
                                row.get("FID") == family_id
                                and row.get("Variant") == variant_key
                                and row.get("Sample") == sample
                            ):
                                validations.append(row)

                validation_history_container.clear()

                with validation_history_container:
                    ui.label("Previous validations:").classes("font-semibold mb-2")
                    if not validations:
                        ui.label("No validations recorded yet").classes(
                            "text-gray-500 text-sm italic"
                        )
                    else:
                        with ui.card().classes("w-full p-2"):
                            with ui.column().classes("gap-2"):
                                for validation in validations:
                                    val_status = validation.get("Validation", "")
                                    user = validation.get("User", "")
                                    timestamp = validation.get("Timestamp", "")
                                    comment = validation.get("Comment", "")
                                    curated_start = validation.get("CuratedStart", "")
                                    curated_end = validation.get("CuratedEnd", "")
                                    is_ignored = validation.get("Ignore", "0") == "1"

                                    # Format timestamp
                                    try:
                                        dt = datetime.fromisoformat(timestamp)
                                        formatted_time = dt.strftime(
                                            "%Y-%m-%d %H:%M:%S"
                                        )
                                    except Exception:
                                        formatted_time = timestamp

                                    # Color based on status
                                    if is_ignored:
                                        badge_color = "bg-gray-400"
                                    elif val_status == "present":
                                        badge_color = "bg-green-600"
                                    elif val_status == "absent":
                                        badge_color = "bg-red-600"
                                    elif val_status == "uncertain":
                                        badge_color = "bg-yellow-600"
                                    elif val_status == "different":
                                        badge_color = "bg-orange-600"
                                    else:
                                        badge_color = "bg-gray-500"

                                    with ui.row().classes(
                                        "items-center gap-2 p-2 w-full"
                                    ):
                                        ui.badge(val_status).classes(
                                            f"{badge_color} text-xs"
                                        )
                                        ui.label(user).classes("text-sm font-semibold")
                                        ui.label(formatted_time).classes(
                                            "text-xs text-gray-500"
                                        )
                                        if curated_start or curated_end:
                                            ui.label(
                                                f"Curated: {curated_start or start}-{curated_end or end}"
                                            ).classes("text-xs text-blue-600")
                                        if comment:
                                            ui.label(f'"{comment}"').classes(
                                                "text-sm text-gray-600 italic"
                                            )
                                        if is_ignored:
                                            ui.badge("IGNORED").classes(
                                                "bg-gray-400 text-xs"
                                            )

            load_sv_validation_history()

        dialog.open()
