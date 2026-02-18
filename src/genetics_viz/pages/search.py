"""Search page for cohort-wide variant search."""

import asyncio
import csv
import io
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import uuid4

import polars as pl
from nicegui import app as nicegui_app
from nicegui import context, ui

from genetics_viz.components.column_selector import build_column_selector
from genetics_viz.components.search_stats import show_stats_dialog
from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.utils.locus import filter_bed_dataframe, filter_dataframe, parse_locus_query
from genetics_viz.utils.wisecondorx import parse_wisecondorx_bed
from genetics_viz.components.header import create_header
from genetics_viz.components.tanstack_table import DataTable
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.sv_dialog import show_sv_dialog
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.utils.view_presets import VIEW_PRESETS, select_preset_for_config
from genetics_viz.utils.column_names import (
    apply_width_constraints,
    get_column_group,
    get_column_sorting,
    get_display_label,
    get_dropped_columns,
    get_schema_overrides,
    reorder_columns_by_group,
)
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.genesets import load_genesets
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.score_colors import get_score_color
from genetics_viz.utils.clinvar import (
    CLINVAR_COLORS,
    format_clinvar_display,
    get_clinvar_color,
)
from genetics_viz.utils.vep import (
    VEP_CONSEQUENCES,
    format_consequence_display,
    get_consequence_color,
    get_consequence_impact,
)





# Keys that are internal/display-only and should never appear in TSV exports
_INTERNAL_KEYS = {
    "_source_type",
    "_cohort_name",
    "_is_snv",
    "_original_locus",
    "_curated_tooltip",
    "IsCurated",
    "ConsequenceBadges",
    "ClinVarBadges",
    "GeneBadges",
    "VEP_Gene_badges",
}


def _is_export_column(key: str) -> bool:
    """Return True if the key should be included in a TSV export."""
    if key in _INTERNAL_KEYS:
        return False
    if key.endswith("_badge") or key.endswith("_badges"):
        return False
    return True


def _generate_tsv(rows: List[Dict[str, Any]], columns: List[str]) -> str:
    """Generate TSV content string from rows using only the specified columns."""
    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t")
    writer.writerow(columns)
    for row in rows:
        writer.writerow([row.get(col, "") for col in columns])
    return output.getvalue()


def _build_gene_badges(symbols: List[str], gene_scorer) -> List[Dict[str, Any]]:
    """Build score-sorted, capped gene badges (svs_tab pattern)."""
    badges = []
    for s in symbols:
        score, _ = gene_scorer.get_gene_score_and_sets(s)
        badges.append(
            {
                "label": s,
                "color": gene_scorer.get_gene_color(s),
                "tooltip": gene_scorer.get_gene_tooltip(s),
                "score": score,
            }
        )
    badges.sort(key=lambda x: x["score"], reverse=True)
    total = len(badges)
    if total > 6:
        badges = badges[:6]
        remaining = total - 6
        badges.append(
            {
                "label": f"+{remaining} genes",
                "color": "#9e9e9e",
                "tooltip": f"{remaining} more genes",
            }
        )
    return badges


def _apply_curated_coordinates(
    row: Dict[str, Any],
    validation_map: Dict,
    sv_variant_key: str,
    sample_id: str,
) -> None:
    """Apply curated coordinate bounds from validation data to a WCX row.

    If a 'present' validation with curated start/end exists, updates
    chr:start-end and Variant with curated values, sets IsCurated=True,
    and adds a tooltip showing original vs curated coordinates.
    """
    row["IsCurated"] = False
    row["_curated_tooltip"] = ""
    map_key = (sv_variant_key, sample_id)
    if map_key not in validation_map:
        return
    validations = validation_map[map_key]
    present_with_curated = [
        v for v in validations if v[0] == "present" and v[3] != "1" and (v[4] or v[5])
    ]
    if not present_with_curated:
        return
    present_with_curated.sort(key=lambda v: v[6], reverse=True)
    most_recent = present_with_curated[0]
    curated_start = most_recent[4]
    curated_end = most_recent[5]

    original_locus = row.get("_original_locus", row.get("chr:start-end", ""))
    parts = original_locus.split(":")
    if len(parts) != 2:
        return
    chrom = parts[0]
    range_parts = parts[1].split("-")
    if len(range_parts) != 2:
        return
    orig_start = range_parts[0]
    orig_end = range_parts[1]
    new_start = curated_start if curated_start else orig_start
    new_end = curated_end if curated_end else orig_end
    row["chr:start-end"] = f"{chrom}:{new_start}-{new_end}"
    row["Variant"] = row["chr:start-end"]
    row["IsCurated"] = True
    row["_curated_tooltip"] = f"Original: {original_locus}\nCurated: {row['chr:start-end']}"
    # Recompute svlen from curated boundaries
    try:
        row["svlen"] = int(new_end) - int(new_start)
    except (ValueError, TypeError):
        pass


def _infer_sv_type(row: Dict[str, Any]) -> str:
    """Infer SV type (dup/del) from WCX row data.

    Checks wisecondorX call label, then call column, then raw type, then ratio.
    Same logic as svs_tab.py's type inference.
    """
    # Check wisecondorX call classification (search page)
    wcx_call = str(row.get("wisecondorX", "")).upper()
    if "GAIN" in wcx_call:
        return "dup"
    if "LOSS" in wcx_call:
        return "del"
    # Check call column (svs_tab compatibility)
    call = str(row.get("call", "")).upper()
    if "GAIN" in call:
        return "dup"
    if "LOSS" in call:
        return "del"
    # Check raw type column
    raw_type = str(row.get("type", "")).lower()
    if raw_type in ("dup", "del"):
        return raw_type
    # Fall back to ratio sign
    try:
        ratio_val = float(row.get("ratio", 0) or 0)
        return "dup" if ratio_val > 0 else "del"
    except (ValueError, TypeError):
        return "del"


# Values that represent "unknown/missing" in pedigree files
_PED_MISSING = {"", "0", "-9"}


def _pedigree_data_from_cohort(cohort_name: str) -> Dict[str, Dict[str, str]]:
    """Build pedigree lookup from the already-parsed Cohort object.

    Uses the DataStore's Cohort (which handles all header formats robustly)
    instead of re-parsing the pedigree file.

    Returns dict mapping sample_id -> {FID, Father, Mother, Sex, Phenotype}.
    """
    store = get_data_store()
    cohort = store.get_cohort(cohort_name)
    if cohort is None:
        return {}

    pedigree_data: Dict[str, Dict[str, str]] = {}
    for family in cohort.families.values():
        for sample in family.samples:
            pedigree_data[sample.sample_id] = {
                "FID": sample.family_id,
                "Father": sample.father_id or "",
                "Mother": sample.mother_id or "",
                "Sex": sample.sex or "",
                "Phenotype": sample.phenotype or "",
            }
    return pedigree_data


@ui.page("/search/{cohort_name}")
def search_cohort_page(cohort_name: str) -> None:
    """Search page for cohort-wide variant search."""
    create_header(cohort_name)

    # Add IGV.js library
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

        # Scan for source files (wombat + SVS)
        wombat_dir = store.data_dir / "cohorts" / cohort_name / "wombat"
        svs_dir = store.data_dir / "cohorts" / cohort_name / "svs"

        with ui.column().classes("w-full px-6 py-6"):
            # Title
            ui.label(f"🔍 Search: {cohort_name}").classes(
                "text-3xl font-bold text-blue-900 mb-6"
            )

            # Scan for wombat files matching pattern
            source_files: List[Dict[str, Any]] = []

            if wombat_dir.exists():
                wombat_pattern = re.compile(
                    rf"{re.escape(cohort_name)}\.rare\.([^.]+)\.(.+?)\.results\.tsv$"
                )
                for tsv_file in wombat_dir.glob("*.tsv"):
                    match = wombat_pattern.match(tsv_file.name)
                    if match:
                        vep_config = match.group(1)
                        wombat_config = match.group(2)
                        source_files.append(
                            {
                                "file_path": tsv_file,
                                "vep_config": vep_config,
                                "wombat_config": wombat_config,
                                "display_name": wombat_config,
                                "source_type": "wombat",
                                "parent_dir": "wombat",
                            }
                        )

            # Scan for SVS subdirectories
            if svs_dir.exists():
                for subdir in sorted(svs_dir.iterdir()):
                    if not subdir.is_dir():
                        continue
                    subdir_name = subdir.name
                    if subdir_name == "wisecondorx":
                        bed_file = subdir / f"{cohort_name}_aberrations.bed"
                        if not bed_file.exists():
                            bed_file = subdir / f"{cohort_name}_aberrations.annotated.bed"
                        if bed_file.exists():
                            source_files.append(
                                {
                                    "file_path": bed_file,
                                    "vep_config": "",
                                    "wombat_config": "",
                                    "display_name": subdir_name,
                                    "source_type": "wisecondorx",
                                    "parent_dir": "svs",
                                }
                            )

            if not source_files:
                ui.label(
                    "No source files found (checked wombat/ and svs/ directories)"
                ).classes("text-gray-500 text-lg italic")
                return

            def _source_key(sf: Dict[str, Any]) -> str:
                return f"{sf['parent_dir']}/{sf['display_name']}"

            # Load pedigree data from the already-parsed Cohort object
            pedigree_data = _pedigree_data_from_cohort(cohort_name)
            sample_to_family = {
                sid: ped["FID"] for sid, ped in pedigree_data.items()
            }

            # Derive unique sex and phenotype values for individual filters
            # Exclude missing-value sentinels ("", "0", "-9") from the option lists
            available_sex_values = sorted(
                {v.get("Sex", "") for v in pedigree_data.values()}
                - _PED_MISSING
            )
            available_phenotype_values = sorted(
                {v.get("Phenotype", "") for v in pedigree_data.values()}
                - _PED_MISSING
            )

            # Load genesets from params/genesets
            available_genesets = load_genesets(store.data_dir)

            # Separate source lists by type for UI
            wombat_sources = [sf for sf in source_files if sf["source_type"] == "wombat"]
            wcx_sources = [sf for sf in source_files if sf["source_type"] == "wisecondorx"]

            # WisecondorX call values (from svs_tab.py pattern)
            ALL_CALL_VALUES = [
                "Robust LOSS",
                "Robust GAIN",
                "Permissive LOSS",
                "Permissive Gain",
                "Below threshold",
            ]

            def _make_wombat_component(source_key: str = "") -> Dict[str, Any]:
                """Create a new wombat source filter component state dict."""
                if not source_key and wombat_sources:
                    source_key = _source_key(wombat_sources[0])
                return {
                    "id": str(uuid4()),
                    "type": "wombat",
                    "source_key": source_key,
                    "locus": {"value": ""},
                    "genesets": {"value": []},
                    "impacts": {"value": list(VEP_CONSEQUENCES.keys())},
                    "validations": {"value": ["present", "absent", "uncertain", "conflicting", "TODO"]},
                    "exclude_lcr": {"value": True},
                    "exclude_gnomad": {"value": True},
                    "exclude_gnomad_wgs": {"value": False},
                }

            def _make_wcx_component() -> Dict[str, Any]:
                """Create a new wisecondorx source filter component state dict."""
                wcx_key = _source_key(wcx_sources[0]) if wcx_sources else ""
                return {
                    "id": str(uuid4()),
                    "type": "wisecondorx",
                    "source_key": wcx_key,
                    "locus": {"value": ""},
                    "genesets": {"value": []},
                    "selected_calls": {"value": [c for c in ALL_CALL_VALUES if c != "Below threshold"]},
                    "exonic_only": {"value": False},
                    "validations": {"value": ["present", "absent", "uncertain", "conflicting", "TODO"]},
                    "ratio_min": {"value": None},
                    "ratio_max": {"value": None},
                }

            # Initialize component list with one default wombat component
            source_components: List[Dict[str, Any]] = []
            if wombat_sources:
                source_components.append(_make_wombat_component())

            # Individual filters state (global, not per-component)
            filter_sex: Dict[str, List[str]] = {"value": []}
            filter_phenotype: Dict[str, List[str]] = {"value": []}
            filter_has_parents: Dict[str, bool] = {"value": False}

            # --- Helper: build geneset menu for a component ---
            def _build_geneset_menu(comp: Dict[str, Any]) -> None:
                """Build an inline genesets dropdown menu for a component."""
                if not available_genesets:
                    return
                comp_genesets = comp["genesets"]
                geneset_btn = (
                    ui.button("Genesets", icon="list")
                    .props(
                        ("outline" if not comp_genesets["value"] else "unelevated color=green")
                        + " dense size=sm"
                    )
                )
                btn_ref = {"button": geneset_btn}
                with geneset_btn:
                    with ui.menu():
                        ui.label("Select Genesets:").classes("px-4 py-2 font-semibold text-sm")
                        ui.separator()
                        with ui.column().classes("p-2"):
                            gs_cbs: Dict[str, Any] = {}

                            def _update_gs_btn():
                                if btn_ref["button"]:
                                    if comp_genesets["value"]:
                                        btn_ref["button"].props(remove="outline", add="unelevated color=green")
                                    else:
                                        btn_ref["button"].props(remove="unelevated color=green", add="outline")
                                    btn_ref["button"].update()

                            with ui.row().classes("gap-2 mb-2"):
                                def _gs_all(cbs=gs_cbs, st=comp_genesets, upd=_update_gs_btn):
                                    st["value"] = list(available_genesets.keys())
                                    for cb in cbs.values():
                                        cb.value = True
                                    upd()

                                def _gs_none(cbs=gs_cbs, st=comp_genesets, upd=_update_gs_btn):
                                    st["value"] = []
                                    for cb in cbs.values():
                                        cb.value = False
                                    upd()

                                ui.button("All", on_click=_gs_all).props("size=sm flat dense").classes("text-xs")
                                ui.button("None", on_click=_gs_none).props("size=sm flat dense").classes("text-xs")

                            ui.separator()
                            for gs_name in sorted(available_genesets.keys()):
                                def make_gs_handler(name, st=comp_genesets, upd=_update_gs_btn):
                                    def handler(e):
                                        if e.value:
                                            if name not in st["value"]:
                                                st["value"].append(name)
                                        else:
                                            if name in st["value"]:
                                                st["value"].remove(name)
                                        upd()
                                    return handler

                                gs_cbs[gs_name] = ui.checkbox(
                                    f"{gs_name} ({len(available_genesets[gs_name])} genes)",
                                    value=gs_name in comp_genesets["value"],
                                    on_change=make_gs_handler(gs_name),
                                ).classes("text-sm")

            # --- Helper: build impact menu for a wombat component ---
            def _build_impact_menu(comp: Dict[str, Any]) -> None:
                """Build an inline impacts dropdown menu for a wombat component."""
                comp_impacts = comp["impacts"]
                impact_btn = ui.button("Impacts", icon="filter_list").props("outline dense size=sm")
                btn_ref = {"button": impact_btn}
                with impact_btn:
                    with ui.menu():
                        ui.label("Select Impact Types:").classes("px-4 py-2 font-semibold text-sm")
                        ui.separator()
                        with ui.column().classes("p-2"):
                            imp_cbs: Dict[str, Any] = {}

                            def _update_imp_btn():
                                if btn_ref["button"]:
                                    if len(comp_impacts["value"]) == len(VEP_CONSEQUENCES):
                                        btn_ref["button"].props(remove="unelevated color=orange", add="outline")
                                    else:
                                        btn_ref["button"].props(remove="outline", add="unelevated color=orange")
                                    btn_ref["button"].update()

                            with ui.row().classes("gap-2 mb-2 flex-wrap"):
                                def _imp_all(cbs=imp_cbs, st=comp_impacts, upd=_update_imp_btn):
                                    st["value"] = list(VEP_CONSEQUENCES.keys())
                                    for cb in cbs.values():
                                        cb.value = True
                                    upd()

                                def _imp_none(cbs=imp_cbs, st=comp_impacts, upd=_update_imp_btn):
                                    st["value"] = []
                                    for cb in cbs.values():
                                        cb.value = False
                                    upd()

                                def _make_level_handler(level, cbs=imp_cbs, st=comp_impacts, upd=_update_imp_btn):
                                    def handler():
                                        selected = [c for c, (imp, _) in VEP_CONSEQUENCES.items() if imp == level]
                                        st["value"] = selected
                                        for name, cb in cbs.items():
                                            cb.value = name in selected
                                        upd()
                                    return handler

                                ui.button("All", on_click=_imp_all).props("size=sm flat dense").classes("text-xs")
                                ui.button("None", on_click=_imp_none).props("size=sm flat dense").classes("text-xs")
                                ui.button("HIGH", on_click=_make_level_handler("HIGH")).props("size=sm flat dense color=red").classes("text-xs")
                                ui.button("MODERATE", on_click=_make_level_handler("MODERATE")).props("size=sm flat dense color=orange").classes("text-xs")
                                ui.button("LOW", on_click=_make_level_handler("LOW")).props("size=sm flat dense color=yellow-8").classes("text-xs")
                                ui.button("MODIFIER", on_click=_make_level_handler("MODIFIER")).props("size=sm flat dense color=grey").classes("text-xs")

                            ui.separator()
                            with ui.column().classes("gap-1"):
                                def make_imp_handler(cons, st=comp_impacts, upd=_update_imp_btn):
                                    def handler(e):
                                        if e.value:
                                            if cons not in st["value"]:
                                                st["value"].append(cons)
                                        else:
                                            if cons in st["value"]:
                                                st["value"].remove(cons)
                                        upd()
                                    return handler

                                for impact_level in ["HIGH", "MODERATE", "LOW", "MODIFIER"]:
                                    consequences = [c for c, (imp, _) in VEP_CONSEQUENCES.items() if imp == impact_level]
                                    if consequences:
                                        ui.label(f"{impact_level}:").classes("text-xs font-bold text-gray-600 mt-2")
                                        for cons in sorted(consequences):
                                            imp_cbs[cons] = ui.checkbox(
                                                format_consequence_display(cons),
                                                value=cons in comp_impacts["value"],
                                                on_change=make_imp_handler(cons),
                                            ).classes("text-sm")

            # --- Helper: build impacts (call filter) menu for a wcx component ---
            def _build_call_filter_menu(comp: Dict[str, Any]) -> None:
                """Build an Impacts dropdown menu for a wisecondorx component."""
                comp_calls = comp["selected_calls"]
                call_btn = ui.button("Impacts", icon="filter_list").props("outline dense size=sm")
                btn_ref = {"button": call_btn}
                with call_btn:
                    with ui.menu():
                        ui.label("Filter by Impact:").classes("px-4 py-2 font-semibold text-sm")
                        ui.separator()
                        with ui.column().classes("p-2"):
                            call_cbs: Dict[str, Any] = {}

                            def _update_call_btn():
                                if btn_ref["button"]:
                                    if len(comp_calls["value"]) == len(ALL_CALL_VALUES):
                                        btn_ref["button"].props(remove="unelevated color=orange", add="outline")
                                    else:
                                        btn_ref["button"].props(remove="outline", add="unelevated color=orange")
                                    btn_ref["button"].update()

                            with ui.row().classes("gap-2 mb-2"):
                                def _call_all(cbs=call_cbs, st=comp_calls, upd=_update_call_btn):
                                    st["value"] = list(ALL_CALL_VALUES)
                                    for cb in cbs.values():
                                        cb.value = True
                                    upd()

                                def _call_none(cbs=call_cbs, st=comp_calls, upd=_update_call_btn):
                                    st["value"] = []
                                    for cb in cbs.values():
                                        cb.value = False
                                    upd()

                                ui.button("All", on_click=_call_all).props("size=sm flat dense").classes("text-xs")
                                ui.button("None", on_click=_call_none).props("size=sm flat dense").classes("text-xs")

                            ui.separator()
                            for call_val in ALL_CALL_VALUES:
                                def make_call_handler(cv, st=comp_calls, upd=_update_call_btn):
                                    def handler(e):
                                        if e.value:
                                            if cv not in st["value"]:
                                                st["value"].append(cv)
                                        else:
                                            if cv in st["value"]:
                                                st["value"].remove(cv)
                                        upd()
                                    return handler

                                call_cbs[call_val] = ui.checkbox(
                                    call_val,
                                    value=call_val in comp_calls["value"],
                                    on_change=make_call_handler(call_val),
                                ).classes("text-sm")

            # --- Forward reference for search handler + locus input tracking ---
            _search_handler: Dict[str, Any] = {"fn": None}
            locus_inputs: List[Any] = []

            def _bind_locus_enter(loc_input) -> None:
                """Bind Enter key on a locus input to trigger search."""
                if _search_handler["fn"]:
                    loc_input.on("keydown.enter", _search_handler["fn"])

            # Search panel
            with ui.card().classes("w-full p-2 mb-2").props("flat bordered"):
                # Side-by-side panels
                with ui.row().classes("w-full gap-4 items-start flex-nowrap"):
                    # --- LEFT PANEL: Variants (2/3 width) ---
                    with ui.column().classes("flex-[2]"):
                        ui.label("Variants").classes("text-sm font-semibold text-gray-600")

                        @ui.refreshable
                        def render_component_list():
                            locus_inputs.clear()
                            for i, comp in enumerate(source_components):
                                if i > 0:
                                    ui.label("OR").classes(
                                        "text-xs font-bold text-center w-full text-gray-400 my-0 py-0"
                                    )
                                if comp["type"] == "wombat":
                                    _render_wombat_card(comp)
                                elif comp["type"] == "wisecondorx":
                                    _render_wcx_card(comp)

                            # Add buttons row
                            with ui.row().classes("gap-2 mt-1"):
                                if wombat_sources:
                                    def _add_wombat():
                                        source_components.append(_make_wombat_component())
                                        render_component_list.refresh()
                                    ui.button("+ Wombat", on_click=_add_wombat).props(
                                        "flat dense size=sm color=blue no-caps"
                                    )
                                if wcx_sources:
                                    def _add_wcx():
                                        source_components.append(_make_wcx_component())
                                        render_component_list.refresh()
                                    ui.button("+ WisecondorX", on_click=_add_wcx).props(
                                        "flat dense size=sm color=teal no-caps"
                                    )

                        def _render_wombat_card(comp: Dict[str, Any]) -> None:
                            """Render a compact wombat source filter card."""
                            with ui.card().classes("w-full p-2 mb-0").props("flat bordered"):
                                # Row 1: label + file selector + locus + remove
                                with ui.row().classes("items-center gap-2 w-full"):
                                    ui.label("Wombat").classes("text-sm font-semibold text-blue-700")
                                    wombat_options = {
                                        _source_key(sf): sf["display_name"]
                                        for sf in wombat_sources
                                    }
                                    ui.select(
                                        options=wombat_options,
                                        value=comp["source_key"],
                                        label="Source",
                                        on_change=lambda e, c=comp: c.update({"source_key": e.value}),
                                    ).props("outlined dense").classes("w-64")

                                    loc_input = (
                                        ui.input(
                                            label="Locus (optional)",
                                            placeholder="chr1:10000-10100, SHANK3, ENSG...",
                                            value=comp["locus"]["value"],
                                            on_change=lambda e, c=comp: c["locus"].update({"value": e.value or ""}),
                                        )
                                        .props("outlined dense")
                                        .classes("flex-grow")
                                    )
                                    locus_inputs.append(loc_input)
                                    _bind_locus_enter(loc_input)

                                    def make_remove_handler(comp_id):
                                        def handler():
                                            source_components[:] = [
                                                c for c in source_components if c["id"] != comp_id
                                            ]
                                            render_component_list.refresh()
                                        return handler

                                    ui.button(
                                        icon="close",
                                        on_click=make_remove_handler(comp["id"]),
                                    ).props("flat round dense size=sm color=grey")

                                # Row 2: genesets + impacts + validation + exclude checkboxes
                                with ui.row().classes("items-center gap-2 w-full flex-wrap"):
                                    _build_geneset_menu(comp)
                                    _build_impact_menu(comp)
                                    create_validation_filter_menu(
                                        all_statuses=["present", "absent", "uncertain", "conflicting", "TODO"],
                                        filter_state=comp["validations"],
                                        on_change=lambda: None,
                                        label="Validation",
                                        button_classes="",
                                        button_props="dense size=sm",
                                    )
                                    ui.checkbox(
                                        "Exclude LCR",
                                        value=comp["exclude_lcr"]["value"],
                                        on_change=lambda e, c=comp: c["exclude_lcr"].update({"value": e.value}),
                                    ).props("dense").classes("text-xs")
                                    ui.checkbox(
                                        "Exclude gnomAD filtered",
                                        value=comp["exclude_gnomad"]["value"],
                                        on_change=lambda e, c=comp: c["exclude_gnomad"].update({"value": e.value}),
                                    ).props("dense").classes("text-xs")
                                    ui.checkbox(
                                        "Exclude gnomAD WGS",
                                        value=comp["exclude_gnomad_wgs"]["value"],
                                        on_change=lambda e, c=comp: c["exclude_gnomad_wgs"].update({"value": e.value}),
                                    ).props("dense").classes("text-xs")

                        def _render_wcx_card(comp: Dict[str, Any]) -> None:
                            """Render a compact wisecondorx source filter card."""
                            with ui.card().classes("w-full p-2 mb-0").props("flat bordered"):
                                # Row 1: label + exonic only + locus + remove
                                with ui.row().classes("items-center gap-2 w-full"):
                                    ui.label("WisecondorX").classes("text-sm font-semibold text-teal-700")

                                    ui.checkbox(
                                        "Exonic only",
                                        value=comp["exonic_only"]["value"],
                                        on_change=lambda e, c=comp: c["exonic_only"].update({"value": e.value}),
                                    ).props("dense").classes("text-xs")

                                    loc_input = (
                                        ui.input(
                                            label="Locus (optional)",
                                            placeholder="chr1:10000-10100, SHANK3, ENSG...",
                                            value=comp["locus"]["value"],
                                            on_change=lambda e, c=comp: c["locus"].update({"value": e.value or ""}),
                                        )
                                        .props("outlined dense")
                                        .classes("flex-grow")
                                    )
                                    locus_inputs.append(loc_input)
                                    _bind_locus_enter(loc_input)

                                    def make_remove_handler(comp_id):
                                        def handler():
                                            source_components[:] = [
                                                c for c in source_components if c["id"] != comp_id
                                            ]
                                            render_component_list.refresh()
                                        return handler

                                    ui.button(
                                        icon="close",
                                        on_click=make_remove_handler(comp["id"]),
                                    ).props("flat round dense size=sm color=grey")

                                # Row 2: genesets + impacts + validation + ratio filter
                                with ui.row().classes("items-center gap-2 w-full flex-wrap"):
                                    _build_geneset_menu(comp)
                                    _build_call_filter_menu(comp)
                                    create_validation_filter_menu(
                                        all_statuses=["present", "absent", "uncertain", "conflicting", "TODO"],
                                        filter_state=comp["validations"],
                                        on_change=lambda: None,
                                        label="Validation",
                                        button_classes="",
                                        button_props="dense size=sm",
                                    )

                                    ui.label("|ratio|:").classes("text-xs text-gray-600")

                                    def _make_ratio_handler(comp_ref, key):
                                        def handler(e):
                                            raw = (e.value or "").replace(",", ".")
                                            if not raw:
                                                comp_ref[key].update({"value": None})
                                                e.sender.props(remove="error")
                                                return
                                            try:
                                                comp_ref[key].update({"value": float(raw)})
                                                e.sender.props(remove="error")
                                            except ValueError:
                                                e.sender.props("error")
                                                ui.notify(
                                                    f"Invalid number: {raw}",
                                                    type="warning",
                                                    position="top",
                                                    timeout=2000,
                                                )
                                        return handler

                                    ui.input(
                                        label="Min",
                                        value=str(comp["ratio_min"]["value"] or ""),
                                        on_change=_make_ratio_handler(comp, "ratio_min"),
                                    ).props("outlined dense").classes("w-20")
                                    ui.input(
                                        label="Max",
                                        value=str(comp["ratio_max"]["value"] or ""),
                                        on_change=_make_ratio_handler(comp, "ratio_max"),
                                    ).props("outlined dense").classes("w-20")

                        render_component_list()

                    ui.separator().props("vertical")

                    # --- RIGHT PANEL: Individuals (1/3 width) ---
                    with ui.column().classes("flex-[1] min-w-0"):
                        ui.label("Individuals").classes("text-sm font-semibold text-gray-600")
                        with ui.column().classes("gap-2"):
                            sex_select = ui.select(
                                options=available_sex_values,
                                label="Sex",
                                value=filter_sex["value"],
                                multiple=True,
                                on_change=lambda e: filter_sex.update({"value": e.value or []}),
                            ).props("outlined dense use-chips").classes("w-full")

                            phenotype_select = ui.select(
                                options=available_phenotype_values,
                                label="Phenotype",
                                value=filter_phenotype["value"],
                                multiple=True,
                                on_change=lambda e: filter_phenotype.update({"value": e.value or []}),
                            ).props("outlined dense use-chips").classes("w-full")

                            has_parents_cb = ui.checkbox(
                                "Only samples with both parents",
                                value=filter_has_parents["value"],
                                on_change=lambda e: filter_has_parents.update({"value": e.value}),
                            ).props("dense")

                # Search + Clear buttons centered below both panels
                with ui.row().classes("justify-center gap-4 mt-2"):
                    search_button = ui.button("Search", icon="search").props("color=blue dense")

                    def _clear_all():
                        source_components.clear()
                        if wombat_sources:
                            source_components.append(_make_wombat_component())
                        filter_sex["value"] = []
                        filter_phenotype["value"] = []
                        filter_has_parents["value"] = False
                        sex_select.value = []
                        phenotype_select.value = []
                        has_parents_cb.value = False
                        render_component_list.refresh()

                    ui.button("Clear", icon="clear_all", on_click=_clear_all).props("outline dense")

            # Results container
            results_container = ui.column().classes("w-full")

            # Capture client context for callbacks
            page_client = context.client

            async def perform_search():
                """Execute the search — iterate over source components, union results."""
                results_container.clear()

                # Show progress indicator while loading
                with results_container:
                    with ui.column().classes("items-center gap-4 justify-center py-8 w-full"):
                        progress = ui.circular_progress(
                            min=0, max=100, value=0, size="xl", color="blue"
                        )
                        status_label = ui.label("Starting search...").classes(
                            "text-lg text-gray-600"
                        )

                # Validate: need at least one component
                if not source_components:
                    results_container.clear()
                    with results_container:
                        ui.label("Please add at least one source component").classes(
                            "text-orange-600"
                        )
                    return

                # Validate: every component must have locus or genesets
                for i, comp in enumerate(source_components):
                    locus_val = comp["locus"]["value"]
                    gs_val = comp["genesets"]["value"]
                    has_locus = bool(locus_val and locus_val.strip())
                    has_genesets = bool(gs_val)
                    if not has_locus and not has_genesets:
                        comp_label = f"Component {i + 1} ({comp['type']})"
                        results_container.clear()
                        with results_container:
                            ui.label(
                                f"{comp_label}: please enter a locus or select at least one geneset"
                            ).classes("text-orange-600")
                        return

                # Update progress
                progress.set_value(10)
                status_label.set_text("Loading data...")
                await asyncio.sleep(0)

                try:
                    # --- Helper: load a single wombat file ---
                    def _load_wombat_file(file_path: Path) -> Optional[pl.DataFrame]:
                        _df = pl.read_csv(
                            file_path,
                            separator="\t",
                            infer_schema_length=10000,
                            schema_overrides=get_schema_overrides(),
                            null_values=[".", ""],
                        )
                        _drop = get_dropped_columns() & set(_df.columns)
                        if _drop:
                            _df = _df.drop(list(_drop))
                        return _df

                    # --- Helper: apply geneset filter ---
                    def _apply_geneset_filter(df, comp, source_type):
                        if not comp["genesets"]["value"]:
                            return df
                        combined_genes = set()
                        for gs_name in comp["genesets"]["value"]:
                            combined_genes.update(available_genesets.get(gs_name, set()))
                        if not combined_genes:
                            return df

                        if source_type == "wombat" and "VEP_SYMBOL" in df.columns:
                            def matches_gs(symbol_str):
                                if not symbol_str:
                                    return False
                                symbols = [s.strip().upper() for s in str(symbol_str).split("&")]
                                return any(s in combined_genes for s in symbols)
                            return df.filter(
                                pl.col("VEP_SYMBOL").map_elements(matches_gs, return_dtype=pl.Boolean)
                            )
                        elif source_type == "wisecondorx":
                            exonic = comp.get("exonic_only", {}).get("value", False)
                            sym_col = "exonic_symbol" if exonic else "genic_symbol"
                            if sym_col in df.columns:
                                def matches_gs_bed(symbol_str):
                                    if not symbol_str:
                                        return False
                                    symbols = [s.strip().upper() for s in str(symbol_str).split(",")]
                                    return any(s in combined_genes for s in symbols)
                                return df.filter(
                                    pl.col(sym_col).map_elements(matches_gs_bed, return_dtype=pl.Boolean)
                                )
                        return df

                    # --- Process each component ---
                    all_wombat_dfs: List[pl.DataFrame] = []
                    all_wcx_dfs: List[pl.DataFrame] = []

                    for comp in source_components:
                        sf = next(
                            (s for s in source_files if _source_key(s) == comp["source_key"]),
                            None,
                        )
                        if sf is None:
                            continue

                        locus_val = comp["locus"]["value"]
                        query_params = None
                        if locus_val and locus_val.strip():
                            query_params = parse_locus_query(locus_val)

                        if comp["type"] == "wombat":
                            df = await asyncio.to_thread(_load_wombat_file, sf["file_path"])
                            if df is None or len(df) == 0:
                                continue
                            # Locus filter
                            if query_params:
                                df = filter_dataframe(df, query_params)
                            # Geneset filter
                            df = _apply_geneset_filter(df, comp, "wombat")
                            if df is not None and len(df) > 0:
                                all_wombat_dfs.append(df)

                        elif comp["type"] == "wisecondorx":
                            df = await asyncio.to_thread(parse_wisecondorx_bed, sf["file_path"])
                            if df is None or len(df) == 0:
                                continue
                            # Locus filter (use exonic columns if exonic_only is checked)
                            exonic = comp.get("exonic_only", {}).get("value", False)
                            if query_params:
                                df = filter_bed_dataframe(df, query_params, exonic=exonic)
                            # Geneset filter
                            df = _apply_geneset_filter(df, comp, "wisecondorx")
                            # Call filter
                            if "wisecondorX" in df.columns and comp.get("selected_calls"):
                                selected_calls = comp["selected_calls"]["value"]
                                if selected_calls:
                                    df = df.filter(pl.col("wisecondorX").is_in(selected_calls))
                            # Ratio filter (on |ratio|)
                            if "ratio" in df.columns:
                                ratio_min = comp.get("ratio_min", {}).get("value")
                                ratio_max = comp.get("ratio_max", {}).get("value")
                                if ratio_min is not None:
                                    df = df.filter(
                                        pl.col("ratio").cast(pl.Float64, strict=False).abs() >= ratio_min
                                    )
                                if ratio_max is not None:
                                    df = df.filter(
                                        pl.col("ratio").cast(pl.Float64, strict=False).abs() <= ratio_max
                                    )
                            if df is not None and len(df) > 0:
                                all_wcx_dfs.append(df)

                    # Update progress: data loaded
                    progress.set_value(35)
                    status_label.set_text("Filtering data...")
                    await asyncio.sleep(0)

                    # Deduplicate wombat rows across components
                    wombat_filtered = None
                    if all_wombat_dfs:
                        combined = pl.concat(all_wombat_dfs, how="diagonal_relaxed")
                        grouping_cols = ["#CHROM", "POS", "REF", "ALT", "sample"]
                        agg_cols = [col for col in combined.columns if col not in grouping_cols]
                        agg_exprs = [pl.len().alias("n_grouped")]
                        for col in agg_cols:
                            agg_exprs.append(
                                pl.col(col)
                                .cast(pl.Utf8)
                                .filter(
                                    (pl.col(col).is_not_null())
                                    & (pl.col(col).cast(pl.Utf8) != "")
                                    & (pl.col(col).cast(pl.Utf8) != ".")
                                )
                                .unique()
                                .str.join(",")
                                .alias(col)
                            )
                        wombat_filtered = combined.group_by(grouping_cols, maintain_order=True).agg(agg_exprs)

                    # Deduplicate wcx rows across components
                    wcx_filtered = None
                    if all_wcx_dfs:
                        combined = pl.concat(all_wcx_dfs, how="diagonal_relaxed")
                        dedup_cols = [c for c in ["chr", "start", "end", "type", "sample"] if c in combined.columns]
                        if dedup_cols:
                            combined = combined.unique(subset=dedup_cols, keep="first")
                        wcx_filtered = combined

                    # Update progress: filtering complete
                    progress.set_value(60)
                    status_label.set_text("Processing variants...")
                    await asyncio.sleep(0)

                    wombat_count = len(wombat_filtered) if wombat_filtered is not None else 0
                    wcx_count = len(wcx_filtered) if wcx_filtered is not None else 0

                    if wombat_count == 0 and wcx_count == 0:
                        results_container.clear()
                        with results_container:
                            ui.label("No results found").classes(
                                "text-gray-500 text-lg italic"
                            )
                        return

                    # Load validation data
                    snv_validation_map = {}
                    sv_validation_map = {}
                    validation_file = store.data_dir / "validations" / "snvs.tsv"
                    sv_validation_file = store.data_dir / "validations" / "svs.tsv"
                    if wombat_filtered is not None:
                        snv_validation_map = load_validation_map(validation_file, None)
                    if wcx_filtered is not None:
                        sv_validation_map = load_validation_map(sv_validation_file, None)

                    # Yield to event loop before badge processing
                    await asyncio.sleep(0)

                    # Track unknown terms for warnings
                    unknown_consequences = set()
                    unknown_clinvar_terms = set()
                    gene_scorer = get_gene_scorer()

                    all_rows: List[Dict[str, Any]] = []

                    # --- Process wombat rows ---
                    if wombat_filtered is not None and wombat_count > 0:
                        wombat_rows = wombat_filtered.to_dicts()

                        # Collect per-component filter state for wombat
                        # Use the first wombat component's filters as baseline
                        # (after dedup, we apply the union of all component filters)
                        wombat_comps = [c for c in source_components if c["type"] == "wombat"]

                        for row in wombat_rows:
                            row["_source_type"] = "wombat"
                            row["IsCurated"] = False
                            row["_curated_tooltip"] = ""
                            chrom = row.get("#CHROM", "")
                            pos = row.get("POS", "")
                            ref = row.get("REF", "")
                            alt = row.get("ALT", "")
                            sample_id = row.get("sample", "")
                            variant_key = f"{chrom}:{pos}:{ref}:{alt}"
                            row["Variant"] = variant_key
                            row["_cohort_name"] = cohort_name

                            ped_info = pedigree_data.get(sample_id, {})
                            row["FID"] = ped_info.get(
                                "FID", sample_to_family.get(sample_id, "")
                            )
                            row["Phenotype"] = ped_info.get("Phenotype", "")

                            # Consequence badges
                            consequence_str = row.get("VEP_Consequence", "")
                            if consequence_str:
                                consequences = []
                                for part in str(consequence_str).split(","):
                                    for cons in part.split("&"):
                                        cons = cons.strip()
                                        if cons:
                                            consequences.append(cons)
                                row["ConsequenceBadges"] = []
                                seen_badges = set()
                                for cons in consequences:
                                    if cons and cons not in VEP_CONSEQUENCES:
                                        unknown_consequences.add(cons)
                                    label = format_consequence_display(cons)
                                    color = get_consequence_color(cons)
                                    badge_key = (label, color)
                                    if badge_key not in seen_badges:
                                        seen_badges.add(badge_key)
                                        row["ConsequenceBadges"].append(
                                            {"label": label, "color": color}
                                        )
                            else:
                                row["ConsequenceBadges"] = []

                            # ClinVar badges
                            clinvar_str = row.get("VEP_CLIN_SIG", "")
                            if clinvar_str:
                                clinvar_sigs = []
                                for part in str(clinvar_str).split(","):
                                    for sig in part.split("&"):
                                        sig = sig.strip()
                                        if sig and sig != ".":
                                            clinvar_sigs.append(sig)
                                row["ClinVarBadges"] = []
                                seen_badges = set()
                                for sig in clinvar_sigs:
                                    sig_lower = sig.lower()
                                    is_known = any(
                                        key.lower() == sig_lower
                                        for key in CLINVAR_COLORS.keys()
                                    )
                                    if sig and not is_known:
                                        unknown_clinvar_terms.add(sig)
                                    label = format_clinvar_display(sig)
                                    color = get_clinvar_color(sig)
                                    badge_key = (label, color)
                                    if badge_key not in seen_badges:
                                        seen_badges.add(badge_key)
                                        row["ClinVarBadges"].append(
                                            {"label": label, "color": color}
                                        )
                            else:
                                row["ClinVarBadges"] = []

                            # Gene badges from VEP_SYMBOL (score-sorted, capped)
                            symbol_str = row.get("VEP_SYMBOL", "")
                            if symbol_str:
                                all_symbols: List[str] = []
                                for s in str(symbol_str).split(","):
                                    for part in s.split("&"):
                                        p = part.strip()
                                        if p and p not in all_symbols:
                                            all_symbols.append(p)
                                row["GeneBadges"] = _build_gene_badges(all_symbols, gene_scorer)
                            else:
                                row["GeneBadges"] = []
                            # Synthetic "gene" column for unified display
                            row["gene"] = ",".join(
                                b["label"] for b in row["GeneBadges"] if not b["label"].startswith("+")
                            )

                            # Gene badges from VEP_Gene (ENSG IDs)
                            gene_str = row.get("VEP_Gene", "")
                            if gene_str:
                                genes = [
                                    g.strip()
                                    for g in str(gene_str).split(",")
                                    if g.strip()
                                ]
                                row["VEP_Gene_badges"] = [
                                    {
                                        "label": g,
                                        "color": gene_scorer.get_gene_color(g),
                                        "tooltip": gene_scorer.get_gene_tooltip(g),
                                    }
                                    for g in genes
                                ]
                            else:
                                row["VEP_Gene_badges"] = []

                            add_validation_status_to_row(
                                row, snv_validation_map, variant_key, sample_id
                            )

                            # Continuous score badges
                            for col_name, value_str in list(row.items()):
                                if value_str and value_str != ".":
                                    try:
                                        value = float(value_str)
                                        badge_info = get_score_color(col_name, value)
                                        if badge_info:
                                            row[f"{col_name}_badge"] = {
                                                "label": f"{value:.3f}",
                                                "color": badge_info["color"],
                                                "tooltip": f"{col_name}: {value:.3f} ({badge_info['label']})",
                                            }
                                    except (ValueError, TypeError):
                                        pass

                        # Apply per-component wombat filters (impact, exclude, validation)
                        # A row passes if ANY wombat component would keep it
                        def _wombat_row_passes(row, comps):
                            """Check if a wombat row passes at least one component's filters."""
                            for wc in comps:
                                # Impact filter
                                impacts_ok = True
                                if set(wc["impacts"]["value"]) != set(VEP_CONSEQUENCES.keys()):
                                    consequence_str = row.get("VEP_Consequence", "")
                                    if consequence_str:
                                        cons_list = []
                                        for part in str(consequence_str).split(","):
                                            for c in part.split("&"):
                                                c = c.strip()
                                                if c:
                                                    cons_list.append(c)
                                        impacts_ok = any(c in wc["impacts"]["value"] for c in cons_list)
                                    else:
                                        impacts_ok = False

                                # Exclude LCR
                                if wc["exclude_lcr"]["value"]:
                                    if row.get("LCR") and "true" in str(row.get("LCR", "")).lower():
                                        continue
                                # Exclude gnomAD filtered
                                if wc["exclude_gnomad"]["value"]:
                                    if row.get("genomes_filters"):
                                        continue
                                # Exclude gnomAD WGS
                                if wc["exclude_gnomad_wgs"]["value"]:
                                    if row.get("fafmax_faf95_max_genomes"):
                                        continue

                                if not impacts_ok:
                                    continue

                                # Validation filter
                                val_status = row.get("Validation", "")
                                if wc["validations"]["value"]:
                                    if val_status in wc["validations"]["value"]:
                                        return True
                                    if "TODO" in wc["validations"]["value"] and not val_status:
                                        return True
                                    continue
                                return True
                            return False

                        wombat_rows = [r for r in wombat_rows if _wombat_row_passes(r, wombat_comps)]
                        all_rows.extend(wombat_rows)

                    # --- Process wisecondorx rows ---
                    if wcx_filtered is not None and wcx_count > 0:
                        wcx_rows = wcx_filtered.to_dicts()
                        for row in wcx_rows:
                            row["_source_type"] = "wisecondorx"
                            sample_id = row.get("sample", "")
                            row["_original_locus"] = row.get("chr:start-end", "")
                            row["Variant"] = row.get("chr:start-end", "")
                            row["_cohort_name"] = cohort_name
                            # Compute svlen from start/end
                            try:
                                row["svlen"] = int(row.get("end", 0)) - int(row.get("start", 0))
                            except (ValueError, TypeError):
                                row["svlen"] = None

                            ped_info = pedigree_data.get(sample_id, {})
                            row["FID"] = ped_info.get(
                                "FID", sample_to_family.get(sample_id, "")
                            )
                            row["Phenotype"] = ped_info.get("Phenotype", "")

                            # No VEP badges for SVS rows
                            row["ConsequenceBadges"] = []
                            row["ClinVarBadges"] = []

                            # Gene badges from genic_symbol (score-sorted, capped)
                            genic_str = row.get("genic_symbol", "")
                            if genic_str and str(genic_str).strip():
                                symbols = [
                                    s.strip()
                                    for s in str(genic_str).split(",")
                                    if s.strip()
                                ]
                                row["GeneBadges"] = _build_gene_badges(symbols, gene_scorer)
                            else:
                                row["GeneBadges"] = []
                            # Synthetic "gene" column for unified display
                            row["gene"] = ",".join(
                                b["label"] for b in row["GeneBadges"] if not b["label"].startswith("+")
                            )

                            row["VEP_Gene_badges"] = []

                            # SVS validation: variant_key is chr:start-end:type
                            sv_type = _infer_sv_type(row)
                            sv_variant_key = f"{row['_original_locus']}:{sv_type}"
                            add_validation_status_to_row(
                                row, sv_validation_map, sv_variant_key, sample_id
                            )
                            _apply_curated_coordinates(
                                row, sv_validation_map, sv_variant_key, sample_id
                            )

                        # Apply per-component WCX filters (validation)
                        wcx_comps = [c for c in source_components if c["type"] == "wisecondorx"]

                        def _wcx_row_passes(row, comps):
                            """Check if a WCX row passes at least one component's validation filter."""
                            for wc in comps:
                                val_status = row.get("Validation", "")
                                allowed = wc.get("validations", {}).get("value", [])
                                if allowed:
                                    if val_status in allowed:
                                        return True
                                    if "TODO" in allowed and not val_status:
                                        return True
                                    continue
                                return True
                            return False

                        wcx_rows = [r for r in wcx_rows if _wcx_row_passes(r, wcx_comps)]
                        all_rows.extend(wcx_rows)

                    # Update progress: badge processing complete
                    progress.set_value(85)
                    status_label.set_text("Rendering table...")
                    await asyncio.sleep(0)

                    # Display warnings for unknown terms
                    if unknown_consequences:
                        ui.notify(
                            f"⚠️ Unknown VEP consequence terms found: {', '.join(sorted(unknown_consequences))}. "
                            "Please add to vep_consequences.yaml",
                            type="warning",
                            timeout=10000,
                            position="top",
                        )

                    if unknown_clinvar_terms:
                        ui.notify(
                            f"⚠️ Unknown ClinVar terms found: {', '.join(sorted(unknown_clinvar_terms))}. "
                            "Please add to clinvar_colors.yaml",
                            type="warning",
                            timeout=10000,
                            position="top",
                        )

                    # Get all columns from both source types
                    all_col_set: set[str] = set()
                    if wombat_filtered is not None:
                        all_col_set.update(wombat_filtered.columns)
                    if wcx_filtered is not None:
                        all_col_set.update(wcx_filtered.columns)
                    # Remove internal/coordinate columns kept for filtering
                    all_col_set -= {"_source_type", "chr", "start", "end"}
                    all_col_set.add("gene")
                    if wcx_filtered is not None:
                        all_col_set.add("svlen")
                    all_columns = list(all_col_set)
                    if "Variant" not in all_columns:
                        all_columns.insert(0, "Variant")
                    for ensure_col in ["FID", "Phenotype", "Validation"]:
                        if ensure_col not in all_columns:
                            all_columns.append(ensure_col)

                    # Group same-group columns together
                    all_columns = reorder_columns_by_group(all_columns)

                    # Default visible columns
                    default_visible = [
                        "Variant",
                        "VEP_Consequence",
                        "gene",
                        "VEP_CLIN_SIG",
                        "fafmax_faf95_max_genomes",
                        "FID",
                        "Phenotype",
                        "sample",
                        "sample_gt",
                        "father_gt",
                        "mother_gt",
                        "Validation",
                    ]
                    # Auto-select preset based on component mix
                    has_wombat_comps = any(
                        c["type"] == "wombat" for c in source_components
                    )
                    has_wcx_comps = any(
                        c["type"] == "wisecondorx" for c in source_components
                    )
                    if has_wcx_comps and not has_wombat_comps:
                        # WCX-only: use SV View preset
                        initial_preset = next(
                            (p for p in VIEW_PRESETS if p["name"] == "SV View"),
                            VIEW_PRESETS[0],
                        )
                    elif has_wcx_comps and has_wombat_comps:
                        # Mixed: use Mix View preset
                        initial_preset = next(
                            (p for p in VIEW_PRESETS if p["name"] == "Mix View"),
                            VIEW_PRESETS[0],
                        )
                    else:
                        # Wombat-only: select based on first wombat config file
                        first_wombat_comp = next(
                            (c for c in source_components if c["type"] == "wombat"),
                            None,
                        )
                        wombat_config = ""
                        if first_wombat_comp:
                            sf = next(
                                (
                                    s
                                    for s in source_files
                                    if _source_key(s) == first_wombat_comp["source_key"]
                                ),
                                None,
                            )
                            if sf:
                                wombat_config = sf.get("wombat_config", "")
                        initial_preset = select_preset_for_config(
                            wombat_config, VIEW_PRESETS
                        )
                    selected_preset = {"name": initial_preset["name"]}

                    # Override with preset columns if available
                    preset_columns = initial_preset.get("columns", [])
                    initial_selected = [col for col in preset_columns if col in all_columns]

                    selected_cols = {
                        "value": initial_selected if initial_selected else [col for col in default_visible if col in all_columns]
                    }

                    # Table state for persistence across refreshes
                    table_state: Dict[str, Any] = {"sorting": [], "page": 0}

                    # Apply individual filters (sex, phenotype, has-parents)
                    if filter_sex["value"]:
                        all_rows = [
                            r for r in all_rows
                            if pedigree_data.get(r.get("sample", ""), {}).get("Sex", "") in filter_sex["value"]
                        ]

                    if filter_phenotype["value"]:
                        all_rows = [
                            r for r in all_rows
                            if pedigree_data.get(r.get("sample", ""), {}).get("Phenotype", "") in filter_phenotype["value"]
                        ]

                    if filter_has_parents["value"]:
                        def _has_parents(ped: Dict[str, str]) -> bool:
                            father = ped.get("Father", "")
                            mother = ped.get("Mother", "")
                            return (
                                father not in _PED_MISSING
                                and mother not in _PED_MISSING
                            )

                        all_rows = [
                            r for r in all_rows
                            if _has_parents(pedigree_data.get(r.get("sample", ""), {}))
                        ]

                    # Impact, exclude, validation filters are now applied per-component
                    # (wombat rows filtered by _wombat_row_passes, wcx rows filtered by call filter)

                    # Update progress: ready to display
                    progress.set_value(100)
                    status_label.set_text("Complete!")
                    await asyncio.sleep(0)

                    # Clear progress indicator and show results
                    results_container.clear()
                    with results_container:

                        @ui.refreshable
                        def render_results_table():
                            # All per-component filters (impact, exclude, validation, calls)
                            # are already applied — use all_rows directly
                            rows = all_rows.copy()



                            def get_columns():
                                cols: List[Dict[str, Any]] = [
                                    {
                                        "id": "actions",
                                        "header": "",
                                        "cellType": "action",
                                        "actionName": "view_variant",
                                        "actionIcon": "visibility",
                                        "actionColor": "#1976d2",
                                        "actionTooltip": "View in IGV",
                                        "sortable": False,
                                    }
                                ]
                                for col in all_columns:
                                    col_def: Dict[str, Any] = {
                                        "id": col,
                                        "header": get_display_label(col),
                                        "group": get_column_group(col),
                                        "sorting": get_column_sorting(col),
                                        "sortable": True,
                                    }
                                    if col == "Variant":
                                        col_def["cellType"] = "curated_locus"
                                        col_def["curatedField"] = "IsCurated"
                                        col_def["tooltipField"] = "_curated_tooltip"
                                    elif col == "Validation":
                                        col_def["cellType"] = "validation"
                                    elif col == "VEP_Consequence":
                                        col_def["cellType"] = "badge_list"
                                        col_def["badgesField"] = "ConsequenceBadges"
                                    elif col == "VEP_CLIN_SIG":
                                        col_def["cellType"] = "badge_list"
                                        col_def["badgesField"] = "ClinVarBadges"
                                    elif col == "gene":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "GeneBadges"
                                    elif col == "VEP_SYMBOL":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "GeneBadges"
                                    elif col == "VEP_Gene":
                                        col_def["cellType"] = "gene_badge"
                                        col_def["badgesField"] = "VEP_Gene_badges"
                                    elif col == "FID":
                                        col_def["cellType"] = "link"
                                        col_def["href"] = "/cohort/{_cohort_name}/family/{FID}"
                                    else:
                                        col_def["cellType"] = "score_badge"
                                    apply_width_constraints(col_def, col)
                                    cols.append(col_def)
                                return cols

                            # Reference to the DataTable for column visibility updates
                            search_dt: Dict[str, Any] = {"ref": None}

                            def _apply_col_visibility():
                                if search_dt["ref"]:
                                    visible = ["actions"] + list(selected_cols["value"])
                                    search_dt["ref"].set_column_visibility(visible)

                            with ui.row().classes("items-center gap-4 mt-4 mb-2 w-full"):
                                ui.label(f"Results ({len(rows)} rows)").classes(
                                    "text-lg font-semibold text-blue-700"
                                )

                                # Preset dropdown
                                preset_select = ui.select(
                                    options={p["name"]: p["name"] for p in VIEW_PRESETS},
                                    value=selected_preset["name"],
                                    label="Preset"
                                ).classes("w-48")

                                ui.space()  # Push column selector to the right

                                # Column selector dialog
                                col_dialog, _sync_col_selector = build_column_selector(
                                    all_columns=all_columns,
                                    selected_cols=selected_cols,
                                    on_visibility_change=_apply_col_visibility,
                                    presets=VIEW_PRESETS,
                                )
                                ui.button(
                                    "Columns", icon="view_column",
                                    on_click=col_dialog.open,
                                ).props("outline color=blue size=sm")

                                # --- Stats button + dialog ---
                                ui.button(
                                    "Stats", icon="bar_chart",
                                    on_click=lambda _e=None, r=rows: show_stats_dialog(r),
                                ).props("outline color=blue size=sm")

                                # --- Save / Download ---
                                def _save_validation_file(filename: str):
                                    """Write TSV to to_validate directory."""
                                    standard_preset = next(
                                        (p for p in VIEW_PRESETS if p["name"] == "Standard View"),
                                        VIEW_PRESETS[0],
                                    )
                                    save_cols = [
                                        col
                                        for col in standard_preset.get("columns", [])
                                        if col != "Validation"
                                        and col in all_columns
                                        and _is_export_column(col)
                                    ]
                                    tsv_content = _generate_tsv(rows, save_cols)
                                    to_validate_dir = store.data_dir / "to_validate"
                                    to_validate_dir.mkdir(parents=True, exist_ok=True)
                                    file_path = to_validate_dir / filename
                                    if file_path.exists():
                                        ui.notify(
                                            f"File '{filename}' already exists!",
                                            type="warning",
                                        )
                                        return
                                    file_path.write_text(tsv_content)
                                    ui.notify(
                                        f"Saved {filename} to to_validate/",
                                        type="positive",
                                        position="top",
                                    )

                                def _download_tsv(filename: str):
                                    """Trigger browser download of TSV with displayed columns."""
                                    download_cols = [
                                        col
                                        for col in selected_cols["value"]
                                        if _is_export_column(col)
                                    ]
                                    tsv_content = _generate_tsv(rows, download_cols)
                                    ui.download(
                                        tsv_content.encode("utf-8"),
                                        filename,
                                        media_type="text/tab-separated-values",
                                    )

                                def _build_save_dialog(title, confirm_label, confirm_icon, on_confirm):
                                    """Pre-build a save/download dialog. Returns (dialog, name_input)."""
                                    with ui.dialog() as dlg, ui.card().classes("p-4 min-w-[400px]"):
                                        ui.label(title).classes("text-lg font-semibold mb-2")

                                        with ui.row().classes("items-center w-full gap-1"):
                                            name_input = ui.input(
                                                label="Filename",
                                                value="",
                                            ).props("outlined dense").classes("flex-grow")
                                            ui.label(".tsv").classes("text-gray-500 text-sm mt-2")

                                        with ui.row().classes("justify-end gap-2 w-full mt-4"):
                                            ui.button("Cancel", on_click=dlg.close).props("flat")

                                            def _on_confirm(inp=name_input, d=dlg, cb=on_confirm):
                                                raw_name = inp.value.strip()
                                                if not raw_name:
                                                    ui.notify("Please enter a filename", type="warning")
                                                    return
                                                if not raw_name.endswith(".tsv"):
                                                    raw_name += ".tsv"
                                                d.close()
                                                cb(raw_name)

                                            ui.button(
                                                confirm_label, icon=confirm_icon, on_click=_on_confirm
                                            ).props("unelevated color=blue")
                                    return dlg, name_input

                                save_dlg, save_input = _build_save_dialog(
                                    "Save as validation file", "Save", "save", _save_validation_file,
                                )
                                dl_dlg, dl_input = _build_save_dialog(
                                    "Download TSV", "Download", "download", _download_tsv,
                                )

                                def _open_dialog(dlg, name_input):
                                    """Reset filename, open dialog, auto-focus and select text."""
                                    name_input.value = (
                                        f"genetics-viz.{datetime.now().strftime('%Y-%m-%d.%H-%M')}"
                                    )
                                    dlg.open()
                                    ui.timer(
                                        0.2,
                                        lambda inp=name_input: ui.run_javascript(
                                            f'var el = getElement({inp.id}).$el.querySelector("input");'
                                            f"if(el){{el.focus();el.select();}}"
                                        ),
                                        once=True,
                                    )

                                with ui.button(
                                    "Save", icon="save",
                                ).props("outline color=blue size=sm"):
                                    with ui.menu():
                                        ui.menu_item(
                                            "Save as validation file",
                                            on_click=lambda: _open_dialog(save_dlg, save_input),
                                        )
                                        ui.menu_item(
                                            "Download TSV",
                                            on_click=lambda: _open_dialog(dl_dlg, dl_input),
                                        )

                            # Preset change handler
                            def on_preset_change(e):
                                """Handle preset selection change."""
                                preset_name = e.value
                                preset = next((p for p in VIEW_PRESETS if p["name"] == preset_name), None)
                                if not preset:
                                    return

                                # Filter columns to only those available in the data
                                available = [col for col in preset.get("columns", []) if col in all_columns]

                                selected_cols["value"] = available
                                selected_preset["name"] = preset_name
                                _apply_col_visibility()
                                _sync_col_selector()

                            # Connect preset change handler
                            preset_select.on_value_change(on_preset_change)

                            # Handle view variant click
                            def on_view_variant(e):
                                row_data = e.get("row", {})
                                variant_str = row_data.get("Variant", "")
                                sample_id = row_data.get("sample", "")
                                source_type = row_data.get("_source_type", "wombat")

                                # Get family from sample
                                family_id = sample_to_family.get(sample_id)

                                if not family_id:
                                    ui.notify(
                                        f"Could not find family for sample {sample_id}",
                                        type="warning",
                                    )
                                    return

                                try:
                                    if source_type == "wisecondorx":
                                        # SV dialog: use original locus (before curated update)
                                        locus_for_dialog = row_data.get(
                                            "_original_locus", variant_str
                                        )
                                        parts = locus_for_dialog.split(":")
                                        if len(parts) == 2:
                                            chrom = parts[0]
                                            range_parts = parts[1].split("-")
                                            if len(range_parts) == 2:
                                                def on_sv_save():
                                                    sv_val_updated = load_validation_map(
                                                        sv_validation_file, None
                                                    )
                                                    for row in all_rows:
                                                        if row.get("_source_type") == "wisecondorx":
                                                            orig = row.get(
                                                                "_original_locus",
                                                                row.get("chr:start-end", ""),
                                                            )
                                                            sv_t = _infer_sv_type(row)
                                                            sv_vk = f"{orig}:{sv_t}"
                                                            # Reset chr:start-end to original before re-applying
                                                            row["chr:start-end"] = orig
                                                            row["Variant"] = orig
                                                            add_validation_status_to_row(
                                                                row, sv_val_updated,
                                                                sv_vk, row.get("sample", ""),
                                                            )
                                                            _apply_curated_coordinates(
                                                                row, sv_val_updated,
                                                                sv_vk, row.get("sample", ""),
                                                            )
                                                    with page_client:
                                                        ui.timer(
                                                            0.1,
                                                            render_results_table.refresh,
                                                            once=True,
                                                        )

                                                show_sv_dialog(
                                                    cohort_name=cohort_name,
                                                    family_id=family_id,
                                                    chrom=chrom,
                                                    start=range_parts[0],
                                                    end=range_parts[1],
                                                    sample=sample_id,
                                                    sv_data={
                                                        **row_data,
                                                        "call": row_data.get(
                                                            "wisecondorX",
                                                            row_data.get("call", ""),
                                                        ),
                                                    },
                                                    on_validation_saved=on_sv_save,
                                                )
                                            else:
                                                ui.notify(
                                                    "Invalid SV format", type="warning"
                                                )
                                        else:
                                            ui.notify(
                                                "Invalid SV format", type="warning"
                                            )
                                    else:
                                        # Wombat variant dialog
                                        parts = variant_str.split(":")
                                        if len(parts) == 4:
                                            chrom, pos, ref, alt = parts
                                            variant_data = dict(row_data)

                                            def on_save(validation_status: str):
                                                snv_val_updated = load_validation_map(
                                                    validation_file, None
                                                )
                                                for row in all_rows:
                                                    if row.get("_source_type") != "wisecondorx":
                                                        v_key = row.get("Variant", "")
                                                        s_id = row.get("sample", "")
                                                        add_validation_status_to_row(
                                                            row, snv_val_updated, v_key, s_id,
                                                        )
                                                with page_client:
                                                    ui.timer(
                                                        0.1,
                                                        render_results_table.refresh,
                                                        once=True,
                                                    )

                                            show_variant_dialog(
                                                cohort_name=cohort_name,
                                                family_id=family_id,
                                                chrom=chrom,
                                                pos=pos,
                                                ref=ref,
                                                alt=alt,
                                                sample=sample_id,
                                                variant_data=variant_data,
                                                on_save_callback=on_save,
                                            )
                                        else:
                                            ui.notify(
                                                "Invalid variant format. Expected chr:pos:ref:alt",
                                                type="warning",
                                            )
                                except Exception as ex:
                                    ui.notify(
                                        f"Error parsing variant: {ex}", type="warning"
                                    )

                            # Restore table state (sorting / page) across refreshes
                            saved_sorting = table_state.get("sorting", [])
                            if saved_sorting:
                                col_id = saved_sorting[0]["id"]
                                desc = saved_sorting[0].get("desc", False)
                                col_def = next(
                                    (c for c in get_columns() if c.get("id") == col_id), {}
                                )
                                sort_field = col_def.get("sortField", col_id)
                                sort_type = col_def.get("sorting", "")
                                if sort_type == "genomic":
                                    from genetics_viz.utils.column_names import genomic_sort_key
                                    rows.sort(
                                        key=lambda r: (
                                            r.get(sort_field) is None,
                                            genomic_sort_key(r.get(sort_field, "")),
                                        ),
                                        reverse=desc,
                                    )
                                elif sort_type == "numerical":
                                    def _num_key(r):
                                        v = r.get(sort_field)
                                        if v is None:
                                            return (True, 0.0)
                                        try:
                                            return (False, float(v))
                                        except (ValueError, TypeError):
                                            return (True, 0.0)
                                    rows.sort(key=_num_key, reverse=desc)
                                else:
                                    rows.sort(
                                        key=lambda r: (
                                            r.get(sort_field) is None,
                                            r.get(sort_field, ""),
                                        ),
                                        reverse=desc,
                                    )

                            # Create table with all columns, initial visibility from preset
                            search_dt["ref"] = DataTable(
                                columns=get_columns(),
                                rows=rows,
                                row_key="Variant",
                                pagination={"rowsPerPage": 50},
                                visible_columns=["actions"] + list(selected_cols["value"]),
                                on_row_action=on_view_variant,
                                initial_sorting=saved_sorting,
                                initial_page=table_state.get("page", 0),
                                state_holder=table_state,
                            )

                        render_results_table()

                except Exception as e:
                    import traceback

                    results_container.clear()
                    with results_container:
                        ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
                        ui.label("Traceback:").classes("text-red-500 font-semibold")
                        ui.label(traceback.format_exc()).classes(
                            "text-red-500 text-xs font-mono whitespace-pre"
                        )

            # Set up handlers after function definition
            search_button.on_click(perform_search)
            _search_handler["fn"] = perform_search
            # Bind Enter key on initial locus inputs
            # (new inputs from add-component get bound in _render_*_card via _bind_locus_enter)
            for loc_inp in locus_inputs:
                loc_inp.on("keydown.enter", perform_search)

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
