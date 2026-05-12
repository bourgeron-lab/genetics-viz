"""Standalone family page — displays family data independently of any cohort/project."""

import asyncio
from datetime import datetime
from typing import Dict, List

from nicegui import ui

from genetics_viz.components.diagnostic_loader import (
    ensure_diagnostic_file,
    load_family_diagnostics,
)
from genetics_viz.components.header import create_header
from genetics_viz.components.notes_loader import (
    delete_note,
    ensure_notes_file,
    load_family_notes,
    save_note,
)
from genetics_viz.components.sample_dialog import show_sample_dialog
from genetics_viz.pages.cohort.components.svs_tab import render_svs_tab
from genetics_viz.pages.cohort.components.wombat_tab import render_wombat_tab
from genetics_viz.pages.cohort.family import (
    _parse_sv_cytoband,
    _render_sv_gene_badges,
)
from genetics_viz.utils.auth import can_write, check_auth, get_current_user
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.data_availability import check_family_availability
from genetics_viz.utils.gene_scoring import get_gene_scorer
from genetics_viz.utils.pedigree import load_family_object, load_family_pedigree
from genetics_viz.utils.sharding import get_family_path


def _render_availability_badges(avail: dict) -> None:
    """Render colored badges for data availability."""
    _ITEMS = [
        ("pedigree", "Pedigree"),
        ("vcfs", "VCFs"),
        ("wombat", "Wombat"),
        ("wisecondorx", "WisecondorX"),
        ("extractor", "Extractor"),
    ]
    for key, label in _ITEMS:
        present = avail.get(key, False)
        color = "green" if present else "grey"
        ui.badge(label, color=color).props("outline" if not present else "")


@ui.page("/family/{family_id}")
async def standalone_family_page(family_id: str) -> None:
    """Render the standalone family detail page (no cohort context)."""
    if redirect := check_auth():
        return redirect
    create_header()

    # Add IGV.js library at page level
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>'
    )

    try:
        store = get_data_store()
        family_path = get_family_path(store.data_dir, family_id)

        if not family_path.is_dir():
            with ui.column().classes("w-full px-6 py-6"):
                ui.label(f"Family not found: {family_id}").classes(
                    "text-xl text-red-500"
                )
                ui.button("← Back to Home", on_click=lambda: ui.navigate.to("/"))
            return

        pedigree_file = family_path / f"{family_id}.pedigree.tsv"

        # Load pedigree data asynchronously
        if pedigree_file.exists():
            members_data = await asyncio.to_thread(load_family_pedigree, pedigree_file)
            family = await asyncio.to_thread(
                load_family_object, pedigree_file, family_id
            )
        else:
            members_data = []
            family = None

        sample_ids = [m["Sample ID"] for m in members_data]

        # Check data availability asynchronously
        avail = await asyncio.to_thread(
            check_family_availability, store.data_dir, family_id, sample_ids
        )

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.label(f"Family: {family_id}").classes("font-semibold")

            # Family header
            with ui.row().classes("items-center gap-4 mb-4"):
                ui.label(f"Family: {family_id}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                if family:
                    ui.badge(f"{family.num_samples} members").props("color=blue")
                    ui.badge(f"{family.num_founders} founders").props("color=teal")

            # Data availability panel
            with ui.card().classes("w-full mb-4"):
                with ui.column().classes("p-4 gap-2"):
                    ui.label("Data Availability").classes(
                        "text-lg font-semibold text-blue-700 mb-2"
                    )
                    # Family-level availability
                    with ui.row().classes("items-center gap-2 mb-2"):
                        ui.label("Family:").classes("font-semibold text-sm w-16")
                        _render_availability_badges(avail)

                    # Per-sample availability
                    sample_avail = avail.get("samples", {})
                    if sample_avail:
                        ui.separator()
                        ui.label("Samples:").classes(
                            "font-semibold text-sm text-gray-600 mt-1"
                        )
                        for sid, s_avail in sample_avail.items():
                            with ui.row().classes("items-center gap-2"):
                                ui.label(sid).classes("font-mono text-sm w-28")
                                _SAMPLE_ITEMS = [
                                    ("cram", "CRAM"),
                                    ("bedgraph", "Bedgraph"),
                                    ("vaf_bedgraph", "VAF"),
                                    ("deepvariant", "DeepVariant"),
                                    ("svs", "SVs"),
                                ]
                                for key, label in _SAMPLE_ITEMS:
                                    present = s_avail.get(key, False)
                                    color = "green" if present else "grey"
                                    ui.badge(label, color=color).props(
                                        "outline" if not present else ""
                                    ).classes("text-xs")

                                def make_viz_handler(sample_id):
                                    def handler():
                                        show_sample_dialog(sample_id)

                                    return handler

                                if s_avail.get("bedgraph") or s_avail.get("cram"):
                                    ui.button(
                                        icon="visibility",
                                        on_click=make_viz_handler(sid),
                                    ).props("flat dense size=sm color=blue").tooltip(
                                        f"Visualize {sid}"
                                    )

            if not members_data:
                ui.label("No pedigree file found for this family.").classes(
                    "text-lg text-amber-700 italic"
                )

            # Track selected members for filtering (default: all selected)
            selected_members: Dict[str, List[str]] = {
                "value": [m["Sample ID"] for m in members_data]
            }
            member_checkboxes: dict = {}
            data_table_refreshers: List = []

            if members_data:
                with ui.row().classes("w-full gap-4 items-start"):
                    with ui.card().classes("flex-1"):
                        # Member selection checkboxes
                        with ui.column().classes("p-4 bg-blue-50"):
                            with ui.row().classes("items-center gap-2 mb-2"):
                                ui.label("Select Members to Display:").classes(
                                    "font-semibold text-blue-800"
                                )

                                def select_all():
                                    selected_members["value"] = [
                                        m["Sample ID"] for m in members_data
                                    ]
                                    for cb in member_checkboxes.values():
                                        cb.value = True
                                    for refresher in data_table_refreshers:
                                        refresher()

                                def select_none():
                                    selected_members["value"] = []
                                    for cb in member_checkboxes.values():
                                        cb.value = False
                                    for refresher in data_table_refreshers:
                                        refresher()

                                ui.button("All", on_click=select_all).props(
                                    "size=sm flat dense"
                                ).classes("text-xs")
                                ui.button("None", on_click=select_none).props(
                                    "size=sm flat dense"
                                ).classes("text-xs")

                            # Native NiceGUI grid layout — no JS DOM manipulation.
                            _GRID_STYLE = (
                                "display: grid;"
                                " grid-template-columns:"
                                " auto auto 1fr 1fr 1fr auto auto;"
                                " gap: 0;"
                                " align-items: center;"
                            )

                            # Header row
                            with (
                                ui.element("div")
                                .classes("w-full bg-blue-100 text-sm")
                                .style(_GRID_STYLE)
                            ):
                                for header in (
                                    "Select",
                                    "",
                                    "Sample ID",
                                    "Father",
                                    "Mother",
                                    "Sex",
                                    "Phenotype",
                                ):
                                    ui.label(header).classes("px-3 py-2 font-semibold")

                            # Data rows
                            for idx, member in enumerate(members_data):
                                sid = member["Sample ID"]
                                bg = "bg-white" if idx % 2 == 0 else "bg-gray-50"

                                def make_change(s):
                                    def handler(e):
                                        if (
                                            e.value
                                            and s not in selected_members["value"]
                                        ):
                                            selected_members["value"].append(s)
                                        elif (
                                            not e.value
                                            and s in selected_members["value"]
                                        ):
                                            selected_members["value"].remove(s)
                                        for refresher in data_table_refreshers:
                                            refresher()

                                    return handler

                                def make_only(s):
                                    def handler():
                                        selected_members["value"] = [s]
                                        for s_id, cb in member_checkboxes.items():
                                            cb.value = s_id == s
                                        for refresher in data_table_refreshers:
                                            refresher()

                                    return handler

                                with (
                                    ui.element("div")
                                    .classes(
                                        f"w-full {bg} border-b border-gray-200 text-sm"
                                    )
                                    .style(_GRID_STYLE)
                                ):
                                    with ui.element("div").classes("px-3 py-2"):
                                        member_checkboxes[sid] = ui.checkbox(
                                            "",
                                            value=True,
                                            on_change=make_change(sid),
                                        )
                                    with ui.element("div").classes("px-3 py-2"):
                                        ui.button(
                                            "only", on_click=make_only(sid)
                                        ).props(
                                            "size=xs flat dense color=blue"
                                        ).classes("text-xs")
                                    ui.label(sid).classes("px-3 py-2 font-medium")
                                    ui.label(member.get("Father", "-")).classes(
                                        "px-3 py-2 text-gray-600"
                                    )
                                    ui.label(member.get("Mother", "-")).classes(
                                        "px-3 py-2 text-gray-600"
                                    )
                                    ui.label(member.get("Sex", "-")).classes(
                                        "px-3 py-2 text-gray-600"
                                    )
                                    ui.label(member.get("Phenotype", "-")).classes(
                                        "px-3 py-2 text-gray-600"
                                    )

                    # ---- Right column: Notes + Diagnostics ----
                    with ui.column().classes("flex-1 gap-4"):
                        # ---- Notes panel ----
                        with ui.card().classes("w-full"):
                            with ui.column().classes("p-4 gap-2"):
                                with ui.row().classes("items-center gap-2"):
                                    ui.label("Notes").classes(
                                        "text-lg font-semibold text-blue-700"
                                    )
                                    ui.icon("info_outline", color="grey").classes(
                                        "text-sm cursor-help"
                                    ).tooltip(
                                        "User-provided information relative to"
                                        " the family or some of its samples"
                                    )

                                @ui.refreshable
                                def render_notes_panel():
                                    notes_file = store.data_dir / "notes" / "notes.tsv"
                                    ensure_notes_file(notes_file)
                                    entries = load_family_notes(notes_file, family_id)

                                    if can_write():
                                        with ui.row().classes("w-full items-end gap-2"):
                                            note_input = (
                                                ui.input(placeholder="Add a note...")
                                                .props("outlined dense")
                                                .classes("flex-1")
                                            )
                                            sample_options = {"": "Family"} | {
                                                m["Sample ID"]: m["Sample ID"]
                                                for m in members_data
                                            }
                                            sample_select = (
                                                ui.select(sample_options, value="")
                                                .props("outlined dense")
                                                .classes("w-36")
                                            )

                                            def add_note():
                                                msg = (
                                                    note_input.value.strip()
                                                    if note_input.value
                                                    else ""
                                                )
                                                if not msg:
                                                    return
                                                save_note(
                                                    notes_file,
                                                    family_id,
                                                    sample_select.value or "",
                                                    msg,
                                                    get_current_user(),
                                                    datetime.now().isoformat(),
                                                )
                                                render_notes_panel.refresh()

                                            ui.button(
                                                icon="add", on_click=add_note
                                            ).props("flat dense color=blue")

                                    if not entries:
                                        ui.label("No notes").classes(
                                            "text-gray-400 text-sm italic"
                                        )
                                    else:
                                        for entry in entries:
                                            with ui.row().classes(
                                                "items-center gap-2 w-full"
                                                " px-2 py-1 border-b"
                                                " border-gray-100"
                                            ):
                                                ui.label(
                                                    entry.get("Message", "")
                                                ).classes("text-sm flex-1")
                                                sample_val = entry.get("Sample", "")
                                                if sample_val:
                                                    ui.badge(sample_val).props(
                                                        "color=blue outline"
                                                    ).classes("text-xs")
                                                ui.label(entry.get("User", "")).classes(
                                                    "text-xs text-gray-400"
                                                )
                                                ts = entry.get("Timestamp", "")
                                                if ts:
                                                    short_ts = ts[:16].replace("T", " ")
                                                    ui.label(short_ts).classes(
                                                        "text-xs text-gray-400"
                                                    )
                                                if can_write():

                                                    def make_del(e_ts, e_user):
                                                        def handler():
                                                            delete_note(
                                                                notes_file,
                                                                family_id,
                                                                e_ts,
                                                                e_user,
                                                            )
                                                            render_notes_panel.refresh()

                                                        return handler

                                                    ui.button(
                                                        icon="close",
                                                        on_click=make_del(
                                                            ts,
                                                            entry.get("User", ""),
                                                        ),
                                                    ).props(
                                                        "flat dense size=xs color=red"
                                                    )

                                render_notes_panel()
                                data_table_refreshers.append(render_notes_panel.refresh)

                        # ---- Diagnostics panel ----
                        with ui.card().classes("w-full"):
                            with ui.column().classes("p-4 gap-2"):
                                ui.label("Diagnostics").classes(
                                    "text-lg font-semibold text-blue-700"
                                )

                            @ui.refreshable
                            def render_diagnostics():
                                snv_file = store.data_dir / "diagnostics" / "snvs.tsv"
                                sv_file = store.data_dir / "diagnostics" / "svs.tsv"
                                ensure_diagnostic_file(snv_file)
                                ensure_diagnostic_file(sv_file)
                                entries = load_family_diagnostics(
                                    snv_file,
                                    sv_file,
                                    family_id,
                                    selected_members["value"],
                                )
                                if not entries:
                                    ui.label(
                                        "No diagnostics recorded for selected members"
                                    ).classes("text-gray-500 text-sm italic")
                                    return

                                _DIAG_COLORS: Dict[str, str] = {
                                    "pathogenic": "red",
                                    "uncertain": "orange",
                                    "benign": "green",
                                    "conflicting": "amber",
                                }
                                gene_scorer = get_gene_scorer()
                                for entry in entries:
                                    diag = entry.get("Diagnostic", "")
                                    color = _DIAG_COLORS.get(diag, "grey")
                                    is_sv = entry.get("_source") == "sv"
                                    with ui.row().classes(
                                        "items-center gap-2 w-full px-2 py-1 border-b border-gray-100"
                                    ):
                                        ui.badge(
                                            diag.upper()[:3] if diag else "?",
                                            color=color,
                                        ).classes("text-xs")
                                        variant = entry.get("Variant", "")
                                        if is_sv and variant:
                                            ui.label(variant).classes(
                                                "text-xs font-mono"
                                            )
                                            cytoband = _parse_sv_cytoband(variant)
                                            if cytoband:
                                                ui.label(cytoband).classes(
                                                    "text-xs text-purple-600 italic"
                                                )
                                        else:
                                            ui.label(variant).classes(
                                                "text-xs font-mono"
                                            )
                                        gene_str = entry.get("Gene", "")
                                        if gene_str and is_sv:
                                            _render_sv_gene_badges(
                                                gene_str, gene_scorer
                                            )
                                        elif gene_str:
                                            ui.label(gene_str).classes(
                                                "text-xs text-blue-700 font-medium"
                                            )
                                        ui.space()
                                        ui.label(entry.get("Sample", "")).classes(
                                            "text-xs text-gray-600"
                                        )
                                        ui.label(entry.get("User", "")).classes(
                                            "text-xs text-gray-400"
                                        )

                            render_diagnostics()
                            data_table_refreshers.append(render_diagnostics.refresh)

            # Analysis tabs
            wombat_state = {"loaded": False}
            svs_state = {"loaded": False}

            with ui.tabs().classes("w-full") as tabs:
                wombat_tab = ui.tab("Wombat")
                svs_tab = ui.tab("SVs")

            with ui.tab_panels(tabs, value=wombat_tab).classes("w-full"):
                with ui.tab_panel(wombat_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def wombat_content():
                        if wombat_state["loaded"]:
                            render_wombat_tab(
                                store=store,
                                family_id=family_id,
                                cohort_name="",
                                selected_members=selected_members,
                                data_table_refreshers=data_table_refreshers,
                            )
                        else:
                            with ui.column().classes(
                                "w-full items-center justify-center py-16"
                            ):
                                ui.spinner(size="xl", color="blue")
                                ui.label("Loading variants...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    wombat_content()

                with ui.tab_panel(svs_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def svs_content():
                        if svs_state["loaded"]:
                            render_svs_tab(
                                store=store,
                                family_id=family_id,
                                cohort_name="",
                                selected_members=selected_members,
                                data_table_refreshers=data_table_refreshers,
                            )
                        else:
                            with ui.column().classes(
                                "w-full items-center justify-center py-16"
                            ):
                                ui.spinner(size="xl", color="blue")
                                ui.label("Loading structural variants...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    svs_content()

            def load_wombat():
                wombat_state["loaded"] = True
                wombat_content.refresh()

            ui.timer(0.1, load_wombat, once=True)

            def on_tab_change(e):
                tab_value = e.args
                if tab_value == "SVs" and not svs_state["loaded"]:

                    def load_svs():
                        svs_state["loaded"] = True
                        svs_content.refresh()

                    ui.timer(0.1, load_svs, once=True)

            tabs.on("update:model-value", on_tab_change)

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
