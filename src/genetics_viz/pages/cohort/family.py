"""Family detail page - displays family members and analysis tabs."""

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
from genetics_viz.pages.cohort.components.svs_tab import probe_svs_data, render_svs_tab
from genetics_viz.utils.auth import can_write, check_auth, get_current_user
from genetics_viz.pages.cohort.components.wombat_tab import (
    probe_wombat_data,
    render_wombat_tab,
)
from genetics_viz.utils.cytobands import get_cytoband_range
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.gene_scoring import GeneScorer, get_gene_scorer


def _parse_sv_cytoband(variant: str) -> str:
    """Extract cytoband range from an SV variant key like 'chr1:1000-2000:dup'."""
    try:
        # Strip trailing :dup/:del if present
        parts = variant.split(":")
        if len(parts) >= 2:
            chrom = parts[0]
            range_part = parts[1]
            if "-" in range_part:
                start_str, end_str = range_part.split("-", 1)
                return get_cytoband_range(chrom, int(start_str), int(end_str))
    except (ValueError, IndexError):
        pass
    return ""


def _render_sv_gene_badges(gene_str: str, gene_scorer: GeneScorer) -> None:
    """Render gene badges for an SV entry, sorted by descending geneset score.

    Shows up to 6 genes as colored badges, then a "+X genes" overflow label.
    Mirrors the pattern used in svs_tab.py for DataTable gene badges.
    """
    genes = []
    for gene_part in str(gene_str).split(","):
        gene_part = gene_part.strip()
        if not gene_part or gene_part == "-":
            continue
        # Handle "SYMBOL:type" format
        symbol = gene_part.split(":")[0].strip()
        if symbol:
            score, _ = gene_scorer.get_gene_score_and_sets(symbol)
            color = gene_scorer.get_gene_color(symbol)
            tooltip = gene_scorer.get_gene_tooltip(symbol)
            genes.append(
                {"symbol": symbol, "score": score, "color": color, "tooltip": tooltip}
            )

    # Sort by score descending
    genes.sort(key=lambda g: g["score"], reverse=True)

    total = len(genes)
    display_genes = genes[:6]

    for g in display_genes:
        ui.badge(g["symbol"], color=g["color"]).classes("text-xs").tooltip(g["tooltip"])

    if total > 6:
        remaining = total - 6
        ui.label(f"+{remaining} genes").classes("text-xs text-gray-500 italic")


@ui.page("/cohort/{cohort_name}/family/{family_id}")
def family_page(cohort_name: str, family_id: str) -> None:
    """Render the family detail page."""
    if redirect := check_auth():
        return redirect
    create_header(cohort_name)

    # Add IGV.js library at page level
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.13/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()

        cohort = store.get_cohort(cohort_name)

        if cohort is None:
            with ui.column().classes("w-full px-6 py-6"):
                ui.label(f"Cohort not found: {cohort_name}").classes(
                    "text-xl text-red-500"
                )
                ui.button("← Back to Home", on_click=lambda: ui.navigate.to("/"))
            return

        family = cohort.families.get(family_id)
        if family is None:
            with ui.column().classes("w-full px-6 py-6"):
                ui.label(f"Family not found: {family_id}").classes(
                    "text-xl text-red-500"
                )
                ui.button(
                    "← Back to Cohort",
                    on_click=lambda: ui.navigate.to(f"/cohort/{cohort_name}"),
                )
            return

        with ui.column().classes("w-full px-6 py-6"):
            # Breadcrumb navigation
            with ui.row().classes("items-center gap-2 mb-4"):
                ui.link("Home", "/").classes("text-blue-600 hover:text-blue-800")
                ui.label("/").classes("text-gray-400")
                ui.link(cohort_name, f"/cohort/{cohort_name}").classes(
                    "text-blue-600 hover:text-blue-800"
                )
                ui.label("/").classes("text-gray-400")
                ui.label(family_id).classes("font-semibold")

            # Family header with prev/next navigation
            sorted_fids = sorted(cohort.families.keys())
            current_idx = (
                sorted_fids.index(family_id) if family_id in sorted_fids else -1
            )

            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"👨‍👩‍👧‍👦 Family: {family_id}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{family.num_samples} members").props("color=blue")
                ui.badge(f"{family.num_founders} founders").props("color=teal")

                if len(sorted_fids) > 1 and current_idx >= 0:
                    ui.space()
                    if current_idx > 0:
                        prev_fid = sorted_fids[current_idx - 1]
                        ui.button(
                            icon="arrow_back",
                            on_click=lambda _, f=prev_fid: ui.navigate.to(
                                f"/cohort/{cohort_name}/family/{f}"
                            ),
                        ).props("flat round color=blue").tooltip(
                            f"Previous: {prev_fid}"
                        )
                    if current_idx < len(sorted_fids) - 1:
                        next_fid = sorted_fids[current_idx + 1]
                        ui.button(
                            icon="arrow_forward",
                            on_click=lambda _, f=next_fid: ui.navigate.to(
                                f"/cohort/{cohort_name}/family/{f}"
                            ),
                        ).props("flat round color=blue").tooltip(f"Next: {next_fid}")

            members_data = cohort.get_family_members(family_id)

            # Track selected members for filtering (default: all selected)
            selected_members = {"value": [m["Sample ID"] for m in members_data]}
            member_checkboxes = {}

            # Store refresh functions for all data tables
            data_table_refreshers: List = []

            with ui.row().classes("w-full gap-4 items-start"):
                with ui.card().classes("flex-1"):
                    # Member selection checkboxes
                    with ui.column().classes("p-4 bg-blue-50"):
                        with ui.row().classes("items-center gap-2 mb-2"):
                            ui.label("Select Members to Display:").classes(
                                "font-semibold text-blue-800"
                            )

                            def select_all_members():
                                selected_members["value"] = [
                                    m["Sample ID"] for m in members_data
                                ]
                                for cb in member_checkboxes.values():
                                    cb.value = True
                                for refresher in data_table_refreshers:
                                    refresher()

                            def select_none_members():
                                selected_members["value"] = []
                                for cb in member_checkboxes.values():
                                    cb.value = False
                                for refresher in data_table_refreshers:
                                    refresher()

                            ui.button("All", on_click=select_all_members).props(
                                "size=sm flat dense"
                            ).classes("text-xs")
                            ui.button("None", on_click=select_none_members).props(
                                "size=sm flat dense"
                            ).classes("text-xs")

                        # Native NiceGUI grid layout — no JS DOM manipulation.
                        # Each row is a CSS grid container; checkboxes and
                        # buttons are placed directly in their cells.
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
                            sample_id = member["Sample ID"]
                            bg_class = "bg-white" if idx % 2 == 0 else "bg-gray-50"

                            def make_change_handler(sid):
                                def handler(e):
                                    if e.value and sid not in selected_members["value"]:
                                        selected_members["value"].append(sid)
                                    elif (
                                        not e.value and sid in selected_members["value"]
                                    ):
                                        selected_members["value"].remove(sid)
                                    for refresher in data_table_refreshers:
                                        refresher()

                                return handler

                            def make_only_handler(sid):
                                def handler():
                                    selected_members["value"] = [sid]
                                    for s_id, checkbox in member_checkboxes.items():
                                        checkbox.value = s_id == sid
                                    for refresher in data_table_refreshers:
                                        refresher()

                                return handler

                            with (
                                ui.element("div")
                                .classes(
                                    f"w-full {bg_class} border-b border-gray-200 text-sm"
                                )
                                .style(_GRID_STYLE)
                            ):
                                # Select column (checkbox)
                                with ui.element("div").classes("px-3 py-2"):
                                    member_checkboxes[sample_id] = ui.checkbox(
                                        "",
                                        value=True,
                                        on_change=make_change_handler(sample_id),
                                    )
                                # "only" button column
                                with ui.element("div").classes("px-3 py-2"):
                                    ui.button(
                                        "only", on_click=make_only_handler(sample_id)
                                    ).props("size=xs flat dense color=blue").classes(
                                        "text-xs"
                                    )
                                # Data columns
                                ui.label(sample_id).classes("px-3 py-2 font-medium")
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

                                        ui.button(icon="add", on_click=add_note).props(
                                            "flat dense color=blue"
                                        )

                                if not entries:
                                    ui.label("No notes").classes(
                                        "text-gray-400 text-sm italic"
                                    )
                                else:
                                    for entry in entries:
                                        with ui.row().classes(
                                            "items-center gap-2 w-full px-2"
                                            " py-1 border-b border-gray-100"
                                        ):
                                            ui.label(entry.get("Message", "")).classes(
                                                "text-sm flex-1"
                                            )
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
                                                ).props("flat dense size=xs color=red")

                            render_notes_panel()
                            data_table_refreshers.append(render_notes_panel.refresh)

                    # ---- Diagnostics panel ----
                    with ui.card().classes("w-full"):
                        with ui.column().classes("p-4 gap-2"):
                            ui.label("Diagnostics").classes(
                                "text-lg font-semibold text-blue-700"
                            )

                        @ui.refreshable
                        def render_diagnostics_panel():
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

                            # Build display rows
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

                                    # Variant + cytoband for SVs
                                    variant = entry.get("Variant", "")
                                    if is_sv and variant:
                                        ui.label(variant).classes("text-xs font-mono")
                                        cytoband = _parse_sv_cytoband(variant)
                                        if cytoband:
                                            ui.label(cytoband).classes(
                                                "text-xs text-purple-600 italic"
                                            )
                                    else:
                                        ui.label(variant).classes("text-xs font-mono")

                                    # Gene display: badges for SVs, plain label for SNVs
                                    gene_str = entry.get("Gene", "")
                                    if gene_str and is_sv:
                                        _render_sv_gene_badges(gene_str, gene_scorer)
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

                        render_diagnostics_panel()
                        data_table_refreshers.append(render_diagnostics_panel.refresh)

            # Pre-check data existence for each analysis tab
            has_wombat = probe_wombat_data(store.data_dir, family_id)
            has_svs = probe_svs_data(store.data_dir, family_id)

            # Track loading state for each tab
            wombat_state = {"loaded": False}
            svs_state = {"loaded": False}

            # Pick the first enabled tab as the default
            default_tab_name = (
                "wombat" if has_wombat else "svs" if has_svs else "wombat"
            )

            # Analysis tabs section
            with ui.tabs().classes("w-full") as tabs:
                wombat_tab = ui.tab("Wombat")
                svs_tab = ui.tab("SVs")
                if not has_wombat:
                    wombat_tab.props("disable")
                if not has_svs:
                    svs_tab.props("disable")

            default_tab = wombat_tab if default_tab_name == "wombat" else svs_tab

            with ui.tab_panels(tabs, value=default_tab).classes("w-full"):
                # Wombat tab panel
                with ui.tab_panel(wombat_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def wombat_content():
                        if not has_wombat:
                            ui.label("No Wombat data available").classes(
                                "text-gray-500 italic"
                            )
                        elif wombat_state["loaded"]:
                            render_wombat_tab(
                                store=store,
                                family_id=family_id,
                                cohort_name=cohort_name,
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

                # SVs tab panel
                with ui.tab_panel(svs_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def svs_content():
                        if not has_svs:
                            ui.label("No SVs data available").classes(
                                "text-gray-500 italic"
                            )
                        elif svs_state["loaded"]:
                            render_svs_tab(
                                store=store,
                                family_id=family_id,
                                cohort_name=cohort_name,
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

            # Load the default tab data asynchronously
            if has_wombat:

                def load_wombat_async():
                    wombat_state["loaded"] = True
                    wombat_content.refresh()

                ui.timer(0.1, load_wombat_async, once=True)
            elif has_svs:
                # SVs is the default tab when Wombat has no data
                def load_svs_default():
                    svs_state["loaded"] = True
                    svs_content.refresh()

                ui.timer(0.1, load_svs_default, once=True)

            # Lazy load tabs when clicked
            def on_tab_change(e):
                tab_value = e.args
                if tab_value == "Wombat" and not wombat_state["loaded"] and has_wombat:

                    def load_wombat_lazy():
                        wombat_state["loaded"] = True
                        wombat_content.refresh()

                    ui.timer(0.1, load_wombat_lazy, once=True)
                elif tab_value == "SVs" and not svs_state["loaded"] and has_svs:

                    def load_svs_lazy():
                        svs_state["loaded"] = True
                        svs_content.refresh()

                    ui.timer(0.1, load_svs_lazy, once=True)

            tabs.on("update:model-value", on_tab_change)

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
