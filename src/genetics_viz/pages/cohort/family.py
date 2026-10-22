"""Family detail page - displays family members and analysis tabs."""

from datetime import datetime
from typing import Any, Dict, List

from fastapi import Request
from nicegui import ui

from genetics_viz.components.diagnostics_panel import render_family_diagnostics
from genetics_viz.components.family_report_dialog import render_report_buttons
from genetics_viz.components.header import create_header
from genetics_viz.components.page_link import PageLink, link_members
from genetics_viz.components.notes_loader import (
    delete_note,
    ensure_notes_file,
    load_family_notes,
    save_note,
)
from genetics_viz.pages.cohort.components.ancestry_tab import (
    probe_ancestry_data,
    render_ancestry_tab,
)
from genetics_viz.pages.cohort.components.member_selector import (
    render_member_selector,
)
from genetics_viz.pages.cohort.components.pgs_tab import (
    probe_pgs_data,
    render_pgs_tab,
)
from genetics_viz.pages.cohort.components.roh_tab import (
    probe_roh_data,
    render_roh_tab,
)
from genetics_viz.pages.cohort.components.sample_qc_tab import (
    family_qc_summary,
    probe_qc_data,
    render_sample_qc_tab,
)
from genetics_viz.pages.cohort.components.str_tab import probe_eh_data, render_str_tab
from genetics_viz.pages.cohort.components.svs_tab import probe_svs_data, render_svs_tab
from genetics_viz.utils.auth import can_write, check_auth, get_current_user
from genetics_viz.pages.cohort.components.wombat_tab import (
    probe_wombat_data,
    render_wombat_tab,
)
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.sample_qc import (
    count_by_severity,
    flags_by_sample,
    severity_info,
    severity_summary,
)


@ui.page("/cohort/{cohort_name}/family/{family_id}")
def family_page(cohort_name: str, family_id: str, request: Request) -> None:
    """Render the family detail page.

    The open tab, the selected members and the Wombat and SV filters are in
    the URL query (``components/page_link.py``), so the address is a link to
    this view.
    """
    if redirect := check_auth():
        return redirect
    link = PageLink.from_request(request)
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

            # Sample QC: the family's own flags, read now (a few KB) for the
            # header badge and the member icons; cross-family ones need the
            # cohort's tables and are left to the tab. Flags a curator
            # acknowledged do not count; a review redraws both.
            qc = {"summary": family_qc_summary(store, family_id)}
            # Filled in once the analysis tabs exist: the badge opens one.
            open_qc_tab: Dict[str, Any] = {"fn": None}

            with ui.row().classes("w-full items-center gap-4 mb-6"):
                ui.label(f"👨‍👩‍👧‍👦 Family: {family_id}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{family.num_samples} members").props("color=blue")
                ui.badge(f"{family.num_founders} founders").props("color=teal")

                @ui.refreshable
                def qc_badge() -> None:
                    qc_data, active, acked = qc["summary"]
                    if qc_data is None:
                        return
                    counts = count_by_severity(active)
                    summary = severity_summary(counts)
                    worst = "high" if counts.get("high") else "medium"
                    color = severity_info(worst)["color"] if summary else "#2e7d32"
                    tip = "Sample QC (somalier): open the Sample QC tab"
                    if acked:
                        tip += f". {len(acked)} acknowledged flag(s) not counted."
                    ui.badge(
                        f"QC: {summary}" if summary else "QC ok", color=color
                    ).classes("cursor-pointer").on(
                        "click",
                        lambda: open_qc_tab["fn"] and open_qc_tab["fn"](),
                    ).tooltip(tip)

                qc_badge()

                ui.space()
                render_report_buttons(store.data_dir, family_id, family, cohort)

                if len(sorted_fids) > 1 and current_idx >= 0:
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

            # Store refresh functions for all data tables
            data_table_refreshers: List = []
            link_members(
                link,
                selected_members,
                [m["Sample ID"] for m in members_data],
                data_table_refreshers,
            )

            with ui.row().classes("w-full gap-4 items-start"):
                _, qc_active, qc_acked = qc["summary"]
                update_qc_icons = render_member_selector(
                    members_data=members_data,
                    selected_members=selected_members,
                    data_table_refreshers=data_table_refreshers,
                    qc_flags=(
                        flags_by_sample(qc_active)
                        if qc["summary"][0] is not None
                        else None
                    ),
                    show_qc_detail=can_write(),
                    qc_acked=flags_by_sample(qc_acked),
                )

                # After a QC review or a UPD validation/diagnostic (ROH and
                # Sample QC tabs): the header badge, the member icons and the
                # Diagnostics panel follow.
                diagnostics_refresh: Dict[str, Any] = {"fn": None}

                def on_acks_changed() -> None:
                    qc["summary"] = family_qc_summary(store, family_id)
                    _, active, acked = qc["summary"]
                    qc_badge.refresh()
                    update_qc_icons(flags_by_sample(active), flags_by_sample(acked))
                    if diagnostics_refresh["fn"] is not None:
                        diagnostics_refresh["fn"]()

                # ---- Right column: tabbed info panel ----
                with ui.column().classes("flex-1").style("gap: 0.5rem; min-width: 0"):
                    with ui.tabs().classes("w-full") as info_tabs:
                        general_tab = ui.tab("General")

                    with ui.tab_panels(info_tabs, value=general_tab).classes("w-full"):
                        with ui.tab_panel(general_tab).classes(
                            "w-full border border-gray-300 rounded-lg p-4"
                        ):
                            # ---- Notes section ----
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

                            # ---- Diagnostics section ----
                            ui.label("Diagnostics").classes(
                                "text-lg font-semibold text-blue-700"
                            )

                            @ui.refreshable
                            def render_diagnostics_panel():
                                render_family_diagnostics(
                                    store.data_dir,
                                    family_id,
                                    selected_members["value"],
                                )

                            render_diagnostics_panel()
                            diagnostics_refresh["fn"] = render_diagnostics_panel.refresh
                            data_table_refreshers.append(
                                render_diagnostics_panel.refresh
                            )

            # Pre-check data existence for each analysis tab
            has_wombat = probe_wombat_data(store.data_dir, family_id)
            has_svs = probe_svs_data(store.data_dir, family_id)
            has_ancestry = probe_ancestry_data(store.data_dir, family_id)
            has_pgs = probe_pgs_data(store.data_dir, family_id)
            has_roh = probe_roh_data(store.data_dir, family_id)
            has_str = probe_eh_data(store.data_dir, family_id)
            has_qc = probe_qc_data(store.data_dir, family_id)

            # Track loading state for each tab
            wombat_state = {"loaded": False}
            svs_state = {"loaded": False}
            ancestry_state = {"loaded": False}
            pgs_state = {"loaded": False}
            roh_state = {"loaded": False}
            str_state = {"loaded": False}
            qc_state = {"loaded": False}

            # Analysis tabs section
            with ui.tabs().classes("w-full").mark("family-tabs") as tabs:
                wombat_tab = ui.tab("Wombat")
                svs_tab = ui.tab("SVs")
                ancestry_tab = ui.tab("Ancestry")
                pgs_tab = ui.tab("Polygenic Scores")
                roh_tab = ui.tab("ROH")
                str_tab = ui.tab("STR").tooltip("Short tandem repeats")
                qc_tab = ui.tab("Sample QC").tooltip(
                    "Relatedness, sex and contamination checks (somalier)"
                )
                if not has_wombat:
                    wombat_tab.props("disable")
                if not has_svs:
                    svs_tab.props("disable")
                if not has_ancestry:
                    ancestry_tab.props("disable")
                if not has_pgs:
                    pgs_tab.props("disable")
                if not has_roh:
                    roh_tab.props("disable")
                if not has_str:
                    str_tab.props("disable")
                if not has_qc:
                    qc_tab.props("disable")

            # Pick the first enabled tab as the default
            if has_wombat:
                default_tab_name = "wombat"
            elif has_svs:
                default_tab_name = "svs"
            elif has_ancestry:
                default_tab_name = "ancestry"
            elif has_pgs:
                default_tab_name = "pgs"
            elif has_roh:
                default_tab_name = "roh"
            elif has_str:
                default_tab_name = "str"
            elif has_qc:
                default_tab_name = "qc"
            else:
                default_tab_name = "wombat"
            first_tab_name = default_tab_name

            # A shared link opens its tab, when the family has its data
            tab_has_data = {
                "wombat": has_wombat,
                "svs": has_svs,
                "ancestry": has_ancestry,
                "pgs": has_pgs,
                "roh": has_roh,
                "str": has_str,
                "qc": has_qc,
            }
            if tab_has_data.get(link.get("tab") or ""):
                default_tab_name = link.get("tab")

            default_tab = {
                "wombat": wombat_tab,
                "svs": svs_tab,
                "ancestry": ancestry_tab,
                "pgs": pgs_tab,
                "roh": roh_tab,
                "str": str_tab,
                "qc": qc_tab,
            }[default_tab_name]

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
                                family=family,
                                on_open_qc=(
                                    (lambda: open_qc_tab["fn"] and open_qc_tab["fn"]())
                                    if has_qc
                                    else None
                                ),
                                link=link,
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
                                link=link,
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

                # Ancestry tab panel
                with ui.tab_panel(ancestry_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def ancestry_content():
                        if not has_ancestry:
                            ui.label("No ancestry data available").classes(
                                "text-gray-500 italic"
                            )
                        elif ancestry_state["loaded"]:
                            render_ancestry_tab(
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
                                ui.label("Loading ancestry...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    ancestry_content()

                # Polygenic scores tab panel
                with ui.tab_panel(pgs_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def pgs_content():
                        if not has_pgs:
                            ui.label("No polygenic score data available").classes(
                                "text-gray-500 italic"
                            )
                        elif pgs_state["loaded"]:
                            render_pgs_tab(
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
                                ui.label("Loading polygenic scores...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    pgs_content()

                # Runs of homozygosity tab panel
                with ui.tab_panel(roh_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def roh_content():
                        if not has_roh:
                            ui.label("No ROH data available").classes(
                                "text-gray-500 italic"
                            )
                        elif roh_state["loaded"]:
                            render_roh_tab(
                                store=store,
                                family_id=family_id,
                                cohort_name=cohort_name,
                                selected_members=selected_members,
                                data_table_refreshers=data_table_refreshers,
                                on_acks_changed=on_acks_changed,
                            )
                        else:
                            with ui.column().classes(
                                "w-full items-center justify-center py-16"
                            ):
                                ui.spinner(size="xl", color="blue")
                                ui.label("Loading runs of homozygosity...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    roh_content()

                # Short tandem repeats tab panel
                with ui.tab_panel(str_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def str_content():
                        if not has_str:
                            ui.label("No STR data available").classes(
                                "text-gray-500 italic"
                            )
                        elif str_state["loaded"]:
                            render_str_tab(
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
                                ui.label("Loading short tandem repeats...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    str_content()

                # Sample QC tab panel
                with ui.tab_panel(qc_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):

                    @ui.refreshable
                    def qc_content():
                        if not has_qc:
                            ui.label("No sample QC data available").classes(
                                "text-gray-500 italic"
                            )
                        elif qc_state["loaded"]:
                            render_sample_qc_tab(
                                store=store,
                                family_id=family_id,
                                cohort_name=cohort_name,
                                selected_members=selected_members,
                                data_table_refreshers=data_table_refreshers,
                                on_acks_changed=on_acks_changed,
                            )
                        else:
                            with ui.column().classes(
                                "w-full items-center justify-center py-16"
                            ):
                                ui.spinner(size="xl", color="blue")
                                ui.label("Loading sample QC...").classes(
                                    "text-lg text-gray-600 mt-4"
                                )

                    qc_content()

            # Each tab loads its data on first view: keyed by the tab name that
            # the "update:model-value" event reports.
            lazy_tabs = {
                "Wombat": (has_wombat, wombat_state, wombat_content),
                "SVs": (has_svs, svs_state, svs_content),
                "Ancestry": (has_ancestry, ancestry_state, ancestry_content),
                "Polygenic Scores": (has_pgs, pgs_state, pgs_content),
                "ROH": (has_roh, roh_state, roh_content),
                "STR": (has_str, str_state, str_content),
                "Sample QC": (has_qc, qc_state, qc_content),
            }
            tab_names = {
                "wombat": "Wombat",
                "svs": "SVs",
                "ancestry": "Ancestry",
                "pgs": "Polygenic Scores",
                "roh": "ROH",
                "str": "STR",
                "qc": "Sample QC",
            }

            def make_loader(state, content):
                def load():
                    state["loaded"] = True
                    content.refresh()

                return load

            # Load the default tab data asynchronously
            default_available, default_state, default_content = lazy_tabs[
                tab_names[default_tab_name]
            ]
            if default_available:
                ui.timer(0.1, make_loader(default_state, default_content), once=True)

            # Lazy load tabs when clicked
            tab_slugs = {label: slug for slug, label in tab_names.items()}

            def link_tab(label: str) -> None:
                link.set("tab", tab_slugs.get(label), default=first_tab_name)
                link.push()

            def on_tab_change(e):
                entry = lazy_tabs.get(e.args)
                if entry is None:
                    return
                link_tab(e.args)
                available, state, content = entry
                if available and not state["loaded"]:
                    ui.timer(0.1, make_loader(state, content), once=True)

            tabs.on("update:model-value", on_tab_change)

            def show_qc_tab() -> None:
                if not has_qc:
                    return
                tabs.set_value(qc_tab)
                link_tab("Sample QC")
                if not qc_state["loaded"]:
                    ui.timer(0.1, make_loader(qc_state, qc_content), once=True)

            open_qc_tab["fn"] = show_qc_tab

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
