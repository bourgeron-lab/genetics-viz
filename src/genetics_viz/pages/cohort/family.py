"""Family detail page - displays family members and analysis tabs."""

from typing import List

from nicegui import app as nicegui_app
from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.pages.cohort.components.dnm_tab import render_dnm_tab
from genetics_viz.pages.cohort.components.wombat_tab import render_wombat_tab
from genetics_viz.utils.data import get_data_store


@ui.page("/cohort/{cohort_name}/family/{family_id}")
def family_page(cohort_name: str, family_id: str) -> None:
    """Render the family detail page."""
    create_header()

    # Add IGV.js library at page level
    ui.add_head_html("""
        <script src="https://cdn.jsdelivr.net/npm/igv@2.15.11/dist/igv.min.js"></script>
    """)

    try:
        store = get_data_store()

        # Serve data files for IGV.js
        nicegui_app.add_static_files("/data", str(store.data_dir))

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

            # Family header
            with ui.row().classes("items-center gap-4 mb-6"):
                ui.label(f"👨‍👩‍👧‍👦 Family: {family_id}").classes(
                    "text-3xl font-bold text-blue-900"
                )
                ui.badge(f"{family.num_samples} members").props("color=blue")
                ui.badge(f"{family.num_founders} founders").props("color=teal")

            members_data = cohort.get_family_members(family_id)

            # Track selected members for filtering (default: all selected)
            selected_members = {"value": [m["Sample ID"] for m in members_data]}
            member_checkboxes = {}

            # Store refresh functions for all data tables
            data_table_refreshers: List = []

            with ui.card().classes("w-full"):
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

                    # Create HTML table with checkboxes and member info
                    table_html = """
                    <table class="w-full text-sm">
                        <thead class="bg-blue-100">
                            <tr>
                                <th class="px-3 py-2 text-left font-semibold">Select</th>
                                <th class="px-3 py-2 text-left font-semibold"></th>
                                <th class="px-3 py-2 text-left font-semibold">Sample ID</th>
                                <th class="px-3 py-2 text-left font-semibold">Father</th>
                                <th class="px-3 py-2 text-left font-semibold">Mother</th>
                                <th class="px-3 py-2 text-left font-semibold">Sex</th>
                                <th class="px-3 py-2 text-left font-semibold">Phenotype</th>
                            </tr>
                        </thead>
                        <tbody>
                    """

                    for idx, member in enumerate(members_data):
                        sample_id = member["Sample ID"]
                        bg_class = "bg-white" if idx % 2 == 0 else "bg-gray-50"
                        table_html += f'''
                            <tr class="{bg_class} border-b border-gray-200">
                                <td class="px-3 py-2" id="checkbox-cell-{idx}"></td>
                                <td class="px-3 py-2" id="only-button-cell-{idx}"></td>
                                <td class="px-3 py-2 font-medium">{sample_id}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Father", "-")}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Mother", "-")}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Sex", "-")}</td>
                                <td class="px-3 py-2 text-gray-600">{member.get("Phenotype", "-")}</td>
                            </tr>
                        '''

                    table_html += """
                        </tbody>
                    </table>
                    """

                    ui.html(table_html, sanitize=False)

                    # Create checkboxes and insert them into the table cells
                    for idx, member in enumerate(members_data):
                        sample_id = member["Sample ID"]

                        def make_change_handler(sid):
                            def handler(e):
                                if e.value and sid not in selected_members["value"]:
                                    selected_members["value"].append(sid)
                                elif not e.value and sid in selected_members["value"]:
                                    selected_members["value"].remove(sid)
                                for refresher in data_table_refreshers:
                                    refresher()

                            return handler

                        with ui.element().classes(f"checkbox-cell-{idx}"):
                            member_checkboxes[sample_id] = ui.checkbox(
                                "",
                                value=True,
                                on_change=make_change_handler(sample_id),
                            )

                        def make_only_handler(sid):
                            def handler():
                                selected_members["value"] = [sid]
                                for s_id, checkbox in member_checkboxes.items():
                                    checkbox.value = s_id == sid
                                for refresher in data_table_refreshers:
                                    refresher()

                            return handler

                        with ui.element().classes(f"only-button-cell-{idx}"):
                            ui.button(
                                "only", on_click=make_only_handler(sample_id)
                            ).props("size=xs flat dense color=blue").classes("text-xs")

                    # Move checkboxes and only buttons into table cells using JavaScript
                    ui.run_javascript(f"""
                        for (let i = 0; i < {len(members_data)}; i++) {{
                            const checkbox = document.querySelector('.checkbox-cell-' + i);
                            const checkboxCell = document.getElementById('checkbox-cell-' + i);
                            if (checkbox && checkboxCell) {{
                                checkboxCell.appendChild(checkbox);
                            }}
                            
                            const onlyButton = document.querySelector('.only-button-cell-' + i);
                            const onlyButtonCell = document.getElementById('only-button-cell-' + i);
                            if (onlyButton && onlyButtonCell) {{
                                onlyButtonCell.appendChild(onlyButton);
                            }}
                        }}
                    """)

            # Analysis tabs section
            with ui.tabs().classes("w-full") as tabs:
                wombat_tab = ui.tab("Wombat")
                snvs_tab = ui.tab("SNVs dnm")

            with ui.tab_panels(tabs, value=wombat_tab).classes("w-full"):
                # Wombat tab panel
                with ui.tab_panel(wombat_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):
                    render_wombat_tab(
                        store=store,
                        family_id=family_id,
                        cohort_name=cohort_name,
                        selected_members=selected_members,
                        data_table_refreshers=data_table_refreshers,
                    )

                # SNVs dnm tab panel
                with ui.tab_panel(snvs_tab).classes(
                    "border border-gray-300 rounded-lg p-4"
                ):
                    render_dnm_tab(
                        store=store,
                        family_id=family_id,
                        cohort_name=cohort_name,
                        selected_members=selected_members,
                        data_table_refreshers=data_table_refreshers,
                    )

    except Exception as e:
        import traceback

        with ui.column().classes("w-full px-6 py-6"):
            ui.label(f"Error: {e}").classes("text-red-500 text-xl mb-4")
            ui.label("Traceback:").classes("text-red-500 font-semibold")
            ui.label(traceback.format_exc()).classes(
                "text-red-500 text-xs font-mono whitespace-pre"
            )
