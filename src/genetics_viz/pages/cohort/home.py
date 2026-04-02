"""Home page - displays available cohorts and quick search."""

import asyncio

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.components.sample_dialog import show_sample_dialog
from genetics_viz.utils.auth import check_auth
from genetics_viz.utils.data import get_data_store
from genetics_viz.utils.data_availability import (
    check_family_availability,
    check_sample_availability,
)
from genetics_viz.utils.pedigree import load_family_pedigree
from genetics_viz.utils.sharding import get_family_path, get_sample_path


@ui.page("/")
async def home_page() -> None:
    """Render the home/welcome page."""
    if redirect := check_auth():
        return redirect
    create_header()

    with ui.column().classes("w-full max-w-6xl mx-auto p-6"):
        ui.label("Welcome to Genetics-Viz").classes(
            "text-3xl font-bold mb-2 text-blue-900"
        )
        ui.label("Select a cohort to explore").classes("text-lg text-gray-600 mb-6")

        try:
            store = get_data_store()

            # Quick Search section
            with ui.card().classes("w-full mb-6 bg-gray-50"):
                with ui.column().classes("p-4 gap-3"):
                    with ui.row().classes("items-center gap-2"):
                        ui.icon("search", color="blue").classes("text-xl")
                        ui.label("Quick Search").classes(
                            "text-lg font-semibold text-blue-800"
                        )

                    search_result = ui.column().classes("w-full")

                    with ui.row().classes("items-end gap-4 w-full"):
                        sample_input = (
                            ui.input(
                                label="Sample Barcode",
                                placeholder="e.g. B00EYRL",
                            )
                            .props("outlined dense")
                            .classes("w-64")
                        )

                        async def search_sample():
                            barcode = (
                                sample_input.value.strip() if sample_input.value else ""
                            )
                            if not barcode:
                                ui.notify("Enter a barcode", type="warning")
                                return
                            search_result.clear()
                            with search_result:
                                with ui.row().classes("items-center gap-2"):
                                    ui.spinner(size="sm")
                                    ui.label("Searching...").classes("text-gray-500")

                            sample_path = get_sample_path(store.data_dir, barcode)
                            exists = await asyncio.to_thread(sample_path.is_dir)

                            search_result.clear()
                            with search_result:
                                if exists:
                                    avail = await asyncio.to_thread(
                                        check_sample_availability,
                                        store.data_dir,
                                        barcode,
                                    )
                                    with ui.card().classes(
                                        "w-full border-l-4 border-green-500"
                                    ):
                                        with ui.row().classes("items-center gap-3 p-3"):
                                            ui.icon(
                                                "check_circle", color="green"
                                            ).classes("text-xl")
                                            ui.label(
                                                f"Sample found: {barcode}"
                                            ).classes("font-semibold text-green-800")
                                            for key, label in [
                                                ("cram", "CRAM"),
                                                ("bedgraph", "Bedgraph"),
                                                ("vaf_bedgraph", "VAF"),
                                                ("deepvariant", "DeepVariant"),
                                            ]:
                                                color = (
                                                    "green"
                                                    if avail.get(key)
                                                    else "grey"
                                                )
                                                ui.badge(label, color=color).props(
                                                    "outline"
                                                    if not avail.get(key)
                                                    else ""
                                                ).classes("text-xs")
                                            ui.space()
                                            if avail.get("bedgraph") or avail.get(
                                                "cram"
                                            ):
                                                ui.button(
                                                    "Visualize",
                                                    icon="visibility",
                                                    on_click=lambda: show_sample_dialog(
                                                        barcode
                                                    ),
                                                ).props("color=blue dense")
                                else:
                                    with ui.card().classes(
                                        "w-full border-l-4 border-amber-500"
                                    ):
                                        with ui.row().classes("items-center gap-2 p-3"):
                                            ui.icon("warning", color="amber").classes(
                                                "text-xl"
                                            )
                                            ui.label(
                                                f"Sample not found: {barcode}"
                                            ).classes("text-amber-800")

                        ui.button(
                            "Search Sample",
                            on_click=search_sample,
                            icon="person_search",
                        ).props("color=blue dense")
                        sample_input.on("keydown.enter", search_sample)

                    with ui.row().classes("items-end gap-4 w-full"):
                        family_input = (
                            ui.input(
                                label="Family ID",
                                placeholder="e.g. C0733-011-068",
                            )
                            .props("outlined dense")
                            .classes("w-64")
                        )

                        async def search_family():
                            fid = (
                                family_input.value.strip() if family_input.value else ""
                            )
                            if not fid:
                                ui.notify("Enter a family ID", type="warning")
                                return
                            search_result.clear()
                            with search_result:
                                with ui.row().classes("items-center gap-2"):
                                    ui.spinner(size="sm")
                                    ui.label("Searching...").classes("text-gray-500")

                            family_path = get_family_path(store.data_dir, fid)
                            exists = await asyncio.to_thread(family_path.is_dir)

                            search_result.clear()
                            with search_result:
                                if exists:
                                    ped_file = family_path / f"{fid}.pedigree.tsv"
                                    member_count = 0
                                    if ped_file.exists():
                                        members = await asyncio.to_thread(
                                            load_family_pedigree, ped_file
                                        )
                                        member_count = len(members)

                                    avail = await asyncio.to_thread(
                                        check_family_availability,
                                        store.data_dir,
                                        fid,
                                        [],
                                    )
                                    with ui.card().classes(
                                        "w-full border-l-4 border-green-500"
                                    ):
                                        with ui.row().classes("items-center gap-3 p-3"):
                                            ui.icon(
                                                "check_circle", color="green"
                                            ).classes("text-xl")
                                            ui.label(f"Family found: {fid}").classes(
                                                "font-semibold text-green-800"
                                            )
                                            if member_count:
                                                ui.badge(
                                                    f"{member_count} members"
                                                ).props("color=blue")
                                            for key, label in [
                                                ("pedigree", "Pedigree"),
                                                ("vcfs", "VCFs"),
                                                ("wombat", "Wombat"),
                                                ("wisecondorx", "WisecondorX"),
                                            ]:
                                                color = (
                                                    "green"
                                                    if avail.get(key)
                                                    else "grey"
                                                )
                                                ui.badge(label, color=color).props(
                                                    "outline"
                                                    if not avail.get(key)
                                                    else ""
                                                ).classes("text-xs")
                                            ui.space()
                                            ui.button(
                                                "Open Family Page",
                                                icon="open_in_new",
                                                on_click=lambda: ui.navigate.to(
                                                    f"/family/{fid}"
                                                ),
                                            ).props("color=blue dense")
                                else:
                                    with ui.card().classes(
                                        "w-full border-l-4 border-amber-500"
                                    ):
                                        with ui.row().classes("items-center gap-2 p-3"):
                                            ui.icon("warning", color="amber").classes(
                                                "text-xl"
                                            )
                                            ui.label(
                                                f"Family not found: {fid}"
                                            ).classes("text-amber-800")

                        ui.button(
                            "Search Family",
                            on_click=search_family,
                            icon="family_restroom",
                        ).props("color=blue dense")
                        family_input.on("keydown.enter", search_family)

            if not store.cohorts:
                with ui.card().classes("w-full p-6 bg-yellow-50"):
                    ui.label("⚠️ No cohorts found").classes(
                        "text-xl font-semibold text-yellow-800"
                    )
                    ui.label(
                        f"No valid cohorts were found in: {store.cohorts_dir}"
                    ).classes("text-gray-600")
                    ui.label(
                        "Make sure each cohort directory contains a .pedigree.tsv file."
                    ).classes("text-gray-500 text-sm")
                return

            # Display cohorts as cards
            ui.label("Available Cohorts").classes(
                "text-2xl font-semibold mb-4 text-blue-800"
            )

            with ui.row().classes("w-full flex-wrap gap-4"):
                for cohort in sorted(store.cohorts.values(), key=lambda c: c.name):
                    cohort_name = cohort.name  # Capture name for lambda
                    with (
                        ui.card()
                        .classes(
                            "cursor-pointer hover:shadow-lg transition-shadow w-80 border-l-4 border-blue-500"
                        )
                        .on(
                            "click",
                            lambda _, n=cohort_name: ui.navigate.to(f"/cohort/{n}"),
                        )
                    ):
                        with ui.card_section():
                            ui.label(cohort.name).classes(
                                "text-xl font-bold text-blue-700"
                            )

                        with ui.card_section():
                            with ui.row().classes("gap-6"):
                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_families)).classes(
                                        "text-3xl font-bold text-blue-600"
                                    )
                                    ui.label("Families").classes(
                                        "text-sm text-gray-500"
                                    )

                                with ui.column().classes("items-center"):
                                    ui.label(str(cohort.num_samples)).classes(
                                        "text-3xl font-bold text-green-600"
                                    )
                                    ui.label("Samples").classes("text-sm text-gray-500")

                        with ui.card_section().classes("bg-gray-50"):
                            ui.label(f"📄 {cohort.pedigree_file.name}").classes(
                                "text-xs text-gray-400"
                            )

        except RuntimeError as e:
            ui.label(f"Error: {e}").classes("text-red-500")
