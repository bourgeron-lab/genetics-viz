"""Header component for genetics-viz."""

from __future__ import annotations

from nicegui import app, context, ui

from genetics_viz import __version__
from genetics_viz.models import ChangeReport
from genetics_viz.utils.auth import get_current_user, is_admin
from genetics_viz.utils.change_monitor import check_now, subscribe
from genetics_viz.utils.data import get_data_dir_options, get_data_store


def create_header(cohort_name: str | None = None) -> None:
    """Create the application header with navigation menu.

    Args:
        cohort_name: The currently active cohort/project name, or None
                     if no cohort is selected (e.g. home, validation pages).
    """
    with ui.header().classes("bg-blue-700 text-white items-center justify-between"):
        with ui.row().classes("items-center gap-4"):
            with ui.column().classes("gap-0"):
                ui.label("🧬 Genetics-Viz").classes("text-xl font-bold leading-tight")
                ui.label(f"v{__version__}").classes(
                    "text-xs text-blue-200 leading-tight"
                )

            with ui.row().classes("gap-2 items-center"):
                # Home button (always visible)
                ui.button(
                    "Home", on_click=lambda: ui.navigate.to("/"), icon="home"
                ).props("flat color=white")

                # Validation dropdown (always visible)
                with ui.button("Validation", icon="verified").props("flat color=white"):
                    with ui.menu():
                        try:
                            store = get_data_store()
                            to_validate_dir = store.data_dir / "to_validate"
                            if to_validate_dir.exists() and to_validate_dir.is_dir():
                                tsv_files = sorted(
                                    [f.stem for f in to_validate_dir.glob("*.tsv")]
                                )
                                for file_name in tsv_files:
                                    ui.menu_item(
                                        file_name,
                                        on_click=lambda fn=file_name: ui.navigate.to(
                                            f"/validation/file/{fn}"
                                        ),
                                    )
                            ui.separator()
                            ui.menu_item(
                                "See All",
                                on_click=lambda: ui.navigate.to("/validation/all"),
                            )
                            ui.menu_item(
                                "Statistics",
                                on_click=lambda: ui.navigate.to(
                                    "/validation/statistics"
                                ),
                            )
                            ui.separator()
                            ui.menu_item(
                                "Waves",
                                on_click=lambda: ui.navigate.to("/validation/waves"),
                            )
                        except RuntimeError:
                            ui.menu_item("Loading...", auto_close=False)

                # Diagnostic dropdown (always visible)
                with ui.button("Diagnostic", icon="medical_services").props(
                    "flat color=white"
                ):
                    with ui.menu():
                        ui.menu_item(
                            "See All",
                            on_click=lambda: ui.navigate.to("/diagnostic/all"),
                        )
                        ui.menu_item(
                            "Statistics",
                            on_click=lambda: ui.navigate.to("/diagnostic/statistics"),
                        )

                # Project selector dropdown (always visible)
                try:
                    store = get_data_store()
                    cohort_names = sorted(store.cohorts.keys())
                except RuntimeError:
                    cohort_names = []

                def on_project_change(e) -> None:
                    if e.value:
                        ui.navigate.to(f"/cohort/{e.value}")

                ui.select(
                    options=cohort_names,
                    value=cohort_name,
                    label="Project",
                    on_change=on_project_change,
                ).props("outlined dense dark color=white label-color=white").classes(
                    "w-48"
                )

                # Cohort button (visible only when a project is selected)
                if cohort_name:
                    ui.button(
                        "Cohort",
                        icon="folder",
                        on_click=lambda n=cohort_name: ui.navigate.to(f"/cohort/{n}"),
                    ).props("flat color=white")

                # Search button (visible only when a project is selected)
                if cohort_name:
                    ui.button(
                        "Search",
                        icon="search",
                        on_click=lambda n=cohort_name: ui.navigate.to(f"/search/{n}"),
                    ).props("flat color=white")

        # Right side — data directory selector + user menu
        with ui.row().classes("items-center gap-2"):
            # Data directory dropdown
            dir_options = get_data_dir_options()
            current_dir = app.storage.user.get("data_dir", "")

            def on_data_dir_change(e) -> None:
                app.storage.user["data_dir"] = e.value
                ui.navigate.to(ui.context.client.page.path)

            if len(dir_options) > 1:
                ui.select(
                    options={opt["value"]: opt["label"] for opt in dir_options},
                    value=current_dir,
                    on_change=on_data_dir_change,
                ).props("outlined dense dark color=white label-color=white").classes(
                    "w-48"
                ).tooltip("Select data directory")
            elif dir_options:
                ui.label(f"📁 {dir_options[0]['label']}").classes("text-sm opacity-75")

            # Refresh button
            async def _manual_refresh() -> None:
                reports = await check_now()
                if any(r.has_changes for r in reports):
                    for r in reports:
                        for line in r.summary_lines():
                            ui.notify(line, type="positive", position="top-right")
                    ui.notify("Refreshing page...", type="info", position="top-right")
                    await ui.run_javascript(
                        "setTimeout(() => window.location.reload(), 500)"
                    )
                else:
                    ui.notify("No changes detected", type="info", position="top-right")

            (
                ui.button(icon="refresh", on_click=_manual_refresh)
                .props("flat color=white round")
                .tooltip("Check for data changes")
            )

            # Help menu
            _MAILTO = (
                "mailto:96e626ca.pasteurfr.onmicrosoft.com@emea.teams.ms"
                "?subject=Bug%20report"
            )
            _TEAMS_URL = (
                "https://teams.microsoft.com/l/channel/"
                "19%3Af6165ca30232476b89a6faf84cbcede7%40thread.tacv2/"
                "genetics-viz?groupId=d78bea08-70bb-4484-a4f8-2597008bd925"
                "&tenantId=096815dc-d9eb-4bc3-a5a3-53c77e7d34e2"
            )
            _GITHUB_ISSUE_URL = (
                "https://github.com/bourgeron-lab/genetics-viz/issues/new"
            )

            with (
                ui.button(icon="help_outline")
                .props("flat color=white round")
                .tooltip("Help")
            ):
                with ui.menu():
                    with ui.menu_item(
                        on_click=lambda: ui.navigate.to(_MAILTO),
                    ):
                        with ui.row().classes("items-center gap-2"):
                            ui.icon("mail")
                            ui.label("Email")
                    with ui.menu_item(
                        on_click=lambda: ui.navigate.to(_TEAMS_URL, new_tab=True),
                    ):
                        with ui.row().classes("items-center gap-2"):
                            ui.icon("chat")
                            ui.label("Teams channel")
                    with ui.menu_item(
                        on_click=lambda: ui.navigate.to(
                            _GITHUB_ISSUE_URL, new_tab=True
                        ),
                    ):
                        with ui.row().classes("items-center gap-2"):
                            ui.icon("bug_report")
                            ui.label("Report bug")

            # User menu
            username = get_current_user()
            with ui.button(icon="person").props("flat color=white round"):
                with ui.menu():
                    ui.menu_item(
                        f"{username}",
                        auto_close=False,
                    ).props("disable").classes("font-bold")
                    ui.separator()
                    ui.menu_item(
                        "Profile",
                        on_click=lambda: ui.navigate.to("/profile"),
                    )
                    if is_admin():
                        ui.separator()
                        ui.menu_item(
                            "Manage Directories",
                            on_click=lambda: ui.navigate.to("/admin/directories"),
                        )
                        ui.menu_item(
                            "Manage Users",
                            on_click=lambda: ui.navigate.to("/admin/users"),
                        )
                    ui.separator()

                    def logout() -> None:
                        app.storage.user.clear()
                        ui.navigate.to("/login")

                    ui.menu_item("Logout", on_click=logout)

    # Per-client change notification subscription
    client = context.client
    if not getattr(client, "_change_monitor_subscribed", False):

        async def _on_change(report: ChangeReport) -> None:
            try:
                store = get_data_store()
            except RuntimeError:
                return
            if str(store.data_dir) != report.data_dir:
                return
            with client:
                for line in report.summary_lines():
                    ui.notify(line, type="info", position="top-right", timeout=10_000)

        unsub = subscribe(_on_change)
        client.on_disconnect(unsub)
        client._change_monitor_subscribed = True  # type: ignore[attr-defined]
