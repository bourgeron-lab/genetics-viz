"""Login page for genetics-viz."""

from nicegui import app, ui

from genetics_viz.config_model import get_config, get_default_data_dir, verify_password


@ui.page("/login")
def login_page() -> None:
    """Login page — no auth check required."""
    # If already authenticated, redirect to home
    if app.storage.user.get("authenticated"):
        ui.navigate.to("/")
        return

    with ui.column().classes("absolute-center items-center gap-4"):
        ui.label("🧬 Genetics-Viz").classes("text-3xl font-bold text-blue-700")
        ui.label("Sign in to continue").classes("text-gray-500")

        with ui.card().classes("w-80 p-6"):
            username_input = (
                ui.input("Username").props("outlined dense").classes("w-full")
            )
            password_input = (
                ui.input("Password", password=True, password_toggle_button=True)
                .props("outlined dense")
                .classes("w-full")
            )
            error_label = ui.label("").classes("text-red-500 text-sm hidden")

            def try_login() -> None:
                username = username_input.value.strip()
                password = password_input.value

                if not username or not password:
                    error_label.text = "Please enter username and password"
                    error_label.classes(remove="hidden")
                    return

                config = get_config()
                for user in config.user_list:
                    if user.username == username and verify_password(
                        password, user.password
                    ):
                        app.storage.user.update(
                            {
                                "authenticated": True,
                                "username": user.username,
                                "role": user.role,
                                "data_dir": get_default_data_dir(config),
                            }
                        )
                        ui.navigate.to("/")
                        return

                error_label.text = "Invalid username or password"
                error_label.classes(remove="hidden")
                password_input.value = ""

            ui.button("Sign in", on_click=try_login).props(
                "color=blue-7 unelevated"
            ).classes("w-full mt-2")
            password_input.on("keydown.enter", try_login)
