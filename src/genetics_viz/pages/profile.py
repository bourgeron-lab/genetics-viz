"""Profile page — view role and change password."""

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.config_model import (
    get_config,
    hash_password,
    save_config,
    verify_password,
)
from genetics_viz.utils.auth import check_auth, get_current_role, get_current_user


@ui.page("/profile")
def profile_page() -> None:
    """Render the user profile page."""
    if redirect := check_auth():
        return redirect
    create_header()

    username = get_current_user()
    role = get_current_role()

    with ui.column().classes("w-full max-w-xl mx-auto p-6"):
        ui.label("Profile").classes("text-3xl font-bold text-blue-900 mb-6")

        # User info card
        with ui.card().classes("w-full p-6 mb-6"):
            with ui.row().classes("items-center gap-4 mb-4"):
                ui.icon("person", size="lg").classes("text-blue-600")
                ui.label(username).classes("text-2xl font-bold")
            with ui.row().classes("items-center gap-2"):
                ui.label("Role:").classes("font-semibold")
                role_colors = {
                    "administrator": "red",
                    "curator": "blue",
                    "reader": "grey",
                }
                ui.badge(role, color=role_colors.get(role, "grey")).classes("text-sm")

        # Change password card
        with ui.card().classes("w-full p-6"):
            ui.label("Change Password").classes("text-xl font-semibold mb-4")

            current_pw = (
                ui.input("Current password", password=True, password_toggle_button=True)
                .props("outlined dense")
                .classes("w-full mb-2")
            )
            new_pw = (
                ui.input("New password", password=True, password_toggle_button=True)
                .props("outlined dense")
                .classes("w-full mb-2")
            )
            confirm_pw = (
                ui.input(
                    "Confirm new password", password=True, password_toggle_button=True
                )
                .props("outlined dense")
                .classes("w-full mb-4")
            )
            error_label = ui.label("").classes("text-red-500 text-sm hidden")

            def change_password() -> None:
                current = current_pw.value
                new = new_pw.value
                confirm = confirm_pw.value

                if not current or not new or not confirm:
                    error_label.text = "All fields are required"
                    error_label.classes(remove="hidden")
                    return

                if len(new) < 8:
                    error_label.text = "New password must be at least 8 characters"
                    error_label.classes(remove="hidden")
                    return

                if new != confirm:
                    error_label.text = "New passwords do not match"
                    error_label.classes(remove="hidden")
                    return

                config = get_config()
                user_entry = next(
                    (u for u in config.user_list if u.username == username), None
                )
                if user_entry is None:
                    error_label.text = "User not found in config"
                    error_label.classes(remove="hidden")
                    return

                if not verify_password(current, user_entry.password):
                    error_label.text = "Current password is incorrect"
                    error_label.classes(remove="hidden")
                    return

                # Update password
                user_entry.password = hash_password(new)
                save_config(config)

                # Clear form
                current_pw.value = ""
                new_pw.value = ""
                confirm_pw.value = ""
                error_label.classes(add="hidden")
                ui.notify("Password changed successfully", type="positive")

            ui.button("Change Password", icon="lock", on_click=change_password).props(
                "color=blue unelevated"
            )
