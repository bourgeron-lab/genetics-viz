"""Admin page for managing users."""

from nicegui import ui

from genetics_viz.components.header import create_header
from genetics_viz.config_model import (
    UserConfig,
    generate_random_password,
    get_config,
    hash_password,
    save_config,
)
from genetics_viz.utils.auth import check_auth, get_current_user, is_admin


@ui.page("/admin/users")
def admin_users_page() -> None:
    """Manage users."""
    if redirect := check_auth():
        return redirect
    if not is_admin():
        ui.navigate.to("/")
        return
    create_header()

    current_username = get_current_user()

    with ui.column().classes("w-full max-w-4xl mx-auto p-6"):
        ui.label("Manage Users").classes("text-3xl font-bold text-blue-900 mb-6")

        # User list
        user_container = ui.column().classes("w-full mb-6")

        def _show_password_dialog(title: str, password: str) -> None:
            """Show a dialog with a generated password to copy."""
            with ui.dialog() as dlg, ui.card().classes("p-6 min-w-[400px]"):
                ui.label(title).classes("text-lg font-semibold mb-4")
                ui.label("Give this password to the user:").classes(
                    "text-sm text-gray-600 mb-2"
                )
                ui.input("Password", value=password).props(
                    "outlined dense readonly"
                ).classes("w-full mb-4 font-mono")

                with ui.row().classes("justify-end gap-2"):
                    ui.button(
                        "Copy",
                        icon="content_copy",
                        on_click=lambda: (
                            ui.run_javascript(
                                f'navigator.clipboard.writeText("{password}")'
                            ),
                            ui.notify("Copied to clipboard", type="positive"),
                        ),
                    ).props("outline color=blue")
                    ui.button("Close", on_click=dlg.close).props(
                        "color=blue unelevated"
                    )
            dlg.open()

        def refresh_users() -> None:
            user_container.clear()
            config = get_config()
            with user_container:
                if not config.user_list:
                    ui.label("No users configured").classes("text-gray-500 italic")
                    return

                for user in config.user_list:
                    role_colors = {
                        "administrator": "red",
                        "curator": "blue",
                        "reader": "grey",
                    }
                    with ui.card().classes("w-full p-4 mb-2"):
                        with ui.row().classes("items-center justify-between w-full"):
                            with ui.row().classes("items-center gap-3"):
                                ui.icon("person").classes("text-blue-600")
                                ui.label(user.username).classes("font-bold")
                                ui.badge(
                                    user.role,
                                    color=role_colors.get(user.role, "grey"),
                                ).classes("text-xs")
                                if user.username == current_username:
                                    ui.badge("you", color="orange").classes("text-xs")

                            with ui.row().classes("items-center gap-2"):
                                # Role change dropdown
                                def make_role_change(uname):
                                    def handler(e):
                                        if uname == current_username:
                                            ui.notify(
                                                "Cannot change your own role",
                                                type="warning",
                                            )
                                            refresh_users()
                                            return
                                        cfg = get_config()
                                        admin_count = sum(
                                            1
                                            for u in cfg.user_list
                                            if u.role == "administrator"
                                        )
                                        target = next(
                                            u
                                            for u in cfg.user_list
                                            if u.username == uname
                                        )
                                        if (
                                            target.role == "administrator"
                                            and e.value != "administrator"
                                            and admin_count <= 1
                                        ):
                                            ui.notify(
                                                "Cannot remove the last administrator",
                                                type="negative",
                                            )
                                            refresh_users()
                                            return
                                        target.role = e.value
                                        save_config(cfg)
                                        ui.notify(
                                            f"Role updated for {uname}",
                                            type="positive",
                                        )
                                        refresh_users()

                                    return handler

                                ui.select(
                                    ["reader", "curator", "administrator"],
                                    value=user.role,
                                    on_change=make_role_change(user.username),
                                ).props("outlined dense").classes("w-36")

                                # Reset password button
                                def make_reset(uname):
                                    def handler():
                                        cfg = get_config()
                                        target = next(
                                            u
                                            for u in cfg.user_list
                                            if u.username == uname
                                        )
                                        new_pw = generate_random_password()
                                        target.password = hash_password(new_pw)
                                        save_config(cfg)
                                        _show_password_dialog(
                                            f"New password for {uname}", new_pw
                                        )

                                    return handler

                                ui.button(
                                    icon="key",
                                    on_click=make_reset(user.username),
                                ).props("flat color=orange size=sm round").tooltip(
                                    "Reset password"
                                )

                                # Remove user button
                                if user.username != current_username:

                                    def make_remove(uname):
                                        def handler():
                                            cfg = get_config()
                                            admin_count = sum(
                                                1
                                                for u in cfg.user_list
                                                if u.role == "administrator"
                                            )
                                            target = next(
                                                u
                                                for u in cfg.user_list
                                                if u.username == uname
                                            )
                                            if (
                                                target.role == "administrator"
                                                and admin_count <= 1
                                            ):
                                                ui.notify(
                                                    "Cannot remove the last administrator",
                                                    type="negative",
                                                )
                                                return
                                            cfg.user_list = [
                                                u
                                                for u in cfg.user_list
                                                if u.username != uname
                                            ]
                                            save_config(cfg)
                                            ui.notify(
                                                f"User {uname} removed",
                                                type="positive",
                                            )
                                            refresh_users()

                                        return handler

                                    ui.button(
                                        icon="delete",
                                        on_click=make_remove(user.username),
                                    ).props("flat color=red size=sm round").tooltip(
                                        "Remove user"
                                    )

        refresh_users()

        # Add user form
        with ui.card().classes("w-full p-6"):
            ui.label("Add User").classes("text-xl font-semibold mb-4")
            with ui.row().classes("items-end gap-4 w-full"):
                username_input = (
                    ui.input("Username").props("outlined dense").classes("w-48")
                )
                role_select = (
                    ui.select(
                        ["reader", "curator", "administrator"],
                        value="reader",
                        label="Role",
                    )
                    .props("outlined dense")
                    .classes("w-40")
                )

                def add_user() -> None:
                    uname = username_input.value.strip()
                    role = role_select.value

                    if not uname:
                        ui.notify("Please enter a username", type="warning")
                        return

                    config = get_config()
                    if any(u.username == uname for u in config.user_list):
                        ui.notify("Username already exists", type="warning")
                        return

                    # Generate random password
                    raw_pw = generate_random_password()
                    config.user_list.append(
                        UserConfig(
                            username=uname,
                            password=hash_password(raw_pw),
                            role=role,
                        )
                    )
                    save_config(config)

                    username_input.value = ""
                    refresh_users()
                    _show_password_dialog(f"Password for {uname}", raw_pw)

                ui.button("Add User", icon="person_add", on_click=add_user).props(
                    "color=blue unelevated"
                )
