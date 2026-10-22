"""Login page for genetics-viz."""

import functools
import math
import secrets

from fastapi import Request
from nicegui import app, ui

from genetics_viz.config_model import (
    get_config,
    get_default_data_dir,
    hash_password,
    needs_rehash,
    save_config,
    verify_password,
)
from genetics_viz.utils import login_throttle


@functools.cache
def _dummy_hash() -> str:
    """A hash checked for unknown usernames, so they take as long as real ones."""
    return hash_password(secrets.token_urlsafe(16))


@ui.page("/login")
def login_page(request: Request) -> None:
    """Login page — no auth check required."""
    # If already authenticated, redirect to home
    if app.storage.user.get("authenticated"):
        ui.navigate.to("/")
        return

    address = request.client.host if request.client else "unknown"

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

                wait = login_throttle.retry_after(address, username)
                if wait > 0:
                    error_label.text = (
                        f"Too many failed attempts. Try again in {math.ceil(wait)} s."
                    )
                    error_label.classes(remove="hidden")
                    password_input.value = ""
                    return

                # Checked on the event loop (about 50 ms): the throttle check,
                # the hash and the recorded failure must not interleave with
                # another attempt, or parallel attempts would all pass the check.
                config = get_config()
                user = next(
                    (u for u in config.user_list if u.username == username), None
                )
                stored = user.password if user is not None else _dummy_hash()
                if verify_password(password, stored) and user is not None:
                    login_throttle.record_success(username)
                    if needs_rehash(user.password):
                        # Legacy SHA-512 (or older scrypt parameters): upgrade
                        # now that the plaintext is known to be right.
                        user.password = hash_password(password)
                        save_config(config)
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

                login_throttle.record_failure(address, username)
                error_label.text = "Invalid username or password"
                error_label.classes(remove="hidden")
                password_input.value = ""

            ui.button("Sign in", on_click=try_login).props(
                "color=blue-7 unelevated"
            ).classes("w-full mt-2")
            password_input.on("keydown.enter", try_login)
