"""Authentication and authorization helpers."""

from starlette.responses import RedirectResponse

from nicegui import app


def check_auth() -> RedirectResponse | None:
    """Return a redirect to /login if the user is not authenticated.

    Usage at the top of every protected page handler::

        if redirect := check_auth():
            return redirect
    """
    if not app.storage.user.get("authenticated"):
        return RedirectResponse("/login")
    return None


def get_current_user() -> str:
    """Return the logged-in username from the session."""
    return app.storage.user.get("username", "unknown")


def get_current_role() -> str:
    """Return the user's role from the session."""
    return app.storage.user.get("role", "reader")


def can_write() -> bool:
    """Return True if the user has curator or administrator role."""
    return get_current_role() in ("curator", "administrator")


def is_admin() -> bool:
    """Return True if the user has administrator role."""
    return get_current_role() == "administrator"
