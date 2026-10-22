"""Authentication and authorization helpers."""

from urllib.parse import quote

from starlette.responses import RedirectResponse

from nicegui import app, context


def check_auth() -> RedirectResponse | None:
    """Return a redirect to /login if the user is not authenticated.

    Usage at the top of every protected page handler::

        if redirect := check_auth():
            return redirect
    """
    if not is_authenticated():
        return RedirectResponse(_login_url())
    return None


def _login_url() -> str:
    """``/login``, with the page asked for as ``next``: a shared link survives."""
    try:
        url = context.client.request.url
    except (AttributeError, RuntimeError):
        return "/login"
    target = url.path + (f"?{url.query}" if url.query else "")
    if target == "/":
        return "/login"
    return f"/login?next={quote(target, safe='')}"


def is_authenticated() -> bool:
    """Return True if the session belongs to a logged-in user."""
    return bool(app.storage.user.get("authenticated"))


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
