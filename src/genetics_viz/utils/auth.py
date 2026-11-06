"""Authentication and authorization helpers.

The session holds who logged in; the role and the data directories granted
are read from the config at each call, so an administrator's change applies at
the user's next request, and a removed user's session stops working.
"""

from typing import Optional
from urllib.parse import quote

from starlette.responses import RedirectResponse

from nicegui import app, context

from genetics_viz.config_model import UserConfig, find_user, is_config_loaded
from genetics_viz.utils.access import may_download


def check_auth() -> RedirectResponse | None:
    """Return a redirect to /login if the user is not authenticated.

    Usage at the top of every protected page handler::

        if redirect := check_auth():
            return redirect

    A user left with no data directory (every one granted was removed) is
    logged out and told so on the login page.
    """
    if not is_authenticated():
        return RedirectResponse(_login_url())
    from genetics_viz.utils.data import get_all_data_stores, user_data_dirs

    if get_all_data_stores() and not user_data_dirs():
        app.storage.user.update(authenticated=False)
        return RedirectResponse("/login?error=no-directory")
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


def current_account() -> Optional[UserConfig]:
    """The config's entry for the logged-in user, or ``None``."""
    if not app.storage.user.get("authenticated"):
        return None
    return find_user(app.storage.user.get("username", ""))


def is_authenticated() -> bool:
    """Return True if the session belongs to a user the config still has."""
    if not app.storage.user.get("authenticated"):
        return False
    # Without a loaded config (unit tests) the session alone decides.
    return not is_config_loaded() or current_account() is not None


def get_current_user() -> str:
    """Return the logged-in username from the session."""
    return app.storage.user.get("username", "unknown")


def get_current_role() -> str:
    """Return the user's role: the config's, else the one stored at login."""
    account = current_account()
    if account is not None:
        return account.role
    return app.storage.user.get("role", "reader")


def granted_data_dirs() -> Optional[list[str]]:
    """The data directory paths the user was granted; ``None`` means all."""
    account = current_account()
    if account is None or not account.data_directories:
        return None
    return list(account.data_directories)


def can_write() -> bool:
    """Return True if the user has curator or administrator role."""
    return get_current_role() in ("curator", "administrator")


def can_download() -> bool:
    """Return True if the user gets the Download buttons (not a student)."""
    return may_download(get_current_role())


def is_admin() -> bool:
    """Return True if the user has administrator role."""
    return get_current_role() == "administrator"
