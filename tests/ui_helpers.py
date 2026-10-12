"""Helpers for the simulated-browser tests."""

from nicegui import ui
from nicegui.testing import User

from conftest import USERS


async def login(user: User, username: str) -> None:
    """Log in through the real form, as a person would."""
    password, _role = USERS[username]
    await user.open("/login")
    user.find("Username").type(username)
    user.find("Password").type(password)
    user.find(kind=ui.button, content="Sign in").click()
    await user.should_see("Welcome to Genetics-Viz", retries=30)


def open_tab(user: User, marker: str, label: str) -> None:
    """Switch a tab bar as a click in the browser would.

    Clicking a ``ui.tab`` does nothing in the simulation; the pages react to the
    tab bar's ``update:model-value`` event, so that is what is sent. NiceGUI
    registers listeners under the camelCase name.
    """
    user.find(marker=marker).trigger("update:modelValue", label)


def capture_tables(monkeypatch) -> list:
    """Record the columns, rows and filter callback of every DataTable built from now on.

    The tables render in JavaScript, so the simulation cannot see their cells;
    what the Python side hands them is what the browser would show. Call it
    after the app has booted: the app re-imports its modules for each test.
    """
    import sys

    table_cls = sys.modules["genetics_viz.components.tanstack_table"].DataTable
    original = table_cls.__init__
    captured: list = []

    def recording_init(self, columns, rows, *args, **kwargs):
        captured.append(
            {
                "columns": [c["id"] for c in columns],
                "rows": rows,
                "on_filter": kwargs.get("on_filter"),
            }
        )
        original(self, columns, rows, *args, **kwargs)

    monkeypatch.setattr(table_cls, "__init__", recording_init)
    return captured


def roh_threshold_labels(user: User) -> list:
    """The labels of the fixed-threshold lines on the ROH strip charts."""
    labels = []
    for chart in user.find(ui.echart).elements:
        for series in chart.options.get("series", []):
            for line in series.get("markLine", {}).get("data", []):
                labels.append(line.get("name", ""))
    return labels
