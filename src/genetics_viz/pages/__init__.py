"""Pages for genetics-viz."""

# Import page modules to register routes when this module is imported
from genetics_viz.pages import (
    admin,
    cohort,
    diagnostic,
    login,
    profile,
    search,
    validation,
)

__all__ = ["admin", "cohort", "diagnostic", "login", "profile", "search", "validation"]
