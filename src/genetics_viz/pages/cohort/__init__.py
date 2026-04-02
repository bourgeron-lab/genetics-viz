"""Cohort pages module."""

# Import all page modules to register their routes
from genetics_viz.pages.cohort import cohort, family, home, standalone_family

# Export page functions for direct access if needed
from genetics_viz.pages.cohort.cohort import cohort_page
from genetics_viz.pages.cohort.family import family_page
from genetics_viz.pages.cohort.home import home_page
from genetics_viz.pages.cohort.standalone_family import standalone_family_page

__all__ = ["home_page", "cohort_page", "family_page", "standalone_family_page"]
