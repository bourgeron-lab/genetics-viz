"""Cohort pages module."""

# Import all page modules to register their routes
from genetics_viz.pages.cohort import cohort, family, home, variant

# Export page functions for direct access if needed
from genetics_viz.pages.cohort.cohort import cohort_page
from genetics_viz.pages.cohort.family import family_page
from genetics_viz.pages.cohort.home import home_page
from genetics_viz.pages.cohort.variant import variant_page

__all__ = ["home_page", "cohort_page", "family_page", "variant_page"]
