"""Cohort pages module."""

# Importing the page functions imports their modules, which registers the routes.
from genetics_viz.pages.cohort.cohort import cohort_page
from genetics_viz.pages.cohort.family import family_page
from genetics_viz.pages.cohort.home import home_page
from genetics_viz.pages.cohort.standalone_family import standalone_family_page

__all__ = ["home_page", "cohort_page", "family_page", "standalone_family_page"]
