"""Validation pages for genetics-viz."""

from genetics_viz.pages.validation.all import validation_all_page
from genetics_viz.pages.validation.file import validation_file_page
from genetics_viz.pages.validation.statistics import validation_statistics_page
from genetics_viz.pages.validation.wave import wave_validation_page
from genetics_viz.pages.validation.waves import waves_validation_page

__all__ = [
    "validation_statistics_page",
    "validation_file_page",
    "validation_all_page",
    "waves_validation_page",
    "wave_validation_page",
]
