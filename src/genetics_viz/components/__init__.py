"""Reusable UI components for genetics-viz."""

from genetics_viz.components.filters import create_validation_filter_menu
from genetics_viz.components.header import create_header
from genetics_viz.components.icons import VALIDATION_STATUS_ICONS, get_validation_icon
from genetics_viz.components.tables import VALIDATION_TABLE_SLOT
from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.components.variant_dialog import show_variant_dialog
from genetics_viz.components.waves_loader import (
    get_wave_category,
    get_wave_color,
    get_wave_score_color,
    load_waves_validations,
    load_waves_validations_full,
    save_wave_validation,
)

__all__ = [
    "create_header",
    "VALIDATION_STATUS_ICONS",
    "get_validation_icon",
    "create_validation_filter_menu",
    "VALIDATION_TABLE_SLOT",
    "load_validation_map",
    "add_validation_status_to_row",
    "show_variant_dialog",
    "load_waves_validations",
    "load_waves_validations_full",
    "save_wave_validation",
    "get_wave_category",
    "get_wave_color",
    "get_wave_score_color",
]
