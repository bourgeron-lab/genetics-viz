"""Validation status consensus and filtering (utils/validation_status.py)."""

import pytest

from genetics_viz.utils.validation_status import (
    FILTER_STATUSES,
    STATUS_COLORS,
    STATUS_ICONS,
    UNVALIDATED,
    VALIDATED_STATUSES,
    consensus,
    filter_by_status,
    row_status,
    status_selected,
)


@pytest.mark.parametrize(
    "validations, expected",
    [
        ([], ("", "")),
        ([("present", "maternal")], ("present", "maternal")),
        ([("present", "maternal"), ("present", "de novo")], ("present", "de novo")),
        ([("present", ""), ("present", "homozygous")], ("present", "homozygous")),
        ([("present", "paternal"), ("in phase MNV", "")], ("in phase MNV", "paternal")),
        ([("present", ""), ("absent", "")], ("conflicting", "")),
        ([("absent", "maternal")], ("absent", "")),
        ([("uncertain", "")], ("uncertain", "")),
        ([("different", "")], ("uncertain", "")),
    ],
)
def test_consensus(validations, expected):
    assert consensus(validations) == expected


def test_unvalidated_rows_count_as_todo():
    assert row_status({"Validation": ""}) == row_status({}) == UNVALIDATED
    assert row_status({"Validation": "absent"}) == "absent"


def test_empty_selection_filters_nothing():
    assert status_selected("absent", []) and status_selected("", [])


def test_todo_selects_unvalidated_rows_only():
    assert status_selected("", [UNVALIDATED]) and status_selected(None, [UNVALIDATED])
    assert not status_selected("present", [UNVALIDATED])


def test_in_phase_mnv_follows_present():
    """The menu has no entry for it: it used to be hidden by every selection."""
    assert status_selected("in phase MNV", ["present"])
    assert not status_selected("in phase MNV", ["absent", UNVALIDATED])
    rows = [{"Validation": "in phase MNV"}, {"Validation": "absent"}, {}]
    assert filter_by_status(rows, FILTER_STATUSES) == rows
    assert filter_by_status(rows, ["present"]) == rows[:1]


def test_every_menu_status_has_an_icon_and_a_colour():
    assert set(FILTER_STATUSES) == set(STATUS_ICONS) == set(STATUS_COLORS)
    assert UNVALIDATED not in VALIDATED_STATUSES
