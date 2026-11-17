"""The member card's Tree view: the pedigree drawn in place of the table.

Families 2001 and 2002 are ``datagen.DUMMY_PEDIGREE``; 1001 is the DEMO trio,
whose Sample QC has open high flags on every member.
"""

from pathlib import Path

import datagen
import pytest
from nicegui import ui
from ui_helpers import login


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    data = datagen.build_data_dir(tmp_path)
    datagen.add_dummy_case(data)
    return data


def _one(user, marker: str):
    (element,) = user.find(marker=marker).elements
    return element


def _tree(user) -> str:
    return _one(user, "member-tree").content


async def _open_tree(user, path: str, title: str) -> None:
    await user.open(path)
    await user.should_see(title, retries=10)
    _one(user, "member-view").set_value("tree")
    await user.should_see(marker="member-tree")


async def test_the_table_is_the_default_and_the_tree_draws_every_member(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/2002")
    await user.should_see("Family: 2002", retries=10)
    await user.should_not_see(marker="member-tree")
    _one(user, "member-view").set_value("tree")
    await user.should_see(marker="member-tree")
    svg = _tree(user)
    assert "<script" not in svg and 'class="pypedigree"' in svg
    for iid in ("S9", "S10", "DUMMY_S9_M", "DUMMY_S9_MF", "DUMMY_S9_MM", "DUMMY_S9_F"):
        assert f'data-id="{iid}"' in svg
    # The four dummies, S9's father (no row) included, are drawn not sampled.
    assert svg.count("Not sampled") == 4
    await user.should_see("Dashed: not sequenced.")
    # Back to the table: the tree is hidden, the grid shown.
    _one(user, "member-view").set_value("table")
    await user.should_not_see(marker="member-tree")


async def test_the_tree_shows_the_selection_and_follows_it(user):
    await login(user, "reader")
    await _open_tree(user, "/cohort/DEMO/family/2002", "Family: 2002")
    assert "node deselected" not in _tree(user)
    user.find(marker="members-none").click()
    svg = _tree(user)
    assert 'class="node deselected" data-id="S9"' in svg
    assert 'class="node deselected" data-id="S10"' in svg
    # A dummy is never selected, so never faded either.
    assert 'deselected" data-id="DUMMY' not in svg
    await user.should_see("Faded: not selected.", retries=5)
    user.find(marker="members-all").click()
    assert "node deselected" not in _tree(user)


async def test_the_view_is_a_saved_preference(user):
    await login(user, "reader")
    await _open_tree(user, "/cohort/DEMO/family/2001", "Family: 2001")
    await user.open("/cohort/DEMO/family/2002")
    await user.should_see("Family: 2002", retries=10)
    assert _one(user, "member-view").value == "tree"
    assert 'data-id="S9"' in _tree(user)


async def test_open_qc_flags_colour_the_outline_and_enlarge_redraws(user):
    await login(user, "curator")
    await _open_tree(user, "/cohort/DEMO/family/1001", "Family: 1001")
    assert 'class="node qc-high"' in _tree(user)
    await user.should_see("Coloured outline: open Sample QC issue (high).")
    user.find(marker="member-tree-enlarge").click()
    big = _one(user, "member-tree-big").content
    assert 'id="member-tree-big-main-group"' in big and 'data-id="S3"' in big
    # Two drawings on one page, ids apart.
    assert 'id="member-tree-main-group"' in _tree(user)


async def test_standalone_family_page_has_the_switch(user):
    await login(user, "reader")
    await user.open("/family/1001")
    await user.should_see("1001", retries=10)
    assert isinstance(_one(user, "member-view"), ui.toggle)
