"""The cohort page's tab row: one level of tabs, the filter bar, the tab badges."""

import asyncio

import datagen
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab


async def wait_for(condition, tries: int = 100) -> None:
    for _ in range(tries):
        if condition():
            return
        await asyncio.sleep(0.05)
    raise AssertionError("condition never met")


def pill_texts(user, marker: str) -> list:
    """The texts of a tab's badges, in order; none drawn is ``[]``."""
    try:
        badges = sorted(user.find(marker=marker).elements, key=lambda b: b.id)
        return [b.text for b in badges]
    except AssertionError:
        return []


def count_label(user) -> str:
    return user.find(marker="cohort-count").elements.pop().text


async def _open(user, monkeypatch, username: str = "reader"):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    return next(t for t in tables if "Sample ID" in t["columns"])


async def test_one_tab_row_with_icons(user, monkeypatch):
    await _open(user, monkeypatch)
    tabs = sorted(
        (t for t in user.find(ui.tab).elements if t.props.get("icon")),
        key=lambda t: t.id,
    )
    assert [t.props["name"] for t in tabs] == [
        "Individuals",
        "Statistics",
        "Ancestry",
        "Polygenic Scores",
        "ROH",
        "UPD",
        "De novo",
        "STR",
        "Sample QC",
        "Workflow",
    ]
    assert user.find(marker="cohort-tabs").elements.pop().value == "Individuals"


async def test_filter_bar_chips_clear_the_table_filters(user, monkeypatch):
    individuals = await _open(user, monkeypatch)
    table = individuals["table"]
    assert count_label(user) == "6 individuals"
    with user.client:
        individuals["on_filter"]({"filters": {"Phenotype": ["2"], "Father": True}})
    await user.should_see("1 of 6 individuals", retries=10)
    chips = {c.text: c for c in user.find(marker="cohort-filter-chip").elements}
    assert set(chips) == {"Phenotype: 2", "Has sequenced father"}

    # A chip's cross clears its filter in the table's filter row too.
    with user.client:
        chips["Has sequenced father"].set_value(False)
    await wait_for(lambda: count_label(user) == "2 of 6 individuals")
    assert {r["Sample ID"] for r in table.rows} == {"S3", "S5"}
    assert table.initial_filters == {"Phenotype": ["2"], "Father": False}

    user.find(marker="cohort-clear-filters").click()
    await wait_for(lambda: count_label(user) == "6 individuals")
    assert len(table.rows) == 6 and table.initial_filters == {}


async def test_edit_filters_goes_back_to_the_table(user, monkeypatch):
    await _open(user, monkeypatch)
    tabs = user.find(marker="cohort-tabs").elements.pop()
    open_tab(user, "cohort-tabs", "ROH")
    await user.should_see(marker="cohort-edit-filters", retries=10)
    user.find(marker="cohort-edit-filters").click()
    await wait_for(lambda: tabs.value == "Individuals")


async def test_workflow_tab_hides_the_filter_bar(user, monkeypatch):
    import sys
    from datetime import datetime

    docs = sys.modules["genetics_viz.utils.ghfc_ngs_docs"]
    # A fresh cached document: the tab reads it instead of fetching GitHub.
    docs._cached = docs.FetchedDoc(
        markdown="# Stub reference", error=None, fetched_at=datetime.now()
    )
    await _open(user, monkeypatch)
    box = user.find(marker="cohort-filter-box").elements.pop()
    assert box.visible
    open_tab(user, "cohort-tabs", "Workflow")
    await wait_for(lambda: not box.visible)
    open_tab(user, "cohort-tabs", "Statistics")
    await wait_for(lambda: box.visible)


async def test_tab_badges_follow_the_filters(user, monkeypatch, data_dir):
    datagen.add_upd_case(data_dir)
    datagen.add_upd_calls(data_dir)
    individuals = await _open(user, monkeypatch, "curator")
    assert pill_texts(user, "individuals-tab-pills") == []
    # High, then medium: S1, S3 and S5 have a high flag; S2's is its UPD
    # candidate, which the UPD tab counts.
    assert pill_texts(user, "qc-tab-pills") == ["3"]
    assert pill_texts(user, "upd-tab-pills") == ["1", "1"]  # S3 chr2, S2 chr15
    # STR is read off-thread: S3's HTT (high), S2's FMR1 and S6's ARX_2.
    await wait_for(lambda: pill_texts(user, "str-tab-pills") == ["1", "2"])
    # The run record says success: an icon, no count.
    await wait_for(lambda: pill_texts(user, "workflow-tab-pills") == [""])
    (run,) = user.find(marker="workflow-tab-pills").elements
    icons = [c for c in run.default_slot.children if isinstance(c, ui.icon)]
    assert [i.props["name"] for i in icons] == ["check_circle"]

    with user.client:
        individuals["on_filter"]({"filters": {"Phenotype": ["2"]}})
    await wait_for(lambda: pill_texts(user, "individuals-tab-pills") == ["2"])
    assert pill_texts(user, "qc-tab-pills") == ["2"]
    assert pill_texts(user, "upd-tab-pills") == ["1"]
    assert pill_texts(user, "str-tab-pills") == ["1"]


async def test_readers_get_no_upd_or_str_badge(user, monkeypatch, data_dir):
    datagen.add_upd_calls(data_dir)
    await _open(user, monkeypatch, "reader")
    await wait_for(lambda: pill_texts(user, "workflow-tab-pills") == [""])
    assert pill_texts(user, "qc-tab-pills")
    assert pill_texts(user, "upd-tab-pills") == []
    assert pill_texts(user, "str-tab-pills") == []
