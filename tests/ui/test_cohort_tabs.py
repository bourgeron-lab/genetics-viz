"""The cohort page's tab row: grouped tabs, the Workflow button, the filter bar, the badges."""

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


async def test_one_row_of_grouped_tabs_without_icons(user, monkeypatch):
    await _open(user, monkeypatch)
    bar = user.find(marker="cohort-tabs").elements.pop()
    tabs = [t for t in bar.descendants() if isinstance(t, ui.tab)]
    assert [t.props["name"] for t in tabs] == [
        "Families",
        "Individuals",
        "Sample QC",
        "De novo",
        "STR",
        "UPD",
        "Ancestry",
        "ROH",
        "Polygenic Scores",
    ]
    captions = [
        c.text
        for c in bar.descendants()
        if isinstance(c, ui.label) and "text-[11px]" in c.classes
    ]
    assert captions == ["Browse", "Checks", "Variants", "Background"]
    assert not any(t.props.get("icon") for t in tabs)
    assert bar.value == "Families"


async def test_statistics_tab_is_gone(user, monkeypatch):
    """An older ``tab=stats`` link opens the default tab, Families."""
    await login(user, "reader")
    await user.open("/cohort/DEMO?tab=stats")
    await user.should_see("6 samples", retries=10)
    assert user.find(marker="cohort-tabs").elements.pop().value == "Families"
    await user.should_not_see("Carriers per geneset")


async def test_filter_bar_chips_clear_the_table_filters(user, monkeypatch):
    individuals = await _open(user, monkeypatch)
    table = individuals["table"]
    assert count_label(user) == "3 families · 6 individuals"
    with user.client:
        individuals["on_filter"]({"filters": {"Phenotype": ["2"], "Father": True}})
    await user.should_see("1 of 3 families · 1 of 6 individuals", retries=10)
    chips = {c.text: c for c in user.find(marker="cohort-filter-chip").elements}
    assert set(chips) == {"Phenotype: 2", "Has sequenced father"}

    # A chip's cross clears its filter in the table's filter row too.
    with user.client:
        chips["Has sequenced father"].set_value(False)
    await wait_for(lambda: count_label(user) == "2 of 3 families · 2 of 6 individuals")
    assert {r["Sample ID"] for r in table.rows} == {"S3", "S5"}
    assert table.initial_filters == {"Phenotype": ["2"], "Father": False}

    user.find(marker="cohort-clear-filters").click()
    await wait_for(lambda: count_label(user) == "3 families · 6 individuals")
    assert len(table.rows) == 6 and table.initial_filters == {}


async def test_edit_filters_goes_back_to_the_table(user, monkeypatch):
    await _open(user, monkeypatch)
    tabs = user.find(marker="cohort-tabs").elements.pop()
    open_tab(user, "cohort-tabs", "ROH")
    await user.should_see(marker="cohort-edit-filters", retries=10)
    user.find(marker="cohort-edit-filters").click()
    await wait_for(lambda: tabs.value == "Families")


async def test_edit_filters_goes_back_to_the_table_last_filtered(user, monkeypatch):
    individuals = await _open(user, monkeypatch)
    tabs = user.find(marker="cohort-tabs").elements.pop()
    with user.client:
        individuals["on_filter"]({"filters": {"Phenotype": ["2"]}})
    open_tab(user, "cohort-tabs", "ROH")
    await user.should_see(marker="cohort-edit-filters", retries=10)
    user.find(marker="cohort-edit-filters").click()
    await wait_for(lambda: tabs.value == "Individuals")


async def _families(user, monkeypatch, username="reader"):
    """Log in and open DEMO: the families and the individuals tables."""
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    return (
        next(t for t in tables if "Structure" in t["columns"]),
        next(t for t in tables if "Sample ID" in t["columns"]),
    )


async def test_families_tab_lists_each_family(user, monkeypatch, data_dir):
    datagen.add_upd_case(data_dir)
    datagen.add_upd_calls(data_dir)
    families, _ = await _families(user, monkeypatch, "curator")
    table = families["table"]
    # The ancestry and FROH columns are read off-thread, then the rows redrawn.
    await wait_for(lambda: all(r.get("Ancestry") for r in table.rows))
    rows = {r["FID"]: r for r in table.rows}
    assert list(rows) == ["1001", "1002", "AU-X-1"]
    assert (rows["1001"]["Structure"], rows["1001"]["Trios"]) == ("Trio", 1)
    assert rows["1002"]["Structure"] == "Singleton"
    assert rows["AU-X-1"]["Affected"] == "1 / 2"
    assert rows["1001"]["Ancestry"] == "EUR"
    # S3's chr2 and S2's chr15 wait for a curator, then S3's HTT, S2's FMR1.
    await wait_for(
        lambda: (
            [b["label"] for b in table.rows[0]["Review_badges"]] == ["UPD 2", "STR 2"]
        )
    )
    assert {"Review", "Conclusion"} <= set(families["columns"])


async def test_readers_get_no_review_columns(user, monkeypatch):
    families, _ = await _families(user, monkeypatch, "reader")
    assert {"Review", "Conclusion"}.isdisjoint(families["columns"])


async def test_family_filter_narrows_every_tab(user, monkeypatch):
    families, individuals = await _families(user, monkeypatch)
    with user.client:
        families["on_filter"]({"filters": {"Structure": ["Trio"]}})
    await wait_for(lambda: count_label(user) == "1 of 3 families · 3 of 6 individuals")
    assert {r["Sample ID"] for r in individuals["table"].rows} == {"S1", "S2", "S3"}
    assert pill_texts(user, "families-tab-pills") == ["1"]
    assert pill_texts(user, "individuals-tab-pills") == ["3"]
    chip = user.find(marker="cohort-filter-chip").elements.pop()
    assert chip.text == "Structure: Trio"
    # The chip's cross clears it in the families table too.
    with user.client:
        chip.set_value(False)
    await wait_for(lambda: count_label(user) == "3 families · 6 individuals")
    assert families["table"].initial_filters == {"Structure": []}


async def test_family_filters_are_in_the_link(user, monkeypatch):
    await login(user, "reader")
    await user.open("/cohort/DEMO?fstruct=Duo")
    await user.should_see("1 of 3 families · 2 of 6 individuals", retries=30)


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
    tabs = user.find(marker="cohort-tabs").elements.pop()
    button = user.find(marker="cohort-workflow").elements.pop()
    assert box.visible and "outline" in button.props
    # Workflow is a button on the title line, not a tab of the row.
    assert "Workflow" not in [
        t.props["name"] for t in tabs.descendants() if isinstance(t, ui.tab)
    ]
    user.find(marker="cohort-workflow").click()
    await wait_for(lambda: tabs.value == "Workflow" and not box.visible)
    assert "unelevated" in button.props and "outline" not in button.props
    await user.should_see("Pipeline status", retries=30)
    open_tab(user, "cohort-tabs", "ROH")
    await wait_for(lambda: box.visible and "outline" in button.props)


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


def _position(user):
    return user.find(marker="family-position").elements.pop()


async def test_family_arrows_follow_the_families_tab(user, monkeypatch):
    from types import SimpleNamespace

    families, _ = await _families(user, monkeypatch)
    table = families["table"]
    with user.client:
        families["on_filter"]({"filters": {"Structure": ["Trio", "Duo"]}})
        table._default_sort_handler(
            SimpleNamespace(args={"sorting": [{"id": "FID", "desc": True}]})
        )
    await wait_for(lambda: [r["FID"] for r in table.rows] == ["AU-X-1", "1001"])

    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    assert _position(user).text == "2 of 2"
    previous = user.find(marker="family-previous").elements.pop()
    assert previous.props.get("icon") == "arrow_back"
    await user.should_see("Previous: AU-X-1")
    await user.should_not_see(marker="family-next")

    # 1002 is not in the list: every family, by ID.
    await user.open("/cohort/DEMO/family/1002")
    await user.should_see("Family: 1002", retries=10)
    assert _position(user).text == "2 of 3"
    await user.should_see("Among all the cohort's families, by ID")

    # A bare cohort page forgets the list.
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    assert _position(user).text == "1 of 3"
