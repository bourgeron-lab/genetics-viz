"""GeneTrek in the pages: profile tiers, table badges, MOI, geneset filter."""

import asyncio
import sys

import yaml
from nicegui import ui
from ui_helpers import capture_tables, login

TIER1 = "#c62828"


async def _snapshot_loaded():
    """The app loads its GeneTrek cache off-thread at startup; wait for it."""
    genetrek = sys.modules["genetics_viz.utils.genetrek"]
    for _ in range(200):
        if genetrek.get_snapshot().ok:
            return genetrek.get_snapshot()
        await asyncio.sleep(0.05)
    raise AssertionError(genetrek.get_snapshot().status_text())


async def _tables(tables, n=1):
    for _ in range(200):
        if len(tables) >= n:
            return tables
        await asyncio.sleep(0.05)
    raise AssertionError("no table was built")


async def test_profile_card_saves_a_tier_change(user, config_file):
    await login(user, "reader")
    await _snapshot_loaded()
    await user.open("/profile")
    await user.should_see("Gene badges", retries=30)
    await user.should_see("Release test-release", retries=30)
    await user.should_see("High-confidence NDD v7")
    await user.should_not_see("Use as default for everyone")

    [select] = user.find(marker="geneset-hc_ndd").elements
    assert select.value == 1  # Tier 1 in the packaged default
    with user.client:
        select.set_value("show")
    user.find(kind=ui.button, content="Save").click()
    await user.should_see("Gene badge settings saved", retries=30)

    prefs = config_file.with_name("genetics_viz.preferences.yaml")
    saved = yaml.safe_load(prefs.read_text())["users"]["reader"]["sets"]
    assert saved["hc_ndd"] == "show"


async def test_admin_can_set_the_default(user):
    await login(user, "admin")
    await _snapshot_loaded()
    await user.open("/profile")
    await user.should_see("Use as default for everyone", retries=30)


async def test_wombat_table_badges_link_and_moi(user, monkeypatch):
    await login(user, "reader")
    await _snapshot_loaded()
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    [table] = await _tables(tables)
    assert "MOI" in table["columns"]
    rows = {r["VEP_SYMBOL"]: r for r in table["rows"]}

    [badge] = rows["SHANK3"]["GeneBadges"]
    assert badge["href"] == "http://genetrek.test/?search=SHANK3"
    assert badge["color"] == TIER1 and badge["textColor"] == "white"
    assert "Tier 1: High-confidence NDD v7" in badge["tooltip"]
    # The old symbol resolves to the gene it names now.
    [old] = rows["KIAA1109"]["GeneBadges"]
    assert old["color"] == TIER1 and old["href"].endswith("?search=BLTP1")

    assert rows["SHANK3"]["MOI"] == "AD ✓"
    assert rows["KIAA1109"]["MOI"] == "AR ✓"
    assert rows["ARONLY"]["MOI"] == "AR ✗"
    assert rows["PLAIN"]["MoiBadges"] == []


async def test_search_filters_on_a_genetrek_tier(user, monkeypatch):
    await login(user, "reader")
    await _snapshot_loaded()
    tables = capture_tables(monkeypatch)
    await user.open("/search/DEMO")
    await user.should_see("Search: DEMO", retries=30)
    tier1 = next(
        cb for cb in user.find(ui.checkbox).elements if cb.text.startswith("Tier 1 (")
    )
    with user.client:
        tier1.set_value(True)
    user.find(marker="search-button").click()
    [table] = await _tables(tables)
    genes = sorted(r["VEP_SYMBOL"] for r in table["rows"])
    # Packaged default Tier 1 includes hc_ndd: SHANK3, BLTP1 (and its old
    # symbol KIAA1109) and XRG.
    assert genes == ["BLTP1", "KIAA1109", "SHANK3", "XRG"]
    assert all("MoiBadges" in r for r in table["rows"])
