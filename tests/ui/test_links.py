"""Shareable links: a page opened with a query string shows that view."""

import asyncio
import json
import sys
from urllib.parse import parse_qs, quote, urlsplit

from nicegui import ui
from ui_helpers import capture_tables, login, open_tab

from conftest import USERS


def record_pushes(monkeypatch) -> list:
    """The URLs the pages put in the address bar, from now on."""
    cls = sys.modules["genetics_viz.components.page_link"].PageLink
    pushed: list = []
    monkeypatch.setattr(cls, "push", lambda self: pushed.append(self.url()))
    return pushed


async def wait_for(condition, tries: int = 100) -> None:
    for _ in range(tries):
        if condition():
            return
        await asyncio.sleep(0.05)
    raise AssertionError("condition never met")


async def test_login_returns_to_the_linked_page(user):
    target = "/cohort/DEMO/family/1001?tab=roh"
    response = await user.http_client.get(target, follow_redirects=False)
    assert response.headers["location"] == f"/login?next={quote(target, safe='')}"

    await user.open(target)
    await user.should_see("Sign in to continue")
    user.find("Username").type("reader")
    user.find("Password").type(USERS["reader"][0])
    user.find(kind=ui.button, content="Sign in").click()
    await user.should_see("Family: 1001", retries=30)
    await user.should_see("Compared with the reference", retries=30)


async def test_login_ignores_a_next_off_this_site(user):
    await user.open("/login?next=//evil.example/")
    user.find("Username").type("reader")
    user.find("Password").type(USERS["reader"][0])
    user.find(kind=ui.button, content="Sign in").click()
    await user.should_see("Welcome to Genetics-Viz", retries=30)


async def test_family_link_restores_members_and_wombat_filters(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/1001?m=S3&wimp=HIGH")
    await user.should_see("Wombat Configuration: default", retries=30)
    await wait_for(lambda: tables)
    rows = tables[-1]["rows"]
    assert rows and {r["sample"] for r in rows} == {"S3"}
    assert all("stop" in r["VEP_Consequence"] for r in rows)
    # The member boxes show the linked selection (S1 unticked), not everyone
    boxes = sorted(user.find(kind=ui.checkbox).elements, key=lambda c: c.id)
    members = [c for c in boxes if c.text == ""]
    assert [c.value for c in members[:3]] == [False, False, True]


async def test_family_link_opens_its_wombat_config(user, monkeypatch):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001?wcfg=strict")
    await user.should_see("Wombat Configuration: strict", retries=30)
    await user.should_not_see("Wombat Configuration: default")


async def test_address_bar_follows_the_view(user, monkeypatch):
    await login(user, "reader")
    pushed = record_pushes(monkeypatch)
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    open_tab(user, "family-tabs", "ROH")
    await wait_for(lambda: pushed and pushed[-1].endswith("?tab=roh"))
    open_tab(user, "family-tabs", "Wombat")
    await wait_for(lambda: pushed[-1] == "/cohort/DEMO/family/1001")
    open_tab(user, "wombat-subtabs", "strict")
    await wait_for(lambda: pushed[-1].endswith("?wcfg=strict"))


async def test_cohort_link_restores_filters_and_panel(user, monkeypatch):
    await login(user, "reader")
    await user.open("/cohort/DEMO?fid=1002&panel=roh")
    await user.should_see("1 individuals (of 6 total)", retries=30)
    await user.should_see("Runs of homozygosity", retries=30)


async def test_search_link_runs_the_search(user, monkeypatch):
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    pushed = record_pushes(monkeypatch)
    q = json.dumps({"c": [{"t": "x", "l": "chr2:165000000-165500000"}]})
    await user.open(f"/search/DEMO?q={quote(q)}")
    await wait_for(lambda: tables)
    assert {r["sample"] for r in tables[-1]["rows"]} == {"S3"}
    # The form was filled from the link: the search it ran writes it back.
    written = json.loads(parse_qs(urlsplit(pushed[-1]).query)["q"][0])
    [comp] = written["c"]
    assert comp["t"] == "x" and comp["l"] == "chr2:165000000-165500000"
    assert "/" not in comp["src"]  # a display name, not a server path
