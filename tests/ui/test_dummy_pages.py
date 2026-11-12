"""Dummy pedigree members on the pages: listed and marked, never looked up.

The data are ``datagen.DUMMY_PEDIGREE``: family 2001 (MZ twins S7/S8 with
dummy parents) and 2002 (a real cousin S9 whose mother is a dummy with dummy
parents, shared with the real aunt S10; S9's dummy father has no row).
"""

import sys
from pathlib import Path
from typing import List

import datagen
import pytest
from ui_helpers import capture_tables, login, open_tab

from genetics_viz.models import is_dummy_id


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    data = datagen.build_data_dir(tmp_path)
    datagen.add_dummy_case(data)
    return data


def spy_sample_lookups(monkeypatch) -> List[str]:
    """Record every sample ID whose files a page resolves.

    Sample files are found through ``utils/sharding`` only; its public
    functions call ``get_entity_path`` / ``get_entity_url_segment`` by global
    name, so patching those catches each one. Call it after the app booted.
    """
    sharding = sys.modules["genetics_viz.utils.sharding"]
    seen: List[str] = []
    for name in ("get_entity_path", "get_entity_url_segment"):
        original = getattr(sharding, name)

        def recording(data_dir, entity_type, entity_id, _original=original):
            if entity_type == "samples":
                seen.append(entity_id)
            return _original(data_dir, entity_type, entity_id)

        monkeypatch.setattr(sharding, name, recording)
    return seen


def _dummies(seen: List[str]) -> List[str]:
    return [s for s in seen if is_dummy_id(s)]


async def test_family_page_marks_dummies_and_never_looks_them_up(user, monkeypatch):
    await login(user, "reader")
    seen = spy_sample_lookups(monkeypatch)
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO/family/2001")
    await user.should_see("Family: 2001", retries=10)
    # Counts are the sequenced twins; the twins' parents are dummies.
    await user.should_see("2 members")
    await user.should_see("2 founders")
    # Both dummies listed in the member selector, badged, without a checkbox.
    marked = user.find(marker="dummy-member").elements
    assert len(marked) == 2
    await user.should_see("DUMMY_S7_F")
    await user.should_see("DUMMY_S7_M")

    # The Wombat tab (the default) shows S7's row; its variant dialog adds
    # the parents' tracks, and the dummy parents have none to look for.
    await user.should_see("Wombat Configuration: default", retries=30)
    wombat = next(t for t in tables if "VEP_SYMBOL" in t["rows"][0])
    assert {r["sample"] for r in wombat["rows"]} == {"S7"}
    with user.client:
        wombat["on_row_action"]({"action": "view_variant", "row": wombat["rows"][0]})
    await user.should_see("S7", retries=30)

    # ROH: the twins are full sibs through their dummy parents.
    open_tab(user, "family-tabs", "ROH")
    # Not the title: the ROH tab's tooltip says it too.
    await user.should_see("Pairwise kinship checked against the pedigree.", retries=30)
    kinship = [r for t in tables for r in t["rows"] if r.get("pair") == "S7 – S8"]
    assert kinship and kinship[0]["declared"] == "Full sibs"

    open_tab(user, "family-tabs", "Ancestry")
    await user.should_see("bundle 1.1.0 / dp10gq20", retries=30)
    open_tab(user, "family-tabs", "Polygenic Scores")
    await user.should_see("2 polygenic scores", retries=30)

    assert "S7" in seen, "the spy sees sample lookups"
    assert _dummies(seen) == []
    for table in tables:
        for row in table["rows"]:
            for key in ("sample", "IID", "Sample ID"):
                assert not is_dummy_id(row.get(key)), (table["columns"], row)


async def test_standalone_family_lists_every_dummy_and_checks_real_members_only(
    user, monkeypatch
):
    await login(user, "reader")
    seen = spy_sample_lookups(monkeypatch)
    await user.open("/family/2002")
    await user.should_see("Family: 2002", retries=10)
    await user.should_see("2 members")
    # Three generations of dummies, including S9's father who has no row.
    assert len(user.find(marker="dummy-member").elements) == 4
    for dummy in ("DUMMY_S9_M", "DUMMY_S9_MF", "DUMMY_S9_MM", "DUMMY_S9_F"):
        await user.should_see(dummy)
    assert set(seen) >= {"S9", "S10"}
    assert _dummies(seen) == []


async def test_cohort_table_counts_samples_and_badges_dummy_parents(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("10 samples", retries=10)
    individuals = next(t for t in tables if "Sample ID" in t["columns"])
    rows = {r["Sample ID"]: r for r in individuals["rows"]}
    assert set(rows) == {f"S{i}" for i in range(1, 11)}
    assert [b["label"] for b in rows["S7"]["_father_badges"]] == ["DUMMY_S7_F"]
    assert [b["label"] for b in rows["S9"]["_mother_badges"]] == ["DUMMY_S9_M"]
    assert rows["S3"]["_father_badges"] == []

    # "Has sequenced father": S3 and S6, not the twins with a dummy father.
    with user.client:
        individuals["on_filter"]({"filters": {"Father": True}})
    assert {r["Sample ID"] for r in individuals["table"].rows} == {"S3", "S6"}
