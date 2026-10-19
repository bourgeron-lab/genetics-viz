"""The SV dialog's copy-number and read-evidence panels, opened from search."""

import asyncio

import datagen
from nicegui import ui
from ui_helpers import capture_tables, login


async def _open_sv(user, monkeypatch, locus: str, sample: str) -> None:
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/search/DEMO")
    await user.should_see("Search: DEMO", retries=30)
    user.find(kind=ui.button, content="+ WisecondorX").click()
    await user.should_see("Exonic only", retries=30)
    user.find("Locus (optional)").type(locus)
    user.find(marker="search-button").click()
    for _ in range(100):
        if tables:
            break
        await asyncio.sleep(0.05)
    [row] = [r for r in tables[-1]["rows"] if r["sample"] == sample]
    with user.client:
        tables[-1]["on_row_action"]({"action": "view_variant", "row": row})


async def _chart(user, marker: str):
    for _ in range(100):
        try:  # find() asserts until the off-thread load has rendered it
            return next(iter(user.find(marker=marker).elements))
        except AssertionError:
            await asyncio.sleep(0.05)
    raise AssertionError(f"no {marker} chart")


async def test_copy_number_and_read_evidence(user, monkeypatch):
    await _open_sv(user, monkeypatch, "chr2:165000000-165500000", "S3")
    await user.should_see("Copy Number View", retries=30)

    cn = await _chart(user, "cn-chart")
    # The carrier, then both parents, added because they have coverage.
    assert len(cn.options["grid"]) == 3
    assert cn.options["title"][0]["text"].startswith("S3 (carrier) · XY")
    carrier = cn.options["series"][0]["data"]
    lo, hi = datagen.CNV_LOSS
    inside = [v for p, v in carrier if lo + 500 < p < hi - 500]
    assert inside and all(abs(v - 1.0) < 0.01 for v in inside)
    father = cn.options["series"][2]["data"]
    assert all(abs(v - 2.0) < 0.01 for _, v in father)

    reads = await _chart(user, "read-chart")
    await user.should_see("No CRAM for: S1, S2", retries=30)
    assert len(reads.options["grid"]) == 4  # split: reads + coverage per breakpoint
    names = {s["name"] for s in reads.options["series"]}
    assert "Deletion (+- large insert)" in names

    # A 200 kb loss forced into one window is too wide for reads.
    toggle = next(iter(user.find(marker="read-mode").elements))
    with user.client:
        toggle.set_value("single")
    await user.should_see("switch to Split", retries=30)


async def test_without_coverage_the_igv_view_stays(user, monkeypatch):
    await _open_sv(user, monkeypatch, "chr2:199000000-201000000", "S6")
    await user.should_see("CNV Coverage View", retries=30)
    await user.should_not_see("Copy Number View")
