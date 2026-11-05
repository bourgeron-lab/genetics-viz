"""The SVs tab's Vedette sub-tab: presets, link, the SV dialog and validations."""

import asyncio
import csv
import io

import datagen
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab


async def _table(tables, n_before: int = 0):
    """The newest DataTable built after ``n_before`` tables."""
    for _ in range(200):
        if len(tables) > n_before:
            return tables[-1]
        await asyncio.sleep(0.05)
    raise AssertionError("no table built")


def _keys(table):
    return [(r["site_id"], r["sample"]) for r in table["rows"]]


async def _open_vedette(user, monkeypatch, role="curator", query=""):
    await login(user, role)
    tables = capture_tables(monkeypatch)
    await user.open(f"/cohort/DEMO/family/1001{query}")
    if "svtab=vedette" not in query:
        open_tab(user, "family-tabs", "SVs")
        await user.should_see("WisecondorX Aberrations", retries=30)
        n = len(tables)
        open_tab(user, "svs-subtabs", "Vedette")
    else:
        n = 0
    await user.should_see("Vedette calls", retries=60)
    await user.should_see(datagen.VEDETTE_CALLER, retries=10)
    for _ in range(200):
        if any(r.get("site_id") for t in tables[n:] for r in t["rows"]):
            break
        await asyncio.sleep(0.05)
    return tables


async def test_vedette_tab_opens_on_rare_calls_and_presets_switch(user, monkeypatch):
    tables = await _open_vedette(user, monkeypatch)
    table = tables[-1]
    # Rare (the opening preset): carrier fraction <= 1 %, GQ >= 20, a missing
    # value filtering nothing (the mosaic row); highest impact first.
    assert _keys(table) == [
        ("vedette_site_1", "S3"),
        ("vedette_site_5", "S2"),
        ("vedette_site_5", "S3"),
        ("mosaic_S3_chr2_180000001", "S3"),
    ]
    assert {"kind", "own_call", "cohort_carrier_frac"} <= set(table["columns"])
    # Vedette's QC flags S3, as it flags family C0733-011-301's affected daughter.
    await user.should_see("Vedette QC flags", retries=10)
    await user.should_see("carrier_load", retries=10)
    assert {"chr:start-end", "family", "inheritance", "Validation"} <= set(
        table["columns"]
    )

    n = len(tables)
    user.find(marker="vedette-preset-denovo").click()
    assert _keys(await _table(tables, n)) == [("vedette_site_1", "S3")]

    n = len(tables)
    user.find(marker="vedette-preset-all").click()
    assert len((await _table(tables, n))["rows"]) == 10


async def test_non_pass_sites_are_shown_on_request(user, monkeypatch):
    tables = await _open_vedette(user, monkeypatch, role="reader")
    n = len(tables)
    user.find(marker="vedette-preset-all").click()
    rows = (await _table(tables, n))["rows"]
    assert len(rows) == 10 and all(r["filter"] != "CohortArtifact" for r in rows)
    n = len(tables)
    switch = next(iter(user.find(marker="vedette-nonpass").elements))
    with user.client:
        switch.set_value(True)
    table = await _table(tables, n)
    assert len(table["rows"]) == 12 and "filter" in table["columns"]
    assert {r["site_id"] for r in table["rows"]} >= {datagen.VEDETTE_NONPASS}


async def test_a_link_opens_the_sub_tab_with_its_preset(user, monkeypatch):
    tables = await _open_vedette(
        user, monkeypatch, role="reader", query="?tab=svs&svtab=vedette&vpre=impact"
    )
    assert _keys(tables[-1]) == [
        ("vedette_site_1", "S3"),
        ("vedette_site_5", "S2"),
        ("vedette_site_5", "S3"),
    ]


async def test_dialog_shows_the_own_call_and_saves_to_vedette_tsv(
    user, monkeypatch, data_dir
):
    tables = await _open_vedette(user, monkeypatch)
    table = tables[-1]
    [row] = [r for r in table["rows"] if r["site_id"] == "vedette_site_1"]
    with user.client:
        table["on_row_action"]({"action": "view_sv", "row": row})
    await user.should_see("Vedette SV", retries=30)
    await user.should_see("vedette_DEL_7", retries=10)
    await user.should_see("precise", retries=10)
    # Read from the Parquet export (v1.9.0), not matched in the sample's VCF.
    await user.should_see("from Vedette's export", retries=10)

    # Read Evidence: the deletion inside the reads, and its breakpoints.
    for _ in range(200):
        charts = user.find(marker="read-chart").elements if _has(user) else set()
        names = {s["name"] for c in charts for s in c.options.get("series", [])}
        if "Deletion in reads (CIGAR)" in names:
            break
        await asyncio.sleep(0.05)
    assert "Deletion in reads (CIGAR)" in names
    # The first deleted base, as Vedette's start (the reads' last base kept + 1).
    await user.should_see(f"Start: {datagen.VEDETTE_START:,}", retries=30)
    user.find(kind=ui.button, content="Use as curated position").click()
    await user.should_see("from the reads", retries=30)
    start_input = next(iter(user.find(kind=ui.input, content="Curated start").elements))
    end_input = next(iter(user.find(kind=ui.input, content="Curated end").elements))
    assert (start_input.value, end_input.value) == (
        str(datagen.VEDETTE_START),
        str(datagen.VEDETTE_DEL[1]),
    )

    user.find(kind=ui.button, content="Save Validation").click()
    saved = data_dir / "validations" / "vedette.tsv"
    for _ in range(100):
        if saved.exists():
            break
        await asyncio.sleep(0.05)
    header, *rows = [line.split("\t") for line in saved.read_text().splitlines()]
    record = dict(zip(header, rows[0]))
    # Keyed by the cohort site's site_key, stable across Vedette releases.
    assert record["Variant"] == datagen.VEDETTE_SITE_KEY
    assert (record["Sample"], record["Validation"], record["User"]) == (
        "S3",
        "present",
        "curator",
    )
    assert (record["CuratedStart"], record["CuratedEnd"]) == (
        str(datagen.VEDETTE_START),
        str(datagen.VEDETTE_DEL[1]),
    )
    assert record["SiteId"] == "vedette_site_1"
    assert record["CallerVersion"] == datagen.VEDETTE_CALLER
    assert record["ModelVersion"] == "v5n"
    # Nothing went to the WisecondorX file.
    assert not (data_dir / "validations" / "svs.tsv").exists()


def _has(user) -> bool:
    try:
        user.find(marker="read-chart")
        return True
    except AssertionError:
        return False


async def test_readers_cannot_save(user, monkeypatch):
    tables = await _open_vedette(user, monkeypatch, role="reader")
    row = tables[-1]["rows"][0]
    with user.client:
        tables[-1]["on_row_action"]({"action": "view_sv", "row": row})
    await user.should_see("Read-only access", retries=30)


async def test_download_is_the_rows_shown(user, monkeypatch):
    await _open_vedette(user, monkeypatch, role="reader")
    user.find(marker="download-vedette").click()
    response = await user.download.next(timeout=5)
    reader = csv.DictReader(io.StringIO(response.content.decode()), delimiter="\t")
    rows = list(reader)
    # The Rare view, highest impact first; badge fields stay out.
    assert [r["sample"] for r in rows] == ["S3", "S2", "S3", "S3"]
    assert "family" in reader.fieldnames and "family_badges" not in reader.fieldnames


async def test_a_vedette_diagnostic_records_the_carrier_call(
    user, monkeypatch, data_dir
):
    tables = await _open_vedette(user, monkeypatch)
    table = tables[-1]
    site4 = next(r for r in table["rows"] if r["id"] == "vedette_site_5|S3")
    with user.client:
        table["on_row_action"](
            {"action": "open_diagnostic", "row": {"id": site4["id"]}}
        )
    await user.should_see("Recorded with the diagnostic", retries=30)
    user.find(marker="diagnostic-save").click()
    await user.should_see("Diagnostic saved", retries=30)

    with open(data_dir / "diagnostics" / "svs.tsv") as f:
        [saved] = csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE)
    assert saved["Variant"] == site4["_variant_key"]
    assert saved["Source"].startswith("Vedette")
    assert saved["CopyNumber"] == str(site4["CN"])
    # S2 carries the same site; the pedigree gives her phenotype.
    assert saved["Carriers"].startswith("S2 (unaffected, ")
    assert saved["Zygosity"] in ("het", "hom", "hemi")
