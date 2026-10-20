"""The Sample QC views, in NiceGUI's user simulation, logged in as each role."""

import asyncio

import datagen
from nicegui import ui
from ui_helpers import capture_tables, login, open_tab


def _tab(user, name):
    return next(t for t in user.find(ui.tab).elements if t.props.get("name") == name)


async def open_family_qc(user, monkeypatch, username, fid="1001"):
    await login(user, username)
    tables = capture_tables(monkeypatch)
    await user.open(f"/cohort/DEMO/family/{fid}")
    await user.should_see(f"Family: {fid}", retries=10)
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    # The cohort's tables come in off-thread and add the cross-family flag.
    await user.should_not_see("Reading the cohort's tables", retries=30)
    return tables


def _issues(tables):
    return [t for t in tables if "severity" in t["columns"]][-1]


async def test_curator_sees_flags_evidence_and_cross_family_relatives(
    user, monkeypatch
):
    tables = await open_family_qc(user, monkeypatch, "curator")
    issues = _issues(tables)
    assert "evidence" in issues["columns"] and "other_family" in issues["columns"]
    # The acknowledged flags (the cross-family second degree, the declared
    # consanguinity) come last, greyed, with their note.
    assert [r["flag"] for r in issues["rows"]] == [
        "Not the declared parent",
        "Possible contamination",
        "Related across families",
        "Related across families",
        "Consanguinity",
    ]
    third = issues["rows"][2]
    assert third["samples"] == "S1 – S5" and third["review"] == ""
    cross = next(r for r in issues["rows"] if r["samples"] == "S2 – S5")
    assert cross["other_family"] == "AU-X-1" and cross["samples"] == "S2 – S5"
    assert cross["review"] == "Explained"
    assert cross["note"] == "Known first cousins: their mothers are sisters"
    assert cross["severity_color"] == "#bdbdbd"
    assert {"actions", "review", "note"} <= set(issues["columns"])

    pairs = [t for t in tables if "expected_kinship" in t["columns"]][-1]
    assert "observed" in pairs["columns"]
    # Mismatches first: the declared father, then the relative in AU-X-1.
    assert [r["pair"] for r in pairs["rows"][:2]] == ["S1 – S3", "S2 – S5"]
    assert {r["pair"] for r in pairs["rows"]} >= {"S1 – S2"}  # the parents' pair

    sex = [t for t in tables if "karyotype" in t["columns"]][-1]
    assert {r["karyotype"] for r in sex["rows"]} == {"XX", "XY"}
    charts = [c.options for c in user.find(ui.echart).elements]
    assert any(s.get("type") == "heatmap" for c in charts for s in c.get("series", []))


async def test_reader_gets_neutral_wording_and_declared_pairs_only(user, monkeypatch):
    tables = await open_family_qc(user, monkeypatch, "reader")
    issues = _issues(tables)
    assert "evidence" not in issues["columns"]
    assert "Pedigree relationship inconsistent" in {r["flag"] for r in issues["rows"]}
    assert all(
        "S1 is not" not in b["tooltip"]
        for r in issues["rows"]
        for b in r["flag_badges"]
    )
    # Readers see the review, not a note that may name a relationship.
    assert "actions" not in issues["columns"] and issues["on_row_action"] is None
    cross = next(r for r in issues["rows"] if r["samples"] == "S2 – S5")
    assert cross["review"] == "Explained" and cross["note"] == "(hidden)"
    assert all("cousins" not in b["tooltip"] for b in cross["review_badges"])
    pairs = [t for t in tables if "expected_kinship" in t["columns"]][-1]
    assert "observed" not in pairs["columns"]
    # The parents' pair (declared unrelated) and the cross-family one are hidden.
    assert {r["pair"] for r in pairs["rows"]} == {"S1 – S3", "S2 – S3"}
    await user.should_see("not declared related in the pedigree are not shown")


async def test_karyotype_is_neutral_for_readers(user, monkeypatch):
    tables = await open_family_qc(user, monkeypatch, "reader", fid="AU-X-1")
    sex = [t for t in tables if "karyotype" in t["columns"]][-1]
    assert {r["sample"]: r["karyotype"] for r in sex["rows"]} == {
        "S5": "Atypical",
        "S6": "XY",
    }
    # A note on a sex flag names no relationship: readers read it.
    s6 = next(r for r in _issues(tables)["rows"] if r["samples"] == "S6")
    assert s6["note"] == "Pedigree entry error, fixed at the next release"


async def test_family_header_badge_and_member_icons(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("QC: 1 high · 1 medium", retries=10)
    icons = [i for i in user.find(ui.icon).elements if i.props.get("name") == "warning"]
    # S1 and S3: NOT_PARENT is about both; S3's contamination shares its icon.
    assert len(icons) == 2
    # Clicking the badge opens the tab.
    user.find("QC: 1 high · 1 medium").trigger("click")
    await user.should_see("Sex and karyotype", retries=30)


async def test_family_without_qc_has_the_tab_disabled(user):
    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1002")
    await user.should_see("Family: 1002", retries=10)
    assert "disable" in _tab(user, "Sample QC").props
    await user.should_not_see("QC ok")


async def test_standalone_family_qc_has_no_cohort_context(user):
    await login(user, "curator")
    await user.open("/family/1001")
    open_tab(user, "family-tabs", "Sample QC")
    await user.should_see("Sex and karyotype", retries=30)
    await user.should_see("open the family from its cohort", retries=10)


async def test_denovo_config_carries_the_qc_warning(user):
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    await user.should_not_see("de novo calls of this family may be unreliable")
    open_tab(user, "wombat-subtabs", "rare_denovo")
    await user.should_see("de novo calls of this family may be unreliable", retries=30)
    await user.should_see("Not the declared parent", retries=10)


async def test_cohort_qc_panel_follows_the_filters(user, monkeypatch):
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-panel-tabs", "Sample QC")
    await user.should_see("Flagged individuals", retries=30)
    await user.should_see("5 of 6 shown individuals checked", retries=30)
    flagged = [t for t in tables if "issues" in t["columns"]][-1]
    # Most severe first, then the most flags: S1's NOT_PARENT and third degree
    # in AU-X-1. S2 and S6 have only acknowledged flags, which the panel
    # leaves out by default.
    assert [r["sample"] for r in flagged["rows"]][0] == "S1"
    assert {r["sample"] for r in flagged["rows"]} == {"S1", "S3", "S5"}
    await user.should_see("3 acknowledged flag(s) not counted")
    await user.should_see("2 related pair(s) across families", retries=10)

    individuals = tables[0]
    with user.client:
        individuals["on_filter"]({"filters": {"FID": "AU"}})
    await user.should_see("2 of 2 shown individuals checked", retries=30)
    flagged = [t for t in tables if "issues" in t["columns"]][-1]
    assert {r["sample"] for r in flagged["rows"]} == {"S5"}


async def test_cohort_table_qc_column_filters(user, monkeypatch):
    await login(user, "reader")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    individuals = tables[0]
    assert "QC" in individuals["columns"]
    by_sample = {r["Sample ID"]: r for r in individuals["rows"]}
    assert by_sample["S4"]["_qc_cats"] == ["Not checked"]
    # The acknowledged cross-family flag is one grey badge, after the others;
    # the third degree with S1 is not acknowledged.
    assert [b["label"] for b in by_sample["S5"]["QC_badges"]] == [
        "Across families",
        "Sex chromosomes",
        "Acknowledged",
    ]
    assert "cousins" not in by_sample["S5"]["QC_badges"][2]["tooltip"]
    assert by_sample["S6"]["_qc_cats"] == ["Acknowledged", "No issue"]
    assert by_sample["S6"]["_qc_rank"] == 0
    with user.client:
        individuals["on_filter"]({"filters": {"QC": ["Sex chromosomes"]}})
    await user.should_see("1 individuals (of 6 total)", retries=10)
    with user.client:
        individuals["on_filter"]({"filters": {"QC": ["Any issue"]}})
    await user.should_see("3 individuals (of 6 total)", retries=10)
    with user.client:
        individuals["on_filter"]({"filters": {"QC": ["Acknowledged"]}})
    await user.should_see("4 individuals (of 6 total)", retries=10)


async def test_home_card_shows_the_qc_figure(user):
    await login(user, "reader")
    await user.should_see("QC issues", retries=30)


# ---------------------------------------------------------------------------
# Acknowledgements (validations/sampleqc.tsv)
# ---------------------------------------------------------------------------


def write_acks(data_dir, rows):
    """Replace the synthetic acknowledgements with ``rows``."""
    datagen.write_tsv(
        data_dir / "validations" / "sampleqc.tsv", datagen.QC_ACK_COLUMNS, rows
    )


async def _wait_for_table(tables, built):
    """A refresh is fire-and-forget: wait for the rebuilt table."""
    for _ in range(100):
        if len(tables) > built:
            return
        await asyncio.sleep(0.05)


def _icons(user, name):
    return [i for i in user.find(ui.icon).elements if i.props.get("name") == name]


async def test_curator_acknowledges_and_withdraws_a_flag(user, monkeypatch, data_dir):
    tables = await open_family_qc(user, monkeypatch, "curator")
    await user.should_see("QC: 1 high · 1 medium")
    issues = _issues(tables)
    row = next(r for r in issues["rows"] if r["flag"] == "Not the declared parent")
    with user.client:
        issues["on_row_action"]({"action": "review_flag", "row": row})
    await user.should_see("Review a sample QC flag", retries=10)
    await user.should_see("not this note: it may name a relationship")
    user.find(marker="qc-ack-note").type("Adopted:\tknown to the clinic")
    built = len(tables)
    user.find(marker="qc-ack-save").click()
    await user.should_see("Acknowledgement saved", retries=30)

    rows = (data_dir / "validations" / "sampleqc.tsv").read_text().splitlines()
    assert rows[-1].split("\t")[:7] == [
        "NOT_PARENT", "1001", "S1", "S3", "explained",
        "Adopted: known to the clinic", "curator",
    ]  # fmt: skip
    await _wait_for_table(tables, built)
    issues = _issues(tables)
    assert [r["flag"] for r in issues["rows"]][0] == "Possible contamination"
    # The header badge and the member icons follow: S1 has nothing left, and
    # S2's only flag, the declared consanguinity, was acknowledged already.
    await user.should_see("QC: 1 medium", retries=10)
    assert len(_icons(user, "warning")) == 1 and len(_icons(user, "task_alt")) == 2

    row = next(r for r in issues["rows"] if r["flag"] == "Not the declared parent")
    with user.client:
        issues["on_row_action"]({"action": "review_flag", "row": row})
    await user.should_see("Review a sample QC flag", retries=10)
    user.find(marker="qc-ack-withdraw").click()
    await user.should_see("Acknowledgement withdrawn", retries=30)
    await user.should_see("QC: 1 high · 1 medium", retries=10)
    assert len(_icons(user, "warning")) == 2


async def test_show_acknowledged_switch_hides_them(user, monkeypatch):
    tables = await open_family_qc(user, monkeypatch, "reader")
    built = len(tables)
    user.find(marker="qc-show-acked").click()
    await user.should_see("2 acknowledged flag(s) hidden", retries=10)
    await _wait_for_table(tables, built)
    assert [r["flag"] for r in _issues(tables)["rows"]] == [
        "Pedigree relationship inconsistent",
        "Possible contamination",
        "Related to another family",
    ]


async def test_reader_cannot_save_an_acknowledgement(user, monkeypatch, data_dir):
    import sys

    await login(user, "reader")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Family: 1001", retries=10)
    acks = sys.modules["genetics_viz.utils.sample_qc_acks"]
    flag = (
        sys.modules["genetics_viz.utils.sample_qc"]
        .load_family_qc(data_dir, "1001")
        .flags[0]
    )
    before = (data_dir / "validations" / "sampleqc.tsv").read_text()
    with user.client:
        try:
            acks.save_flag_ack(data_dir, flag, user="reader", status="explained")
        except PermissionError:
            pass
        else:
            raise AssertionError("a reader saved an acknowledgement")
    assert (data_dir / "validations" / "sampleqc.tsv").read_text() == before


async def test_cohort_panel_switch_counts_acknowledged_flags(user, monkeypatch):
    await login(user, "curator")
    tables = capture_tables(monkeypatch)
    await user.open("/cohort/DEMO")
    await user.should_see("6 samples", retries=10)
    open_tab(user, "cohort-panel-tabs", "Sample QC")
    await user.should_see("Flagged individuals", retries=30)
    built = len(tables)
    user.find(marker="cohort-qc-include-acked").click()
    await user.should_see("3 acknowledged flag(s) counted", retries=30)
    await _wait_for_table(tables, built)
    flagged = [t for t in tables if "issues" in t["columns"]][-1]
    assert {r["sample"] for r in flagged["rows"]} == {"S1", "S2", "S3", "S5", "S6"}
    s6 = next(r for r in flagged["rows"] if r["sample"] == "S6")
    assert [b["color"] for b in s6["issues_badges"]] == ["#bdbdbd"]


async def test_denovo_warning_drops_explained_flags_only(user, data_dir):
    write_acks(
        data_dir,
        [
            ["NOT_PARENT", "1001", "S1", "S3", "explained", "Adopted", "curator",
             "2026-10-01T10:00:00"],
            ["POSSIBLE_CONTAMINATION", "1001", "S3", "", "confirmed", "", "curator",
             "2026-10-01T10:00:00"],
        ],
    )  # fmt: skip
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("Wombat Configuration: default", retries=30)
    open_tab(user, "wombat-subtabs", "rare_denovo")
    await user.should_see("de novo calls of this family may be unreliable", retries=30)
    await user.should_see("Confirmed issue by curator", retries=10)
    warning = user.find(marker="denovo-qc-warning").elements.pop()
    lines = [e.text for e in warning.descendants() if isinstance(e, ui.label)]
    assert not any(line.startswith("Not the declared parent") for line in lines)
    assert any(line.startswith("Possible contamination") for line in lines)


async def test_denovo_warning_gone_when_every_flag_is_explained(user, data_dir):
    write_acks(
        data_dir,
        [
            [code, "1001", a, b, "explained", "", "curator", "2026-10-01T10:00:00"]
            for code, a, b in (
                ("NOT_PARENT", "S1", "S3"),
                ("POSSIBLE_CONTAMINATION", "S3", ""),
            )
        ],
    )
    await login(user, "curator")
    await user.open("/cohort/DEMO/family/1001")
    await user.should_see("QC ok", retries=10)
    await user.should_see("Wombat Configuration: default", retries=30)
    open_tab(user, "wombat-subtabs", "rare_denovo")
    await user.should_see("Wombat Configuration: rare_denovo", retries=30)
    await user.should_not_see("de novo calls of this family may be unreliable")


async def test_home_card_leaves_out_acknowledged_flags(user, data_dir):
    write_acks(
        data_dir,
        [
            [code, fid, a, b, "explained", "", "curator", "2026-10-01T10:00:00"]
            for code, fid, a, b in (
                ("NOT_PARENT", "1001", "S1", "S3"),
                ("KARYOTYPE_XXY", "AU-X-1", "S5", ""),
                ("SEX_MISMATCH", "AU-X-1", "S6", ""),
                ("CRYPTIC_RELATED", "1001,AU-X-1", "S2", "S5"),
            )
        ],
    )
    await login(user, "reader")
    await user.should_see("QC ok", retries=30)
    await user.should_not_see("QC issues")
