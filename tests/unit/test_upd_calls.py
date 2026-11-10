"""The upd step's tables (utils/upd_calls.py), the candidates merged from every
UPD source (utils/upd.py), their plots, and what a validation records."""

import datagen
import pytest

from genetics_viz.utils import roh, upd, upd_calls, upd_plots
from genetics_viz.utils import sample_qc as qc
from genetics_viz.utils.validation_files import validation_path


@pytest.fixture
def calls(data_dir):
    datagen.add_upd_calls(data_dir)
    return upd_calls.load_family_upd(data_dir, "1001")


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(upd, "can_write", lambda: True)


def _flag(chrom="2"):
    evidence = datagen.UPD_FLAG[4].replace("chr2", f"chr{chrom}")
    return qc.FlagRow("family", "1001", "S2", "S3", "UPD_CANDIDATE", "high", evidence)


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------


def test_the_newest_written_label_is_shown(data_dir, calls):
    assert calls.labels == ["upd_v4", "upd_v3"]
    assert calls.label == "upd_v4" and calls.has_windows
    older = upd_calls.load_family_upd(data_dir, "1001", "upd_v3")
    assert older.label == "upd_v3" and not older.has_windows
    assert upd_calls.load_family_upd(data_dir, "1001", "upd_v9").label == "upd_v4"
    assert upd_calls.probe_upd_data(data_dir, "1001")
    assert not upd_calls.probe_upd_data(data_dir, "1002")


def test_events_keep_ids_as_strings_and_older_columns_as_nulls(tmp_path):
    path = tmp_path / "88171200009.upd_v1.events.tsv"
    old_columns = [
        c
        for c in datagen.UPD_EVENT_COLUMNS
        if c not in ("support_fraction", "shared_with")
    ]
    row = dict(
        zip(datagen.UPD_EVENT_COLUMNS, datagen._upd_event_rows("88171200009", False)[0])
    )
    row["chrom"] = "chr2"
    datagen.write_tsv(path, old_columns, [[row[c] for c in old_columns]])
    data = upd_calls.load_upd(tmp_path, "88171200009")
    event = data.events[0]
    assert event["family_id"] == "88171200009"
    assert event["chrom"] == "2"
    assert isinstance(event["start"], int)
    assert event["support_fraction"] is None and event["shared_with"] is None
    assert data.segments == [] and data.members == []


def test_windows_of_a_member_and_chromosome(calls):
    windows = upd_calls.windows_of(upd_calls.load_windows(calls), "S3", "2")
    assert windows[0]["start"] == 1 and windows[0]["end"] == 500_000
    assert [w["start"] for w in windows] == sorted(w["start"] for w in windows)
    inside = windows[0]
    assert (inside["n_het"], inside["n_pat"], inside["n_lacks_mother"]) == (0, 0, 90)
    # Solo members have no parent's columns.
    solo = upd_calls.windows_of(upd_calls.load_windows(calls), "S1", "2")[0]
    assert solo["n_mat"] is None and solo["n_lacks_father"] is None


def test_member_status_is_worded_for_the_role():
    status = "parentage: S1 excluded; no_roh_tables"
    assert upd_calls.member_status_label(status, show_detail=False) == "No ROH tables"
    assert "S1 excluded" in upd_calls.member_status_label(status, show_detail=True)
    kinship = "not_isolated: parents' kinship 0.0625"
    assert upd_calls.member_status_label(kinship, show_detail=False) == "Analysed"
    assert upd_calls.member_status_label("something_new", False) == "something_new"


def test_member_export_hides_the_parents_used_from_readers(calls):
    rows, columns = upd_calls.table_export(calls.members, "members", show_detail=False)
    assert "father_used" not in columns and columns[:2] == ["FID", "IID"]
    rows, columns = upd_calls.table_export(calls.members, "members", show_detail=True)
    assert "father_used" in columns
    events, columns = upd_calls.table_export(calls.events, "events", show_detail=False)
    assert "recurrent" not in columns  # the family table has none


def test_event_wording():
    event = {"composition": "iso", "extent": "segmental", "terminal": "yes"}
    assert upd_calls.event_type(event) == "isodisomy · segmental · terminal"
    assert upd_calls.span_text(79_974_693, 101_859_834) == "80.0–101.9 Mb"
    assert upd_calls.span_text(None, 5) == ""


# ---------------------------------------------------------------------------
# Candidates
# ---------------------------------------------------------------------------


def test_candidates_merge_every_source(calls):
    froh = [{"IID": "S3", "upd_candidate_chrom": "2"}]
    validations = [
        {"FID": "1001", "Variant": "chr7:UPD", "Sample": "S1", "Validation": "absent"},
    ]
    found = upd.upd_candidates(
        calls.events,
        froh,
        [_flag()],
        validations,
        ["high", "medium"],
        ["S1", "S2", "S3"],
    )
    assert [(c.sample, c.chrom) for c in found] == [
        ("S1", "7"),
        ("S2", "15"),
        ("S3", "2"),
    ]
    s3 = found[-1]
    assert s3.roh and len(s3.flags) == 1 and s3.confidence == "high"
    assert s3.span == (100_000, 39_700_000)
    assert found[0].reviewed and not found[0].events  # a validation alone keeps it
    # S3's low chr16 event is only listed when Low is shown.
    assert ("S3", "16") not in {(c.sample, c.chrom) for c in found}
    every = upd.upd_candidates(calls.events, [], [], [], None, ["S3"])
    assert [c.chrom for c in every] == ["2", "16"]


def test_a_candidate_keeps_its_low_events_once_listed(calls):
    found = upd.upd_candidates(calls.events, [], [_flag("16")], [], ["high"], ["S3"])
    chr16 = next(c for c in found if c.chrom == "16")
    assert chr16.flags and chr16.events[0]["confidence"] == "low"


def test_standing_validation_and_review_count():
    rows = [
        {"FID": "1001", "Variant": "chr2:UPD", "Sample": "S3", "Validation": "present"},
        {"FID": "1001", "Variant": "chr2:UPD", "Sample": "S3", "Validation": "absent",
         "Ignore": "1"},
        {"FID": "1001", "Variant": "chr15:UPD", "Sample": "S2", "Validation": "uncertain"},
    ]  # fmt: skip
    standing = upd.standing_validations(rows)
    assert standing[("S3", "2")]["Validation"] == "present"
    found = upd.upd_candidates([], [], [], rows, None, None)
    assert upd.count_to_review(found) == 1  # uncertain still needs a look


def test_family_helpers(data_dir, calls):
    assert upd.probe_family_upd(data_dir, "1001")
    assert not upd.probe_family_upd(data_dir, "1002")
    found = upd.family_candidates(data_dir, "1001")
    assert {(c.sample, c.chrom) for c in found} == {("S2", "15"), ("S3", "2")}
    evidence = upd.gather_evidence(data_dir, "1001", "S3", "2")
    assert evidence.upd_label == "upd_v4" and evidence.source == "UPD step"
    assert upd.default_origin([], "S1", "S2", evidence.events) == "paternal"


def test_parse_upd_link():
    assert upd.parse_upd_link("S3:2") == ("S3", "2")
    assert upd.parse_upd_link("S3:chrX") == ("S3", "X")
    assert upd.parse_upd_link("S3") is None
    assert upd.parse_upd_link(None) is None


def test_span_from_qc_evidence():
    assert upd.span_from_evidence(_flag()) == (100_000, 39_700_000)


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def test_a_validation_records_span_type_and_label(tmp_path, curator):
    upd.save_upd_validation(
        tmp_path, family_id="1001", sample="S3", chrom="2", user="curator",
        validation="present", origin="paternal", start=100_000, end=39_700_000,
        composition="iso", extent="segmental", upd_label="upd_v4",
    )  # fmt: skip
    row = upd.read_upd_validations(tmp_path, "1001")[0]
    assert (row["Start"], row["End"], row["Composition"], row["Extent"]) == (
        "100000", "39700000", "iso", "segmental",
    )  # fmt: skip
    assert row["UpdLabel"] == "upd_v4"
    assert upd.validated_span(row) == (100_000, 39_700_000)
    for bad in (
        dict(composition="tri"),
        dict(extent="arm"),
        dict(start=5),
        dict(start=9, end=3),
    ):
        with pytest.raises(ValueError):
            upd.save_upd_validation(
                tmp_path, family_id="1001", sample="S3", chrom="2", user="c",
                validation="present", **bad,
            )  # fmt: skip


def test_an_older_file_gains_the_new_columns(tmp_path, curator):
    path = validation_path(tmp_path, "UPD")
    path.parent.mkdir(parents=True)
    path.write_text(
        "\t".join(upd.HEADER[:9]) + "\n"
        "1001\tchr2:UPD\tS3\tcurator\tpaternal\tpresent\tSNP array\t0"
        "\t2026-10-01T10:00:00\n"
    )
    upd.save_upd_validation(
        tmp_path, family_id="1001", sample="S3", chrom="2", user="curator",
        validation="absent", composition="het",
    )  # fmt: skip
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == upd.HEADER
    rows = upd.read_upd_validations(tmp_path, "1001")
    assert rows[0]["Comment"] == "SNP array" and rows[0].get("Start") in (None, "")
    assert rows[1]["Composition"] == "het"


# ---------------------------------------------------------------------------
# ROH badges follow the validation
# ---------------------------------------------------------------------------


def test_roh_badge_follows_the_validation():
    row = {"upd_candidate_chrom": "2", "qc_pass": True}
    assert roh.flag_badges(row)[0]["label"] == "UPD? chr2"
    present = {"Validation": "present", "Inheritance": "paternal", "User": "curator"}
    assert roh.flag_badges(row, present)[0]["label"] == "UPD chr2 · paternal"
    assert roh.flag_badges(row, {"Validation": "absent"}) == []
    assert roh.flag_badges(row, {"Validation": "uncertain"})[0]["label"] == "UPD? chr2"
    assert roh.is_flagged(row)
    assert not roh.is_flagged(row, {"Validation": "absent"})


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def test_chromosome_view_draws_parents_up_and_down(calls):
    windows = upd_calls.windows_of(upd_calls.load_windows(calls), "S3", "2")
    option = upd_plots.chromosome_option(
        "2", windows, segments=calls.segments_of("S3", "2"),
        qc_span=(100_000, 39_700_000), events=calls.events_of("S3", "2"),
    )  # fmt: skip
    series = {s["name"]: s for s in option["series"]}
    lacks_mother = [p for p in series["Lacks the mother's allele"]["data"] if p]
    assert lacks_mother[0][1] == pytest.approx(-90 / 400)
    assert all(p[1] >= 0 for p in series["Allele only from the mother"]["data"] if p)
    assert series["Het rate"]["data"][0][1] == 0
    assert "Segments: isodisomy paternal" in option["legend"]["data"]
    assert "Sample QC opposite homozygotes" in option["legend"]["data"]
    assert "gneg" not in option["legend"]["data"]
    # chr2 has no windows past its end, so no gap point is drawn there.
    assert series["Depth ratio"]["data"][-1] is not None


def test_chromosome_view_without_windows_draws_lanes_only(calls):
    option = upd_plots.chromosome_option("2", [], segments=calls.segments_of("S3", "2"))
    names = {s["name"] for s in option["series"]}
    assert "Het rate" not in names and "Segments: isodisomy paternal" in names
    assert upd_plots.chromosome_height(False) < upd_plots.chromosome_height(True)


def test_a_solo_member_has_no_parent_bars(calls):
    windows = upd_calls.windows_of(upd_calls.load_windows(calls), "S1", "2")
    option = upd_plots.chromosome_option("2", windows)
    names = {s["name"] for s in option["series"]}
    assert "Allele only from the mother" not in names and "Het rate" in names


def test_window_colours():
    balanced = {"n_sites": 400, "n_mat": 50, "n_pat": 50, "n_lacks_father": 0,
                "n_lacks_mother": 0, "n_het": 120, "n_hom_alt": 80}  # fmt: skip
    assert upd_plots.window_color(balanced) == "#e5e7eb"
    paternal = dict(balanced, n_mat=0, n_pat=0, n_lacks_mother=90)
    assert upd_plots.window_color(paternal) == "#0d47a1"
    solo = {"n_sites": 400, "n_het": 0, "n_hom_alt": 80, "n_mat": None, "n_pat": None,
            "n_lacks_father": None, "n_lacks_mother": None}  # fmt: skip
    assert upd_plots.window_color(solo) == "#4a148c"


def test_genome_overview_outlines_events(calls):
    windows = upd_calls.load_windows(calls)
    by_member = {
        m: {c: upd_calls.windows_of(windows, m, c) for c in ("2", "15", "16")}
        for m in ("S1", "S2", "S3")
    }
    svg = upd_plots.genome_svg(
        ["S1", "S2", "S3"], by_member, {"S3": calls.events_of("S3")}, {"S3": "trio"}
    )
    assert svg.startswith("<svg") and "S3 · chr2 0.1–39.7 Mb" in svg


def test_cohort_chart_stacks_confidences(calls):
    option = upd_plots.cohort_chromosome_option(calls.events)
    by_name = {s["name"]: s["data"] for s in option["series"]}
    assert by_name["High"][1] == 1  # chr2
    assert by_name["Medium"][14] == 1  # chr15
