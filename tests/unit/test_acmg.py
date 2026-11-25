"""Five-tier verdicts and ACMG criteria (utils/acmg.py, utils/diagnostic_badges.py)."""

import pytest

from genetics_viz.components import diagnostic_loader as dl
from genetics_viz.utils import acmg
from genetics_viz.utils import diagnostic_badges as db
from genetics_viz.utils import family_status as fs


def test_the_five_tiers_strongest_first():
    assert db.VERDICTS == (
        "pathogenic",
        "likely pathogenic",
        "uncertain",
        "likely benign",
        "benign",
    )
    assert db.ALL_STATUSES[-1] == db.CONFLICTING
    assert [db.verdict_symbol(v) for v in db.ALL_STATUSES] == [
        "P", "LP", "VUS", "LB", "B", "!"
    ]  # fmt: skip
    assert db.strongest(["benign", "likely pathogenic", "uncertain"]) == (
        "likely pathogenic"
    )
    assert db.strongest(["conflicting", "other"]) == ""
    # A known verdict's badge is its symbol alone; an unknown value is shown raw.
    assert db.build_diagnostic_badge("likely pathogenic", [])["badge"]["text"] == ""
    assert db.build_diagnostic_badge("weird", [])["badge"]["text"] == "WEI"


def test_likely_pathogenic_solves_a_family():
    assert set(fs.SOLVING_VERDICTS) == {"pathogenic", "likely pathogenic"}


@pytest.mark.parametrize(
    "criteria, points, verdict",
    [
        # PVS1 + PS2: 12 points.
        ("PVS1;PS2", 12, "pathogenic"),
        # PVS1 + PM2_Supporting: 9, likely pathogenic (Tavtigian 2020).
        ("PVS1;PM2_Supporting", 9, "likely pathogenic"),
        ("PS2;PM2_Supporting;PP3", 6, "likely pathogenic"),
        ("PM2_Supporting;PP3", 2, "uncertain"),
        # Conflicting evidence adds up: 4 - 4 = 0.
        ("PS3;BS3", 0, "uncertain"),
        ("BP4;BP7", -2, "likely benign"),
        ("BS1;BS2", -8, "benign"),
        # BA1 is stand-alone, whatever else applies.
        ("BA1;PVS1;PS2", 12, "benign"),
    ],
)
def test_points_and_class(criteria, points, verdict):
    applied = acmg.parse_criteria(criteria)
    assert acmg.total_points(applied) == points
    assert acmg.classify(applied) == verdict


def test_names_order_and_origin():
    applied = [
        acmg.apply("PP3"),
        acmg.apply("PM2"),  # proposed at Supporting (ClinGen SVI 2020)
        acmg.apply("PVS1", "Strong"),
        acmg.apply("BP4", "Moderate", origin="computed"),
    ]
    assert acmg.format_criteria(applied) == (
        "PVS1_Strong;PM2_Supporting;PP3;BP4_Moderate:computed"
    )
    parsed = acmg.parse_criteria("PP3; PVS1_Strong ;BP4_Moderate:computed;PM2")
    assert [str(a) for a in parsed] == [
        "PVS1_Strong", "PM2", "PP3", "BP4_Moderate:computed"
    ]  # fmt: skip
    assert parsed[1].points == 2  # PM2 written bare is at its own strength
    assert acmg.evidence_columns([]) == ["", "", ""]
    assert acmg.evidence_columns(parsed) == [
        "PVS1_Strong;PM2;PP3;BP4_Moderate:computed",
        "5",
        acmg.RULESET,
    ]


@pytest.mark.parametrize(
    "bad",
    ["PX9", "PM2_Huge", "BA1_Strong", "PS2:guessed"],
)
def test_unknown_criteria_are_refused(bad):
    with pytest.raises(ValueError):
        acmg.parse_criteria(bad)


def test_record_evidence_and_override():
    record = {
        "Diagnostic": "pathogenic",
        "Criteria": "PVS1;PM2_Supporting",
        "Points": "9",
        "Ruleset": acmg.RULESET,
    }
    evidence = acmg.record_evidence(record)
    assert (evidence.points, evidence.computed_class) == (9, "likely pathogenic")
    assert acmg.override(record) == "likely pathogenic"
    assert acmg.summary(record, db.verdict_label) == (
        "PVS1, PM2_Supporting · 9 points (criteria give likely pathogenic)"
    )
    record["Diagnostic"] = "likely pathogenic"
    assert acmg.override(record) == ""
    # A legacy record, or one this version cannot read: the verdict alone.
    assert acmg.summary({"Diagnostic": "pathogenic"}) == ""
    assert not acmg.record_evidence({"Criteria": "PZZ"}).recorded


def test_writer_records_criteria_and_refuses_unknown_verdicts(tmp_path, monkeypatch):
    monkeypatch.setattr(dl, "can_write", lambda: True)
    path = tmp_path / "snvs.tsv"
    args = dict(
        diagnostic_file=path,
        family_id="1001",
        variant_key="chr1:100:A:T",
        gene="SHANK3",
        impact="stop_gained",
        sample="S3",
        user="curator",
        timestamp="2026-11-25T10:00:00",
        comment="",
    )
    with pytest.raises(ValueError):
        dl.save_diagnostic_entry(**args, diagnostic="probably bad")
    criteria = acmg.parse_criteria("PVS1;PS2;PM2_Supporting")
    dl.save_diagnostic_entry(**args, diagnostic="pathogenic", criteria=criteria)
    dl.save_diagnostic_entry(**args, diagnostic="likely benign")
    header, first, second = path.read_text().splitlines()
    assert header.split("\t")[-3:] == ["Criteria", "Points", "Ruleset"]
    assert first.split("\t")[-3:] == ["PVS1;PS2;PM2_Supporting", "13", acmg.RULESET]
    assert second.split("\t")[-3:] == ["", "", ""]
    assert second.split("\t")[9] == "likely benign"


def test_an_older_file_gains_the_columns(tmp_path, monkeypatch):
    """A file written before v0.46 keeps its records and reads as before."""
    monkeypatch.setattr(dl, "can_write", lambda: True)
    path = tmp_path / "snvs.tsv"
    old_header = dl.DIAGNOSTIC_HEADER[: -len(dl.ACMG_COLUMNS)]
    old_row = ["1001", "chr1:100:A:T", "G", "x", "S3", "bob", "2026-01-01"]
    old_row += ["", "0", "uncertain"] + [""] * (len(old_header) - 10)
    path.write_text("\t".join(old_header) + "\n" + "\t".join(old_row) + "\n")
    dl.save_diagnostic_entry(
        diagnostic_file=path,
        family_id="1001",
        variant_key="chr1:100:A:T",
        gene="G",
        impact="x",
        sample="S3",
        user="alice",
        timestamp="2026-11-25T10:00:00",
        comment="",
        diagnostic="likely pathogenic",
        criteria=acmg.parse_criteria("PS2;PM2_Supporting;PP3"),
    )
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == dl.DIAGNOSTIC_HEADER
    assert lines[1] == "\t".join(old_row)
    vmap = dl.load_diagnostic_map(path, "1001")
    assert [d[0] for d in vmap[("chr1:100:A:T", "S3")]] == [
        "uncertain",
        "likely pathogenic",
    ]
    row = {}
    dl.add_diagnostic_status_to_row(row, vmap, "chr1:100:A:T", "S3")
    # Two curators, two tiers: they disagree.
    assert row["Diagnostic"] == db.CONFLICTING


def test_report_lists_likely_pathogenic_with_its_criteria(data_dir, monkeypatch):
    from genetics_viz.utils import family_report

    monkeypatch.setattr(dl, "can_write", lambda: True)
    path = data_dir / "diagnostics" / "snvs.tsv"
    common = dict(
        diagnostic_file=path,
        family_id="1001",
        gene="SHANK3",
        impact="missense_variant",
        sample="S3",
        user="curator",
        timestamp="2026-11-25T10:00:00",
        comment="",
    )
    dl.save_diagnostic_entry(
        **common,
        variant_key="chr22:50721500:G:A",
        diagnostic="likely pathogenic",
        criteria=acmg.parse_criteria("PS2;PM2_Supporting;PP3"),
    )
    dl.save_diagnostic_entry(
        **common, variant_key="chr1:100:A:T", diagnostic="likely benign"
    )
    entries = dl.load_family_diagnostics(data_dir, "1001", ["S1", "S2", "S3"])
    section = family_report._diagnostics(data_dir, "1001", entries)
    # Likely benign is left out, as benign is.
    (row,) = section.rows
    assert row["verdict"] == "Likely pathogenic"
    assert row["_verdict_note"] == "PS2, PM2_Supporting, PP3 · 6 points"
    assert section.attention
    assert any("Tavtigian" in f for f in section.footnotes)
    assert section.glance == "1 likely pathogenic"
