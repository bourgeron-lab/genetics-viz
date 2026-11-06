"""Curation files written to the data directory: diagnostics, notes, validations, waves."""

import pytest

from genetics_viz.components import diagnostic_loader as dl
from genetics_viz.components import notes_loader as nl
from genetics_viz.components import validation_loader as vl
from genetics_viz.components import waves_loader as wl


#: Every curation writer here checks the role itself.
_WRITERS = (dl, nl, vl, wl)


@pytest.fixture(autouse=True)
def curator(monkeypatch):
    """Curation writes check the role; these tests write as a curator."""
    for module in _WRITERS:
        monkeypatch.setattr(module, "can_write", lambda: True)


def as_reader(monkeypatch):
    for module in _WRITERS:
        monkeypatch.setattr(module, "can_write", lambda: False)


def save_diag(
    path,
    sample="S1",
    user="u1",
    ts="2026-01-01T10:00:00",
    verdict="pathogenic",
    comment="ok",
    variant="chr1:1:A:T",
):
    dl.save_diagnostic_entry(
        path, "F1", variant, "GENE", "HIGH", sample, user, ts, comment, verdict
    )


def test_diagnostic_file_is_created_with_header(tmp_path):
    path = tmp_path / "diagnostics" / "snvs.tsv"
    save_diag(path)
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == dl.DIAGNOSTIC_HEADER and len(lines) == 2


def test_multiline_diagnostic_comment_keeps_one_record_per_line(tmp_path):
    path = tmp_path / "snvs.tsv"
    save_diag(path, comment="line one\nline\ttwo\r\n")
    entries = dl.load_diagnostic_map(path)[("chr1:1:A:T", "S1")]
    assert len(path.read_text().splitlines()) == 2
    assert (
        entries[0][0] == "pathogenic"
        and "\n" not in entries[0][3]
        and "\t" not in entries[0][3]
    )


def test_conflicting_and_ignored_diagnostics(tmp_path):
    path = tmp_path / "snvs.tsv"
    save_diag(path, user="u1", verdict="pathogenic")
    save_diag(path, user="u2", ts="2026-01-02T10:00:00", verdict="benign")
    row = {}
    dl.add_diagnostic_status_to_row(
        row, dl.load_diagnostic_map(path, "F1"), "chr1:1:A:T", "S1"
    )
    assert row["Diagnostic"] == dl.CONFLICTING_DIAGNOSTIC

    assert dl.update_diagnostic_ignore_status(
        path, "F1", "chr1:1:A:T", "S1", "2026-01-02T10:00:00", "1"
    )
    row = {}
    dl.add_diagnostic_status_to_row(
        row, dl.load_diagnostic_map(path, "F1"), "chr1:1:A:T", "S1"
    )
    assert row["Diagnostic"] == "pathogenic" and row["Diagnostic_badge"]
    assert not dl.update_diagnostic_ignore_status(path, "F1", "nope", "S1", "x", "1")


def test_viewers_cannot_write_diagnostics(tmp_path, monkeypatch):
    path = tmp_path / "snvs.tsv"
    save_diag(path)
    monkeypatch.setattr(dl, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        save_diag(path, user="u2")
    with pytest.raises(PermissionError):
        dl.update_diagnostic_ignore_status(
            path, "F1", "chr1:1:A:T", "S1", "2026-01-01T10:00:00", "1"
        )
    assert len(path.read_text().splitlines()) == 2


def test_family_diagnostics_fold_same_verdict_samples(tmp_path):
    snv = tmp_path / "diagnostics" / "snvs.tsv"
    save_diag(snv, sample="S2", user="u1")
    save_diag(snv, sample="S1", user="u2", ts="2026-02-01T00:00:00")
    save_diag(snv, sample="S3", verdict="benign")
    entries = dl.load_family_diagnostics(tmp_path, "F1", ["S1", "S2", "S3"])
    by_verdict = {e["Diagnostic"]: e for e in entries}
    assert by_verdict["pathogenic"]["Sample"] == "S1, S2"  # pedigree order
    assert by_verdict["pathogenic"]["User"] == "u2, u1"  # newest first
    assert by_verdict["benign"]["Sample"] == "S3"
    assert "on S1" in dl.format_diagnostic_contributors(by_verdict["pathogenic"])


def test_notes_round_trip(tmp_path):
    path = tmp_path / "notes" / "notes.tsv"
    nl.save_note(path, "F1", "S1", "first\tline\nsecond", "u", "2026-01-01T00:00:00")
    nl.save_note(path, "F1", "S1", "later", "u", "2026-02-01T00:00:00")
    nl.save_note(path, "F2", "S9", "other family", "u", "2026-03-01T00:00:00")
    notes = nl.load_family_notes(path, "F1")
    assert [n["Message"] for n in notes] == ["later", "first line second"]
    assert nl.delete_note(path, "F1", "2026-01-01T00:00:00", "u")
    assert not nl.delete_note(path, "F1", "2026-01-01T00:00:00", "u")
    assert [n["Message"] for n in nl.load_family_notes(path, "F1")] == ["later"]


def test_note_starting_with_a_quote_survives(tmp_path):
    path = tmp_path / "notes.tsv"
    nl.save_note(path, "F1", "S1", '"quoted" start', "u", "t")
    assert nl.load_family_notes(path, "F1")[0]["Message"] == '"quoted" start'


def test_deleting_a_note_leaves_the_others_unquoted(tmp_path):
    path = tmp_path / "notes.tsv"
    nl.save_note(path, "F1", "S1", 'say "hi"', "u", "t1")
    nl.save_note(path, "F1", "S1", "gone", "u", "t2")
    assert nl.delete_note(path, "F1", "t2", "u")
    assert [n["Message"] for n in nl.load_family_notes(path, "F1")] == ['say "hi"']
    assert '"""' not in path.read_text()


def test_quoted_diagnostic_comment_survives_an_ignore_toggle(tmp_path):
    path = tmp_path / "snvs.tsv"
    save_diag(path, comment='"likely" pathogenic')
    assert dl.update_diagnostic_ignore_status(
        path, "F1", "chr1:1:A:T", "S1", "2026-01-01T10:00:00", "1"
    )
    [entry] = dl.load_diagnostic_map(path)[("chr1:1:A:T", "S1")]
    assert entry[3] == '"likely" pathogenic'


PARENTS = {"father": "P1", "mother": "P2"}


@pytest.mark.parametrize(
    "inheritance, modes, expected",
    [
        ("maternal", ("either",), ["P2"]),
        ("paternal", ("either",), ["P1"]),
        ("either", ("either",), ["P2", "P1"]),
        ("homozygous", ("either",), []),
        ("homozygous", ("either", "homozygous"), ["P2", "P1"]),
        ("de novo", ("either",), []),
    ],
)
def test_parents_to_validate(inheritance, modes, expected):
    assert vl.parents_to_validate(inheritance, PARENTS, set(), modes) == expected


def test_parents_already_validated_or_missing_are_skipped():
    assert vl.parents_to_validate("either", PARENTS, {"P2"}) == ["P1"]
    assert (
        vl.parents_to_validate("maternal", {"father": "P1", "mother": None}, set())
        == []
    )


def test_validated_samples(tmp_path):
    path = tmp_path / "snvs.tsv"
    assert vl.validated_samples(path, "F1", "v") == set()
    path.write_text(
        "\t".join(vl.SNV_VALIDATION_HEADER)
        + "\n"
        + "".join(
            "\t".join([fid, var, sample, "u", "", "present", "", "0", "t"]) + "\n"
            for fid, var, sample in [
                ("F1", "v", "S1"),
                ("F1", "w", "S2"),
                ("F2", "v", "S3"),
            ]
        )
    )
    assert vl.validated_samples(path, "F1", "v") == {"S1"}


VALIDATION_HEADER = (
    "FID\tVariant\tSample\tUser\tInheritance\tValidation\tComment\tIgnore\tTimestamp\n"
)


def validations(tmp_path, *rows):
    path = tmp_path / "snvs.tsv"
    path.write_text(VALIDATION_HEADER + "".join("\t".join(r) + "\n" for r in rows))
    row = {}
    vl.add_validation_status_to_row(row, vl.load_validation_map(path, "F1"), "v", "S1")
    return row


def vrow(status, inheritance="", ignore="0", user="u"):
    return ("F1", "v", "S1", user, inheritance, status, "", ignore, "2026-01-01")


def test_validation_statuses(tmp_path):
    assert (
        validations(tmp_path, vrow("present"), vrow("absent"))["Validation"]
        == "conflicting"
    )
    assert (
        validations(tmp_path, vrow("present"), vrow("in phase MNV"))["Validation"]
        == "in phase MNV"
    )
    assert validations(tmp_path, vrow("uncertain"))["Validation"] == "uncertain"
    ignored = validations(tmp_path, vrow("present", ignore="1"))
    assert ignored["Validation"] == "" and ignored["Validation_badge"] is None


def test_validation_inheritance_priority(tmp_path):
    row = validations(
        tmp_path,
        vrow("present", "paternal"),
        vrow("present", "homozygous"),
        vrow("present", "de novo"),
    )
    assert row["ValidationInheritance"] == "de novo"
    row = validations(
        tmp_path, vrow("present", "paternal"), vrow("present", "homozygous")
    )
    assert row["ValidationInheritance"] == "homozygous"


def test_waves(tmp_path):
    path = tmp_path / "validations" / "waves.tsv"
    wl.save_wave_validation(path, "S1", "u", 1, "t1")
    wl.save_wave_validation(path, "S1", "u", 3, "t2")
    assert path.read_text().count("sample\tuser") == 1
    assert sorted(wl.load_waves_validations(path)["S1"]) == [1, 3]
    assert wl.get_wave_category([1, 3]) == "high wave"
    assert wl.get_wave_category(None) == "TODO" and wl.get_wave_category([7]) == "TODO"


def test_readers_cannot_write_notes_or_waves(tmp_path, monkeypatch):
    """#57: the writers refuse a reader, not only the hidden controls."""
    notes = tmp_path / "notes.tsv"
    waves = tmp_path / "waves.tsv"
    nl.save_note(notes, "F1", "S1", "kept", "u", "t1")
    as_reader(monkeypatch)
    with pytest.raises(PermissionError):
        nl.save_note(notes, "F1", "S1", "refused", "u", "t2")
    with pytest.raises(PermissionError):
        nl.delete_note(notes, "F1", "t1", "u")
    with pytest.raises(PermissionError):
        wl.save_wave_validation(waves, "S1", "u", 2, "t1")
    assert [n["Message"] for n in nl.load_family_notes(notes, "F1")] == ["kept"]
    assert not waves.exists()


def test_validation_ignore_is_curator_only(tmp_path, monkeypatch):
    """#54: ignoring or restoring an SNV validation checks the role."""
    path = tmp_path / "snvs.tsv"
    vl_row = ["F1", "chr1:1:A:T", "S1", "u", "", "present", "", "0", "t1"]
    from genetics_viz.utils.tsv import append_tsv_row

    append_tsv_row(path, vl.SNV_VALIDATION_HEADER, vl_row)
    assert vl.update_validation_ignore_status(path, "F1", "chr1:1:A:T", "S1", "t1", "1")
    as_reader(monkeypatch)
    with pytest.raises(PermissionError):
        vl.update_validation_ignore_status(path, "F1", "chr1:1:A:T", "S1", "t1", "0")
    assert path.read_text().splitlines()[1].split("\t")[7] == "1"


def test_validation_guide_ships_with_the_package():
    """#56: the guide is read from the package, so an installed copy has it."""
    from genetics_viz.components import variant_dialog

    assert "genetics_viz" in variant_dialog.VALIDATION_GUIDE_PATH.parts
    assert variant_dialog._load_validation_guide() != "Validation guide not found."
