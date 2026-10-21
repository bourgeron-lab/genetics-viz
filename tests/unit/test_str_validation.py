"""STR validations: the strs.tsv writer and its readers (utils/str_validation.py)."""

import pytest

from genetics_viz.components.validation_loader import (
    add_validation_status_to_row,
    load_validation_map,
)
from genetics_viz.utils import str_eh, str_validation
from genetics_viz.utils.validation_files import VALIDATION_TYPES, validation_path

KEY = "chr4:3074876-3074933:HTT"


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(str_validation, "can_write", lambda: True)


def _save(data_dir, **kw):
    args = dict(
        family_id="1001",
        variant_key=KEY,
        sample="S3",
        user="curator",
        validation="present",
    )
    args.update(kw)
    str_validation.save_str_validation(data_dir, **args)


def test_kinds_include_str_and_upd(tmp_path):
    assert VALIDATION_TYPES == ["SNV", "SV", "STR", "UPD"]
    assert validation_path(tmp_path, "STR") == tmp_path / "validations" / "strs.tsv"


def test_viewers_are_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(str_validation, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        _save(tmp_path)
    assert not str_validation.str_validation_file(tmp_path).exists()


@pytest.mark.parametrize(
    "bad",
    [
        {"validation": "maybe"},
        {"inheritance": "sideways"},
        {"curated_repeats": "36-40"},
        {"user": ""},
    ],
)
def test_bad_values_are_refused(tmp_path, curator, bad):
    with pytest.raises(ValueError):
        _save(tmp_path, **bad)


def test_one_header_and_clean_fields(tmp_path, curator):
    _save(tmp_path, comment="long-read\tconfirms\nexpansion", curated_repeats=" 41 ")
    _save(tmp_path, sample="S1", validation="absent")
    lines = str_validation.str_validation_file(tmp_path).read_text().splitlines()
    assert lines[0].split("\t") == str_validation.HEADER
    assert len(lines) == 3
    first = dict(zip(str_validation.HEADER, lines[1].split("\t")))
    assert first["Comment"] == "long-read confirms expansion"
    assert first["CuratedRepeats"] == "41" and first["Ignore"] == "0"


def test_read_filters_family_and_locus(tmp_path, curator):
    _save(tmp_path)
    _save(tmp_path, variant_key="chr9:69037286-69037304:FXN", sample="S2")
    _save(tmp_path, family_id="1002", sample="S4")
    assert len(str_validation.read_str_validations(tmp_path, "1001")) == 2
    rows = str_validation.read_str_validations(tmp_path, "1001", KEY)
    assert [r["Sample"] for r in rows] == ["S3"]
    assert str_validation.read_str_validations(tmp_path / "none", "1001") == []


def test_shared_loader_aggregates_str_rows(tmp_path, curator):
    _save(tmp_path, inheritance="de novo")
    _save(tmp_path, sample="S2", validation="absent")
    _save(tmp_path, sample="S2", validation="present")
    vmap = load_validation_map(str_validation.str_validation_file(tmp_path), "1001")
    row: dict = {}
    add_validation_status_to_row(row, vmap, KEY, "S3")
    assert (row["Validation"], row["ValidationInheritance"]) == ("present", "de novo")
    add_validation_status_to_row(row, vmap, KEY, "S2")
    assert row["Validation"] == "conflicting"


def test_validation_badge_in_the_eh_table(data_dir):
    d = str_eh.get_family_eh_dir(data_dir, "1001")
    fam = str_eh.read_eh_table(str_eh.eh_path(d, "1001", "jmsa50", "genotypes"))
    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(
        fam, members, ["S1", "S2", "S3"], validations={("HTT", "S3"): "present"}
    )
    htt = next(r for r in rows if r["locus"] == "HTT")
    assert htt["m2_badges"][-1]["label"] == "✓"
    assert htt["m2_badges"][-1]["tooltip"] == "Curator validation: present"
    assert all(b["label"] != "✓" for b in htt["m0_badges"])
    assert str_validation.str_variant_key("chr4", 3074876, 3074933, "HTT") == KEY
