"""Sample QC acknowledgements: validations/sampleqc.tsv (utils/sample_qc_acks.py)."""

import datagen
import pytest

from genetics_viz.utils import sample_qc as qc
from genetics_viz.utils import sample_qc_acks as acks


class Cohort:
    def __init__(self, data_dir):
        self.name = "DEMO"
        self.path = data_dir / "cohorts" / "DEMO"


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(acks, "can_write", lambda: True)


def _flag(code="NOT_PARENT", fid="1001", a="S1", b="S3", severity="high"):
    return qc.FlagRow("family", fid, a, b, code, severity, "evidence")


def _cohort_flags(data_dir):
    return {f.ack_key: f for f in qc.load_cohort_qc(Cohort(data_dir)).flags}


def test_ack_key_names_the_flag_and_its_samples(data_dir):
    flags = _cohort_flags(data_dir)
    assert ("CRYPTIC_RELATED", "1001,AU-X-1", "S2", "S5") in flags
    assert ("POSSIBLE_CONTAMINATION", "1001", "S3", "") in flags
    # The family file writes a family flag as the cohort file does.
    family = qc.load_family_qc(data_dir, "1001").flags
    assert {f.ack_key for f in family} <= set(flags)


def test_last_row_wins_and_withdrawn_removes(data_dir):
    standing = acks.read_acks(data_dir)
    assert set(standing) == {
        ("SEX_MISMATCH", "AU-X-1", "S6", ""),
        ("CRYPTIC_RELATED", "1001,AU-X-1", "S2", "S5"),
    }
    cross = standing[("CRYPTIC_RELATED", "1001,AU-X-1", "S2", "S5")]
    assert cross.status == "explained" and cross.label == "Explained"
    assert cross.user == "curator" and cross.date == "2026-09-03"


def test_no_file_means_no_acknowledgement(tmp_path):
    assert acks.read_acks(tmp_path) == {}
    assert acks.load_acks(tmp_path) == {}


def test_save_writes_the_header_once_and_cleans_the_note(tmp_path, curator):
    flag = _flag()
    acks.save_flag_ack(
        tmp_path,
        flag,
        user="alice",
        status="confirmed",
        note='"Adopted"\tnot\nthe father ',
        timestamp="2026-10-20T09:00:00",
    )
    acks.save_flag_ack(tmp_path, _flag("SEX_MISMATCH", b=None), user="bob",
                       status="explained")  # fmt: skip
    lines = acks.ack_file(tmp_path).read_text().splitlines()
    assert lines[0].split("\t") == acks.HEADER
    assert len(lines) == 3
    assert lines[1].split("\t") == [
        "NOT_PARENT", "1001", "S1", "S3", "confirmed", '"Adopted" not the father',
        "alice", "2026-10-20T09:00:00",
    ]  # fmt: skip
    # A note starting with a quote reads back whole.
    ack = acks.read_acks(tmp_path)[flag.ack_key]
    assert ack.note == '"Adopted" not the father' and ack.label == "Confirmed issue"
    assert ("SEX_MISMATCH", "1001", "S1", "") in acks.read_acks(tmp_path)


def test_withdraw_then_acknowledge_again(tmp_path, curator):
    flag = _flag()
    acks.save_flag_ack(tmp_path, flag, user="alice", status="explained")
    acks.withdraw_flag_ack(tmp_path, flag, user="alice")
    assert acks.read_acks(tmp_path) == {}
    acks.save_flag_ack(tmp_path, flag, user="bob", status="confirmed", note="ok")
    assert acks.read_acks(tmp_path)[flag.ack_key].user == "bob"


def test_viewers_and_bad_values_are_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(acks, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        acks.save_flag_ack(tmp_path, _flag(), user="r", status="explained")
    with pytest.raises(PermissionError):
        acks.withdraw_flag_ack(tmp_path, _flag(), user="r")
    monkeypatch.setattr(acks, "can_write", lambda: True)
    with pytest.raises(ValueError):
        acks.save_flag_ack(tmp_path, _flag(), user="c", status="withdrawn")
    with pytest.raises(ValueError):
        acks.save_flag_ack(tmp_path, _flag(), user="", status="explained")
    assert not acks.ack_file(tmp_path).exists()


def test_cache_picks_up_an_append(tmp_path, curator):
    assert acks.load_acks(tmp_path) == {}
    acks.save_flag_ack(tmp_path, _flag(), user="a", status="explained")
    assert len(acks.load_acks(tmp_path)) == 1
    acks.withdraw_flag_ack(tmp_path, _flag(), user="a")
    assert acks.load_acks(tmp_path) == {}


def test_split_keeps_the_order(data_dir):
    flags = qc.load_cohort_qc(Cohort(data_dir)).flags
    active, acked = acks.split_acknowledged(flags, acks.load_acks(data_dir))
    assert [f.code for f in acked] == ["SEX_MISMATCH", "CRYPTIC_RELATED"]
    assert active == [f for f in flags if f not in acked]
    assert acks.active_flags(flags, {}) == flags


def test_note_wording_follows_the_role(data_dir):
    standing = acks.load_acks(data_dir)
    flags = _cohort_flags(data_dir)
    cross = flags[("CRYPTIC_RELATED", "1001,AU-X-1", "S2", "S5")]
    sex = flags[("SEX_MISMATCH", "AU-X-1", "S6", "")]
    assert acks.note_visible(cross, True) and not acks.note_visible(cross, False)
    assert acks.note_visible(sex, False)
    assert not acks.note_visible(_flag("KARYOTYPE_XXY", b=None), False)
    assert (
        acks.ack_summary(cross, standing[cross.ack_key], True)
        == "Explained by curator, 2026-09-03: Known first cousins: their mothers"
        " are sisters"
    )
    assert acks.ack_summary(cross, standing[cross.ack_key], False) == (
        "Explained by curator, 2026-09-03 (note)"
    )


def test_acknowledged_badges_are_grey(data_dir):
    standing = acks.load_acks(data_dir)
    flags = list(_cohort_flags(data_dir).values())
    badges = acks.acknowledged_badges(flags, standing, show_detail=False)
    assert [b["label"] for b in badges] == [
        "Sex inconsistent with pedigree",
        "Related to another family",
    ]
    assert {b["color"] for b in badges} == {acks.muted_color()}
    assert all("Acknowledged — Explained" in b["tooltip"] for b in badges)
    assert "cousins" not in badges[1]["tooltip"]
    assert acks.acknowledged_titles(flags, standing, False) == [
        "Sex inconsistent with pedigree — explained",
        "Related to another family — explained",
    ]


def test_cohort_headline_leaves_out_acknowledged_flags(data_dir):
    cohort = Cohort(data_dir)
    assert qc.cohort_headline(cohort)["high"] == 5
    counts = qc.cohort_headline(cohort, set(acks.load_acks(data_dir)))
    # S1, S3 (NOT_PARENT) and S5 (XXY); S2 and S6 had only acknowledged flags.
    assert counts["high"] == 3


def test_datagen_rows_match_the_header():
    assert datagen.QC_ACK_COLUMNS == acks.HEADER
    assert all(len(r) == len(acks.HEADER) for r in datagen.QC_ACKS)
