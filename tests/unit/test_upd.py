"""UPD validations (utils/upd.py) and the acknowledgements they give QC flags."""

import datagen
import pytest

from genetics_viz.components.validation_loader import load_validation_map
from genetics_viz.utils import sample_qc as qc
from genetics_viz.utils import sample_qc_acks as acks
from genetics_viz.utils import upd
from genetics_viz.utils.validation_files import validation_path


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(upd, "can_write", lambda: True)


def _save(data_dir, **kw):
    args = dict(
        family_id="1001",
        sample="S3",
        chrom="2",
        user="curator",
        validation="present",
    )
    args.update(kw)
    upd.save_upd_validation(data_dir, **args)


def _flag(evidence=datagen.UPD_FLAG[4], a="S2", b="S3", fid="1001"):
    return qc.FlagRow("family", fid, a, b, "UPD_CANDIDATE", "high", evidence)


def test_key_and_chromosome():
    assert upd.upd_key("2") == upd.upd_key("chr2") == "chr2:UPD"
    assert upd.chrom_of_key("chrX:UPD") == "X"
    assert upd.chrom_from_evidence(_flag()) == "2"
    assert upd.chrom_from_evidence(_flag("no span given")) is None
    assert upd.child_of(_flag()) == "S3" and upd.is_upd_flag(_flag())


def test_default_origin_is_the_other_parent():
    assert upd.default_origin([_flag(a="S2")], "S1", "S2") == "paternal"
    assert upd.default_origin([_flag(a="S1")], "S1", "S2") == "maternal"
    assert upd.default_origin([_flag(a="S1"), _flag(a="S2")], "S1", "S2") == "unknown"
    assert upd.default_origin([], "S1", "S2") == "unknown"


def test_viewers_and_bad_values_are_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(upd, "can_write", lambda: False)
    with pytest.raises(PermissionError):
        _save(tmp_path)
    monkeypatch.setattr(upd, "can_write", lambda: True)
    for bad in (
        dict(validation="maybe"),
        dict(origin="grandmaternal"),
        dict(chrom=""),
        dict(user=""),
    ):
        with pytest.raises(ValueError):
            _save(tmp_path, **bad)
    assert not upd.upd_validation_file(tmp_path).exists()


def test_shared_loader_reads_the_file(tmp_path, curator):
    _save(tmp_path, origin="paternal", comment="SNP array")
    assert upd.upd_validation_file(tmp_path) == validation_path(tmp_path, "UPD")
    vmap = load_validation_map(upd.upd_validation_file(tmp_path), "1001")
    validation, origin, comment, *_ = vmap[("chr2:UPD", "S3")][0]
    assert (validation, origin, comment) == ("present", "paternal", "SNP array")


def test_latest_validation_skips_ignored_rows():
    rows = [
        {"Sample": "S3", "Validation": "present", "Ignore": "0"},
        {"Sample": "S3", "Validation": "absent", "Ignore": "1"},
        {"Sample": "S4", "Validation": "absent", "Ignore": "0"},
    ]
    assert upd.latest_validation(rows, "S3")["Validation"] == "present"
    assert upd.latest_validation(rows, "S5") is None


@pytest.mark.parametrize(
    "validation, status",
    [("present", "confirmed"), ("absent", "explained"), ("uncertain", None)],
)
def test_validation_acknowledges_the_flag(tmp_path, curator, validation, status):
    _save(tmp_path, validation=validation, origin="paternal", comment="array")
    ack = acks.load_acks(tmp_path).get(_flag())
    if status is None:
        assert ack is None
        return
    assert (ack.status, ack.source) == (status, "upd")
    assert ack.note == f"UPD validated {validation} (paternal): array"


def test_latest_validation_wins(tmp_path, curator):
    _save(tmp_path, validation="present", timestamp="2026-10-01T10:00:00")
    _save(tmp_path, validation="uncertain", timestamp="2026-10-02T10:00:00")
    assert _flag() not in acks.load_acks(tmp_path)


def test_flag_on_another_chromosome_or_child_is_not_acknowledged(tmp_path, curator):
    _save(tmp_path, chrom="8")
    standing = acks.load_acks(tmp_path)
    assert _flag() not in standing
    assert _flag(datagen.UPD_FLAG[4].replace("chr2", "chr8")) in standing
    assert _flag(datagen.UPD_FLAG[4].replace("chr2", "chr8"), b="S4") not in standing
    # No chromosome in the evidence: nothing says it is about another one.
    assert _flag("no span given") in standing


def test_newer_of_an_ack_and_a_validation_stands(tmp_path, curator, monkeypatch):
    monkeypatch.setattr(acks, "can_write", lambda: True)
    acks.save_flag_ack(
        tmp_path, _flag(), user="a", status="explained", timestamp="2026-10-01"
    )
    _save(tmp_path, validation="present", timestamp="2026-10-02")
    assert acks.load_acks(tmp_path).get(_flag()).status == "confirmed"


def test_evidence_of_the_upd_case(data_dir):
    datagen.add_upd_case(data_dir)
    ev = upd.gather_evidence(data_dir, "1001", "S3", "chr2")
    assert ev.roh and [f.sample_a for f in ev.flags] == ["S2"]
    assert (ev.father, ev.mother) == ("S1", "S2")
    assert ev.source == "ROH + Sample QC"
    assert upd.gather_evidence(data_dir, "1001", "S3", "8").source == "Curator"


def test_datagen_flag_matches_the_parser():
    assert upd.chrom_from_evidence(_flag()) == "2"
