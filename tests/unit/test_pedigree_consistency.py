"""Pedigree consistency across cohorts: the findings, aliases, acknowledgements."""

from pathlib import Path

import datagen
import pytest

from genetics_viz.models import DataStore
from genetics_viz.utils import pedigree_consistency as pc


def _ped(*rows):
    """``cohort -> FID -> member -> Row`` from (FID, IID, father, mother, sex, pheno)."""
    out = {}
    for fid, iid, fa, mo, sex, ph in rows:
        out.setdefault(fid, {})[iid] = (fa, mo, sex, ph)
    return out


TRIO = [("F1", "C", "P", "M", "1", "2"), ("F1", "P", "0", "0", "1", "1"),
        ("F1", "M", "0", "0", "2", "1")]  # fmt: skip


def codes(pedigrees):
    return [f.code for f in pc.scan_pedigrees(pedigrees)]


def test_agreeing_cohorts_give_nothing():
    assert codes({"A": _ped(*TRIO), "B": _ped(*TRIO)}) == []


def test_unlinked_parents_are_rows_that_differ():
    # PMS's C0733-011-021: the parents listed, the child's links at 0.
    unlinked = [("F1", "C", "0", "0", "1", "2")] + TRIO[1:]
    (f,) = pc.scan_pedigrees({"A": _ped(*TRIO), "B": _ped(*unlinked)})
    assert f.code == "SAME_FID_ROWS_DIFFER" and f.severity == "high"
    assert f.places == (("A", "F1"), ("B", "F1"))
    assert f.details == ("C: father P (A) / 0 (B); mother M (A) / 0 (B)",)


def test_dummy_parents_against_none_differ():
    # 420498: no parents in one cohort, shared dummy parents in the other.
    plain = [("F", "X", "0", "0", "1", "1"), ("F", "Y", "0", "0", "1", "1")]
    dummies = [("F", "X", "DUMMY_X_F", "DUMMY_X_M", "1", "1"),
               ("F", "Y", "DUMMY_X_F", "DUMMY_X_M", "1", "1"),
               ("F", "DUMMY_X_F", "0", "0", "1", "0"),
               ("F", "DUMMY_X_M", "0", "0", "2", "0")]  # fmt: skip
    (f,) = pc.scan_pedigrees({"A": _ped(*plain), "B": _ped(*dummies)})
    assert f.code == "SAME_FID_ROWS_DIFFER"
    assert "DUMMY_X_F: only in B" in f.details


def test_codes_are_normalised():
    floats = [("F1", "C", "P", "M", "1.0", "2.0"), ("F1", "P", "-9", "", "1.0", "1"),
              ("F1", "M", "0", ".", "2.0", "1.0")]  # fmt: skip
    assert (
        codes({"A": _ped(*TRIO), "B": pc.rows_from_lines("\t".join(r) for r in floats)})
        == []
    )


def test_members_differ_and_collision():
    # C0733-011-126: one cohort holds a single member of the other's family.
    lone = [("F1", "C", "0", "0", "1", "1")]
    (f,) = pc.scan_pedigrees({"A": _ped(*TRIO), "B": _ped(*lone)})
    assert f.code == "SAME_FID_MEMBERS_DIFFER"
    assert f.details[0] == "B: no M, P"
    assert (
        "C: father P (A) / 0 (B); mother M (A) / 0 (B); phenotype 2 (A) / 1 (B)"
        in f.details
    )
    other = [("F1", "Z", "0", "0", "1", "2")]
    assert codes({"A": _ped(*TRIO), "B": _ped(*other)}) == ["FID_COLLISION"]


def test_phenotype_alone_is_info():
    adhd = [("F1", "C", "P", "M", "1", "1")] + TRIO[1:]
    (f,) = pc.scan_pedigrees({"A": _ped(*TRIO), "B": _ped(*adhd)})
    assert f.code == "SAME_FID_PHENOTYPE_DIFFERS" and f.severity == "info"
    assert not f.acknowledgeable


def test_aliases_by_shared_barcode():
    alias = [(r[0].replace("F1", "X9"), *r[1:]) for r in TRIO]
    assert codes({"A": _ped(*TRIO), "B": _ped(*alias)}) == ["ALIAS_SAME_MEMBERS"]
    # 820910: the child alone under another ID, parents not there.
    lone = [("K", "C", "0", "0", "1", "2")]
    assert codes({"A": _ped(*TRIO), "B": _ped(*lone)}) == ["ALIAS_SUBSET"]
    # Each family has a member the other lacks.
    cross = [("K", "C", "P", "0", "1", "2"), ("K", "P", "0", "0", "1", "1"),
             ("K", "Q", "0", "0", "2", "1")]  # fmt: skip
    assert "ALIAS_OVERLAP" in codes({"A": _ped(*TRIO), "B": _ped(*cross)})


def test_a_sample_with_two_sexes_or_fathers():
    swapped = [("K", "C", "M", "P", "1", "2"), ("K", "P", "0", "0", "1", "1"),
               ("K", "M", "0", "0", "2", "1")]  # fmt: skip
    found = pc.scan_pedigrees({"A": _ped(*TRIO), "B": _ped(*swapped)})
    assert [f.code for f in found] == ["BARCODE_ROWS_DIFFER", "ALIAS_SAME_MEMBERS"]
    assert "C: father P (F1) / M (K)" in found[0].details


def test_dummies_link_no_families():
    a = [("F", "X", "DUMMY_X_F", "0", "1", "2"), ("F", "DUMMY_X_F", "0", "0", "1", "0")]
    b = [("G", "Y", "DUMMY_X_F", "0", "1", "2"), ("G", "DUMMY_X_F", "0", "0", "1", "0")]
    assert codes({"A": _ped(*a), "B": _ped(*b)}) == []


def test_alias_groups_join_same_and_subset_only():
    alias = [(r[0].replace("F1", "X9"), *r[1:]) for r in TRIO]
    lone = [("K", "C", "0", "0", "1", "2")]
    found = pc.Consistency(
        pc.scan_pedigrees({"A": _ped(*TRIO), "B": _ped(*alias), "C": _ped(*lone)})
    )
    assert found.alias_groups()["F1"] == {"F1", "X9", "K"}
    # K is part of F1 and of X9, its alias.
    assert [f.code for f in found.for_family("K")] == ["ALIAS_SUBSET"] * 2


@pytest.fixture
def data_dir(tmp_path: Path) -> Path:
    return datagen.build_data_dir(tmp_path)


def _store(data_dir: Path) -> DataStore:
    store = DataStore(data_dir)
    store.load()
    return store


def test_the_synthetic_data_agree_by_default(data_dir):
    assert pc.scan(_store(data_dir)).findings == []


def test_scan_reads_the_cohort_files(data_dir):
    datagen.add_cohort_disagreements(data_dir)
    found = pc.scan(_store(data_dir))
    assert sorted((f.code, sorted(f.family_ids)) for f in found.findings) == [
        ("ALIAS_SAME_MEMBERS", ["AU-77", "AU-X-1"]),
        ("ALIAS_SUBSET", ["1001", "K3"]),
        ("SAME_FID_ROWS_DIFFER", ["1001"]),
    ]
    (rows,) = [f for f in found.findings if f.code == "SAME_FID_ROWS_DIFFER"]
    assert rows.cohorts == ("DEMO", "NOHEAD")


def test_acknowledgement_and_a_change_of_rows(data_dir, monkeypatch):
    datagen.add_cohort_disagreements(data_dir)
    monkeypatch.setattr(pc, "can_write", lambda: True)
    finding = pc.scan(_store(data_dir)).for_family("AU-X-1")[0]
    pc.save_ack(data_dir, finding, "centre code", user="alice")
    acks = pc.load_acks(data_dir)
    assert acks[finding.key].note == "centre code"
    assert pc.open_findings([finding], acks) == []
    # The alias gains a member: the explanation no longer covers it.
    path = data_dir / "cohorts" / "ALIAS" / "ALIAS.pedigree.tsv"
    path.write_text(path.read_text() + "AU-77\tS4\t0\t0\t2\t1\n")
    changed = pc.scan(_store(data_dir)).for_family("AU-X-1")[0]
    assert changed.code == "ALIAS_OVERLAP" or changed.key != finding.key
    assert pc.open_findings([changed], pc.load_acks(data_dir)) == [changed]
    pc.save_ack(data_dir, finding, "", user="alice", status=pc.WITHDRAWN)
    assert finding.key not in pc.load_acks(data_dir)


def test_acknowledging_is_for_curators(data_dir, monkeypatch):
    monkeypatch.setattr(pc, "can_write", lambda: False)
    finding = pc.Finding("ALIAS_SAME_MEMBERS", (("A", "F"), ("B", "G")))
    with pytest.raises(PermissionError):
        pc.save_ack(data_dir, finding, "", user="eve")
