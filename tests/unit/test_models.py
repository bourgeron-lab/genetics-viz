"""Pedigree loading and cohort discovery (models.py)."""

import os
from pathlib import Path

import pytest

from genetics_viz.models import Cohort, DataStore, Sample

HEADER = "FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n"
ROWS = "F1\tA\t0\t0\t1\t1\nF1\tB\t0\t0\t2\t1\nF1\tC\tA\tB\t1\t2\n"


def make_cohort(tmp_path: Path, name: str, text: str) -> Cohort:
    d = tmp_path / name
    d.mkdir()
    (d / f"{name}.pedigree.tsv").write_text(text)
    return Cohort.from_directory(d)


def test_header_and_headerless_pedigrees_load_identically(tmp_path):
    with_header = make_cohort(tmp_path, "H", HEADER + ROWS)
    without = make_cohort(tmp_path, "N", ROWS)
    assert with_header.dataframe.equals(without.dataframe)
    assert with_header.dataframe.columns == [
        "FID",
        "IID",
        "PAT",
        "MAT",
        "SEX",
        "PHENOTYPE",
    ]
    assert with_header.num_samples == 3


def test_extra_columns_are_named_by_position(tmp_path):
    cohort = make_cohort(
        tmp_path, "X", HEADER.rstrip("\n") + "\tEXTRA\n" + "F1\tA\t0\t0\t1\t1\tz\n"
    )
    assert cohort.dataframe.columns[-1] == "COL7"


@pytest.mark.parametrize("missing", ["", "0", "-9", "0.0", "-9.0"])
def test_missing_parent_sentinels_mean_no_parent(tmp_path, missing):
    cohort = make_cohort(tmp_path, "M", f"F1\tA\t{missing}\t{missing}\t1\t1\n")
    sample = cohort.families["F1"].samples[0]
    assert sample.father_id is None and sample.mother_id is None
    assert sample.is_founder


def test_sex_and_phenotype_zero_are_kept(tmp_path):
    cohort = make_cohort(tmp_path, "Z", "F1\tA\t0\t0\t0\t0\n")
    sample = cohort.families["F1"].samples[0]
    assert (sample.sex, sample.phenotype) == ("0", "0")


def test_ids_stay_strings(tmp_path):
    cohort = make_cohort(tmp_path, "S", "1001\t007\t0\t0\t1\t1\n")
    assert list(cohort.families) == ["1001"]
    assert cohort.families["1001"].samples[0].sample_id == "007"


def test_family_members_render_missing_as_dash(tmp_path):
    cohort = make_cohort(tmp_path, "D", HEADER + ROWS)
    members = cohort.get_family_members("F1")
    assert members[0] == {
        "Sample ID": "A",
        "Father": "-",
        "Mother": "-",
        "Sex": "1",
        "Phenotype": "1",
    }
    assert members[2]["Father"] == "A"
    assert cohort.get_family_members("nope") == []
    assert cohort.families["F1"].num_founders == 2


def test_founder_definition():
    assert Sample("A", "F").is_founder
    assert not Sample("A", "F", father_id="B").is_founder


def test_one_column_pedigree_is_rejected(tmp_path):
    with pytest.raises(ValueError):
        make_cohort(tmp_path, "B", "only\n")


def test_missing_pedigree_raises(tmp_path):
    (tmp_path / "E").mkdir()
    with pytest.raises(FileNotFoundError):
        Cohort.from_directory(tmp_path / "E")


def test_headerless_pedigree_whose_first_family_starts_with_fid(tmp_path):
    cohort = make_cohort(
        tmp_path, "P", "FID001\tA\t0\t0\t1\t1\nFID001\tB\t0\t0\t2\t1\n"
    )
    assert cohort.num_samples == 2


# --- DataStore -----------------------------------------------------------------


def test_datastore_sorts_cohorts_into_explorable_and_incomplete(data_dir):
    store = DataStore(data_dir=data_dir)
    store.load()
    assert sorted(store.cohorts) == ["DEMO", "NOHEAD"]
    stubs = {s.name: s for s in store.incomplete_cohorts}
    assert stubs["EMPTY"].reason == "missing"
    assert stubs["BROKEN"].reason == "unreadable" and stubs["BROKEN"].error
    demo = store.cohorts["DEMO"]
    assert demo.params_file is not None and demo.params_file.name == "DEMO.params.yml"
    assert set(demo.families) == {"1001", "1002", "AU-X-1"}


def test_datastore_without_cohorts_dir_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        DataStore(data_dir=tmp_path).load()


def bump(path: Path) -> None:
    st = path.stat()
    os.utime(path, (st.st_atime, st.st_mtime + 10))


def test_snapshot_detects_changes(data_dir):
    store = DataStore(data_dir=data_dir)
    store.load()
    before = store.take_snapshot()

    bump(data_dir / "cohorts" / "DEMO" / "DEMO.pedigree.tsv")
    (data_dir / "cohorts" / "NEW").mkdir()
    (data_dir / "cohorts" / "EMPTY").rmdir()
    report = DataStore.compare_snapshot(before, store.take_snapshot())

    assert report.modified_cohorts == ["DEMO"]
    assert report.added_cohorts == ["NEW"]
    assert report.removed_cohorts == ["EMPTY"]
    assert report.has_changes and len(report.summary_lines()) == 3


def test_snapshot_watches_params_but_not_the_run_record(data_dir):
    store = DataStore(data_dir=data_dir)
    store.load()
    before = store.take_snapshot()
    (data_dir / "cohorts" / "DEMO" / ".ghfc-ngs.state.json").write_text("{}")
    assert not DataStore.compare_snapshot(before, store.take_snapshot()).has_changes
    bump(data_dir / "cohorts" / "DEMO" / "DEMO.params.yml")
    assert DataStore.compare_snapshot(
        before, store.take_snapshot()
    ).modified_cohorts == ["DEMO"]


def test_gaining_a_pedigree_makes_a_cohort_explorable(data_dir):
    store = DataStore(data_dir=data_dir)
    store.load()
    (data_dir / "cohorts" / "EMPTY" / "EMPTY.pedigree.tsv").write_text(ROWS)
    store.reload()
    assert "EMPTY" in store.cohorts
    assert "EMPTY" not in {s.name for s in store.incomplete_cohorts}


def test_hidden_directories_are_not_cohorts(data_dir):
    (data_dir / "cohorts" / ".snapshot").mkdir()
    store = DataStore(data_dir=data_dir)
    store.load()
    assert ".snapshot" not in {s.name for s in store.incomplete_cohorts}
    assert ".snapshot" not in store.cohorts
