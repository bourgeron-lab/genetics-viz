"""What the family and cohort pages' tab badges count (utils/review_counts.py)."""

import dataclasses

import datagen
import polars as pl
import pytest

from genetics_viz.utils import review_counts as rc
from genetics_viz.utils import str_validation, upd
from genetics_viz.utils.cohort_state import read_cohort_state
from genetics_viz.utils.sample_qc import load_cohort_qc, load_family_qc
from genetics_viz.utils.sample_qc_acks import load_acks, split_acknowledged
from genetics_viz.utils.str_eh import (
    eh_path,
    find_eh_labels,
    flagged_expr,
    get_cohort_eh_dir,
    get_family_eh_dir,
    read_eh_table,
)


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(upd, "can_write", lambda: True)
    monkeypatch.setattr(str_validation, "can_write", lambda: True)


def _active(data_dir, fid="1001"):
    qc = load_family_qc(data_dir, fid)
    return split_acknowledged(qc.flags, load_acks(data_dir))[0]


def _explain_not_parent(data_dir):
    """Acknowledge 1001's NOT_PARENT as explained: the de novo table is read."""
    with open(data_dir / "validations" / "sampleqc.tsv", "a") as f:
        f.write(
            "\t".join(["NOT_PARENT", "1001", "S1", "S3", "explained", "Adopted",
                       "curator", "2026-11-10T10:00:00"]) + "\n"
        )  # fmt: skip


def _add_denovo(data_dir, pos, consequence):
    """One more de novo call of S3, on its own position."""
    path = datagen.family_dir(data_dir, "1001") / "wombat"
    path = path / "1001.rare.vep1.annotated.rare_denovo.tsv"
    row = list(datagen.WOMBAT_ROWS[0])
    row[1] = pos
    row[datagen.WOMBAT_COLUMNS.index("VEP_Consequence")] = consequence
    with open(path, "a") as f:
        f.write("\t".join(map(str, row)) + "\n")


def test_qc_counts_high_and_medium(data_dir):
    count = rc.qc_to_review(_active(data_dir), upd_tab=False)
    assert (count.n, count.worst) == (2, "high")
    assert count.tooltip == "Sample QC: 1 high · 1 medium not acknowledged"
    assert rc.qc_family_to_review(data_dir, "1001", upd_tab=False).n == 2
    assert not rc.qc_family_to_review(data_dir, "1002", upd_tab=False).shown


def test_qc_leaves_upd_flags_to_the_upd_tab(data_dir):
    datagen.add_upd_case(data_dir)
    assert rc.qc_to_review(_active(data_dir), upd_tab=False).n == 3
    count = rc.qc_to_review(_active(data_dir), upd_tab=True)
    assert count.n == 2 and "counted on the UPD tab" in count.tooltip


def test_denovo_held_back_by_a_wrong_relationship(data_dir):
    count = rc.denovo_to_review(data_dir, "1001")
    assert count.held and count.shown and count.n == 0
    assert "held back" in count.tooltip


def test_denovo_counts_coding_calls_never_validated(data_dir, monkeypatch):
    _explain_not_parent(data_dir)
    _add_denovo(data_dir, 50674100, "synonymous_variant")  # LOW: not counted
    _add_denovo(data_dir, 50674200, "stop_gained")  # HIGH
    count = rc.denovo_to_review(data_dir, "1001")
    assert not count.held
    assert (count.n, count.worst) == (2, "high")  # SHANK3 missense + stop
    assert "(1 HIGH)" in count.tooltip
    # A validation, whatever its status, settles a call.
    validations = data_dir / "validations" / "snvs.tsv"
    validations.parent.mkdir(exist_ok=True)
    datagen.write_tsv(
        validations,
        ["FID", "Variant", "Sample", "User", "Inheritance", "Validation",
         "Comment", "Ignore", "Timestamp"],
        [["1001", "chr22:50674200:G:A", "S3", "curator", "de novo", "uncertain",
          "", "0", "2026-11-10T10:00:00"]],
    )  # fmt: skip
    count = rc.denovo_to_review(data_dir, "1001")
    assert (count.n, count.worst) == (1, "medium")


def test_denovo_nothing_without_a_denovo_config(data_dir):
    assert not rc.denovo_to_review(data_dir, "1002").shown


def test_str_counts_flagged_genotypes_never_validated(data_dir, curator):
    count = rc.str_to_review(data_dir, "1001")
    assert (count.n, count.worst) == (2, "high")
    directory = get_family_eh_dir(data_dir, "1001")
    label = find_eh_labels(directory, "1001")[0]
    fam = read_eh_table(eh_path(directory, "1001", label, "genotypes"))
    high = fam.filter(flagged_expr(fam) & (pl.col("severity") == "high")).to_dicts()[0]
    str_validation.save_str_validation(
        data_dir,
        family_id="1001",
        variant_key=str_validation.str_variant_key(
            high["chrom"], high["pos"], high["end"], high["locus"]
        ),
        sample=high["sample"],
        user="curator",
        validation="present",
    )
    count = rc.str_to_review(data_dir, "1001")
    assert (count.n, count.worst) == (1, "medium")


def test_upd_count(data_dir):
    assert not rc.upd_to_review(data_dir, "1001").shown
    datagen.add_upd_calls(data_dir)
    assert rc.upd_to_review(data_dir, "1001").n == 2


# ---------------------------------------------------------------------------
# Cohort page badges
# ---------------------------------------------------------------------------


class DemoCohort:
    name = "DEMO"

    def __init__(self, data_dir):
        self.path = data_dir / "cohorts" / "DEMO"


ALL_DEMO = {
    ("1001", "S1"),
    ("1001", "S2"),
    ("1001", "S3"),
    ("1002", "S4"),
    ("AU-X-1", "S5"),
    ("AU-X-1", "S6"),
}


def _cohort_active(data_dir):
    qc = load_cohort_qc(DemoCohort(data_dir))
    return split_acknowledged(qc.flags, load_acks(data_dir))[0]


def test_filter_pill_only_when_filtered():
    assert rc.filter_pills(6, 6, 0) == []
    (pill,) = rc.filter_pills(2, 6, 1)
    assert pill.text == "2" and "2 of 6 individuals shown (1 filter)" in pill.tooltip


def test_cohort_qc_counts_individuals_by_worst_flag(data_dir):
    datagen.add_upd_case(data_dir)
    worst = rc.qc_worst_by_individual(_cohort_active(data_dir), upd_tab=False)
    # S2's only flag is its UPD candidate: high without a UPD tab, left out with one
    assert worst[("1001", "S2")] == "high"
    with_tab = rc.qc_worst_by_individual(_cohort_active(data_dir), upd_tab=True)
    assert ("1001", "S2") not in with_tab
    pills = rc.qc_cohort_pills(with_tab, ALL_DEMO, upd_tab=True)
    assert [p.text for p in pills] == [str(sum(v == "high" for v in with_tab.values()))]
    assert "counted on the UPD tab" in pills[0].tooltip
    # Only the shown individuals count
    assert rc.qc_cohort_pills(with_tab, {("1002", "S4")}, upd_tab=True) == []


def test_cohort_upd_pills_split_by_confidence():
    pills = rc.upd_cohort_pills([["high"], ["medium", "medium"], []], 3)
    assert [p.text for p in pills] == ["1", "2"]
    assert pills[0].color != pills[1].color
    assert "high-confidence UPD candidate no curator settled" in pills[0].tooltip
    assert rc.upd_cohort_pills([[], []], 2) == []


def test_cohort_str_flags_drop_validated_genotypes(data_dir, curator):
    cohort = DemoCohort(data_dir)
    flags = rc.str_cohort_open_flags(cohort, data_dir)
    assert {(key, locus, sev) for key, locus, sev in flags} == {
        (("1001", "S3"), "HTT", "high"),
        (("1001", "S2"), "FMR1", "medium"),
        (("AU-X-1", "S6"), "ARX_2", "medium"),
    }
    assert [p.text for p in rc.str_cohort_pills(flags, ALL_DEMO)] == ["1", "2"]
    assert [p.text for p in rc.str_cohort_pills(flags, {("1001", "S2")})] == ["1"]

    directory = get_cohort_eh_dir(cohort)
    label = find_eh_labels(directory, "DEMO")[0]
    df = read_eh_table(eh_path(directory, "DEMO", label, "flags"))
    htt = df.filter(pl.col("locus") == "HTT").to_dicts()[0]
    str_validation.save_str_validation(
        data_dir,
        family_id="1001",
        variant_key=str_validation.str_variant_key(
            htt["chrom"], htt["pos"], htt["end"], htt["locus"]
        ),
        sample="S3",
        user="curator",
        validation="present",
    )
    flags = rc.str_cohort_open_flags(cohort, data_dir)
    assert [p.text for p in rc.str_cohort_pills(flags, ALL_DEMO)] == ["2"]


def test_workflow_pill_from_the_run_record(data_dir):
    state = read_cohort_state(DemoCohort(data_dir).path)
    (pill,) = rc.workflow_pills(state, "2/2 steps complete")
    assert (pill.text, pill.icon, pill.color) == ("", "check_circle", "green")
    assert pill.tooltip == "Last run: success (2/2 steps complete)"
    assert rc.workflow_pills(None) == []
    unknown = dataclasses.replace(
        state, last_run=dataclasses.replace(state.last_run, status="not-a-status")
    )
    assert rc.workflow_pills(unknown) == []
