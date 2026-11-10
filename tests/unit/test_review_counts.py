"""What the family page's tab badges count (utils/review_counts.py)."""

import datagen
import polars as pl
import pytest

from genetics_viz.utils import review_counts as rc
from genetics_viz.utils import str_validation, upd
from genetics_viz.utils.sample_qc import load_family_qc
from genetics_viz.utils.sample_qc_acks import load_acks, split_acknowledged
from genetics_viz.utils.str_eh import (
    eh_path,
    find_eh_labels,
    flagged_expr,
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
