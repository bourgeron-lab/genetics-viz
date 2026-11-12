"""Diagnostic kinds (utils/diagnostic_files.py) and the STR kind's readers."""

import pytest

from genetics_viz.components import diagnostic_loader as dl
from genetics_viz.utils import str_eh
from genetics_viz.utils.diagnostic_files import (
    DIAGNOSTIC_TYPES,
    diagnostic_path,
    diagnostic_paths,
)

KEY = "chr4:3074876-3074933:HTT"


@pytest.fixture(autouse=True)
def curator(monkeypatch):
    monkeypatch.setattr(dl, "can_write", lambda: True)


def _save(data_dir, dtype, variant, sample="S3", verdict="pathogenic", ts="t1"):
    dl.save_diagnostic_entry(
        diagnostic_path(data_dir, dtype), "1001", variant, "HTT", "pathogenic",
        sample, "curator", ts, "", verdict,
    )  # fmt: skip


def test_kinds_include_str_and_upd(tmp_path):
    assert DIAGNOSTIC_TYPES == ["SNV", "SV", "STR", "UPD"]
    assert diagnostic_path(tmp_path, "STR") == tmp_path / "diagnostics" / "strs.tsv"
    assert [t for t, _ in diagnostic_paths(tmp_path)] == DIAGNOSTIC_TYPES


def test_family_diagnostics_read_every_kind(tmp_path):
    _save(tmp_path, "SNV", "chr1:1:A:T")
    _save(tmp_path, "STR", KEY)
    entries = dl.load_family_diagnostics(tmp_path, "1001", ["S3"])
    assert {(e["_source"], e["Variant"]) for e in entries} == {
        ("snv", "chr1:1:A:T"),
        ("str", KEY),
    }


def test_diagnostic_badge_in_the_eh_table(data_dir):
    from genetics_viz.pages.cohort.components import str_tab

    fam = str_eh.read_eh_table(
        str_eh.eh_path(str_eh.get_family_eh_dir(data_dir, "1001"), "1001",
                       "jmsa50", "genotypes")
    )  # fmt: skip
    _save(data_dir, "STR", KEY)
    statuses = str_tab._diagnostic_statuses(data_dir, "1001", fam)
    assert statuses == {("HTT", "S3"): "pathogenic"}

    members = str_eh.family_members(fam)
    rows = str_eh.family_locus_rows(
        fam, members, [m.sample for m in members], diagnostics=statuses
    )
    htt = next(r for r in rows if r["locus"] == "HTT")
    s3 = next(i for i, m in enumerate(members) if m.sample == "S3")
    labels = [b["label"] for b in htt[f"{str_eh.member_key(s3)}_badges"]]
    assert "Dx P" in labels


def test_locus_gene():
    assert [str_eh.locus_gene(x) for x in ("HTT", "FXN_A", "NUTM2B-AS1")] == [
        "HTT",
        "FXN",
        "NUTM2B-AS1",
    ]


# A transmitted duplication, called with other breakpoints in each carrier.
CHILD_CALL = "chr2:1020001-1100000:dup"
FATHER_CALL = "chr2:1050001-1090000:dup"
CURATED = "chr2:1004558-1100312:dup"


def _sv_diag(data_dir, variant, sample, genes, ts):
    dl.save_diagnostic_entry(
        diagnostic_path(data_dir, "SV"), "1001", variant, genes, "GAIN",
        sample, "curator", ts, "", "uncertain",
    )  # fmt: skip


def _sv_validation(data_dir, variant, sample, start, end, status="present"):
    from genetics_viz.components.validation_loader import SV_VALIDATION_HEADER
    from genetics_viz.utils.tsv import append_tsv_rows
    from genetics_viz.utils.validation_files import validation_path

    record = {
        "FID": "1001", "Variant": variant, "Sample": sample, "User": "curator",
        "Inheritance": "paternal", "Validation": status, "CuratedStart": start,
        "CuratedEnd": end, "Comment": "", "Ignore": "0", "Timestamp": "t0",
    }  # fmt: skip
    append_tsv_rows(
        validation_path(data_dir, "SV"),
        SV_VALIDATION_HEADER,
        [[record[c] for c in SV_VALIDATION_HEADER]],
    )


def test_svs_curated_to_one_locus_share_a_row(tmp_path):
    _sv_diag(tmp_path, CHILD_CALL, "S3", "GENE1:exonic,GENE2:exonic", "t2")
    _sv_diag(tmp_path, FATHER_CALL, "S1", "GENE1:intronic", "t1")
    _sv_validation(tmp_path, CHILD_CALL, "S3", "1004558", "1100312")
    _sv_validation(tmp_path, FATHER_CALL, "S1", "1004558", "1100312")

    [entry] = dl.load_family_diagnostics(tmp_path, "1001", ["S1", "S3"])
    assert entry["_shown_variant"] == CURATED and entry["_curated"]
    assert entry["Sample"] == "S1, S3"
    assert dict(entry["_sample_variants"]) == {"S3": CHILD_CALL, "S1": FATHER_CALL}
    # The newest record's overlap wins for a gene both calls list.
    assert entry["Gene"] == "GENE1:exonic,GENE2:exonic"


def test_uncurated_or_absent_sv_keeps_its_call(tmp_path):
    _sv_diag(tmp_path, CHILD_CALL, "S3", "GENE1", "t2")
    _sv_diag(tmp_path, FATHER_CALL, "S1", "GENE1", "t1")
    _sv_validation(tmp_path, CHILD_CALL, "S3", "1004558", "")
    _sv_validation(tmp_path, FATHER_CALL, "S1", "1004558", "1100312", "absent")

    entries = dl.load_family_diagnostics(tmp_path, "1001", ["S1", "S3"])
    shown = {e["Sample"]: (e["_shown_variant"], e["_curated"]) for e in entries}
    # A curated start alone keeps the call's end.
    assert shown == {
        "S3": ("chr2:1004558-1100000:dup", True),
        "S1": (FATHER_CALL, False),
    }


def test_report_shows_the_curated_sv_and_each_carriers_validation(tmp_path):
    from genetics_viz.utils import family_report

    _sv_diag(tmp_path, CHILD_CALL, "S3", "GENE1", "t2")
    _sv_diag(tmp_path, FATHER_CALL, "S1", "GENE1", "t1")
    _sv_validation(tmp_path, CHILD_CALL, "S3", "1004558", "1100312")
    _sv_validation(tmp_path, FATHER_CALL, "S1", "1004558", "1100312")

    entries = dl.load_family_diagnostics(tmp_path, "1001", ["S1", "S3"])
    [row] = family_report._diagnostics(tmp_path, "1001", entries).rows
    assert row["finding"].startswith("DUP chr2:1004558-1100312")
    assert row["samples"] == "S1, S3"
    assert row["validation"] == "S1: present (paternal); S3: present (paternal)"
    assert "curated" in row["_finding_note"]
