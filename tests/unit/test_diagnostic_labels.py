"""What the family Diagnostics panel says about STR and UPD records."""

import pytest

from genetics_viz.utils import diagnostic_labels as dlabels
from genetics_viz.utils import str_validation, upd

KEY = "chr4:3074876-3074933:HTT"


@pytest.fixture
def curator(monkeypatch):
    monkeypatch.setattr(str_validation, "can_write", lambda: True)
    monkeypatch.setattr(upd, "can_write", lambda: True)


def _entry(variant=KEY, sample="S3", impact="pathogenic"):
    return {"Variant": variant, "Sample": sample, "Impact": impact}


def test_str_count_is_the_longer_allele(data_dir):
    text, tooltip = dlabels.describe_str(data_dir, "1001", _entry())
    assert text == "~40 repeats"
    assert "Motif: CAG" in tooltip and "genotype 18/40" in tooltip


def test_curated_repeats_beat_the_estimate(data_dir, curator):
    str_validation.save_str_validation(
        data_dir, family_id="1001", variant_key=KEY, sample="S3", user="c",
        validation="present", curated_repeats="43",
    )  # fmt: skip
    text, tooltip = dlabels.describe_str(data_dir, "1001", _entry())
    assert text == "43 repeats (curated)"
    assert "S3: 43 repeats (curated" in tooltip


def test_several_members_and_unknown_locus(data_dir):
    text, _ = dlabels.describe_str(data_dir, "1001", _entry(sample="S1, S3"))
    assert text.endswith("S3 ~40 repeats") and text.startswith("S1 ~")
    text, tooltip = dlabels.describe_str(data_dir, "1001", _entry("chr1:1-2:NOPE"))
    assert text == "" and "no genotype" in tooltip


def test_upd_label_names_the_origin(data_dir, curator):
    entry = _entry("chr2:UPD", impact="ROH + Sample QC")
    assert dlabels.describe_upd(data_dir, "1001", entry)[0] == "chr2"
    upd.save_upd_validation(
        data_dir, family_id="1001", sample="S3", chrom="2", user="c",
        validation="present", origin="paternal",
    )  # fmt: skip
    text, tooltip = dlabels.describe_upd(data_dir, "1001", entry)
    assert text == "chr2 · paternal"
    assert "S3: validated present, paternal" in tooltip
    assert "Found by: ROH + Sample QC" in tooltip
