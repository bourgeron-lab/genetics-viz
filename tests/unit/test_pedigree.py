"""Per-family pedigree loading (utils/pedigree.py)."""

from pathlib import Path

import pytest

from genetics_viz.utils.pedigree import load_family_object, load_family_pedigree

#: As the ghfc-ngs pipeline writes <FID>.pedigree.tsv: no header, numeric FID.
HEADERLESS = (
    "88181200007\tIP00FM8\t0\t0\t1\t1\n"
    "88181200007\tIP00FM9\t-9.0\t\t2.0\t-9.0\n"
    "88181200007\tIP00FMA\tIP00FM8\tIP00FM9\t1\t2\n"
)

EXPECTED = [
    {"Sample ID": "IP00FM8", "Father": "-", "Mother": "-", "Sex": "1", "Phenotype": "1"},
    {"Sample ID": "IP00FM9", "Father": "-", "Mother": "-", "Sex": "2.0", "Phenotype": "-9.0"},
    {"Sample ID": "IP00FMA", "Father": "IP00FM8", "Mother": "IP00FM9", "Sex": "1", "Phenotype": "2"},
]  # fmt: skip


def write(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "F.pedigree.tsv"
    path.write_text(text)
    return path


@pytest.mark.parametrize(
    "text",
    [
        pytest.param(HEADERLESS, id="headerless"),
        pytest.param("FID\tIID\tPAT\tMAT\tSEX\tPHENOTYPE\n" + HEADERLESS, id="plink"),
        pytest.param(
            "#family_id\tsample_id\tfather_id\tmother_id\tsex\tphenotype\n"
            + HEADERLESS,
            id="hash-prefixed",
        ),
    ],
)
def test_family_pedigree_with_and_without_header(tmp_path, text):
    assert load_family_pedigree(write(tmp_path, text)) == EXPECTED


def test_per_family_barcode_header_is_read_by_name(tmp_path):
    # Columns out of the positional order, named the per-family way.
    text = (
        "sample_id\tFatherBarcode\tMotherBarcode\tSex\tPheno_ASD\n"
        "IP00FM8\t0\t0\t1\t1\n"
        "IP00FM9\t0\t0\t2.0\t-9.0\n"
        "IP00FMA\tIP00FM8\tIP00FM9\t1\t2\n"
    )
    members = load_family_pedigree(write(tmp_path, text))
    assert [m["Sample ID"] for m in members] == ["IP00FM8", "IP00FM9", "IP00FMA"]
    assert members[2]["Father"] == "IP00FM8" and members[2]["Phenotype"] == "2"


def test_headerless_family_whose_id_starts_with_fid(tmp_path):
    members = load_family_pedigree(write(tmp_path, "FID001\tA\t0\t0\t1\t1\n"))
    assert [m["Sample ID"] for m in members] == ["A"]


def test_family_object_from_headerless_pedigree(tmp_path):
    family = load_family_object(write(tmp_path, HEADERLESS), "88181200007")
    assert family.num_samples == 3
    assert family.num_founders == 2  # "-9.0" and "" are no parent
    child = family.get_sample("IP00FMA")
    assert child is not None and child.family_id == "88181200007"
    assert (child.father_id, child.mother_id) == ("IP00FM8", "IP00FM9")


def test_header_without_a_sample_column_yields_no_members(tmp_path):
    text = "FID\tPAT\tMAT\n88181200007\t0\t0\n"
    assert load_family_pedigree(write(tmp_path, text)) == []
    assert load_family_object(write(tmp_path, text), "F").num_samples == 0
