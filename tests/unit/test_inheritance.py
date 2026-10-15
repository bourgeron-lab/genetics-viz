"""Inheritance: transmission from genotypes, and fit with the gene's modes."""

import datagen
import pytest
from genetics_viz.models import Sample
from genetics_viz.utils import inheritance as inh


@pytest.mark.parametrize(
    "child, father, mother, expected",
    [
        ("0/1", "0/0", "0/0", "de novo"),
        ("0/1", "0/1", "0/0", "paternal"),
        ("0|1", "0/0", "1|0", "maternal"),  # phased
        ("0/1", "0/1", "0/1", "either"),
        ("1/1", "0/1", "0/1", "homozygous"),
        ("1/1", "0/0", "0/1", "unknown"),
        ("0/1", "0/0", "./.", "not paternal"),
        ("0/1", "./.", "0/0", "not maternal"),
        ("0/1", "0/1", "", "paternal"),
        ("0/0", "0/1", "0/1", "unknown"),
        ("0/1,1/1", "0/0", "0/0", "unknown"),  # aggregated, disagreeing
        ("0/1,0/1", "0/0", "0/0", "de novo"),  # aggregated, agreeing
    ],
)
def test_infer_inheritance(child, father, mother, expected):
    assert inh.infer_inheritance(child, father, mother) == expected


def obs(**kw):
    base = dict(
        chrom="chr1", pos=100, sample_gt="0/1", father_gt="0/0", mother_gt="0/0"
    )
    base.update(kw)
    return inh.observe(**base)


def test_observe_x_in_a_male():
    o = obs(chrom="chrX", pos=70_000_000, sample_gt="0/1", mother_gt="0/1", sex="1")
    assert (o.zygosity, o.origin) == ("hemi", "maternal")
    par = obs(chrom="X", pos=100_000, sample_gt="0/1", sex="1.0")
    assert par.zygosity == "het"  # PAR1: diploid in a male
    female = obs(chrom="chrX", pos=70_000_000, sex="2")
    assert female.zygosity == "het"


def test_observe_de_novo_caller_rows():
    assert obs(sample_gt=None, de_novo=True).origin == "de novo"
    o = obs(chrom="chrX", pos=70_000_000, sample_gt=None, de_novo=True, sex="1")
    assert o.zygosity == "hemi"


def verdict(modes, o, hits=()):
    v = inh.assess_moi(modes, o, list(hits))
    return v and (v.status, v.mode)


def test_dominant():
    assert verdict(["AD"], obs()) == ("consistent", "AD")
    unaffected = obs(father_gt="0/1", father_pheno="1")
    assert verdict(["AD"], unaffected) == ("possible", "AD")
    affected = obs(father_gt="0/1", father_pheno="2")
    assert verdict(["AD"], affected) == ("consistent", "AD")


def test_recessive():
    hom = obs(sample_gt="1/1", father_gt="0/1", mother_gt="0/1")
    assert verdict(["AR"], hom) == ("consistent", "AR")
    assert verdict(["AR"], obs(sample_gt="1/1")) == ("possible", "AR")
    paternal = obs(father_gt="0/1")
    assert verdict(["AR"], paternal) == ("inconsistent", "AR")
    trans = [inh.SecondHit("chr1:200", "maternal")]
    assert verdict(["AR"], paternal, trans) == ("consistent", "AR")
    cis = [inh.SecondHit("chr1:200", "paternal")]
    assert verdict(["AR"], paternal, cis) == ("inconsistent", "AR")
    unphased = [inh.SecondHit("chr1:200", "unknown")]
    assert verdict(["AR"], paternal, unphased) == ("possible", "AR")


def test_x_linked():
    male = obs(chrom="chrX", pos=70_000_000, sample_gt="1/1", mother_gt="0/1", sex="1")
    assert verdict(["XLR"], male) == ("consistent", "XLR")
    female = obs(chrom="chrX", pos=70_000_000, mother_gt="0/1", sex="2")
    assert verdict(["XLR"], female) == ("possible", "XLR")
    assert verdict(["XLD"], female) == ("consistent", "XLD")


def test_union_of_modes_takes_the_best():
    paternal = obs(father_gt="0/1")
    assert verdict(["AR", "AD"], paternal) == ("consistent", "AD")
    assert verdict(["MT"], paternal) is None
    assert verdict([], paternal) is None
    assert verdict(["AD"], obs(sample_gt="0/0")) is None


def test_gene_modes_pool_the_sources(snapshot):
    codes, sources = inh.gene_modes(snapshot, "BOTH")
    assert codes == ["AR", "AD"] and len(sources) == 2
    assert inh.gene_modes(snapshot, "PLAIN") == ([], [])


def test_annotate_a_trio_table(snapshot):
    rows = [dict(zip(datagen.WOMBAT_COLUMNS, r)) for r in datagen.WOMBAT_ROWS]
    samples = [
        Sample(
            iid,
            fid,
            father if father != "0" else None,
            mother if mother != "0" else None,
            sex,
            pheno,
        )
        for fid, iid, father, mother, sex, pheno in datagen.DEMO_PEDIGREE
    ]
    inh.annotate_moi(rows, snapshot, pedigree=inh.family_pedigree(samples))
    by_gene = {r["VEP_SYMBOL"]: r for r in rows}
    assert by_gene["SHANK3"]["MOI"] == "AD ✓"
    assert by_gene["KIAA1109"]["MOI"] == "AR ✓"  # trans with the BLTP1 row
    assert "in trans" in by_gene["KIAA1109"]["MoiBadges"][0]["tooltip"]
    assert by_gene["BLTP1"]["MOI"] == "AR ✓"
    assert by_gene["XRG"]["MOI"] == "XLR ✓"
    assert by_gene["ARONLY"]["MOI"] == "AR ✗"
    assert by_gene["PLAIN"]["MoiBadges"] == [] and by_gene["PLAIN"]["_moi_rank"] == 0
    assert by_gene["SHANK3"]["_moi_rank"] > by_gene["ARONLY"]["_moi_rank"] > 0


def test_every_genetrek_mode_is_mapped(snapshot):
    for info in snapshot.info.values():
        for _, raw in info.inheritance:
            assert raw in inh.INHERITANCE_CONFIG["modes"], raw
