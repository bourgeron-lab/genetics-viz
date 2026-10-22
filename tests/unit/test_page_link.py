"""Shareable links (components/page_link.py and the pages' codecs)."""

import pytest

from genetics_viz.components import page_link as pl
from genetics_viz.utils.vep import VEP_CONSEQUENCES


def test_defaults_stay_out_of_the_link():
    link = pl.PageLink("/p", {})
    link.set("tab", "wombat", default="wombat")
    link.set("m", ["S2", "S1"], default=["S1", "S2"])
    link.set("lcr", True, default=True)
    assert link.url() == "/p"
    link.set("tab", "roh", default="wombat")
    link.set("lcr", False, default=True)
    link.set("m", ["S3"], default=["S1", "S2"])
    assert link.url() == "/p?tab=roh&lcr=0&m=S3"


def test_an_emptied_list_differs_from_an_absent_one():
    link = pl.PageLink("/p", {})
    link.set("gs", [], default=["a"])
    assert link.params["gs"] == "-"
    back = pl.PageLink("/p", link.params)
    assert back.get_list("gs") == [] and back.get_list("other") is None


def test_reading():
    link = pl.PageLink("/p", {"a": "x,y", "b": "0", "c": "", "d": "1"})
    assert link.get_list("a") == ["x", "y"]
    assert link.get_bool("b", True) is False and link.get_bool("d", False) is True
    assert link.get_bool("c", True) is True  # empty values are dropped
    assert link.get("zz", "dflt") == "dflt"


def test_consequences_shorten_to_levels():
    high = [t for t, (imp, _) in VEP_CONSEQUENCES.items() if imp == "HIGH"]
    low = [t for t, (imp, _) in VEP_CONSEQUENCES.items() if imp == "LOW"]
    encoded = pl.encode_consequences(high + low[:1])
    assert encoded == ["HIGH", low[0]]
    assert set(pl.decode_consequences(encoded)) == set(high + low[:1])
    assert pl.decode_consequences(["nonsense"]) == []


def test_sorting_round_trip():
    sorting = [{"id": "POS", "desc": False}, {"id": "CADD", "desc": True}]
    raw = pl.encode_sorting(sorting)
    assert raw == "POS,-CADD" and pl.decode_sorting(raw) == sorting
    assert pl.encode_sorting([]) is None and pl.decode_sorting(None) == []


@pytest.mark.parametrize(
    "target, expected",
    [
        ("/cohort/DEMO?tab=roh", "/cohort/DEMO?tab=roh"),
        ("//evil.example/", None),
        ("/\\evil.example", None),
        ("https://evil.example/", None),
        ("/login?next=/x", None),
        (None, None),
    ],
)
def test_safe_next(target, expected):
    assert pl.safe_next(target) == expected


def test_wombat_filters_round_trip():
    from genetics_viz.pages.cohort.components import wombat_filters as wf

    filters = wf.WombatFilters()
    filters.exclude_mnv["value"] = False
    filters.genesets["value"] = ["hc_ndd"]
    filters.impacts["value"] = [
        t for t, (imp, _) in VEP_CONSEQUENCES.items() if imp == "HIGH"
    ]
    filters.validations["value"] = ["present", "TODO"]
    link = pl.PageLink("/p", {})
    wf.filters_to_link(filters, link)
    assert link.params == {
        "wmnv": "0",
        "wgs": "hc_ndd",
        "wimp": "HIGH",
        "wval": "present,TODO",
    }
    back = wf.filters_from_link(pl.PageLink("/p", link.params))
    assert back == filters
    wf.filters_to_link(wf.WombatFilters(), link)
    assert link.params == {}


def test_search_form_round_trip():
    from genetics_viz.utils import search_rows as sr

    wombat = [{"parent_dir": "/srv/x", "display_name": "strict"}]
    wombat.insert(0, {"parent_dir": "/srv/x", "display_name": "default"})
    wcx = [{"parent_dir": "/srv/y", "display_name": "DEMO_aberrations"}]
    comp_w = sr.new_wombat_component(sr.source_key(wombat[1]))
    comp_w["locus"]["value"] = "SHANK3"
    comp_w["exclude_lcr"]["value"] = False
    comp_x = sr.new_wcx_component(sr.source_key(wcx[0]))
    comp_x["ratio_min"]["value"] = 0.3
    ind = sr.IndividualFilters()
    ind.exclude_samples["value"] = ["S2"]

    raw = sr.form_to_link([comp_w, comp_x], ind, wombat + wcx)
    assert "/srv" not in raw  # no server paths in a link
    components, back = sr.form_from_link(raw, wombat, wcx)
    assert components[0]["source_key"] == sr.source_key(wombat[1])
    assert components[0]["locus"]["value"] == "SHANK3"
    assert components[0]["exclude_lcr"]["value"] is False
    assert components[1]["ratio_min"]["value"] == 0.3
    assert back.exclude_samples["value"] == ["S2"]
    assert sr.form_from_link("not json", wombat, wcx) is None
    assert sr.form_from_link(None, wombat, wcx) is None


def test_cohort_filters_round_trip():
    from genetics_viz.pages.cohort import cohort

    link = pl.PageLink("/cohort/DEMO", {})
    filters = {"FID": "1002", "Sex": ["Female"], "Father": True, "Mother": False}
    cohort._filters_to_link(filters, link)
    assert link.params == {"fid": "1002", "sex": "Female", "father": "1"}
    assert cohort._filters_from_link(link) == {
        "FID": "1002",
        "Sex": ["Female"],
        "Father": True,
    }
