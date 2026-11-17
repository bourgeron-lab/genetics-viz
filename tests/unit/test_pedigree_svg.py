"""Pedigree drawings through pypedigree: structure only, never pixels."""

import re
import warnings
import xml.etree.ElementTree as ET

from genetics_viz.components.pedigree_tree import (
    tree_html,
    tree_key,
    tree_people,
    worst_severity,
)
from genetics_viz.utils.pedigree_svg import (
    PedigreePerson as P,
    pedigree_svg,
    pedigree_text,
)
from genetics_viz.utils.sample_qc import FlagRow

SVG = "{http://www.w3.org/2000/svg}"

TRIO = [P("S1", sex="1", phenotype="1"), P("S2", sex="2", phenotype="1"),
        P("S3", "S1", "S2", "1", "2")]  # fmt: skip


def _nodes(svg: str) -> dict:
    """``data-id`` -> the node's attributes."""
    root = ET.fromstring(svg.split("?>", 1)[-1])
    return {
        g.get("data-id"): g.attrib
        for g in root.iter(f"{SVG}g")
        if "node" in (g.get("class") or "").split()
    }


def _titles(svg: str) -> dict:
    root = ET.fromstring(svg.split("?>", 1)[-1])
    return {
        g.get("data-id"): g.find(f"{SVG}title").text
        for g in root.iter(f"{SVG}g")
        if g.get("data-id")
    }


def test_trio_is_a_static_svg_under_its_prefix():
    drawing = pedigree_svg("1001", TRIO, id_prefix="pedigree-")
    assert drawing is not None
    svg = drawing.svg
    assert "<script" not in svg
    assert 'class="pypedigree"' in svg and 'width="100%"' in svg
    ids = re.findall(r'\sid="([^"]+)"', svg)
    assert ids and all(i.startswith("pedigree-") for i in ids)
    nodes = _nodes(svg)
    assert set(nodes) == {"S1", "S2", "S3"}
    assert nodes["S3"]["data-parents"] == "S1,S2"
    assert drawing.width > 0 and drawing.height > 0
    assert not drawing.has_unsampled and not drawing.has_notes


def test_one_parent_gets_a_dashed_partner_and_stays_connected():
    # pypedigree's builder left this child unconnected; the parser does not.
    drawing = pedigree_svg("AU-X-1", [P("S5", sex="1"), P("S6", father="S5", sex="2")])
    parents = _nodes(drawing.svg)["S6"]["data-parents"].split(",")
    assert parents[0] == "S5" and len(parents) == 2
    assert drawing.has_unsampled


def test_dummies_and_absent_parents_are_drawn_not_sampled_without_warning():
    people = [
        P("S9", "DUMMY_S9_F", "DUMMY_S9_M", "2", "2"),
        P("DUMMY_S9_M", sex="2", phenotype="2", dummy=True),
    ]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        drawing = pedigree_svg("2002", people)
    titles = _titles(drawing.svg)
    assert set(titles) == {"S9", "DUMMY_S9_M", "DUMMY_S9_F"}
    assert "Not sampled" in titles["DUMMY_S9_F"]  # named, but no row
    assert "Not sampled" in titles["DUMMY_S9_M"]
    # A dummy has no status, whatever its row says.
    assert "Phenotype: Unknown" in titles["DUMMY_S9_M"]
    assert "Not sampled" not in titles["S9"]
    assert drawing.has_unsampled


def test_notes_and_ids_are_escaped():
    drawing = pedigree_svg(
        "F", [P("<b>&x", sex="1", note='XXY<i>"'), P("kid", father="<b>&x")]
    )
    _nodes(drawing.svg)  # parses
    assert "<b>" not in drawing.svg and "&lt;b&gt;&amp;x" in drawing.svg
    assert drawing.has_notes


def test_whitespace_in_a_cell_cannot_shift_the_columns():
    text = pedigree_text("F A", [P("a\tb", note="two words")])
    header, row = text.splitlines()
    assert len(row.split("\t")) == len(header.split("\t"))
    cells = dict(zip(header.split("\t"), row.split("\t")))
    assert (cells["FID"], cells["IID"], cells["a_Note"]) == ("F_A", "a_b", "two_words")


def test_long_ids_widen_the_layout():
    short = pedigree_svg("F", [P("A", sex="1"), P("B", sex="2"), P("C", "A", "B")])
    long = pedigree_svg(
        "F",
        [P("DUMMY_B00IDL0_MF", sex="1", dummy=True),
         P("DUMMY_B00IDL0_MM", sex="2", dummy=True),
         P("C", "DUMMY_B00IDL0_MF", "DUMMY_B00IDL0_MM")],
    )  # fmt: skip

    def span(drawing):
        xs = [float(n["data-x"]) for n in _nodes(drawing.svg).values()]
        return max(xs) - min(xs)

    # The couple sits further apart than a 16-character label is wide (about
    # 2.1 units at 10 pt), so the two labels clear each other.
    assert span(short) < 2.1 < span(long)


def test_nobody_gives_no_drawing():
    assert pedigree_svg("F", []) is None


# --- The member card's Tree view (components/pedigree_tree.py) ----------------


def _flag(sample: str, severity: str) -> FlagRow:
    return FlagRow("family", "F", sample, None, "CODE", severity, "")


def test_classes_reach_the_nodes_and_bad_names_are_dropped():
    drawing = pedigree_svg(
        "F", [P("A", sex="1", classes=("deselected", "qc-high", "bad name", "<x>"))]
    )
    assert _nodes(drawing.svg)["A"]["class"] == "node deselected qc-high"


def test_tree_people_fade_the_unselected_and_outline_open_flags():
    members = [
        {"Sample ID": "S1", "Father": "-", "Mother": "-", "Sex": "1", "Phenotype": "1"},
        {"Sample ID": "S2", "Father": "-", "Mother": "-", "Sex": "2.0", "Phenotype": "1"},
        {"Sample ID": "S3", "Father": "S1", "Mother": "S2", "Sex": "1", "Phenotype": "2"},
        {"Sample ID": "DUMMY_S3_F", "Father": "0", "Mother": "-9", "Sex": "1",
         "Phenotype": "0", "Dummy": True},
    ]  # fmt: skip
    flags = {"S1": [_flag("S1", "info")], "S2": [_flag("S2", "medium"),
             _flag("S2", "high")], "DUMMY_S3_F": [_flag("DUMMY_S3_F", "high")]}  # fmt: skip
    people = {p.iid: p for p in tree_people(members, ["S1", "S3"], flags)}
    assert people["S1"].classes == ()  # info is not outlined
    assert people["S2"].classes == ("deselected", "qc-high")  # the worst
    assert people["S3"].classes == ()
    assert people["DUMMY_S3_F"].classes == () and people["DUMMY_S3_F"].dummy
    assert (people["S2"].sex, people["S3"].father, people["DUMMY_S3_F"].mother) == (
        "2", "S1", "",
    )  # fmt: skip
    drawing = pedigree_svg("F", list(people.values()))
    assert tree_key(list(people.values()), drawing) == [
        "Faded: not selected",
        "Dashed: not sequenced",
        "Coloured outline: open Sample QC issue (high)",
    ]
    assert worst_severity([]) is None


def test_tree_box_is_capped_by_card_natural_size_and_height():
    drawing = pedigree_svg("1001", TRIO)
    html = tree_html(drawing, 55, "300px")
    assert "<?xml" not in html and html.count("<svg") == 1
    natural = f"{drawing.width * 55:.0f}px"
    assert f"width: min(100%, {natural}, calc(300px * " in html
