"""The family report as one self-contained HTML document.

The dialog shows it in an ``<iframe srcdoc>``, and the PDF button prints the
same string from a hidden iframe. Isolating it in an iframe keeps Quasar and
Tailwind styles out, so screen and print look alike. Everything it needs is
inline: the styles, and the plots as SVG (:mod:`report_svg`).
"""

from pathlib import Path

from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup

from .family_report import FamilyReport

_TEMPLATES = Path(__file__).parent.parent / "templates"

_ENV = Environment(
    loader=FileSystemLoader(_TEMPLATES),
    autoescape=select_autoescape(["html"]),
    trim_blocks=True,
    lstrip_blocks=True,
)


def render_report_html(report: FamilyReport) -> str:
    """The whole document, ``<!DOCTYPE html>`` included."""
    template = _ENV.get_template("family_report.html")
    # The SVGs are built from escaped text by report_svg; they are marked safe
    # here rather than with ``|safe`` in the template, so no other field can be.
    return template.render(report=report, svg_of=lambda plot: Markup(plot.svg))
