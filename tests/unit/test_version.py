"""The version shown in the header is the one released."""

import tomllib
from pathlib import Path

from genetics_viz import __version__


def test_version_is_pyprojects():
    pyproject = Path(__file__).parents[2] / "pyproject.toml"
    with pyproject.open("rb") as f:
        assert __version__ == tomllib.load(f)["project"]["version"]
