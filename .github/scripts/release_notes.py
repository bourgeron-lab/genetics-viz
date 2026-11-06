"""Print one version's CHANGELOG section, the body of its GitHub release.

    python .github/scripts/release_notes.py 0.31.0 > notes.md

Exits non-zero when CHANGELOG.md has no entry for the version: the publish
workflow checks that before anything reaches PyPI.
"""

import re
import sys
from pathlib import Path

CHANGELOG = Path(__file__).resolve().parents[2] / "CHANGELOG.md"


def section(version: str, text: str) -> str:
    match = re.search(
        rf"^## \[{re.escape(version)}\][^\n]*\n(.*?)(?=^## \[|\Z)",
        text,
        re.S | re.M,
    )
    if match is None or not match.group(1).strip():
        raise SystemExit(f"CHANGELOG.md has no entry for {version}")
    return match.group(1).strip()


if __name__ == "__main__":
    print(section(sys.argv[1].removeprefix("v"), CHANGELOG.read_text()))
