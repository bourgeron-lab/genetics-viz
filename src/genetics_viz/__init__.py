"""
Genetics-Viz: A web-based visualization tool for genetics cohort data.
"""

from importlib.metadata import PackageNotFoundError, version

# The version lives in pyproject.toml only: a release bumps nothing else
# (DEPLOY.md). A hard-coded copy here stayed at 0.36.0 for three releases.
try:
    __version__ = version("genetics-viz")
except PackageNotFoundError:
    __version__ = "unknown"
