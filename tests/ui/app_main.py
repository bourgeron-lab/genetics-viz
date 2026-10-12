"""Main file run by NiceGUI's user simulation, once per test.

NiceGUI resets its globals between tests: it removes every route and drops the
modules that defined pages from sys.modules. genetics_viz registers its pages
and /health on import, and its packages stay imported, so a plain import would
register nothing from the second test on. Purging the whole package first makes
each test boot the app from scratch, module-level state included.
"""

import os
import sys
from pathlib import Path

for name in [
    m for m in sys.modules if m == "genetics_viz" or m.startswith("genetics_viz.")
]:
    del sys.modules[name]

from nicegui import ui  # noqa: E402

from genetics_viz.app import _init_from_config  # noqa: E402
from genetics_viz.config_model import get_config  # noqa: E402

_init_from_config(Path(os.environ["GENETICS_VIZ_TEST_CONFIG"]))
ui.run(storage_secret=get_config().storage_secret)
