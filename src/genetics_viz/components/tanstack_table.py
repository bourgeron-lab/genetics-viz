"""TanStack Table component for NiceGUI.

Provides a reusable DataTable class that bridges NiceGUI (Python)
with a TanStack Table (vanilla JS) rendered in the browser.
"""

import json
import uuid
from pathlib import Path
from typing import Any, Callable

from nicegui import ui

TANSTACK_CDN = "https://cdn.jsdelivr.net/npm/@tanstack/table-core/+esm"

_STATIC_DIR = Path(__file__).parent.parent / "static"
_JS_PATH = _STATIC_DIR / "js" / "data_table.js"
_CSS_PATH = _STATIC_DIR / "css" / "data_table.css"

# Track whether scripts have been injected for the current page builder
_injected_pages: set[int] = set()


def _inject_once() -> None:
    """Inject TanStack CDN + DataTable JS + CSS once per page context."""
    page_id = id(ui.context.client)
    if page_id in _injected_pages:
        return
    _injected_pages.add(page_id)

    css_source = _CSS_PATH.read_text()
    js_source = _JS_PATH.read_text()

    ui.add_head_html(f"<style>{css_source}</style>")
    # Use run_javascript instead of a <script> tag so that injection works
    # even when called from timer callbacks (lazy-loaded tabs).
    ui.run_javascript(f"""
        (function() {{
            if (window.TanStackTable && window.DataTable) return;
            import('{TANSTACK_CDN}').then(function(mod) {{
                window.TanStackTable = {{
                    createTable: mod.createTable,
                    getCoreRowModel: mod.getCoreRowModel,
                }};
                {js_source}
            }}).catch(function(err) {{
                console.error('Failed to load TanStack Table:', err);
            }});
        }})();
    """)


class DataTable:
    """A TanStack Table rendered in the browser, controlled from Python.

    Args:
        columns: List of column definition dicts. Each must have at minimum
            ``id`` and ``header``. Optional keys: ``cellType``, ``sortable``,
            and cell-type-specific metadata (see data_table.js).
        rows: List of row dicts.
        row_key: Field name used as unique row identifier.
        pagination: Dict with ``rowsPerPage`` (default 50).
        dense: Use dense (compact) styling.
        selection: Row selection mode - ``"single"``, ``"multi"``, or ``None``.
        sticky_first_column: Make the first column sticky on horizontal scroll.
        on_row_action: Callback for action button clicks.
            Receives dict ``{"action": str, "row": dict}``.
        on_selection: Callback for row selection changes.
            Receives dict ``{"selected": list, "row": dict}``.
        on_sort: Callback for sort changes.
            Receives dict ``{"sorting": [{"id": str, "desc": bool}]}``.
            When provided, the caller is responsible for re-sorting data
            and calling ``update_data()`` with the sorted result.
        instance_id: Optional explicit instance ID. Auto-generated if omitted.
    """

    def __init__(
        self,
        columns: list[dict[str, Any]],
        rows: list[dict[str, Any]],
        *,
        row_key: str = "id",
        pagination: dict[str, Any] | None = None,
        dense: bool = True,
        selection: str | None = None,
        sticky_first_column: bool = False,
        visible_columns: list[str] | None = None,
        on_row_action: Callable | None = None,
        on_selection: Callable | None = None,
        on_sort: Callable | None = None,
        instance_id: str | None = None,
    ):
        self.instance_id = instance_id or f"dt-{uuid.uuid4().hex[:8]}"
        self.columns = columns
        self.rows = rows
        self.row_key = row_key
        self.pagination = pagination or {"rowsPerPage": 50}
        self.dense = dense
        self.selection = selection
        self.sticky_first_column = sticky_first_column
        self.visible_columns = visible_columns
        self._on_row_action = on_row_action
        self._on_selection = on_selection
        self._on_sort = on_sort

        _inject_once()
        self._create_container()
        self._register_events()
        self._init_table()

    def _create_container(self) -> None:
        self.container = ui.html(
            f'<div id="{self.instance_id}" class="data-table-root"></div>',
            sanitize=False,
        )

    def _event_prefix(self) -> str:
        return self.instance_id.replace("-", "_")

    def _register_events(self) -> None:
        prefix = self._event_prefix()

        if self._on_sort:
            ui.on(f"{prefix}_sort-change", lambda e: self._on_sort(e.args))
        else:
            ui.on(f"{prefix}_sort-change", self._default_sort_handler)

        if self._on_row_action:
            ui.on(f"{prefix}_row-action", lambda e: self._on_row_action(e.args))

        if self._on_selection:
            ui.on(f"{prefix}_selection", lambda e: self._on_selection(e.args))

    def _init_table(self) -> None:
        cfg_dict: dict[str, Any] = {
            "containerId": self.instance_id,
            "columns": self.columns,
            "pagination": self.pagination,
            "dense": self.dense,
            "selection": self.selection,
            "stickyFirstColumn": self.sticky_first_column,
            "rowKey": self.row_key,
        }
        if self.visible_columns is not None:
            cfg_dict["visibleColumns"] = self.visible_columns
        config = json.dumps(cfg_dict, default=str)
        data = json.dumps(self.rows, default=str)
        inst_id = self.instance_id

        ui.timer(
            0.1,
            lambda: ui.run_javascript(f"""
                (function tryInit(n) {{
                    if (!window.TanStackTable || !window.DataTable) {{
                        if (n < 50) setTimeout(function() {{ tryInit(n+1) }}, 200);
                        else console.error('DataTable init failed: TanStack not loaded');
                        return;
                    }}
                    var el = document.getElementById('{inst_id}');
                    if (!el) {{
                        if (n < 50) setTimeout(function() {{ tryInit(n+1) }}, 200);
                        return;
                    }}
                    if (window['_dt_{inst_id}']) return;
                    window['_dt_{inst_id}'] = new DataTable({config}, {data});
                }})(1);
            """),
            once=True,
        )

    def update_data(self, rows: list[dict[str, Any]]) -> None:
        """Push new row data to the browser table."""
        self.rows = rows
        data = json.dumps(rows, default=str)
        inst_id = self.instance_id
        ui.run_javascript(f"""
            (function pushData() {{
                if (!window['_dt_{inst_id}']) {{
                    setTimeout(pushData, 100);
                    return;
                }}
                window['_dt_{inst_id}'].setData({data});
            }})();
        """)

    def set_column_visibility(self, visible_ids: list[str]) -> None:
        """Show only the columns whose IDs are in *visible_ids*."""
        vis = json.dumps(visible_ids)
        inst_id = self.instance_id
        ui.run_javascript(f"""
            (function tryVis(n) {{
                if (!window['_dt_{inst_id}']) {{
                    if (n < 50) setTimeout(function() {{ tryVis(n+1) }}, 200);
                    return;
                }}
                window['_dt_{inst_id}'].setColumnVisibility({vis});
            }})(1);
        """)

    def destroy(self) -> None:
        """Clean up the JS instance."""
        inst_id = self.instance_id
        ui.run_javascript(f"""
            if (window['_dt_{inst_id}']) {{
                window['_dt_{inst_id}'].destroy();
                delete window['_dt_{inst_id}'];
            }}
        """)

    # ---- Default sort handler (server-side in Python) ----

    def _default_sort_handler(self, e: Any) -> None:
        data = e.args
        sorting = data.get("sorting", [])
        if sorting:
            col_id = sorting[0]["id"]
            desc = sorting[0].get("desc", False)
            self.rows.sort(
                key=lambda r: (r.get(col_id) is None, r.get(col_id, "")),
                reverse=desc,
            )
        self.update_data(self.rows)
