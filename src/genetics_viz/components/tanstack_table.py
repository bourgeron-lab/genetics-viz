"""TanStack Table component for NiceGUI.

Provides a reusable DataTable class that bridges NiceGUI (Python)
with a TanStack Table (vanilla JS) rendered in the browser.
"""

import json
import uuid
from pathlib import Path
from typing import Any, Callable

from nicegui import ui

# Pin the major version: the v9 API dropped `getCoreRowModel` in favour of
# feature objects, so an unpinned URL breaks as soon as jsDelivr's `latest`
# moves on.
TANSTACK_CDN = "https://cdn.jsdelivr.net/npm/@tanstack/table-core@8.21.3/+esm"

_STATIC_DIR = Path(__file__).parent.parent / "static"
_JS_PATH = _STATIC_DIR / "js" / "data_table.js"
_MOUNT_JS_PATH = _STATIC_DIR / "js" / "data_table_mount.js"
_CSS_PATH = _STATIC_DIR / "css" / "data_table.css"


def _inject_once() -> None:
    """Inject TanStack CDN + DataTable JS + CSS once per page context."""
    client = ui.context.client
    if getattr(client, "_dt_scripts_injected", False):
        return
    client._dt_scripts_injected = True

    css_source = _CSS_PATH.read_text()
    js_source = _JS_PATH.read_text()

    ui.add_head_html(f"<style>{css_source}</style>")
    # Use run_javascript instead of a <script> tag so that injection works
    # even when called from timer callbacks (lazy-loaded tabs). The mount
    # registry runs first, synchronously: every table registers with it.
    ui.run_javascript(f"""
        {_MOUNT_JS_PATH.read_text()}
        (function() {{
            if (window.TanStackTable && window.DataTable) return;
            import('{TANSTACK_CDN}').then(function(mod) {{
                window.TanStackTable = {{
                    createTable: mod.createTable,
                    getCoreRowModel: mod.getCoreRowModel,
                }};
                {js_source}
                window._dtMount.flush();
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
        on_filter: Callback for filter changes.
            Receives dict ``{"filters": {col_id: value, ...}}``.
            The caller should filter the full dataset and call
            ``update_data()`` with the result.
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
        on_filter: Callable | None = None,
        on_page_change: Callable | None = None,
        initial_sorting: list[dict[str, Any]] | None = None,
        initial_page: int = 0,
        state_holder: dict[str, Any] | None = None,
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
        self._on_filter = on_filter
        self._on_page_change = on_page_change
        self.initial_sorting = initial_sorting or []
        self.initial_page = initial_page
        self._state_holder = state_holder

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

        if self._on_filter:
            ui.on(f"{prefix}_filterChange", lambda e: self._on_filter(e.args))

        if self._on_page_change or self._state_holder is not None:
            ui.on(f"{prefix}_page-change", self._handle_page_change)

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
        if self.initial_sorting:
            cfg_dict["initialSorting"] = self.initial_sorting
        if self.initial_page:
            cfg_dict["initialPage"] = self.initial_page
        config = json.dumps(cfg_dict, default=str)
        data = json.dumps(self.rows, default=str)
        inst_id = self.instance_id

        # The table is created when its container enters the document, however
        # late (data_table_mount.js): a table built inside a tab panel that is
        # not shown has no container until the panel is opened.
        ui.run_javascript(f"window._dtMount.add('{inst_id}', {config}, {data});")

    def update_data(self, rows: list[dict[str, Any]]) -> None:
        """Push new row data to the browser table."""
        self.rows = rows
        data = json.dumps(rows, default=str)
        inst_id = self.instance_id
        ui.run_javascript(
            f"window._dtMount.whenReady('{inst_id}', "
            f"function(t) {{ t.setData({data}); }});"
        )

    def set_column_visibility(self, visible_ids: list[str]) -> None:
        """Show only the columns whose IDs are in *visible_ids*."""
        vis = json.dumps(visible_ids)
        inst_id = self.instance_id
        ui.run_javascript(
            f"window._dtMount.whenReady('{inst_id}', "
            f"function(t) {{ t.setColumnVisibility({vis}); }});"
        )

    def destroy(self) -> None:
        """Clean up the JS instance."""
        inst_id = self.instance_id
        ui.run_javascript(f"window._dtMount.remove('{inst_id}');")

    # ---- Page change handler ----

    def _handle_page_change(self, e: Any) -> None:
        page = e.args.get("page", 0)
        if self._state_holder is not None:
            self._state_holder["page"] = page
        if self._on_page_change:
            self._on_page_change(e.args)

    # ---- Default sort handler (server-side in Python) ----

    def _default_sort_handler(self, e: Any) -> None:
        data = e.args
        sorting = data.get("sorting", [])
        if sorting:
            col_id = sorting[0]["id"]
            desc = sorting[0].get("desc", False)

            # Look up column def for sort metadata
            col_def = next((c for c in self.columns if c.get("id") == col_id), {})
            sort_field = col_def.get("sortField", col_id)
            sort_type = col_def.get("sorting", "")

            if sort_type == "genomic":
                from genetics_viz.utils.column_names import genomic_sort_key

                self.rows.sort(
                    key=lambda r: (
                        r.get(sort_field) is None,
                        genomic_sort_key(r.get(sort_field, "")),
                    ),
                    reverse=desc,
                )
            elif sort_type == "numerical":

                def _num_key(r: dict) -> tuple:
                    v = r.get(sort_field)
                    if v is None:
                        return (True, 0.0)
                    try:
                        return (False, float(v))
                    except (ValueError, TypeError):
                        return (True, 0.0)

                self.rows.sort(key=_num_key, reverse=desc)
            else:
                self.rows.sort(
                    key=lambda r: (
                        r.get(sort_field) is None,
                        r.get(sort_field, ""),
                    ),
                    reverse=desc,
                )
        # Save sorting state for persistence across refreshes
        if self._state_holder is not None:
            self._state_holder["sorting"] = sorting
            self._state_holder["page"] = 0  # Sort resets to page 0
        self.update_data(self.rows)
