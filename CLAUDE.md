# AG Grid Migration Notes

## Branch: test_aggrid

This branch migrates the search page (`src/genetics_viz/pages/search.py`) from NiceGUI's `ui.table()` (Quasar QTable) to `ui.aggrid()` (AG Grid Community).

## Key Architecture

### AG Grid in NiceGUI
- `ui.aggrid(options_dict)` — options dict is serialized to JSON
- **`:` prefix** on property keys (e.g., `":cellRenderer"`) tells NiceGUI to treat the value as a JavaScript expression instead of a string literal
- Cell renderers use inline JS functions: `":cellRenderer": "(params) => { ... }"`
- `agHtmlCellRenderer` and custom component registration do NOT work in NiceGUI's AG Grid — always use `:cellRenderer` with inline JS

### HTML Rendering in AG Grid Cells
- Badge columns (Consequence, ClinVar, Gene Symbol, Gene ID) have pre-rendered HTML stored in `_*_html` fields on each row
- The cell renderer creates a div and sets `innerHTML`: `"(params) => { const div = document.createElement('div'); div.innerHTML = params.value || ''; return div; }"`
- Score badge columns also use this pattern (e.g., `_REVEL_html`, `_MPC_html`)
- FID column renders clickable links to `/cohort/{cohort_name}/family/{FID}`

### Column Configuration (column_names.yaml)
- `src/genetics_viz/config/column_names.yaml` — defines display `name` and `group` for columns
- `src/genetics_viz/utils/column_config.py` — loads the YAML once at import time, exposes:
  - `get_column_display_name(col)` — returns display name (falls back to column name)
  - `get_column_group(col)` — returns group name or None
- Both `search.py` and `wombat_tab.py` delegate display labels to this config
- Column groups in AG Grid are created dynamically: columns sharing the same `group` value are wrapped in `{"headerName": group_name, "children": [...]}`

### CSS for AG Grid
- `.search-grid` class on the aggrid element
- `.search-grid .ag-cell { display: flex !important; align-items: center !important; }` for vertical centering
- `rowHeight: 32` for compact rows
- `autoSizeStrategy: {"type": "fitCellContents"}` for column auto-width
- Badge styling: `padding: 1px 6px; border-radius: 3px; line-height: 1.2; font-size: 0.75em`

### Wombat Tab Stats Dialog
- Stats button next to "+ column" in wombat_tab.py
- Shows chromosome bar chart, consequence pie, validation pie (ECharts)
- Ideogram view toggle with SVG cytobands from `config/cytobands_hg38.tsv`
- All stats based on unique variants (deduplicated by #CHROM/POS/REF/ALT)

## Files Changed (vs main)

- `src/genetics_viz/pages/search.py` — AG Grid migration, column grouping, autofit, HTML cell renderers
- `src/genetics_viz/pages/cohort/components/wombat_tab.py` — Stats dialog, uses shared column config
- `src/genetics_viz/config/column_names.yaml` — NEW: column display names and groups
- `src/genetics_viz/utils/column_config.py` — NEW: shared column config loader

## Known Patterns / Gotchas
- NiceGUI AG Grid does NOT support `components` registration — no way to register custom cell renderer classes
- `:cellRenderer` must be an inline function string, not a reference to a named component
- `autoHeight` on columns conflicts with `rowHeight` — don't use both
- Column groups in AG Grid require `children` array — the group itself has no `field`
- When inserting groups into the flat column list, calculate position based on non-grouped columns before the group's first member
