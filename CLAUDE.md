# CLAUDE.md

Project-level instructions for Claude Code.

## Project Overview

NiceGUI-based web application for genetics cohort data visualization. Uses TanStack tables for dataframes, Polars for data processing, and ECharts for plots. Managed with `uv`.

## Commands

- **Run app**: `uv run genetics-viz /path/to/data`
- **Run tests**: `uv run pytest`
- **Lint**: `uv run ruff check .`
- **Format**: `uv run ruff format .`
- **Compile check**: `uv run python -m py_compile src/genetics_viz/path/to/file.py`
- **Always use `uv run`** to execute Python commands, never bare `python`.

## Project Structure

```
src/genetics_viz/
  cli.py                          # Typer CLI entry point
  models.py                       # Cohort, Family, Sample dataclasses + DataStore
  components/
    tanstack_table.py             # DataTable: TanStack table bridge (JS/CSS injection)
    column_selector.py            # Column visibility dialog with presets
    filters.py                    # Reusable filter menu components
    header.py                     # Shared page header with config reload
    validation_loader.py          # Validation TSV loading and badge logic
    variant_dialog.py             # IGV variant inspection dialog
  config/                         # YAML config files (loaded at module level)
    clinvar_colors.yaml           # ClinVar significance term -> color
    column_names.yaml             # Column display names, groups, sorting, drop flags
    cytobands_hg38.tsv            # Cytoband data for ideogram rendering
    score_colors.yaml             # Continuous score color ranges
    vep_consequences.yaml         # VEP consequence terms, impacts, colors
    view_presets.yaml             # Column visibility presets
  pages/
    search.py                     # Cohort-wide variant search with individual filters
    cohort/
      __init__.py                 # Cohort page routing
      cohort.py                   # Cohort overview page
      family.py                   # Family detail page with tabs
      components/
        wombat_tab.py             # WOMBAT analysis tab (per-family)
        dnm_tab.py                # DNM analysis tab (per-family)
        svs_tab.py                # SV analysis tab (per-family)
        stats_panel.py            # Carrier stats box/bar plots
    validation/
      file.py                     # Per-file validation page
      all.py                      # Aggregated validation page
  utils/
    clinvar.py                    # ClinVar color/display utilities (from YAML)
    cytobands.py                  # Cytoband, chromosome, ideogram constants
    vep.py                        # VEP consequence utilities (from YAML)
    view_presets.py               # View preset loading with reload support
    column_names.py               # Column name/group/sorting/schema utilities
    data.py                       # DataStore singleton
    gene_scoring.py               # Gene scoring and color coding
    score_colors.py               # Continuous score color ranges
  static/
    css/data_table.css            # TanStack table styles
    js/data_table.js              # TanStack table JS bridge
```

## Key Patterns

### State Management
- Filter state uses mutable dicts for pass-by-reference in closures: `{"value": [...]}` for lists, `{"value": False}` for booleans.
- `@ui.refreshable` functions support forward references (callbacks defined after the UI that calls `.refresh()`).

### Closure Factory
Use `make_handler(param)` returning inner `handler(e)` for loop variable capture:
```python
for name in names:
    def make_handler(n):
        def handler(e):
            # n is captured correctly
        return handler
    ui.checkbox(name, on_change=make_handler(name))
```

### Button Visual State
Toggle button appearance with Quasar props:
```python
button.props(remove="outline", add="unelevated color=green")  # active
button.props(remove="unelevated color=green", add="outline")  # inactive
```

### Config Utilities
Config files (YAML) are loaded once at module level. Each utility module exposes a `reload_*()` function called from `header.py` to hot-reload configs without restarting the app.

### Pedigree Missing Values
The sentinel set `{"", "0", "-9"}` represents unknown/missing in pedigree fields (parent IDs, sex, phenotype). Defined as `_PED_MISSING` in `search.py` and handled in `models.py` via `treat_missing_as_null`.

### TanStack Table
- `DataTable` in `tanstack_table.py` bridges Python to vanilla JS TanStack Table.
- Per-client JS/CSS injection via `client._dt_scripts_injected` attribute (avoids global set memory leak).
- MutationObserver with 30-second timeout for container detection.

### Async
- Use `asyncio.to_thread()` to offload blocking I/O (e.g., `pl.read_csv()`) in NiceGUI page handlers.
- NiceGUI runs a single-threaded event loop; blocking calls freeze the UI for all users.

## Conventions

- Utility code shared across pages lives in `utils/`, not duplicated in page components.
- YAML configs live in `config/`, loaded by corresponding `utils/*.py` modules.
- All new ClinVar terms must be added to `config/clinvar_colors.yaml`.
- All new VEP consequences must be added to `config/vep_consequences.yaml`.
- Column display/group/sort/drop config goes in `config/column_names.yaml`.
