/**
 * DataTable - TanStack Table wrapper for NiceGUI.
 *
 * Provides server-side sorting with client-side pagination,
 * custom cell rendering driven by column metadata, and
 * a JS<->Python bridge via emitEvent / ui.run_javascript.
 */
class DataTable {
    constructor(config, data) {
        this.containerId = config.containerId;
        this.container = document.getElementById(config.containerId);
        this.config = config;
        this.data = data || [];
        this.sorting = [];
        this.currentPage = 0;
        this.rowsPerPage = (config.pagination && config.pagination.rowsPerPage) || 50;
        this.dense = config.dense !== false;
        this.selection = config.selection || null;
        this.selectedRowKeys = [];
        this.stickyFirstColumn = config.stickyFirstColumn || false;
        this.rowKey = config.rowKey || 'id';

        this._allColumnDefs = config.columns;
        this._visibleColumnIds = config.visibleColumns
            ? config.visibleColumns
            : config.columns.map(function(c) { return c.id; });

        this._initTable();
        this.render();
    }

    _initTable() {
        var TanStack = window.TanStackTable;
        var self = this;

        this.table = TanStack.createTable({
            data: this._getPageData(),
            columns: this._buildTanStackColumns(),
            state: { sorting: this.sorting },
            onStateChange: function(updater) {
                self.table.setOptions(function(prev) {
                    return Object.assign({}, prev, {
                        state: updater(self.table.getState()),
                    });
                });
            },
            onSortingChange: function(updater) {
                self.sorting = typeof updater === 'function'
                    ? updater(self.sorting) : updater;
                self.currentPage = 0;
                self._updateTableState();
                self.render();
                self._emitEvent('sortChange', {
                    sorting: self.sorting.map(function(s) {
                        return { id: s.id, desc: s.desc };
                    })
                });
            },
            getCoreRowModel: TanStack.getCoreRowModel(),
            enableSorting: true,
            enableSortingRemoval: true,
        });

        // CRITICAL: merge initialState into active state
        this.table.setOptions(function(prev) {
            return Object.assign({}, prev, {
                state: Object.assign({}, self.table.initialState, prev.state),
            });
        });
    }

    // ---- Column building ----

    _getVisibleColumnDefs() {
        var ids = this._visibleColumnIds;
        return this._allColumnDefs.filter(function(c) {
            return ids.indexOf(c.id) !== -1;
        });
    }

    _buildTanStackColumns() {
        var self = this;
        return this._getVisibleColumnDefs().map(function(def) {
            return {
                accessorKey: def.id,
                header: def.header || def.id,
                enableSorting: def.sortable !== false,
                meta: def,
            };
        });
    }

    // ---- Pagination ----

    _getPageData() {
        var start = this.currentPage * this.rowsPerPage;
        return this.data.slice(start, start + this.rowsPerPage);
    }

    _getTotalPages() {
        return Math.max(1, Math.ceil(this.data.length / this.rowsPerPage));
    }

    // ---- State updates ----

    _updateTableState() {
        var self = this;
        this.table.setOptions(function(prev) {
            return Object.assign({}, prev, {
                data: self._getPageData(),
                state: Object.assign({}, self.table.getState(), {
                    sorting: self.sorting,
                }),
            });
        });
    }

    _updateTableData() {
        var self = this;
        this.table.setOptions(function(prev) {
            return Object.assign({}, prev, {
                data: self._getPageData(),
                state: self.table.getState(),
            });
        });
    }

    // ---- Public API ----

    setData(data) {
        this.data = data;
        this.currentPage = 0;
        this._updateTableData();
        this.render();
    }

    setColumns(columns) {
        this._allColumnDefs = columns;
        this._visibleColumnIds = columns.map(function(c) { return c.id; });
        this._rebuildTable();
    }

    setColumnVisibility(visibleIds) {
        this._visibleColumnIds = visibleIds;
        this._rebuildTable();
    }

    _rebuildTable() {
        var self = this;
        this.table.setOptions(function(prev) {
            return Object.assign({}, prev, {
                columns: self._buildTanStackColumns(),
                data: self._getPageData(),
                state: self.table.getState(),
            });
        });
        this.render();
    }

    destroy() {
        if (this.container) this.container.innerHTML = '';
        this.table = null;
        this.data = null;
    }

    // ---- Rendering ----

    render() {
        if (!this.container) return;
        this.container.innerHTML = '';

        var wrapper = document.createElement('div');
        wrapper.className = 'dt-wrapper';

        var scrollContainer = document.createElement('div');
        scrollContainer.className = 'dt-scroll-container';

        var tableContainer = document.createElement('div');
        tableContainer.className = 'dt-table-container';

        var tableClasses = 'dt-table';
        if (this.dense) tableClasses += ' dt-dense';
        if (this.stickyFirstColumn) tableClasses += ' dt-sticky-col';

        var table = document.createElement('table');
        table.className = tableClasses;

        table.appendChild(this._renderHeader());
        table.appendChild(this._renderBody());
        tableContainer.appendChild(table);
        scrollContainer.appendChild(tableContainer);
        wrapper.appendChild(scrollContainer);

        if (this.data.length > this.rowsPerPage) {
            wrapper.appendChild(this._renderPagination());
        } else {
            wrapper.appendChild(this._renderRowCount());
        }

        this.container.appendChild(wrapper);
    }

    _renderHeader() {
        var thead = document.createElement('thead');
        var headerGroups = this.table.getHeaderGroups();
        if (headerGroups.length === 0) return thead;

        var headers = headerGroups[0].headers;

        // Collect group info from column metadata
        var colInfos = [];
        var hasGroups = false;
        for (var i = 0; i < headers.length; i++) {
            var meta = headers[i].column.columnDef.meta || {};
            var group = meta.group || '';
            if (group) hasGroups = true;
            colInfos.push({ header: headers[i], group: group });
        }

        if (!hasGroups) {
            // Single-row header (no groups)
            var tr = document.createElement('tr');
            for (var i = 0; i < colInfos.length; i++) {
                tr.appendChild(this._renderHeaderCell(colInfos[i].header));
            }
            thead.appendChild(tr);
        } else {
            // Two-row header: group row + column row
            var groupRow = document.createElement('tr');
            groupRow.className = 'dt-group-header-row';
            var columnRow = document.createElement('tr');

            var i = 0;
            while (i < colInfos.length) {
                var info = colInfos[i];
                if (info.group === '') {
                    // Ungrouped: rowspan=2 in group row, skip column row
                    var th = this._renderHeaderCell(info.header);
                    th.rowSpan = 2;
                    groupRow.appendChild(th);
                    i++;
                } else {
                    // Count consecutive columns with same group
                    var groupName = info.group;
                    var startIdx = i;
                    while (i < colInfos.length && colInfos[i].group === groupName) {
                        i++;
                    }
                    var span = i - startIdx;

                    // Group row: one th with colspan
                    var groupTh = document.createElement('th');
                    groupTh.className = 'dt-group-header';
                    groupTh.colSpan = span;
                    groupTh.textContent = groupName;
                    groupRow.appendChild(groupTh);

                    // Column row: individual th cells
                    for (var j = startIdx; j < startIdx + span; j++) {
                        columnRow.appendChild(this._renderHeaderCell(colInfos[j].header));
                    }
                }
            }

            thead.appendChild(groupRow);
            thead.appendChild(columnRow);
        }

        return thead;
    }

    _renderHeaderCell(header) {
        var th = document.createElement('th');
        var canSort = header.column.getCanSort();

        var label = document.createElement('span');
        label.className = 'dt-header-label';
        label.textContent = header.column.columnDef.header;

        if (canSort) {
            th.className = 'dt-sortable';
            var sortDir = header.column.getIsSorted();
            if (sortDir === 'asc') {
                var arrow = document.createElement('span');
                arrow.className = 'dt-sort-indicator';
                arrow.textContent = ' \u25B2';
                label.appendChild(arrow);
            } else if (sortDir === 'desc') {
                var arrow = document.createElement('span');
                arrow.className = 'dt-sort-indicator';
                arrow.textContent = ' \u25BC';
                label.appendChild(arrow);
            }

            var handler = header.column.getToggleSortingHandler();
            if (handler) {
                (function(h) {
                    th.addEventListener('click', function(e) { h(e); });
                })(handler);
            }
        }

        th.appendChild(label);
        return th;
    }

    _renderBody() {
        var tbody = document.createElement('tbody');
        var rows = this.table.getRowModel().rows;

        if (rows.length === 0) {
            var tr = document.createElement('tr');
            var td = document.createElement('td');
            td.colSpan = this._getVisibleColumnDefs().length || 1;
            td.className = 'dt-empty';
            td.textContent = 'No data';
            tr.appendChild(td);
            tbody.appendChild(tr);
            return tbody;
        }

        for (var ri = 0; ri < rows.length; ri++) {
            var row = rows[ri];
            var tr = document.createElement('tr');
            var isSelected = this.selectedRowKeys.indexOf(
                row.original[this.rowKey]
            ) !== -1;

            tr.className = ri % 2 === 0 ? 'dt-row-even' : 'dt-row-odd';
            if (isSelected) tr.className += ' dt-row-selected';

            if (this.selection) {
                var self = this;
                (function(rowData) {
                    tr.addEventListener('click', function() {
                        self._handleRowSelection(rowData);
                    });
                    tr.style.cursor = 'pointer';
                })(row.original);
            }

            var cells = row.getVisibleCells();
            for (var ci = 0; ci < cells.length; ci++) {
                var td = document.createElement('td');
                td.innerHTML = this._renderCell(cells[ci], row.original);
                tr.appendChild(td);
            }
            tbody.appendChild(tr);
        }
        return tbody;
    }

    // ---- Cell rendering ----

    _renderCell(cell, rowData) {
        var value = cell.getValue();
        var meta = cell.column.columnDef.meta || {};
        var cellType = meta.cellType || 'text';

        switch (cellType) {
            case 'action': return this._renderActionCell(meta, rowData);
            case 'validation': return this._renderValidationCell(value, rowData);
            case 'badge': return this._renderBadgeCell(value, meta, rowData);
            case 'badge_list': return this._renderBadgeListCell(value, meta, rowData);
            case 'gene_badge': return this._renderGeneBadgeCell(value, meta, rowData);
            case 'link': return this._renderLinkCell(value, meta, rowData);
            case 'color_scale': return this._renderColorScaleCell(value, meta, rowData);
            case 'score_badge': return this._renderScoreBadgeCell(value, meta, rowData);
            case 'cnv_call': return this._renderCnvCallCell(value, meta);
            case 'number': return this._renderNumberCell(value);
            default: return this._renderTextCell(value);
        }
    }

    _renderTextCell(value) {
        if (value === null || value === undefined || value === '') return '';
        return this._escapeHtml(String(value));
    }

    _renderNumberCell(value) {
        if (value === null || value === undefined || value === '') return '';
        var n = Number(value);
        if (isNaN(n)) return this._escapeHtml(String(value));
        return Number.isInteger(n) ? String(n) : n.toFixed(4);
    }

    _renderActionCell(meta, rowData) {
        var actionName = meta.actionName || 'row_action';
        var icon = meta.actionIcon || 'visibility';
        var color = meta.actionColor || '#1976d2';
        var tooltip = meta.actionTooltip || '';
        var instId = this.containerId;

        var html = '<button class="dt-action-btn" style="color: ' + color + ';"'
            + ' onclick="DataTable._emitRowAction(\'' + instId + '\', \'' + actionName + '\', this)"'
            + ' data-row=\'' + this._escapeAttr(JSON.stringify(rowData)) + '\'';
        if (tooltip) html += ' data-tooltip="' + this._escapeAttr(tooltip) + '"';
        html += '><span class="material-icons" style="font-size: 20px;">' + icon + '</span></button>';

        // Show grouped count badge if configured
        if (meta.showGroupBadge && rowData.n_grouped && rowData.n_grouped > 1) {
            html += '<span class="dt-badge dt-group-badge">' + rowData.n_grouped + '</span>';
        }

        return html;
    }

    _renderValidationCell(value, rowData) {
        if (!value) return '';
        var v = String(value).toLowerCase();

        if (v === 'present' || v === 'in phase mnv') {
            var html = '<span class="dt-validation-wrap">';
            html += '<span class="material-icons dt-validation-icon" style="color: #4caf50;" data-tooltip="Validated as ' + this._escapeAttr(value) + '">check_circle</span>';
            var inh = rowData.ValidationInheritance || rowData.Inheritance || '';
            if (inh === 'de novo') html += '<strong class="dt-inh-label">dnm</strong>';
            else if (inh === 'homozygous') html += '<strong class="dt-inh-label">hom</strong>';
            if (v === 'in phase mnv') html += '<span class="dt-mnv-label">MNV</span>';
            html += '</span>';
            return html;
        }
        if (v === 'absent') {
            return '<span class="material-icons dt-validation-icon" style="color: #f44336;" data-tooltip="Validated as absent">cancel</span>';
        }
        if (v === 'uncertain' || v === 'different') {
            return '<span class="material-icons dt-validation-icon" style="color: #ff9800;" data-tooltip="Validation uncertain or different">help</span>';
        }
        if (v === 'conflicting') {
            return '<span class="material-icons dt-validation-icon" style="color: #ff8f00;" data-tooltip="Conflicting validations">bolt</span>';
        }
        return this._escapeHtml(String(value));
    }

    _renderBadgeCell(value, meta, rowData) {
        if (value === null || value === undefined || value === '') return '';
        var color = (meta.colorField && rowData[meta.colorField]) || meta.badgeColor || '#9e9e9e';
        var textColor = meta.textColor || 'white';
        return '<span class="dt-badge" style="background-color: ' + color + '; color: ' + textColor + ';">'
            + this._escapeHtml(String(value)) + '</span>';
    }

    _renderBadgeListCell(value, meta, rowData) {
        var field = meta.badgesField || (meta.id + '_badges');
        var badges = rowData[field];
        if (!badges || !badges.length) {
            return value ? this._escapeHtml(String(value)) : '';
        }
        var html = '<div class="dt-badge-wrap">';
        for (var i = 0; i < badges.length; i++) {
            var b = badges[i];
            var textColor = b.color === '#ffffff' ? 'black' : 'white';
            html += '<span class="dt-badge" style="background-color: ' + b.color + '; color: ' + textColor + ';"';
            if (b.tooltip) html += ' data-tooltip="' + this._escapeAttr(b.tooltip) + '"';
            html += '>' + this._escapeHtml(b.label) + '</span>';
        }
        html += '</div>';
        return html;
    }

    _renderGeneBadgeCell(value, meta, rowData) {
        var field = meta.badgesField || (meta.id + '_badges');
        var badges = rowData[field];
        if (!badges || !badges.length) {
            return value ? this._escapeHtml(String(value)) : '-';
        }
        var html = '<div class="dt-badge-wrap">';
        for (var i = 0; i < badges.length; i++) {
            var b = badges[i];
            var textColor = b.color === '#ffffff' ? 'black' : 'white';
            var borderStyle = '';
            if (b.borderColor) borderStyle = 'border: 2px solid ' + b.borderColor + ';';
            html += '<span class="dt-badge dt-gene-badge" style="background-color: ' + b.color + '; color: ' + textColor + '; ' + borderStyle + '"';
            if (b.tooltip) html += ' data-tooltip="' + this._escapeAttr(b.tooltip) + '"';
            html += '>' + this._escapeHtml(b.label) + '</span>';
        }
        html += '</div>';
        return html;
    }

    _renderLinkCell(value, meta, rowData) {
        if (value === null || value === undefined || value === '') return '';
        var href = meta.href || '#';
        // Interpolate {field} placeholders from rowData
        href = href.replace(/\{([^}]+)\}/g, function(match, field) {
            return rowData[field] !== undefined ? encodeURIComponent(rowData[field]) : match;
        });
        return '<a href="' + this._escapeAttr(href) + '" class="dt-link">' + this._escapeHtml(String(value)) + '</a>';
    }

    _renderColorScaleCell(value, meta, rowData) {
        if (value === null || value === undefined || value === '') return '';
        var n = parseFloat(value);
        if (isNaN(n)) return this._escapeHtml(String(value));

        var color = 'inherit';
        var weight = 'normal';
        var thresholds = meta.thresholds || [];

        for (var i = 0; i < thresholds.length; i++) {
            var t = thresholds[i];
            var match = false;
            if (t.op === '<=' && n <= t.value) match = true;
            else if (t.op === '>=' && n >= t.value) match = true;
            else if (t.op === '<' && n < t.value) match = true;
            else if (t.op === '>' && n > t.value) match = true;

            if (match) {
                color = t.color;
                weight = t.weight || 'normal';
                break;
            }
        }

        var display = Number.isInteger(n) ? String(n) : n.toFixed(4);
        return '<span style="color: ' + color + '; font-weight: ' + weight + ';">' + display + '</span>';
    }

    _renderScoreBadgeCell(value, meta, rowData) {
        var badgeField = meta.id + '_badge';
        var badge = rowData[badgeField];
        if (badge) {
            var textColor = badge.color === '#ffffff' ? 'black' : 'white';
            var html = '<span class="dt-badge dt-score-badge" style="background-color: ' + badge.color + '; color: ' + textColor + ';"';
            if (badge.tooltip) html += ' data-tooltip="' + this._escapeAttr(badge.tooltip) + '"';
            html += '>' + this._escapeHtml(String(badge.label)) + '</span>';
            return html;
        }
        if (value === null || value === undefined || value === '') return '';
        return this._escapeHtml(String(value));
    }

    _renderCnvCallCell(value, meta) {
        if (value === null || value === undefined || value === '') return '';
        var callColors = meta.callColors || {};
        var color = callColors[String(value)];
        if (color) {
            return '<span class="dt-badge" style="background-color: ' + color + '; color: white;">'
                + this._escapeHtml(String(value)) + '</span>';
        }
        return '<span class="text-grey-6">' + this._escapeHtml(String(value)) + '</span>';
    }

    // ---- Pagination ----

    _renderPagination() {
        var self = this;
        var totalPages = this._getTotalPages();
        var totalRows = this.data.length;
        var start = this.currentPage * this.rowsPerPage + 1;
        var end = Math.min((this.currentPage + 1) * this.rowsPerPage, totalRows);

        var bar = document.createElement('div');
        bar.className = 'dt-pagination';

        // Row count
        var info = document.createElement('span');
        info.className = 'dt-pagination-info';
        info.textContent = start + '-' + end + ' of ' + totalRows;
        bar.appendChild(info);

        // Rows per page selector
        var rppWrap = document.createElement('span');
        rppWrap.className = 'dt-rpp-wrap';
        var rppLabel = document.createElement('span');
        rppLabel.textContent = 'Rows per page: ';
        rppWrap.appendChild(rppLabel);
        var select = document.createElement('select');
        select.className = 'dt-rpp-select';
        var options = [10, 25, 50, 100];
        for (var i = 0; i < options.length; i++) {
            var opt = document.createElement('option');
            opt.value = options[i];
            opt.textContent = options[i];
            if (options[i] === self.rowsPerPage) opt.selected = true;
            select.appendChild(opt);
        }
        select.addEventListener('change', function() {
            self.rowsPerPage = parseInt(this.value);
            self.currentPage = 0;
            self._updateTableData();
            self.render();
        });
        rppWrap.appendChild(select);
        bar.appendChild(rppWrap);

        // Navigation buttons
        var nav = document.createElement('span');
        nav.className = 'dt-pagination-nav';

        var prevBtn = document.createElement('button');
        prevBtn.className = 'dt-page-btn';
        prevBtn.innerHTML = '<span class="material-icons">chevron_left</span>';
        prevBtn.disabled = self.currentPage === 0;
        prevBtn.addEventListener('click', function() {
            if (self.currentPage > 0) {
                self.currentPage--;
                self._updateTableData();
                self.render();
            }
        });
        nav.appendChild(prevBtn);

        var pageInfo = document.createElement('span');
        pageInfo.className = 'dt-page-info';
        pageInfo.textContent = (self.currentPage + 1) + ' / ' + totalPages;
        nav.appendChild(pageInfo);

        var nextBtn = document.createElement('button');
        nextBtn.className = 'dt-page-btn';
        nextBtn.innerHTML = '<span class="material-icons">chevron_right</span>';
        nextBtn.disabled = self.currentPage >= totalPages - 1;
        nextBtn.addEventListener('click', function() {
            if (self.currentPage < totalPages - 1) {
                self.currentPage++;
                self._updateTableData();
                self.render();
            }
        });
        nav.appendChild(nextBtn);

        bar.appendChild(nav);
        return bar;
    }

    _renderRowCount() {
        var bar = document.createElement('div');
        bar.className = 'dt-pagination';
        var info = document.createElement('span');
        info.className = 'dt-pagination-info';
        info.textContent = this.data.length + ' row' + (this.data.length !== 1 ? 's' : '');
        bar.appendChild(info);
        return bar;
    }

    // ---- Selection ----

    _handleRowSelection(rowData) {
        var key = rowData[this.rowKey];
        if (this.selection === 'single') {
            if (this.selectedRowKeys.length === 1 && this.selectedRowKeys[0] === key) {
                this.selectedRowKeys = [];
            } else {
                this.selectedRowKeys = [key];
            }
        } else if (this.selection === 'multi') {
            var idx = this.selectedRowKeys.indexOf(key);
            if (idx === -1) {
                this.selectedRowKeys.push(key);
            } else {
                this.selectedRowKeys.splice(idx, 1);
            }
        }
        this.render();
        this._emitEvent('selection', {
            selected: this.selectedRowKeys,
            row: rowData,
        });
    }

    // ---- Events ----

    _emitEvent(name, data) {
        if (typeof emitEvent !== 'function') return;
        var prefix = this.containerId.replace(/-/g, '_');
        emitEvent(prefix + '_' + name, data);
    }

    static _emitRowAction(containerId, actionName, btnEl) {
        var rowData = JSON.parse(btnEl.getAttribute('data-row'));
        var prefix = containerId.replace(/-/g, '_');
        if (typeof emitEvent === 'function') {
            emitEvent(prefix + '_rowAction', { action: actionName, row: rowData });
        }
    }

    // ---- Utilities ----

    _escapeHtml(str) {
        var div = document.createElement('div');
        div.appendChild(document.createTextNode(str));
        return div.innerHTML;
    }

    _escapeAttr(str) {
        return String(str)
            .replace(/&/g, '&amp;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#39;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;');
    }
}

window.DataTable = DataTable;
