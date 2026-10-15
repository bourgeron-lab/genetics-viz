/**
 * DataTable mount registry.
 *
 * A table is built in Python before its container is necessarily in the
 * document: Quasar mounts only the active tab panel, and a keep-alive panel
 * that has been left sits detached until it is shown again. That can be any
 * time later, so tables wait here, with no timeout, and are created when a
 * single MutationObserver sees their container appear.
 *
 * Runs before the TanStack library has loaded; data_table.js calls flush()
 * once it has.
 */
if (!window._dtMount) {
    window._dtMount = (function() {
        var pending = {};
        var observer = null;

        function flush() {
            if (!window.TanStackTable || !window.DataTable) return;
            Object.keys(pending).forEach(function(id) {
                if (!document.getElementById(id)) return;
                var entry = pending[id];
                delete pending[id];
                // Already initialized (re-render guard)
                if (window['_dt_' + id]) return;
                var inst = new DataTable(entry.config, entry.data);
                window['_dt_' + id] = inst;
                entry.calls.forEach(function(fn) { fn(inst); });
            });
            if (observer && Object.keys(pending).length === 0) {
                observer.disconnect();
                observer = null;
            }
        }

        function add(id, config, data) {
            pending[id] = { config: config, data: data, calls: [] };
            if (!observer) {
                observer = new MutationObserver(flush);
                observer.observe(document.body, { childList: true, subtree: true });
            }
            flush();
        }

        // Run fn on the table now, or once it is created.
        function whenReady(id, fn) {
            var inst = window['_dt_' + id];
            if (inst) fn(inst);
            else if (pending[id]) pending[id].calls.push(fn);
        }

        function remove(id) {
            delete pending[id];
            var inst = window['_dt_' + id];
            if (inst) {
                inst.destroy();
                delete window['_dt_' + id];
            }
        }

        return { add: add, flush: flush, whenReady: whenReady, remove: remove };
    })();
}
