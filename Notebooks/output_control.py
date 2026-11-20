"""Route optional notebook products to temporary storage until the run ends.

Authoritative master catalogs and expensive external-query caches are deliberately
written by their original code. Only explicitly routed derived products use this
policy; disabling exports must never delete pre-existing scientific data.
"""
from pathlib import Path
import atexit
import tempfile


class NotebookOutputs:
    """Keep file-based intermediate dependencies working without permanent exports."""

    def __init__(self, *, save_analysis=False, save_diagnostics=False):
        self.save_analysis = save_analysis
        self.save_diagnostics = save_diagnostics
        self._temporary = None
        self._paths = {}
        self._shell = None
        atexit.register(self.cleanup)

    def enabled(self, path):
        """Classify by directory component, with analysis as the default category."""
        return self.save_diagnostics if 'diagnostics' in Path(path).parts else self.save_analysis

    def destination(self, path, *, diagnostic=None):
        """Return a persistent or temporary destination, creating parents on demand."""
        path = Path(path).resolve()
        if diagnostic is None:
            enabled = self.enabled(path)
        else:
            enabled = self.save_diagnostics if diagnostic else self.save_analysis
        if enabled:
            destination = path
        else:
            if self._temporary is None:
                self._temporary = tempfile.TemporaryDirectory(prefix='7dt-notebook-')
            # Retain the full relative tree to avoid collisions between fields.
            destination = Path(self._temporary.name) / path.relative_to(path.anchor)
        destination.parent.mkdir(parents=True, exist_ok=True)
        self._paths[path] = destination
        return destination

    def resolve(self, path):
        """Prefer this run's intermediate over a potentially stale saved product."""
        path = Path(path).resolve()
        return self._paths.get(path, path)

    def exists(self, path):
        return self.resolve(path).exists()

    def read_table(self, path, *args, **kwargs):
        """Read both existing inputs and temporary products through Astropy."""
        from astropy.table import Table
        return Table.read(self.resolve(path), *args, **kwargs)

    def write_table(self, table, path, *args, **kwargs):
        """Serialize intermediate tables temporarily when their export is disabled."""
        return table.write(self.destination(path), *args, **kwargs)

    def write_text(self, path, *args, **kwargs):
        return self.destination(path).write_text(*args, **kwargs)

    def install_error_cleanup(self):
        """Clean up also when Jupyter stops a run because a cell raises an error."""
        try:
            from IPython import get_ipython
            self._shell = get_ipython()
            if self._shell is not None:
                self._shell.events.register('post_run_cell', self._after_cell)
        except ImportError:
            pass  # Plain Python execution still has explicit and atexit cleanup.

    def _after_cell(self, result):
        if result.error_before_exec is not None or result.error_in_exec is not None:
            self.cleanup()

    def cleanup(self):
        """Remove only this run's temporary files; safe to call repeatedly.

        Called by the final notebook cell and again at kernel exit as a fallback.
        After manually stopping between cells, call OUTPUTS.cleanup() or restart
        the kernel; cell errors are handled by install_error_cleanup(). A forcibly
        killed process cannot run Python cleanup handlers.
        """
        if self._shell is not None:
            self._shell.events.unregister('post_run_cell', self._after_cell)
            self._shell = None
        if self._temporary is not None:
            self._temporary.cleanup()
            self._temporary = None
        self._paths.clear()
