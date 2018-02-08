"""
Pytest plugin to discover doctests from Numpy ufuncs
"""


from _pytest.doctest import (_is_doctest, _is_setup_py, _get_checker,
                             get_optionflags, DoctestTextfile, DoctestItem)
import pytest


### Copied from pytest
def pytest_collect_file(path, parent):
    config = parent.config
    if path.ext == ".py":
        if config.option.doctestmodules and not _is_setup_py(config, path, parent):
            return DoctestModule(path, parent)
    elif _is_doctest(config, path, parent):
        return DoctestTextfile(path, parent)
### End copied from pytest


class DoctestModule(pytest.Module):
    def collect(self):
        import numpy as np
        ### Copied from pytest
        import doctest
        if self.fspath.basename == "conftest.py":
            module = self.config.pluginmanager._importconftest(self.fspath)
        else:
            try:
                module = self.fspath.pyimport()
            except ImportError:
                if self.config.getvalue('doctest_ignore_import_errors'):
                    pytest.skip('unable to import module %r' % self.fspath)
                else:
                    raise
        # uses internal doctest module parsing mechanism
        finder = doctest.DocTestFinder()
        optionflags = get_optionflags(self)
        runner = doctest.DebugRunner(verbose=0, optionflags=optionflags,
                                     checker=_get_checker())

        for test in finder.find(module, module.__name__):
            if test.examples:  # skip empty doctests
                yield DoctestItem(test.name, self, runner, test)
        ### End copied from pytest

        for method in module.__dict__.values():
            if isinstance(method, np.ufunc):
                for test in finder.find(method, module=module):
                    if test.examples:  # skip empty doctests
                        yield DoctestItem(test.name, self, runner, test)
