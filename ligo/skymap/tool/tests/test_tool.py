from os.path import dirname
from importlib import import_module
from pkgutil import iter_modules

import pytest

from ... import tool

entry_points = [getattr(import_module('.' + m.name, tool.__name__), 'main')
                for m in iter_modules([dirname(tool.__file__)]) if not m.ispkg]


@pytest.mark.parametrize('entry_point', entry_points)
def test_help(entry_point):
    """Check that --help works."""
    with pytest.raises(SystemExit) as exc_info:
        entry_point(['--help'])
    assert exc_info.value.code == 0


@pytest.mark.parametrize('entry_point', entry_points)
def test_version(entry_point, capsys):
    """Check that --version works."""
    with pytest.raises(SystemExit) as exc_info:
        entry_point(['--version'])
    assert exc_info.value.code == 0

    out, err = capsys.readouterr()

    from ... import version
    name, _, _ = version.__name__.rpartition('.')
    assert out == name + ' ' + version.version + '\n'
    assert err == ''
