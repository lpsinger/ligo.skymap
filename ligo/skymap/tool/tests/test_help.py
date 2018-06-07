import pytest

from ... import version
from . import entry_points, run_entry_point


@pytest.mark.parametrize('entry_point', entry_points)
def test_help(entry_point):
    """Check that --help works."""
    run_entry_point(entry_point, '--help')


@pytest.mark.parametrize('entry_point', entry_points)
def test_version(entry_point, capsys):
    """Check that --version works."""
    run_entry_point(entry_point, '--version')
    out, err = capsys.readouterr()
    assert out.strip().endswith(version.version)
    assert err == ''
