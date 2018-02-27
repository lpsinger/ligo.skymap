import pkg_resources
import pytest

dist = 'ligo.skymap'
group = 'console_scripts'

entry_points = [
    pkg_resources.load_entry_point(dist, group, name)
    for name in pkg_resources.get_entry_map(dist, group).keys()]


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
