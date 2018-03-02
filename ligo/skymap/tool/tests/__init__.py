import pkg_resources
import pytest

dist = 'ligo.skymap'
group = 'console_scripts'

__all__ = ('entry_points', 'run_entry_point')

entry_points = list(pkg_resources.get_entry_map(dist, group).keys())


def run_entry_point(name, *args):
    main = pkg_resources.load_entry_point(dist, group, name)
    with pytest.raises(SystemExit) as exc_info:
        main(args)
    assert exc_info.value.code == 0
