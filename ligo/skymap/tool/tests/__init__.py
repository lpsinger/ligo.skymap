import pkg_resources

dist = 'ligo.skymap'
group = 'console_scripts'

__all__ = ('entry_points', 'run_entry_point')

entry_points = list(pkg_resources.get_entry_map(dist, group).keys())


def run_entry_point(name, *args):
    main = pkg_resources.load_entry_point(dist, group, name)
    try:
        main(args)
    except SystemExit as e:
        assert e.code == 0
