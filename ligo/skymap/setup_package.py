from astropy_helpers.openmp_helpers import add_openmp_flags_if_available
from astropy_helpers.setup_helpers import pkg_config
from distutils.core import Extension


def _get_build_options():
    kwargs = pkg_config(['chealpix', 'gsl'], [])
    kwargs['include_dirs'].extend(['ligo/skymap/src', 'numpy'])
    kwargs['language'] = 'c'
    return kwargs


def get_extensions():
    kwargs = _get_build_options()
    kwargs['libraries'].append('bayestar')

    extensions = [
        Extension(
            name='ligo.skymap.bayestar._sky_map',
            sources=['ligo/skymap/bayestar/_sky_map.c'],
            **kwargs),
        Extension(
            name='ligo.skymap._distance',
            sources=['ligo/skymap/_distance.c'],
            **kwargs),
        Extension(
            name='ligo.skymap._moc',
            sources=['ligo/skymap/_moc.c'],
            **kwargs)
        ]

    for extension in extensions:
        add_openmp_flags_if_available(extension)

    return extensions


def get_package_info():
    kwargs = _get_build_options()
    library = 'bayestar', dict(sources=[
        'ligo/skymap/src/bayestar_distance.c',
        'ligo/skymap/src/bayestar_moc.c',
        'ligo/skymap/src/bayestar_sky_map.c',
        'ligo/skymap/src/cubic_interp.c'
    ], **kwargs)
    return libraries
