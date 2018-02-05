from astropy_helpers.openmp_helpers import add_openmp_flags_if_available
from astropy_helpers.setup_helpers import pkg_config
from distutils.core import Extension

def get_extensions():
    include_dirs = ['numpy']
    libraries = pkg_config(['chealpix', 'gsl'], [])

    extensions = [
        Extension(
            name='ligo.skymap.bayestar._sky_map',
            sources=['ligo/skymap/bayestar/_sky_map.c'],
            language='c',
            libraries=libraries,
            include_dirs=include_dirs),
        Extension(
            name='ligo.skymap._distance',
            sources=['ligo/skymap/_distance.c'],
            language='c',
            libraries=libraries,
            include_dirs=include_dirs),
        Extension(
            name='ligo.skymap._moc',
            sources=['ligo/skymap/_moc.c'],
            language='c',
            libraries=libraries,
            include_dirs=include_dirs)
        ]

    for extension in extensions:
        add_openmp_flags_if_available(extension)

    return extensions
