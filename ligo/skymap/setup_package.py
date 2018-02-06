from astropy_helpers import openmp_helpers, setup_helpers
from distutils.core import Extension


def get_extensions():
    pkg_config_packages = ['gsl']

    sources = [
        'src/cubic_interp.c',
        'src/bayestar_moc.c',
        'src/core.c',
        'src/bayestar_distance.c',
        'src/bayestar_sky_map.c']

    include_dirs = ['numpy']

    if setup_helpers.use_system_library('chealix'):
        pkg_config_packages.append('chealpix')
    else:
        include_dirs.append('cextern/chealpix')
        sources.append('cextern/chealpix/chealpix.c')

    kwargs = setup_helpers.pkg_config(pkg_config_packages, [])
    kwargs['include_dirs'].extend(include_dirs)

    extension = Extension(
        name='ligo.skymap.core', language='c', sources=sources, **kwargs)

    openmp_helpers.add_openmp_flags_if_available(extension)

    return [extension]


def get_external_libraries():
    return ['chealpix']
