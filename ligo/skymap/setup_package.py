from astropy_helpers.openmp_helpers import add_openmp_flags_if_available
from astropy_helpers.setup_helpers import pkg_config
from distutils.core import Extension


def get_extensions():
    kwargs = pkg_config(['chealpix', 'gsl'], [])
    kwargs['include_dirs'].extend(['numpy'])

    extension = Extension(
        name='ligo.skymap.core',
        language='c',
        sources=[
            'src/cubic_interp.c',
            'src/bayestar_moc.c',
            'src/core.c',
            'src/bayestar_distance.c',
            'src/bayestar_sky_map.c'],
        **kwargs)

    add_openmp_flags_if_available(extension)

    return [extension]
