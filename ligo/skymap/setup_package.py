def get_extensions():
    from astropy_helpers import (distutils_helpers, openmp_helpers,
                                 setup_helpers)
    from distutils.core import Extension

    pkg_config_packages = ['gsl']

    sources = [
        'src/bayestar_distance.c',
        'src/bayestar_moc.c',
        'src/bayestar_sky_map.c',
        'src/core.c',
        'src/cubic_interp.c',
        'src/cubic_interp_test.c',
    ]

    include_dirs = ['cextern/numpy', 'numpy']

    if setup_helpers.use_system_library('chealpix'):
        pkg_config_packages.append('chealpix')
    else:
        include_dirs.append('cextern/chealpix')
        sources.append('cextern/chealpix/chealpix.c')

    kwargs = setup_helpers.pkg_config(pkg_config_packages, [])
    kwargs['include_dirs'].extend(include_dirs)
    kwargs['extra_compile_args'].extend(['-std=gnu99',
                                         '-DGSL_RANGE_CHECK_OFF'])

    if distutils_helpers.get_distutils_build_option('with_ittnotify'):
        kwargs.setdefault('define_macros', []).append(('WITH_ITTNOTIFY', 1))
        kwargs.setdefault('libraries', []).append('ittnotify')

    extension = Extension(
        name='ligo.skymap.core', language='c', sources=sources, **kwargs)

    openmp_helpers.add_openmp_flags_if_available(extension)

    return [extension]


def get_external_libraries():
    return ['chealpix']


def get_build_options():
    return [
        ('with-ittnotify', "Instrument code with Intel's ittnotify API.", True)
     ]
