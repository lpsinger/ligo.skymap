#!/usr/bin/env python
# This file is adapted from the Astropy package template, which is licensed
# under a 3-clause BSD style license - see licenses/TEMPLATE_LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file.

import sys

from setuptools import setup


def get_extensions():
    from distutils.core import Extension
    import os

    from extension_helpers import add_openmp_flags_if_available, pkg_config
    import numpy as np

    pkg_config_packages = ['gsl']

    sources = [
        'src/bayestar_distance.c',
        'src/bayestar_moc.c',
        'src/bayestar_sky_map.c',
        'src/core.c',
        'src/cubic_interp.c',
        'src/cubic_interp_test.c',
    ]

    include_dirs = [np.get_include()]

    if os.environ.get('LIGO_SKYMAP_USE_SYSTEM_CHEALPIX'):
        pkg_config_packages.append('chealpix')
    else:
        include_dirs.append('cextern/chealpix')
        sources.append('cextern/chealpix/chealpix.c')

    if os.environ.get('LIGO_SKYMAP_USE_ITTNOTIFY'):
        pkg_config_packages.append('ittnotify')

    kwargs = pkg_config(pkg_config_packages, [])
    kwargs['include_dirs'].extend(include_dirs)
    kwargs['extra_compile_args'].extend(['-std=gnu11',
                                         '-DGSL_RANGE_CHECK_OFF',
                                         '-DHAVE_INLINE'])

    if os.environ.get('LIGO_SKYMAP_USE_ITTNOTIFY'):
        kwargs.setdefault('define_macros', []).append(('WITH_ITTNOTIFY', 1))

    extension = Extension(name='ligo.skymap.core', language='c',
                          py_limited_api=True, sources=sources, **kwargs)

    if not os.environ.get('LIGO_SKYMAP_DISABLE_OPENMP'):
        add_openmp_flags_if_available(extension)

    return [extension]


# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    tox -e test

If you don't already have tox installed, you can install it with:

    pip install tox

If you only want to run part of the test suite, you can also use pytest
directly with::

    pip install -e .[test]
    pytest

For more information, see:

  http://docs.astropy.org/en/latest/development/testguide.html#running-tests
"""

if 'test' in sys.argv:
    print(TEST_HELP)
    sys.exit(1)

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py build_docs'. Instead you will need to run:

    tox -e build_docs

If you don't already have tox installed, you can install it with:

    pip install tox

You can also build the documentation with Sphinx directly using::

    pip install -e .[docs]
    cd docs
    make html

For more information, see:

  http://docs.astropy.org/en/latest/install.html#builddocs
"""

if 'build_docs' in sys.argv or 'build_sphinx' in sys.argv:
    print(DOCS_HELP)
    sys.exit(1)

VERSION_TEMPLATE = """
# Note that we need to fall back to the hard-coded version if either
# setuptools_scm can't be imported or setuptools_scm can't determine the
# version, so we catch the generic 'Exception'.
try:
    from setuptools_scm import get_version
    version = get_version(root='..', relative_to=__file__)
except Exception:
    version = '{version}'
""".lstrip()

setup(use_scm_version={'write_to': 'ligo/skymap/version.py',
                       'write_to_template': VERSION_TEMPLATE},
      ext_modules=get_extensions())
