#!/usr/bin/env python
# This file is adapted from the Astropy package template, which is licensed
# under a 3-clause BSD style license - see licenses/TEMPLATE_LICENSE.rst

# NOTE: The configuration for the package, including the name, version, and
# other information are set in the setup.cfg file.

import os
import sys

from setuptools import setup


def get_extensions():
    from setuptools import Extension
    import os

    from extension_helpers import add_openmp_flags_if_available, pkg_config
    import numpy as np

    pkg_config_packages = ['gsl']

    orig_kwargs = {
        'sources': [
            'src/bayestar_distance.c',
            'src/bayestar_moc.c',
            'src/bayestar_sky_map.c',
            'src/core.c',
            'src/cubic_interp.c',
            'src/cubic_interp_test.c',
            'src/find_floor.c',
        ],
        'include_dirs': [
            np.get_include()
        ],
        'define_macros': [
            ('GSL_RANGE_CHECK_OFF', None),
            ('HAVE_INLINE', None),
            ('Py_LIMITED_API', 0x030A0000),
            ('NPY_TARGET_VERSION', 'NPY_1_23_API_VERSION'),
            ('NPY_NO_DEPRECATED_API', 'NPY_2_0_API_VERSION'),
        ],
        'extra_compile_args': [
            '-std=gnu11',
            '-fvisibility=hidden'
        ],
    }

    if os.environ.get('LIGO_SKYMAP_USE_SYSTEM_CHEALPIX'):
        pkg_config_packages.append('chealpix')
    else:
        orig_kwargs['include_dirs'].append('cextern/chealpix')
        orig_kwargs['sources'].append('cextern/chealpix/chealpix.c')

    if os.environ.get('LIGO_SKYMAP_USE_ITTNOTIFY'):
        pkg_config_packages.append('ittnotify')
        orig_kwargs['define_macros'].append(('WITH_ITTNOTIFY', 1))

    kwargs = pkg_config(pkg_config_packages, [])
    for key, orig_value in orig_kwargs.items():
        kwargs.setdefault(key, []).extend(orig_value)

    extension = Extension(name='ligo.skymap.core', language='c',
                          py_limited_api=True, **kwargs)

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


use_scm_version = {'write_to': 'ligo/skymap/version.py',
                   'write_to_template': VERSION_TEMPLATE}

# If we are building under the GitLab CI pipeline and we are building on the
# default branch, then disable the local part of the version
# (+g<short commit hash>) so that we can upload nightly builds to PyPI.
if (
    os.environ.get('CI') == 'true' and
    os.environ.get('CI_COMMIT_BRANCH') == os.environ['CI_DEFAULT_BRANCH']
):
    use_scm_version['local_scheme'] = 'no-local-version'

setup(use_scm_version=use_scm_version, ext_modules=get_extensions())
