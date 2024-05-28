# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure. It needs to live inside the package in order for it to
# get picked up when running the tests inside an interpreter using
# packagename.test

import warnings

from astropy.version import version as astropy_version
import pytest

try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


def pytest_configure(config):

    if ASTROPY_HEADER:

        config.option.astropy_header = True

        # Customize the following lines to add/remove entries from the list of
        # packages for which version numbers are displayed when running the tests.
        PYTEST_HEADER_MODULES.pop('Pandas', None)
        PYTEST_HEADER_MODULES['astropy'] = 'astropy'
        PYTEST_HEADER_MODULES['astropy-healpix'] = 'astropy_healpix'
        PYTEST_HEADER_MODULES['healpy'] = 'healpy'
        PYTEST_HEADER_MODULES['reproject'] = 'reproject'

        from . import __version__
        packagename = 'ligo.skymap'
        TESTED_VERSIONS[packagename] = __version__


# Uncomment the last two lines in this block to treat all DeprecationWarnings as
# exceptions. For Astropy v2.0 or later, there are 2 additional keywords,
# as follow (although default should work for most cases).
# To ignore some packages that produce deprecation warnings on import
# (in addition to 'compiler', 'scipy', 'pygments', 'ipykernel', and
# 'setuptools'), add:
#     modules_to_ignore_on_import=['module_1', 'module_2']
# To ignore some specific deprecation warning messages for Python version
# MAJOR.MINOR or later, add:
#     warnings_to_ignore_by_pyver={(MAJOR, MINOR): ['Message to ignore']}
# from astropy.tests.helper import enable_deprecations_as_exceptions  # noqa
# enable_deprecations_as_exceptions()


@pytest.fixture(autouse=True, scope='session')
def treat_unclosed_files_as_errors():
    """Treat warnings abnout unclosed files as errors

    Many of the command-line tools in :mod:`ligo.skymap.tool` use
    :class:`arparse.FileType` and therefore might leave files opened. Treat
    this as an error.

    """
    warnings.filterwarnings('error', 'unclosed file .*', ResourceWarning)
    warnings.filterwarnings(
        'error', category=pytest.PytestUnraisableExceptionWarning)
