# This file is used to configure the behavior of pytest when using the Astropy
# test infrastructure.
import os

from astropy.tests.plugins.display import (PYTEST_HEADER_MODULES,
                                           TESTED_VERSIONS)

from astropy.tests.helper import enable_deprecations_as_exceptions

## Uncomment the following line to treat all DeprecationWarnings as
## exceptions. For Astropy v2.0 or later, there are 2 additional keywords,
## as follow (although default should work for most cases).
## To ignore some packages that produce deprecation warnings on import
## (in addition to 'compiler', 'scipy', 'pygments', 'ipykernel', and
## 'setuptools'), add:
##     modules_to_ignore_on_import=['module_1', 'module_2']
## To ignore some specific deprecation warning messages for Python version
## MAJOR.MINOR or later, add:
##     warnings_to_ignore_by_pyver={(MAJOR, MINOR): ['Message to ignore']}
# enable_deprecations_as_exceptions()

# Customize the following lines to add/remove entries from
# the list of packages for which version numbers are displayed when running
# the tests.
PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
PYTEST_HEADER_MODULES.pop('h5py', None)

# This is to figure out the package version, rather than
# using Astropy's
from .version import version, astropy_helpers_version

packagename = 'ligo.skymap'
TESTED_VERSIONS[packagename] = version
TESTED_VERSIONS['astropy_helpers'] = astropy_helpers_version
