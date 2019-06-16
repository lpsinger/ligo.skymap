# This file is adapted from the Astropy affiliated package template, which is
# licensed under a 3-clause BSD style license - see astropy_helpers/LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

# Enforce Python version check during package import.
# This is the same check as the one at the top of setup.py
import sys
from distutils.version import LooseVersion

__minimum_python_version__ = "3.6"

__all__ = []


class UnsupportedPythonError(Exception):
    pass


if LooseVersion(sys.version) < LooseVersion(__minimum_python_version__):
    raise UnsupportedPythonError("ligo.skymap does not support Python < {}"
                                 .format(__minimum_python_version__))

if not _ASTROPY_SETUP_:   # noqa
    # For egg_info test builds to pass, put package imports here.
    from .core import omp   # noqa
    # Then you can be explicit to control what ends up in the namespace,
    __all__ += ['omp']   # noqa
    # or you can keep everything from the subpackage with the following instead
    # __all__ += example_mod.__all__
