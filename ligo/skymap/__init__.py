# This file is adapted from the Astropy affiliated package template, which is
# licensed under a 3-clause BSD style license - see astropy_helpers/LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

if not _ASTROPY_SETUP_:
    # For egg_info test builds to pass, put package imports here.

    from .core import omp
