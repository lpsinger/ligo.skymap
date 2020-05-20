# This file is adapted from the Astropy package template, which is licensed
# under a 3-clause BSD style license - see licenses/TEMPLATE_LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

__all__ = []
if not _ASTROPY_SETUP_:   # noqa
    from .core import omp   # noqa
    __all__ += ['omp']
