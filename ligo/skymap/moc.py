#
# Copyright (C) 2017  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
Multi-order coverage (MOC) HEALPix indexing.
"""


import numpy as np
from .core import nest2uniq, uniq2nest, uniq2order, uniq2pixarea
from .core import rasterize as _rasterize

__all__ = ('nest2uniq', 'uniq2nest', 'uniq2order', 'uniq2pixarea', 'rasterize')


def _add_newdoc_ufunc(func, doc):
    # The function `np.lib.add_newdoc_ufunc` can only change a ufunc's
    # docstring if it is `NULL`. This workaround avoids an exception
    # when the user tries to `reload()` this module.
    try:
        np.lib.add_newdoc_ufunc(func, doc)
    except ValueError as e:
        msg = 'Cannot change docstring of ufunc with non-NULL docstring'
        if e.args[0] == msg:
            pass


_add_newdoc_ufunc(nest2uniq, """\
Convert a pixel index from NESTED to NUNIQ ordering.

Parameters
----------
order : `numpy.ndarray`
    HEALPix resolution order, the logarithm base 2 of `nside`
ipix : `numpy.ndarray`
    NESTED pixel index

Returns
-------
uniq : `numpy.ndarray`
    NUNIQ pixel index
""")


_add_newdoc_ufunc(uniq2order, """\
Determine the HEALPix resolution order of a HEALPix NESTED index.

Parameters
----------
uniq : `numpy.ndarray`
    NUNIQ pixel index

Returns
-------
order : `numpy.ndarray`
    HEALPix resolution order, the logarithm base 2 of `nside`
""")


_add_newdoc_ufunc(uniq2pixarea, """\
Determine the area of a HEALPix NESTED index.

Parameters
----------
uniq : `numpy.ndarray`
    NUNIQ pixel index

Returns
-------
area : `numpy.ndarray`
    The pixel's area in steradians
""")


_add_newdoc_ufunc(uniq2nest, """\
Convert a pixel index from NUNIQ to NESTED ordering.

Parameters
----------
uniq : `numpy.ndarray`
    NUNIQ pixel index

Returns
-------
order : `numpy.ndarray`
    HEALPix resolution order (logarithm base 2 of `nside`)
ipix : `numpy.ndarray`
    NESTED pixel index
""")


def rasterize(moc_data):
    """Convert a multi-order HEALPix dataset to fixed-order NESTED ordering.

    Parameters
    ----------
    moc_data : `numpy.ndarray`
    A multi-order HEALPix dataset stored as a Numpy record array whose first
    column is called UNIQ and contains the NUNIQ pixel index. Every point on
    the unit sphere must be contained in exactly one pixel in the dataset.

    Returns
    -------
    nested_data : `numpy.ndarray`
        A fixed-order, NESTED-ordering HEALPix dataset with all of the columns
        that were in moc_data, with the exception of the UNIQ column.
    """
    return _rasterize(moc_data)
