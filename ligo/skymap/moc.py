#
# Copyright (C) 2017-2018  Leo Singer
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
Support for HEALPix UNIQ pixel indexing [1]_ and multi-order coverage (MOC)
maps [2]_.

References
----------

.. [1] Reinecke & Hivon, 2015. "Efficient data structures for masks on 2D
       grids." AA 580, A132 <https://doi.org/10.1051/0004-6361/201526549>.
.. [2] Boch et al., 2014. "MOC - HEALPix Multi-Order Coverage map." IVOA
       Recommendation <http://ivoa.net/documents/MOC/>.
"""


from .core import nest2uniq, uniq2nest, uniq2order, uniq2pixarea, uniq2ang
from .core import rasterize as _rasterize
from .util.numpy import add_newdoc_ufunc

__all__ = ('nest2uniq', 'uniq2nest', 'uniq2order', 'uniq2pixarea',
           'uniq2ang', 'rasterize')


add_newdoc_ufunc(nest2uniq, """\
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


add_newdoc_ufunc(uniq2order, """\
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


add_newdoc_ufunc(uniq2pixarea, """\
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


add_newdoc_ufunc(uniq2nest, """\
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
        A multi-order HEALPix dataset stored as a Numpy record array whose
        first column is called UNIQ and contains the NUNIQ pixel index. Every
        point on the unit sphere must be contained in exactly one pixel in the
        dataset.

    Returns
    -------
    nested_data : `numpy.ndarray`
        A fixed-order, NESTED-ordering HEALPix dataset with all of the columns
        that were in moc_data, with the exception of the UNIQ column.
    """
    return _rasterize(moc_data)


del add_newdoc_ufunc
