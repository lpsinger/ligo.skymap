#
# Copyright (C) 2013-2016  Leo Singer
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
Multiresolution HEALPix trees
"""
import numpy as np
import healpy as hp
import collections
import itertools

__all__ = ('HEALPIX_MACHINE_ORDER', 'HEALPIX_MACHINE_NSIDE', 'HEALPixTree',
           'adaptive_healpix_histogram', 'interpolate_nested',
           'reconstruct_nested')


# Maximum 64-bit HEALPix resolution.
HEALPIX_MACHINE_ORDER = 29
HEALPIX_MACHINE_NSIDE = hp.order2nside(HEALPIX_MACHINE_ORDER)


_HEALPixTreeVisitExtra = collections.namedtuple(
    'HEALPixTreeVisit', 'nside full_nside ipix ipix0 ipix1 value')


_HEALPixTreeVisit = collections.namedtuple(
    'HEALPixTreeVisit', 'nside ipix')


class HEALPixTree(object):
    """Data structure used internally by the function
    adaptive_healpix_histogram()."""

    def __init__(
            self, samples, max_samples_per_pixel, max_order,
            order=0, needs_sort=True):
        if needs_sort:
            samples = np.sort(samples)
        if len(samples) >= max_samples_per_pixel and order < max_order:
            # All nodes have 4 children, except for the root node,
            # which has 12.
            nchildren = 12 if order == 0 else 4
            self.samples = None
            self.children = [
                HEALPixTree(
                    [], max_samples_per_pixel, max_order, order=order + 1)
                for i in range(nchildren)]
            for ipix, samples in itertools.groupby(
                    samples, self.key_for_order(order)):
                self.children[np.uint64(ipix % nchildren)] = HEALPixTree(
                    list(samples), max_samples_per_pixel, max_order,
                    order=order + 1, needs_sort=False)
        else:
            # There are few enough samples that we can make this cell a leaf.
            self.samples = list(samples)
            self.children = None

    @staticmethod
    def key_for_order(order):
        """Create a function that downsamples full-resolution pixel indices."""
        return lambda ipix: ipix >> np.uint64(
            2 * (HEALPIX_MACHINE_ORDER - order))

    @property
    def order(self):
        """Return the maximum HEALPix order required to represent this tree,
        which is the same as the tree depth."""
        if self.children is None:
            return 0
        else:
            return 1 + max(child.order for child in self.children)

    def _visit(self, order, full_order, ipix, extra):
        if self.children is None:
            nside = 1 << order
            full_nside = 1 << order
            ipix0 = ipix << 2 * (full_order - order)
            ipix1 = (ipix + 1) << 2 * (full_order - order)
            if extra:
                yield _HEALPixTreeVisitExtra(
                    nside, full_nside, ipix, ipix0, ipix1, self.samples)
            else:
                yield _HEALPixTreeVisit(nside, ipix)
        else:
            for i, child in enumerate(self.children):
                yield from child._visit(
                    order + 1, full_order, (ipix << 2) + i, extra)

    def _visit_depthfirst(self, extra):
        order = self.order
        for ipix, child in enumerate(self.children):
            yield from child._visit(0, order, ipix, extra)

    def _visit_breadthfirst(self, extra):
        return sorted(
            self._visit_depthfirst(extra), lambda _: (_.nside, _.ipix))

    def visit(self, order='depthfirst', extra=True):
        """Traverse the leaves of the HEALPix tree.

        Parameters
        ----------
        order : string, optional
            Traversal order: 'depthfirst' (the default) or 'breadthfirst'.

        extra : bool
            Whether to output extra information about the pixel
            (default is True).

        Yields
        ------
        nside : int
            The HEALPix resolution of the node.

        full_nside : int, present if extra=True
            The HEALPix resolution of the deepest node in the tree.

        ipix : int
            The nested HEALPix index of the node.

        ipix0 : int, present if extra=True
            The start index of the range of pixels spanned by the node at the
            resolution `full_nside`.

        ipix1 : int, present if extra=True
            The end index of the range of pixels spanned by the node at the
            resolution `full_nside`.

        samples : list, present if extra=True
            The list of samples contained in the node.

        Examples
        --------

        >>> ipix = np.arange(12, dtype=np.uint64) * HEALPIX_MACHINE_NSIDE**2
        >>> tree = HEALPixTree(ipix, max_samples_per_pixel=1, max_order=1)
        >>> [tuple(_) for _ in tree.visit(extra=False)]
        [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11)]
        """
        funcs = {'depthfirst': self._visit_depthfirst,
                 'breadthfirst': self._visit_breadthfirst}
        func = funcs[order]
        yield from func(extra)

    @property
    def flat_bitmap(self):
        """Return flattened HEALPix representation."""
        m = np.empty(hp.nside2npix(hp.order2nside(self.order)))
        for nside, full_nside, ipix, ipix0, ipix1, samples in self.visit():
            m[ipix0:ipix1] = len(samples) / hp.nside2pixarea(nside)
        return m


def adaptive_healpix_histogram(
        theta, phi, max_samples_per_pixel, nside=-1, max_nside=-1, nest=False):
    """Adaptively histogram the posterior samples represented by the
    (theta, phi) points using a recursively subdivided HEALPix tree. Nodes are
    subdivided until each leaf contains no more than max_samples_per_pixel
    samples. Finally, the tree is flattened to a fixed-resolution HEALPix image
    with a resolution appropriate for the depth of the tree. If nside is
    specified, the result is resampled to another desired HEALPix resolution.
    """
    # Calculate pixel index of every sample, at the maximum 64-bit resolution.
    #
    # At this resolution, each pixel is only 0.2 mas across; we'll use the
    # 64-bit pixel indices as a proxy for the true sample coordinates so that
    # we don't have to do any trigonometry (aside from the initial hp.ang2pix
    # call).
    #
    # FIXME: Cast to uint64 needed because Healpy returns signed indices.
    ipix = hp.ang2pix(
        HEALPIX_MACHINE_NSIDE, theta, phi, nest=True).astype(np.uint64)

    # Build tree structure.
    if nside == -1 and max_nside == -1:
        max_order = HEALPIX_MACHINE_ORDER
    elif nside == -1:
        max_order = hp.nside2order(max_nside)
    elif max_nside == -1:
        max_order = hp.nside2order(nside)
    else:
        max_order = hp.nside2order(min(nside, max_nside))
    tree = HEALPixTree(ipix, max_samples_per_pixel, max_order)

    # Compute a flattened bitmap representation of the tree.
    p = tree.flat_bitmap

    # If requested, resample the tree to the output resolution.
    if nside != -1:
        p = hp.ud_grade(p, nside, order_in='NESTED', order_out='NESTED')

    # Normalize.
    p /= np.sum(p)

    if not nest:
        p = hp.reorder(p, n2r=True)

    # Done!
    return p


def _interpolate_level(m):
    """Recursive multi-resolution interpolation. Modifies `m` in place."""
    # Determine resolution.
    npix = len(m)

    if npix > 12:
        # Determine which pixels comprise multi-pixel tiles.
        ipix = np.flatnonzero(
            (m[0::4] == m[1::4]) &
            (m[0::4] == m[2::4]) &
            (m[0::4] == m[3::4]))

        if len(ipix):
            ipix = 4 * ipix + np.expand_dims(np.arange(4, dtype=np.intp), 1)
            ipix = ipix.T.ravel()

            nside = hp.npix2nside(npix)

            # Downsample.
            m_lores = hp.ud_grade(
                m, nside // 2, order_in='NESTED', order_out='NESTED')

            # Interpolate recursively.
            _interpolate_level(m_lores)

            # Record interpolated multi-pixel tiles.
            m[ipix] = hp.get_interp_val(
                m_lores, *hp.pix2ang(nside, ipix, nest=True), nest=True)


def interpolate_nested(m, nest=False):
    """
    Apply bilinear interpolation to a multiresolution HEALPix map, assuming
    that runs of pixels containing identical values are nodes of the tree. This
    smooths out the stair-step effect that may be noticeable in contour plots.

    Here is how it works. Consider a coarse tile surrounded by base tiles, like
    this::

                +---+---+
                |   |   |
                +-------+
                |   |   |
        +---+---+---+---+---+---+
        |   |   |       |   |   |
        +-------+       +-------+
        |   |   |       |   |   |
        +---+---+---+---+---+---+
                |   |   |
                +-------+
                |   |   |
                +---+---+

    The value within the central coarse tile is computed by downsampling the
    sky map (averaging the fine tiles), upsampling again (with bilinear
    interpolation), and then finally copying the interpolated values within the
    coarse tile back to the full-resolution sky map. This process is applied
    recursively at all successive HEALPix resolutions.

    Note that this method suffers from a minor discontinuity artifact at the
    edges of regions of coarse tiles, because it temporarily treats the
    bordering fine tiles as constant. However, this artifact seems to have only
    a minor effect on generating contour plots.

    Parameters
    ----------

    m: `~numpy.ndarray`
        a HEALPix array

    nest: bool, default: False
        Whether the input array is stored in the `NESTED` indexing scheme
        (True) or the `RING` indexing scheme (False).

    """
    # Convert to nest indexing if necessary, and make sure that we are working
    # on a copy.
    if nest:
        m = m.copy()
    else:
        m = hp.reorder(m, r2n=True)

    _interpolate_level(m)

    # Convert to back ring indexing if necessary
    if not nest:
        m = hp.reorder(m, n2r=True)

    # Done!
    return m


def _reconstruct_nested_breadthfirst(m, extra):
    m = np.asarray(m)
    max_npix = len(m)
    max_nside = hp.npix2nside(max_npix)
    max_order = hp.nside2order(max_nside)
    seen = np.zeros(max_npix, dtype=bool)

    for order in range(max_order + 1):
        nside = hp.order2nside(order)
        npix = hp.nside2npix(nside)
        skip = max_npix // npix
        if skip > 1:
            b = m.reshape(-1, skip)
            a = b[:, 0].reshape(-1, 1)
            b = b[:, 1:]
            aseen = seen.reshape(-1, skip)
            eq = ((a == b) | ((a != a) & (b != b))).all(1) & (~aseen).all(1)
        else:
            eq = ~seen
        for ipix in np.flatnonzero(eq):
            ipix0 = ipix * skip
            ipix1 = (ipix + 1) * skip
            seen[ipix0:ipix1] = True
            if extra:
                yield _HEALPixTreeVisitExtra(
                    nside, max_nside, ipix, ipix0, ipix1, m[ipix0])
            else:
                yield _HEALPixTreeVisit(nside, ipix)


def _reconstruct_nested_depthfirst(m, extra):
    result = sorted(
        _reconstruct_nested_breadthfirst(m, True),
        key=lambda _: _.ipix0)
    if not extra:
        result = (_HEALPixTreeVisit(_.nside, _.ipix) for _ in result)
    return result


def reconstruct_nested(m, order='depthfirst', extra=True):
    """Reconstruct the leaves of a multiresolution tree.

    Parameters
    ----------
    m : `~numpy.ndarray`
        A HEALPix array in the NESTED ordering scheme.

    order : {'depthfirst', 'breadthfirst'}, optional
        Traversal order: 'depthfirst' (the default) or 'breadthfirst'.

    extra : bool
        Whether to output extra information about the pixel (default is True).

    Yields
    ------
    nside : int
        The HEALPix resolution of the node.

    full_nside : int, present if extra=True
        The HEALPix resolution of the deepest node in the tree.

    ipix : int
        The nested HEALPix index of the node.

    ipix0 : int, present if extra=True
        The start index of the range of pixels spanned by the node at the
        resolution `full_nside`.

    ipix1 : int, present if extra=True
        The end index of the range of pixels spanned by the node at the
        resolution `full_nside`.

    value : list, present if extra=True
        The value of the map at the node.

    Examples
    --------

    An nside=1 array of all zeros:

    >>> m = np.zeros(12)
    >>> result = reconstruct_nested(m, order='breadthfirst', extra=False)
    >>> [tuple(_) for _ in result]
    [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11)]

    An nside=1 array of distinct values:

    >>> m = range(12)
    >>> result = reconstruct_nested(m, order='breadthfirst', extra=False)
    >>> [tuple(_) for _ in result]
    [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11)]

    An nside=8 array of zeros:

    >>> m = np.zeros(768)
    >>> result = reconstruct_nested(m, order='breadthfirst', extra=False)
    >>> [tuple(_) for _ in result]
    [(1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11)]

    An nside=2 array, all zeros except for four consecutive distinct elements:

    >>> m = np.zeros(48); m[:4] = range(4)
    >>> result = reconstruct_nested(m, order='breadthfirst', extra=False)
    >>> [tuple(_) for _ in result]
    [(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11), (2, 0), (2, 1), (2, 2), (2, 3)]

    Same, but in depthfirst order:

    >>> result = reconstruct_nested(m, order='depthfirst', extra=False)
    >>> [tuple(_) for _ in result]
    [(2, 0), (2, 1), (2, 2), (2, 3), (1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6), (1, 7), (1, 8), (1, 9), (1, 10), (1, 11)]

    An nside=2 array, all elements distinct except for four consecutive zeros:

    >>> m = np.arange(48); m[:4] = 0
    >>> result = reconstruct_nested(m, order='breadthfirst', extra=False)
    >>> [tuple(_) for _ in result]
    [(1, 0), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (2, 9), (2, 10), (2, 11), (2, 12), (2, 13), (2, 14), (2, 15), (2, 16), (2, 17), (2, 18), (2, 19), (2, 20), (2, 21), (2, 22), (2, 23), (2, 24), (2, 25), (2, 26), (2, 27), (2, 28), (2, 29), (2, 30), (2, 31), (2, 32), (2, 33), (2, 34), (2, 35), (2, 36), (2, 37), (2, 38), (2, 39), (2, 40), (2, 41), (2, 42), (2, 43), (2, 44), (2, 45), (2, 46), (2, 47)]
    """
    funcs = {'depthfirst': _reconstruct_nested_depthfirst,
             'breadthfirst': _reconstruct_nested_breadthfirst}
    func = funcs[order]
    yield from func(m, extra)
