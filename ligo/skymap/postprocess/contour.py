#
# Copyright (C) 2013-2018  Leo Singer
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

import healpy as hp
import networkx as nx
import numpy as np

__all__ = ('contour',)


def _norm(vertices):
    return np.sqrt(np.sum(np.square(vertices), -1))


def _adjacent_triangle_areas(vertices):
    return 0.5 * _norm(np.cross(
        np.roll(vertices, -1, axis=0) - vertices,
        np.roll(vertices, +1, axis=0) - vertices))


def _simplify(vertices, min_area):
    """Visvalingam's algorithm (see http://bost.ocks.org/mike/simplify/)
    for linear rings on a sphere. This is a naive, slow implementation."""
    area = _adjacent_triangle_areas(vertices)

    while True:
        i_min_area = np.argmin(area)
        if area[i_min_area] > min_area:
            break

        vertices = np.delete(vertices, i_min_area, axis=0)
        area = np.delete(area, i_min_area)
        new_area = _adjacent_triangle_areas(vertices)
        area = np.maximum(area, new_area)

    return vertices


def _vec2radec(vertices, degrees=False):
    theta, phi = hp.vec2ang(np.asarray(vertices))
    ret = np.column_stack((phi % (2 * np.pi), 0.5 * np.pi - theta))
    if degrees:
        ret = np.rad2deg(ret)
    return ret


def contour(m, levels, nest=False, degrees=False, simplify=True):
    """Calculate contours from a HEALPix dataset.

    Parameters
    ----------
    m : `numpy.ndarray`
        The HEALPix dataset.
    levels : list
        The list of contour values.
    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).
    degrees : bool, default=False
        Whether the contours are in degrees instead of radians.
    simplify : bool, default=True
        Whether to simplify the paths.

    Returns
    -------
    list
        A list with the same length as `levels`.
        Each item is a list of disjoint polygons, of which each item is a
        list of points, of which each is a list consisting of the right
        ascension and declination.

    Examples
    --------

    A very simply example sky map...

    >>> nside = 32
    >>> npix = hp.nside2npix(nside)
    >>> ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True)
    >>> m = dec
    >>> contour(m, [10, 20, 30], degrees=True)
    [[[[..., ...], ...], ...], ...]

    Output above was rounded for shorter output.
    """
    # Determine HEALPix resolution.
    npix = len(m)
    nside = hp.npix2nside(npix)
    min_area = 0.4 * hp.nside2pixarea(nside)

    neighbors = hp.get_all_neighbours(nside, np.arange(npix), nest=nest).T

    # Loop over the requested contours.
    paths = []
    for level in levels:

        # Find credible region.
        indicator = (m >= level)

        # Find all faces that lie on the boundary.
        # This speeds up the doubly nested ``for`` loop below by allowing us to
        # skip the vast majority of faces that are on the interior or the
        # exterior of the contour.
        tovisit = np.flatnonzero(
            np.any(indicator.reshape(-1, 1) !=
                   indicator[neighbors[:, ::2]], axis=1))

        # Construct a graph of the edges of the contour.
        graph = nx.Graph()
        face_pairs = set()
        for ipix1 in tovisit:
            neighborhood = neighbors[ipix1]
            for _ in range(4):
                neighborhood = np.roll(neighborhood, 2)
                ipix2 = neighborhood[4]

                # Skip this pair of faces if we have already examined it.
                new_face_pair = frozenset((ipix1, ipix2))
                if new_face_pair in face_pairs:
                    continue
                face_pairs.add(new_face_pair)

                # Determine if this pair of faces are on a boundary of the
                # credible level.
                if indicator[ipix1] == indicator[ipix2]:
                    continue

                # Add the common edge of this pair of faces.
                # Label each vertex with the set of faces that they share.
                graph.add_edge(
                    frozenset((ipix1, *neighborhood[2:5])),
                    frozenset((ipix1, *neighborhood[4:7])))
        graph = nx.freeze(graph)

        # Find contours by detecting cycles in the graph.
        cycles = nx.cycle_basis(graph)

        # Construct the coordinates of the vertices by averaging the
        # coordinates of the connected faces.
        cycles = [[
            np.sum(hp.pix2vec(nside, [i for i in v if i != -1], nest=nest), 1)
            for v in cycle] for cycle in cycles]

        # Simplify paths if requested.
        if simplify:
            cycles = [_simplify(cycle, min_area) for cycle in cycles]
            cycles = [cycle for cycle in cycles if len(cycle) > 2]

        # Convert to angles.
        cycles = [
            _vec2radec(cycle, degrees=degrees).tolist() for cycle in cycles]

        # Add to output paths.
        paths.append([cycle + [cycle[0]] for cycle in cycles])

    return paths
