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

    Output above wa rounded for shorter output.
    """
    # Determine HEALPix resolution
    npix = len(m)
    nside = hp.npix2nside(npix)
    min_area = 0.4 * hp.nside2pixarea(nside)

    # Compute faces, vertices, and neighbors.
    # vertices is an N X 3 array of the distinct vertices of the HEALPix faces.
    # faces is an npix X 4 array mapping HEALPix faces to their vertices.
    # neighbors is an npix X 4 array mapping faces to their nearest neighbors.
    faces = np.ascontiguousarray(
        np.rollaxis(hp.boundaries(nside, np.arange(npix), nest=nest), 2, 1))
    dtype = faces.dtype
    faces = faces.view(np.dtype((np.void, dtype.itemsize * 3)))
    vertices, faces = np.unique(faces.ravel(), return_inverse=True)
    faces = faces.reshape(-1, 4)
    vertices = vertices.view(dtype).reshape(-1, 3)
    neighbors = hp.get_all_neighbours(nside, np.arange(npix), nest=nest)[::2].T

    # Loop over the requested contours.
    paths = []
    for level in levels:

        # Find credible region
        indicator = (m >= level)

        # Construct a graph of the edges of the contour.
        graph = nx.Graph()
        face_pairs = set()
        for ipix1, ipix2 in enumerate(neighbors):
            for ipix2 in ipix2:
                # Determine if we have already considered this pair of faces.
                new_face_pair = frozenset((ipix1, ipix2))
                if new_face_pair in face_pairs:
                    continue
                face_pairs.add(new_face_pair)

                # Determine if this pair of faces are on a boundary of the
                # credible level.
                if indicator[ipix1] == indicator[ipix2]:
                    continue

                # Add all common edges of this pair of faces.
                i1 = np.concatenate((faces[ipix1], [faces[ipix1][0]]))
                i2 = np.concatenate((faces[ipix2], [faces[ipix2][0]]))
                edges1 = frozenset(frozenset(_) for _ in zip(i1[:-1], i1[1:]))
                edges2 = frozenset(frozenset(_) for _ in zip(i2[:-1], i2[1:]))
                for edge in edges1 & edges2:
                    graph.add_edge(*edge)
        graph = nx.freeze(graph)

        # Record a closed path for each cycle in the graph.
        cycles = [
            np.take(vertices, cycle, axis=0)
            for cycle in nx.cycle_basis(graph)]

        # Simplify paths if requested
        if simplify:
            cycles = [_simplify(cycle, min_area) for cycle in cycles]
            cycles = [cycle for cycle in cycles if len(cycle) > 2]

        # Convert to lists
        cycles = [
            _vec2radec(cycle, degrees=degrees).tolist() for cycle in cycles]

        # Add to output paths
        paths.append([cycle + [cycle[0]] for cycle in cycles])

    return paths
