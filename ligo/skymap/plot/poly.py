# -*- coding: utf-8 -*-
#
# Copyright (C) 2012-2018  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Plotting tools for drawing polygons
"""
from __future__ import division

from shapely import geometry
import numpy as np
import healpy as hp

# FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
import matplotlib
import distutils.version
mpl_version = distutils.version.LooseVersion(matplotlib.__version__)

from .angle import *

__all__ = ('subdivide_vertices', 'cut_dateline',
           'cut_prime_meridian', 'make_rect_poly')


def subdivide_vertices(vertices, subdivisions):
    """Subdivide a list of vertices by inserting subdivisions additional vertices
    between each original pair of vertices using linear interpolation."""
    subvertices = np.empty((subdivisions * len(vertices), vertices.shape[1]))
    frac = np.atleast_2d(
        np.arange(subdivisions + 1, dtype=float) / subdivisions).T.repeat(
            vertices.shape[1], 1)
    for i in range(len(vertices)):
        subvertices[i*subdivisions:(i+1)*subdivisions] = frac[:0:-1, :] * \
            np.expand_dims(vertices[i-1, :], 0).repeat(subdivisions, 0) + \
            frac[:-1, :] * \
            np.expand_dims(vertices[i, :], 0).repeat(subdivisions, 0)
    return subvertices


def cut_dateline(vertices):
    """Cut a polygon across the dateline, possibly splitting it into multiple
    polygons. Vertices consist of (longitude, latitude) pairs where longitude
    is always given in terms of a reference angle (between -π and π).

    This routine is not meant to cover all possible cases; it will only work
    for convex polygons that extend over less than a hemisphere."""
    vertices = vertices.copy()
    vertices[:, 0] += np.pi
    vertices = cut_prime_meridian(vertices)
    for v in vertices:
        v[:, 0] -= np.pi
    return vertices


def cut_prime_meridian(vertices):
    """Cut a polygon across the prime meridian, possibly splitting it into
    multiple polygons. Vertices consist of (longitude, latitude) pairs where
    longitude is always given in terms of a wrapped angle (between 0 and 2π).

    This routine is not meant to cover all possible cases; it will only work
    for convex polygons that extend over less than a hemisphere."""

    # Ensure that the list of vertices does not contain a repeated endpoint.
    if (vertices[0] == vertices[-1]).all():
        vertices = vertices[:-1]

    # Ensure that the longitudes are wrapped from 0 to 2π.
    vertices = np.column_stack((wrapped_angle(vertices[:, 0]), vertices[:, 1]))

    # Test if the segment consisting of points i-1 and i croses the meridian.
    #
    # If the two longitudes are in [0, 2π), then the shortest arc connecting
    # them crosses the meridian if the difference of the angles is greater
    # than π.
    phis = vertices[:, 0]
    phi0, phi1 = np.sort(np.row_stack((np.roll(phis, 1), phis)), axis=0)
    crosses_meridian = (phi1 - phi0 > np.pi)

    # Count the number of times that the polygon crosses the meridian.
    meridian_crossings = np.sum(crosses_meridian)

    if meridian_crossings == 0:
        # There were zero meridian crossings, so we can use the
        # original vertices as is.
        out_vertices = [vertices]
    elif meridian_crossings == 1:
        # There was one meridian crossing, so the polygon encloses the pole.
        # Any meridian-crossing edge has to be extended
        # into a curve following the nearest polar edge of the map.
        i, = np.flatnonzero(crosses_meridian)
        v0 = vertices[i - 1]
        v1 = vertices[i]

        # Find the latitude at which the meridian crossing occurs by
        # linear interpolation.
        delta_lon = abs(reference_angle(v1[0] - v0[0]))
        lat = (abs(reference_angle(v0[0])) / delta_lon * v0[1] +
               abs(reference_angle(v1[0])) / delta_lon * v1[1])

        # FIXME: Use this simple heuristic to decide which pole to enclose.
        sign_lat = np.sign(np.sum(vertices[:, 1]))

        # Find the closer of the left or the right map boundary for
        # each vertex in the line segment.
        lon_0 = 0. if v0[0] < np.pi else 2*np.pi
        lon_1 = 0. if v1[0] < np.pi else 2*np.pi

        # Set the output vertices to the polar cap plus the original
        # vertices.
        out_vertices = [
            np.vstack((
                vertices[:i],
                [[lon_0, lat],
                 [lon_0, sign_lat * np.pi / 2],
                 [lon_1, sign_lat * np.pi / 2],
                 [lon_1, lat]],
                vertices[i:]))]
    elif meridian_crossings == 2:
        # Since the polygon is assumed to be convex, if there is an even number
        # of meridian crossings, we know that the polygon does not enclose
        # either pole. Then we can use ordinary Euclidean polygon intersection
        # algorithms.

        out_vertices = []

        # Construct polygon representing map boundaries.
        frame_poly = geometry.Polygon(np.asarray([
            [0., np.pi/2],
            [0., -np.pi/2],
            [2*np.pi, -np.pi/2],
            [2*np.pi, np.pi/2]]))

        # Intersect with polygon re-wrapped to lie in [-π, π) or [π, 3π).
        for shift in [0, 2 * np.pi]:
            poly = geometry.Polygon(np.column_stack((
                reference_angle(vertices[:, 0]) + shift, vertices[:, 1])))
            intersection = poly.intersection(frame_poly)
            if intersection:
                assert isinstance(intersection, geometry.Polygon)
                assert intersection.is_simple
                out_vertices += [np.asarray(intersection.exterior)]
    else:
        # There were more than two intersections. Not implemented!
        raise NotImplemented('The polygon intersected the map boundaries two '
                             'or more times, so it is probably not simple and '
                             'convex.')

    # Done!
    return out_vertices


def make_rect_poly(width, height, theta, phi, subdivisions=10):
    """Create a Polygon patch representing a rectangle with half-angles width
    and height rotated from the north pole to (theta, phi)."""

    # Convert width and height to radians, then to Cartesian coordinates.
    w = np.sin(np.deg2rad(width))
    h = np.sin(np.deg2rad(height))

    # Generate vertices of rectangle.
    v = np.asarray([[-w, -h], [w, -h], [w, h], [-w, h]])

    # Subdivide.
    v = subdivide_vertices(v, subdivisions)

    # Project onto sphere by calculating z-coord from normalization condition.
    v = np.hstack((v, np.sqrt(1. - np.expand_dims(np.square(v).sum(1), 1))))

    # Transform vertices.
    v = np.dot(v, hp.rotator.euler_matrix_new(phi, theta, 0, Y=True))

    # Convert to spherical polar coordinates.
    thetas, phis = hp.vec2ang(v)

    # FIXME: Remove this after all Matplotlib monkeypatches are obsolete.
    if mpl_version < '1.2.0':
        # Return list of vertices as longitude, latitude pairs.
        return np.column_stack((reference_angle(phis), 0.5 * np.pi - thetas))
    else:
        # Return list of vertices as longitude, latitude pairs.
        return np.column_stack((wrapped_angle(phis), 0.5 * np.pi - thetas))
