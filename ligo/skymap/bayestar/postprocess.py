# -*- coding: utf-8 -*-
#
# Copyright (C) 2013-2017  Leo Singer
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
Postprocessing utilities for HEALPix sky maps
"""
from __future__ import division
from __future__ import print_function
import collections
import pkg_resources

from astropy.coordinates import (CartesianRepresentation, SkyCoord,
                                 UnitSphericalRepresentation)
from astropy import units as u
from astropy.wcs import WCS
import healpy as hp
import lal
import lalsimulation
import numpy as np
from scipy.interpolate import interp1d

from .. import distance
from .. import moc
from ..healpix_tree import *

try:
    pkg_resources.require('numpy >= 1.10.0')
except pkg_resources.VersionConflict:
    """FIXME: in numpy < 1.10.0, the digitize() function only works on 1D
    arrays. Remove this workaround once we require numpy >= 1.10.0.

    Example:
    >>> digitize([3], np.arange(5))
    array([4])
    >>> digitize(3, np.arange(5))
    array(4)
    >>> digitize([[3]], np.arange(5))
    array([[4]])
    >>> digitize([], np.arange(5))
    array([], dtype=int64)
    >>> digitize([[], []], np.arange(5))
    array([], shape=(2, 0), dtype=int64)
    """
    def digitize(x, *args, **kwargs):
        x_flat = np.ravel(x)
        x_shape = np.shape(x)
        if len(x_flat) == 0:
            return np.zeros(x_shape, dtype=np.intp)
        else:
            return np.digitize(x_flat, *args, **kwargs).reshape(x_shape)

    # FIXME: np.cov() got the aweights argument in version 1.10.
    # Until we require numpy >= 1.10, this is copied from the Numpy source.
    def cov(m, y=None, rowvar=True, bias=False, ddof=None, fweights=None,
            aweights=None):
        from numpy import array, average, dot, sum
        import warnings

        # Check inputs
        if ddof is not None and ddof != int(ddof):
            raise ValueError(
                "ddof must be integer")

        # Handles complex arrays too
        m = np.asarray(m)
        if m.ndim > 2:
            raise ValueError("m has more than 2 dimensions")

        if y is None:
            dtype = np.result_type(m, np.float64)
        else:
            y = np.asarray(y)
            if y.ndim > 2:
                raise ValueError("y has more than 2 dimensions")
            dtype = np.result_type(m, y, np.float64)

        X = array(m, ndmin=2, dtype=dtype)
        if not rowvar and X.shape[0] != 1:
            X = X.T
        if X.shape[0] == 0:
            return np.array([]).reshape(0, 0)
        if y is not None:
            y = array(y, copy=False, ndmin=2, dtype=dtype)
            if not rowvar and y.shape[0] != 1:
                y = y.T
            X = np.concatenate((X, y), axis=0)

        if ddof is None:
            if bias == 0:
                ddof = 1
            else:
                ddof = 0

        # Get the product of frequencies and weights
        w = None
        if fweights is not None:
            fweights = np.asarray(fweights, dtype=float)
            if not np.all(fweights == np.around(fweights)):
                raise TypeError(
                    "fweights must be integer")
            if fweights.ndim > 1:
                raise RuntimeError(
                    "cannot handle multidimensional fweights")
            if fweights.shape[0] != X.shape[1]:
                raise RuntimeError(
                    "incompatible numbers of samples and fweights")
            if any(fweights < 0):
                raise ValueError(
                    "fweights cannot be negative")
            w = fweights
        if aweights is not None:
            aweights = np.asarray(aweights, dtype=float)
            if aweights.ndim > 1:
                raise RuntimeError(
                    "cannot handle multidimensional aweights")
            if aweights.shape[0] != X.shape[1]:
                raise RuntimeError(
                    "incompatible numbers of samples and aweights")
            if any(aweights < 0):
                raise ValueError(
                    "aweights cannot be negative")
            if w is None:
                w = aweights
            else:
                w *= aweights

        avg, w_sum = average(X, axis=1, weights=w, returned=True)
        w_sum = w_sum[0]

        # Determine the normalization
        if w is None:
            fact = X.shape[1] - ddof
        elif ddof == 0:
            fact = w_sum
        elif aweights is None:
            fact = w_sum - ddof
        else:
            fact = w_sum - ddof*sum(w*aweights)/w_sum

        if fact <= 0:
            warnings.warn("Degrees of freedom <= 0 for slice",
                          RuntimeWarning, stacklevel=2)
            fact = 0.0

        X -= avg[:, None]
        if w is None:
            X_T = X.T
        else:
            X_T = (X*w).T
        c = dot(X, X_T.conj())
        c *= 1. / np.float64(fact)
        return c.squeeze()
else:
    cov = np.cov
    digitize = np.digitize


def flood_fill(nside, ipix, m, nest=False):
    """Stack-based flood fill algorithm in HEALPix coordinates.
    Based on <http://en.wikipedia.org/w/index.php?title=Flood_fill&oldid=566525693#Alternative_implementations>.
    """
    # Initialize stack with starting pixel index.
    stack = [ipix]
    while stack:
        # Pop last pixel off of the stack.
        ipix = stack.pop()
        # Is this pixel in need of filling?
        if m[ipix]:
            # Fill in this pixel.
            m[ipix] = False
            # Find the pixels neighbors.
            neighbors = hp.get_all_neighbours(nside, ipix, nest=nest)
            # All pixels have up to 8 neighbors. If a pixel has less than 8
            # neighbors, then some entries of the array are set to -1. We
            # have to skip those.
            neighbors = neighbors[neighbors != -1]
            # Push neighboring pixels onto the stack.
            stack.extend(neighbors)


def count_modes(m, nest=False):
    """Count the number of modes in a binary HEALPix image by repeatedly
    applying the flood-fill algorithm.

    WARNING: The input array is clobbered in the process."""
    npix = len(m)
    nside = hp.npix2nside(npix)
    for nmodes in range(npix):
        nonzeroipix = np.flatnonzero(m)
        if len(nonzeroipix):
            flood_fill(nside, nonzeroipix[0], m, nest=nest)
        else:
            break
    return nmodes


def count_modes_moc(uniq, i):
    n = len(uniq)
    mask = np.concatenate((np.ones(i + 1, dtype=bool),
                           np.zeros(n - i - 1, dtype=bool)))
    sky_map = np.rec.fromarrays((uniq, mask), names=('UNIQ', 'MASK'))
    sky_map = moc.rasterize(sky_map)['MASK']
    return count_modes(sky_map, nest=True)


def indicator(n, i):
    """Create a binary array of length n that is True for every index that is in
    i and False for every other index. Named after the indicator function."""
    m = np.zeros(n, dtype=np.bool)
    np.put(m, i, True)
    return m


def cos_angle_distance(theta0, phi0, theta1, phi1):
    """Cosine of angular separation in radians between two points on the
    unit sphere."""
    cos_angle_distance = (
        np.cos(phi1 - phi0) * np.sin(theta0) * np.sin(theta1)
        + np.cos(theta0) * np.cos(theta1))
    return np.clip(cos_angle_distance, -1, 1)


def angle_distance(theta0, phi0, theta1, phi1):
    """Angular separation in radians between two points on the unit sphere."""
    return np.arccos(cos_angle_distance(theta0, phi0, theta1, phi1))


# Class to hold return value of find_injection method
FoundInjection = collections.namedtuple(
    'FoundInjection',
    'searched_area searched_prob offset searched_modes contour_areas '
    'area_probs contour_modes searched_prob_dist contour_dists '
    'searched_vol searched_prob_vol contour_vols')


def find_injection_moc(sky_map, true_ra=None, true_dec=None, true_dist=None,
                       contours=(), areas=(), modes=False, nest=False):
    """
    Given a sky map and the true right ascension and declination (in radians),
    find the smallest area in deg^2 that would have to be searched to find the
    source, the smallest posterior mass, and the angular offset in degrees from
    the true location to the maximum (mode) of the posterior. Optionally, also
    compute the areas of and numbers of modes within the smallest contours
    containing a given total probability.
    """

    if (true_ra is None) ^ (true_dec is None):
        raise ValueError('Both true_ra and true_dec must be provided or None')

    contours = np.asarray(contours)

    distmean = sky_map.meta.get('distmean', np.nan)

    # Sort the pixels by descending posterior probability.
    sky_map = np.flipud(np.sort(sky_map, order='PROBDENSITY'))

    # Find the pixel that contains the injection.
    order, ipix = moc.uniq2nest(sky_map['UNIQ'])
    max_order = np.max(order)
    max_nside = hp.order2nside(max_order)
    max_ipix = ipix << np.uint64(2 * (max_order - order))
    ipix = ipix.astype(np.int64)
    max_ipix = max_ipix.astype(np.int64)
    if true_ra is not None:
        true_theta = 0.5 * np.pi - true_dec
        true_phi = true_ra
        true_pix = hp.ang2pix(max_nside, true_theta, true_phi, nest=True)
        i = np.argsort(max_ipix)
        true_idx = i[digitize(true_pix, max_ipix[i]) - 1]

    # Find the angular offset between the mode and true locations.
    mode_theta, mode_phi = hp.pix2ang(
        hp.order2nside(order[0]), ipix[0].astype(np.int64), nest=True)
    if true_ra is None:
        offset = np.nan
    else:
        offset = np.rad2deg(
            angle_distance(true_theta, true_phi, mode_theta, mode_phi))

    # Calculate the cumulative area in deg2 and the cumulative probability.
    dA = moc.uniq2pixarea(sky_map['UNIQ'])
    dP = sky_map['PROBDENSITY'] * dA
    prob = np.cumsum(dP)
    area = np.cumsum(dA) * np.square(180 / np.pi)

    # Construct linear interpolants to map between probability and area.
    # This allows us to compute more accurate contour areas and probabilities
    # under the approximation that the pixels have constant probability
    # density.
    prob_padded = np.concatenate(([0], prob))
    area_padded = np.concatenate(([0], area))
    # FIXME: we should use the assume_sorted=True argument below, but
    # it was added in Scipy 0.14.0, and our Scientific Linux 7 clusters
    # only have Scipy 0.12.1.
    prob_for_area = interp1d(area_padded, prob_padded)
    area_for_prob = interp1d(prob_padded, area_padded)

    if true_ra is None:
        searched_area = searched_prob = np.nan
    else:
        # Find the smallest area that would have to be searched to find
        # the true location.
        searched_area = area[true_idx]

        # Find the smallest posterior mass that would have to be searched to find
        # the true location.
        searched_prob = prob[true_idx]

    # Find the contours of the given credible levels.
    contour_idxs = digitize(contours, prob) - 1

    # For each of the given confidence levels, compute the area of the
    # smallest region containing that probability.
    contour_areas = area_for_prob(contours).tolist()

    # For each listed area, find the probability contained within the
    # smallest credible region of that area.
    area_probs = prob_for_area(areas).tolist()

    if modes:
        if true_ra is None:
            searched_modes = np.nan
        else:
            # Count up the number of modes in each of the given contours.
            searched_modes = count_modes_moc(sky_map['UNIQ'], true_idx)
        contour_modes = [
            count_modes_moc(sky_map['UNIQ'], i) for i in contour_idxs]
    else:
        searched_modes = np.nan
        contour_modes = np.nan

    # Distance stats now...
    if 'DISTMU' in sky_map.dtype.names:
        probdensity = sky_map['PROBDENSITY']
        mu = sky_map['DISTMU']
        sigma = sky_map['DISTSIGMA']
        norm = sky_map['DISTNORM']
        args = (dP, mu, sigma, norm)
        if true_dist is None:
            searched_prob_dist = np.nan
        else:
            searched_prob_dist = distance.marginal_cdf(true_dist, *args)
        # FIXME: old verisons of Numpy can't handle passing zero-length
        # arrays to generalized ufuncs. Remove this workaround once LIGO
        # Data Grid clusters provide a more modern version of Numpy.
        if len(contours) == 0:
            contour_dists = []
        else:
            lo, hi = distance.marginal_ppf(
                np.row_stack((
                    0.5 * (1 - contours),
                    0.5 * (1 + contours)
                )), *args)
            contour_dists = (hi - lo).tolist()

        # Set up distance grid.
        n_r = 1000
        max_r = max(6 * distmean, np.max(true_dist))
        d_r = max_r / n_r
        r = d_r * np.arange(1, n_r)

        # Calculate volume of frustum-shaped voxels with distance centered on r
        # and extending from (r - d_r) to (r + d_r).
        dV = (np.square(r) + np.square(d_r) / 12) * d_r * dA.reshape(-1, 1)

        # Calculate probability within each voxel.
        dP = probdensity.reshape(-1, 1) * dV * np.exp(
            -0.5 * np.square(
                (r.reshape(1, -1) - mu.reshape(-1, 1)) / sigma.reshape(-1, 1)
            )
        ) * (norm / sigma).reshape(-1, 1) / np.sqrt(2 * np.pi)
        dP[np.isnan(dP)] = 0  # Suppress invalid values

        # Calculate probability density per unit volume.
        dP_dV = dP / dV
        i = np.flipud(np.argsort(dP_dV.ravel()))

        P_flat = np.cumsum(dP.ravel()[i])
        V_flat = np.cumsum(dV.ravel()[i])

        contour_vols = interp1d(
            P_flat, V_flat, bounds_error=False)(contours).tolist()
        P = np.empty_like(P_flat)
        V = np.empty_like(V_flat)
        P[i] = P_flat
        V[i] = V_flat
        P = P.reshape(dP.shape)
        V = V.reshape(dV.shape)
        if true_dist is None:
            searched_vol = searched_prob_vol = np.nan
        else:
            i_radec = true_idx
            i_dist = digitize(true_dist, r) - 1
            searched_prob_vol = P[i_radec, i_dist]
            searched_vol = V[i_radec, i_dist]
    else:
        searched_vol = searched_prob_vol = searched_prob_dist = np.nan
        contour_dists = [np.nan] * len(contours)
        contour_vols = [np.nan] * len(contours)

    # Done.
    return FoundInjection(
        searched_area, searched_prob, offset, searched_modes, contour_areas,
        area_probs, contour_modes, searched_prob_dist, contour_dists,
        searched_vol, searched_prob_vol, contour_vols)


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
    try:
        import networkx as nx
    except:
        raise RuntimeError('This function requires the networkx package.')

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

        # Construct a graph of the eges of the contour.
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


def find_greedy_credible_levels(p, ranking=None):
    """Find the greedy credible levels of a (possibly multi-dimensional) array.

    Parameters
    ----------

    p : np.ndarray
        The input array, typically a HEALPix image.

    ranking : np.ndarray, optional
        The array to rank in order to determine the greedy order.
        The default is `p` itself.

    Returns
    -------

    cls : np.ndarray
        An array with the same shape as `p`, with values ranging from `0`
        to `p.sum()`, representing the greedy credible level to which each
        entry in the array belongs.
    """
    p = np.asarray(p)
    pflat = p.ravel()
    if ranking is None:
        ranking = pflat
    else:
        ranking = np.ravel(ranking)
    i = np.flipud(np.argsort(ranking))
    cs = np.cumsum(pflat[i])
    cls = np.empty_like(pflat)
    cls[i] = cs
    return cls.reshape(p.shape)


def get_detector_pair_axis(ifo1, ifo2, gmst):
    """Find the sky position where the line between two detectors pierces the
    celestial sphere.

    Parameters
    ----------

    ifo1 : str or `~lal.Detector` or `~np.ndarray`
        The first detector; either the name of the detector (e.g. `'H1'`), or a
        `lal.Detector` object (e.g., as returned by
        `lalsimulation.DetectorPrefixToLALDetector('H1')` or the geocentric
        Cartesian position of the detection in meters.

    ifo2 : str or `~lal.Detector` or `~np.ndarray`
        The second detector; same as described above.

    gmst : float
        The Greenwich mean sidereal time in radians, as returned by
        `lal.GreenwichMeanSiderealTime`.

    Returns
    -------

    pole_ra : float
        The right ascension in radians at which a ray from `ifo1` to `ifo2`
        would pierce the celestial sphere.

    pole_dec : float
        The declination in radians at which a ray from `ifo1` to `ifo2` would
        pierce the celestial sphere.

    light_travel_time : float
        The light travel time from `ifo1` to `ifo2` in seconds.
    """

    # Get location of detectors if ifo1, ifo2 are LAL detector structs
    try:
        ifo1 = lalsimulation.DetectorPrefixToLALDetector(ifo1)
    except TypeError:
        pass
    try:
        ifo1 = ifo1.location
    except AttributeError:
        pass
    try:
        ifo2 = lalsimulation.DetectorPrefixToLALDetector(ifo2)
    except TypeError:
        pass
    try:
        ifo2 = ifo2.location
    except AttributeError:
        pass

    n = ifo2 - ifo1
    light_travel_time = np.sqrt(np.sum(np.square(n))) / lal.C_SI
    (theta,), (phi,) = hp.vec2ang(n)
    pole_ra = (gmst + phi) % (2 * np.pi)
    pole_dec = 0.5 * np.pi - theta
    return pole_ra, pole_dec, light_travel_time


def rotate_map_to_axis(m, ra, dec, nest=False, method='direct'):
    """Rotate a sky map to place a given line of sight on the +z axis.

    Parameters
    ----------

    m : np.ndarray
        The input HEALPix array.

    ra : float
        Right ascension of axis in radians.

        To specify the axis in geocentric coordinates, supply ra=(lon + gmst),
        where lon is the geocentric longitude and gmst is the Greenwich mean
        sidereal time in radians.

    dec : float
        Declination of axis in radians.

        To specify the axis in geocentric coordinates, supply dec=lat,
        where lat is the geocentric latitude in radians.

    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).

    method : 'direct' or 'fft'
        Select whether to use spherical harmonic transformation ('fft') or
        direct coordinate transformation ('direct')

    Returns
    -------

    m_rotated : np.ndarray
        The rotated HEALPix array.
    """
    npix = len(m)
    nside = hp.npix2nside(npix)

    theta = 0.5 * np.pi - dec
    phi = ra

    if method == 'fft':
        if nest:
            m = hp.reorder(m, n2r=True)
        alm = hp.map2alm(m)
        hp.rotate_alm(alm, -phi, -theta, 0.0)
        ret = hp.alm2map(alm, nside, verbose=False)
        if nest:
            ret = hp.reorder(ret, r2n=True)
    elif method == 'direct':
        R = hp.Rotator(
            rot=np.asarray([0, theta, -phi]),
            deg=False, inv=False, eulertype='Y')
        theta, phi = hp.pix2ang(nside, np.arange(npix), nest=nest)
        ipix = hp.ang2pix(nside, *R(theta, phi), nest=nest)
        ret = m[ipix]
    else:
        raise ValueError('Unrecognized method: {0}'.format(method))

    return ret


def polar_profile(m, nest=False):
    """Obtain the marginalized polar profile of sky map.

    Parameters
    ----------

    m : np.ndarray
        The input HEALPix array.

    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).

    Returns
    -------

    theta : np.ndarray
        The polar angles (i.e., the colatitudes) of the isolatitude rings.

    m_int : np.ndarray
        The normalized probability density, such that `np.trapz(m_int, theta)`
        is approximately `np.sum(m)`.
    """
    npix = len(m)
    nside = hp.npix2nside(npix)
    nrings, = hp.pix2ring(nside, np.asarray([npix]))
    startpix, ringpix, costheta, sintheta, _ = hp.ringinfo(
        nside, np.arange(1, nrings))

    if nest:
        m = hp.reorder(m, n2r=True)

    theta = np.arccos(costheta)
    m_int = np.asarray(
        [m[i:i+j].sum() * stheta * 0.5 * npix / j
         for i, j, stheta in zip(startpix, ringpix, sintheta)])

    return theta, m_int


def smooth_ud_grade(m, nside, nest=False):
    """Resample a sky map to a new resolution using bilinear interpolation.

    Parameters
    ----------

    m : np.ndarray
        The input HEALPix array.

    nest : bool, default=False
        Indicates whether the input sky map is in nested rather than
        ring-indexed HEALPix coordinates (default: ring).

    Returns
    -------

    new_m : np.ndarray
        The resampled HEALPix array. The sum of `m` is approximately preserved.
    """
    npix = hp.nside2npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix), nest=nest)
    new_m = hp.get_interp_val(m, theta, phi, nest=nest)
    return new_m * len(m) / len(new_m)


def posterior_mean(prob, nest=False):
    npix = len(prob)
    nside = hp.npix2nside(npix)
    xyz = hp.pix2vec(nside, np.arange(npix), nest=nest)
    mean_xyz = np.average(xyz, axis=1, weights=prob)
    pos = SkyCoord(*mean_xyz, representation=CartesianRepresentation)
    pos.representation = UnitSphericalRepresentation
    return pos


def posterior_max(prob, nest=False):
    npix = len(prob)
    nside = hp.npix2nside(npix)
    i = np.argmax(prob)
    return SkyCoord(
        *hp.pix2ang(nside, i, nest=nest, lonlat=True), unit=u.deg)


def find_ellipse(prob, cl=90, projection='ARC', nest=False):
    """For a HEALPix map, find an ellipse that contains a given probability.

    Parameters
    ----------
    prob : np.ndarray
        The HEALPix probability map.
    cl : float
        The desired credible level (default: 90).
    projection : str
        The WCS projection (default: 'ARC', or zenithal equidistant).
        For a list of possible values, see:
        http://docs.astropy.org/en/stable/wcs/index.html#supported-projections
    nest : bool
        HEALPix pixel ordering (default: False, or ring ordering).

    Returns
    -------
    ra : float
        The ellipse center right ascension in degrees.
    dec : float
        The ellipse center right ascension in degrees.
    pa : float
        The position angle of the semimajor axis in degrees.
    a : float
        The lenth of the semimajor axis in degrees.
    b : float
        The length o the semiminor axis in degrees.

    Examples
    --------

    I'm not showing the `ra` or `pa` output from the examples below because
    the right ascension is arbitary when dec=90° and the position angle is
    arbitrary when a=b; their arbitrary values may vary depending on your math
    library. Also, I add 0.0 to the outputs because on some platforms you tend
    to get values of dec or pa that get rounded to -0.0, which is within
    numerical precision but would break the doctests (see
    https://stackoverflow.com/questions/11010683).

    Example 1
    ~~~~~~~~~

    This is an example sky map that is uniform in sin(theta) out to a given
    radius in degrees. The 90% credible radius should be 0.9 * radius. (There
    will be deviations for small radius due to finite resolution.)

    >>> def make_uniform_in_sin_theta(radius, nside=512):
    ...     npix = hp.nside2npix(nside)
    ...     theta, phi = hp.pix2ang(nside, np.arange(npix))
    ...     theta_max = np.deg2rad(radius)
    ...     prob = np.where(theta <= theta_max, 1 / np.sin(theta), 0)
    ...     return prob / prob.sum()
    ...

    >>> prob = make_uniform_in_sin_theta(1)
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b)
    90.0 0.82241 0.82241

    >>> prob = make_uniform_in_sin_theta(10)
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b)
    90.0 9.05512 9.05512

    >>> prob = make_uniform_in_sin_theta(120)
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b)
    90.0 107.9745 107.9745

    Example 2
    ~~~~~~~~~

    These are approximately Gaussian distributions.

    >>> from scipy import stats
    >>> def make_gaussian(mean, cov, nside=512):
    ...     npix = hp.nside2npix(nside)
    ...     xyz = np.transpose(hp.pix2vec(nside, np.arange(npix)))
    ...     # FIXME: stats.multivariate_normal was added in scipy 0.14,
    ...     # but we still need to support an older version on our
    ...     # Scientific Linux 7 clusters.
    ...     #
    ...     # dist = stats.multivariate_normal(mean, cov)
    ...     # prob = dist.pdf(xyz)
    ...     if np.ndim(cov) == 0:
    ...         cov = [cov] * 3
    ...     if np.ndim(cov) == 1:
    ...         cov = np.diag(cov)
    ...     d = xyz - mean
    ...     prob = np.exp(-0.5 * (np.linalg.solve(cov, d.T) * d.T).sum(0))
    ...     return prob / prob.sum()
    ...

    This one is centered at RA=45°, Dec=0° and has a standard deviation of ~1°.

    >>> prob = make_gaussian(
    ...     [1/np.sqrt(2), 1/np.sqrt(2), 0],
    ...     np.square(np.deg2rad(1)))
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(ra, dec, a, b)
    45.0 0.0 2.14209 2.14209

    This one is centered at RA=45°, Dec=0°, and is elongated in the north-south
    direction.

    >>> prob = make_gaussian(
    ...     [1/np.sqrt(2), 1/np.sqrt(2), 0],
    ...     np.diag(np.square(np.deg2rad([1, 1, 10]))))
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(ra, dec, pa, a, b)
    45.0 0.0 0.0 13.44746 2.1082

    This one is centered at RA=0°, Dec=0°, and is elongated in the east-west
    direction.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     np.diag(np.square(np.deg2rad([1, 10, 1]))))
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, pa, a, b)
    0.0 90.0 13.4194 2.1038

    This one is centered at RA=0°, Dec=0°, and is tilted about 10° to the west
    of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 0.1, -0.15],
    ...      [0, -0.15, 1]])
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, pa, a, b)
    0.0 170.78253 63.82809 34.00824

    This one is centered at RA=0°, Dec=0°, and is tilted about 10° to the east
    of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 0.1, 0.15],
    ...      [0, 0.15, 1]])
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, pa, a, b)
    0.0 9.21747 63.82809 34.00824

    This one is centered at RA=0°, Dec=0°, and is tilted about 80° to the east
    of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 1, 0.15],
    ...      [0, 0.15, 0.1]])
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, pa, a, b)
    0.0 80.78252 63.82533 34.00677

    This one is centered at RA=0°, Dec=0°, and is tilted about 80° to the west
    of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 1, -0.15],
    ...      [0, -0.15, 0.1]])
    ...
    >>> ra, dec, pa, a, b = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, pa, a, b)
    0.0 99.21748 63.82533 34.00677
    """
    npix = len(prob)
    nside = hp.npix2nside(npix)

    # Find mean right ascension and declination.
    xyz0 = (hp.pix2vec(nside, np.arange(npix), nest=nest) * prob).sum(axis=1)
    (ra,), (dec,) = hp.vec2ang(xyz0, lonlat=True)

    # Construct WCS with the specified projection
    # and centered on mean direction.
    w = WCS()
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ['RA---' + projection, 'DEC--' + projection]

    # Transform HEALPix to zenithal equidistant coordinates.
    xy = w.wcs_world2pix(
        np.transpose(
            hp.pix2ang(
                nside, np.arange(npix), nest=nest, lonlat=True)), 1)

    # Keep only values that were inside the projection.
    keep = np.logical_and.reduce(np.isfinite(xy), axis=1)
    xy = xy[keep]
    prob = prob[keep]

    # Find covariance matrix.
    c = cov(xy, aweights=prob, rowvar=False)

    # If each point is n-sigma from the center, find n.
    nsigmas = np.sqrt(np.sum(xy.T * np.linalg.solve(c, xy.T), axis=0))

    # Find the number of sigma that enclose the cl% credible level.
    i = np.argsort(nsigmas)
    nsigmas = nsigmas[i]
    cls = np.cumsum(prob[i])
    nsigma = np.interp(1e-2 * cl, cls, nsigmas)

    # If the credible level is not within the projection,
    # then stop here and return all nans.
    if 1e-2 * cl > cls[-1]:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # Find the eigendecomposition of the covariance matrix.
    w, v = np.linalg.eigh(c)

    # Find the semi-minor and semi-major axes.
    b, a = nsigma * np.sqrt(w)

    # Find the position angle.
    pa = np.rad2deg(np.arctan2(*v[:, 1]))

    # An ellipse is symmetric under rotations of 180°.
    # Return the smallest possible positive position angle.
    pa %= 180

    # Done!
    return ra, dec, pa, a, b
