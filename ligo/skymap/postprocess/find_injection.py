# -*- coding: utf-8 -*-
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
"""
Injection finding for HEALPix sky maps
"""
from collections import namedtuple

import healpy as hp
import numpy as np
from scipy.interpolate import interp1d

from .. import distance
from .. import moc

from .cosmology import dVC_dVL_for_DL

__all__ = ('find_injection_moc', 'FoundInjection')


def flood_fill(nside, ipix, m, nest=False):
    """Stack-based flood fill algorithm in HEALPix coordinates.
    Based on <http://en.wikipedia.org/w/index.php?title=Flood_fill&oldid=566525693#Alternative_implementations>.
    """  # noqa: E501
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


def cos_angle_distance(theta0, phi0, theta1, phi1):
    """Cosine of angular separation in radians between two points on the
    unit sphere."""
    cos_angle_distance = (
        np.cos(phi1 - phi0) * np.sin(theta0) * np.sin(theta1) +
        np.cos(theta0) * np.cos(theta1))
    return np.clip(cos_angle_distance, -1, 1)


def angle_distance(theta0, phi0, theta1, phi1):
    """Angular separation in radians between two points on the unit sphere."""
    return np.arccos(cos_angle_distance(theta0, phi0, theta1, phi1))


# Class to hold return value of find_injection method
FoundInjection = namedtuple(
    'FoundInjection',
    'searched_area searched_prob offset searched_modes contour_areas '
    'area_probs contour_modes searched_prob_dist contour_dists '
    'searched_vol searched_prob_vol contour_vols')


def find_injection_moc(sky_map, true_ra=None, true_dec=None, true_dist=None,
                       contours=(), areas=(), modes=False, nest=False,
                       cosmology=False):
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
        true_idx = i[np.digitize(true_pix, max_ipix[i]) - 1]

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
    prob_for_area = interp1d(area_padded, prob_padded, assume_sorted=True)
    area_for_prob = interp1d(prob_padded, area_padded, assume_sorted=True)

    if true_ra is None:
        searched_area = searched_prob = np.nan
    else:
        # Find the smallest area that would have to be searched to find
        # the true location.
        searched_area = area[true_idx]

        # Find the smallest posterior mass that would have to be searched to
        # find the true location.
        searched_prob = prob[true_idx]

    # Find the contours of the given credible levels.
    contour_idxs = np.digitize(contours, prob) - 1

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
        max_r = 6 * distmean
        if true_dist is not None and np.max(true_dist) > max_r:
            max_r = np.max(true_dist)
        d_r = max_r / n_r
        r = d_r * np.arange(1, n_r)

        # Calculate volume of each voxel, defined as the region within the
        # HEALPix pixel and contained within the two centric spherical shells
        # with radii (r - d_r / 2) and (r + d_r / 2).
        dV = (np.square(r) + np.square(d_r) / 12) * d_r * dA.reshape(-1, 1)

        # Calculate probability within each voxel.
        dP = probdensity.reshape(-1, 1) * dV * np.exp(
            -0.5 * np.square(
                (r.reshape(1, -1) - mu.reshape(-1, 1)) / sigma.reshape(-1, 1)
            )
        ) * (norm / sigma).reshape(-1, 1) / np.sqrt(2 * np.pi)
        dP[np.isnan(dP)] = 0  # Suppress invalid values

        # Calculate probability density per unit volume.

        if cosmology:
            dV *= dVC_dVL_for_DL(r)
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
            i_dist = np.digitize(true_dist, r) - 1
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
