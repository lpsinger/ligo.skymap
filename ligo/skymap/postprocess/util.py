#
# Copyright (C) 2013-2024  Leo Singer
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""Postprocessing utilities for HEALPix sky maps."""

import astropy_healpix as ah
from astropy.coordinates import (CartesianRepresentation, SkyCoord,
                                 UnitSphericalRepresentation)
from astropy import units as u
import healpy as hp
import numpy as np

__all__ = ('find_greedy_credible_levels', 'interp_greedy_credible_levels',
           'smooth_ud_grade', 'posterior_mean', 'posterior_max')


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
    npix = ah.nside_to_npix(nside)
    theta, phi = hp.pix2ang(nside, np.arange(npix), nest=nest)
    new_m = hp.get_interp_val(m, theta, phi, nest=nest)
    return new_m * len(m) / len(new_m)


def interp_greedy_credible_levels(x, xp, fp, right=None):
    """Perform linear interpolation suitable for finding credible levels.

    The linear interpolation is performed with the boundary condition that
    :math:`f(x) = 0`.

    Examples
    --------
    >>> xp = [1, 2, 3, 4, 5]
    >>> fp = [0.2, 0.4, 0.6, 0.8, 1.0]
    >>> interp_greedy_credible_levels([0, 0.5, 1.0, 1.5, 2.0], xp, fp)
    array([0. , 0.1, 0.2, 0.3, 0.4])
    """
    return np.interp(x, np.pad(xp, (1, 0)), np.pad(fp, (1, 0)), right=right)


def posterior_mean(prob, nest=False):
    npix = len(prob)
    nside = ah.npix_to_nside(npix)
    xyz = hp.pix2vec(nside, np.arange(npix), nest=nest)
    mean_xyz = np.average(xyz, axis=1, weights=prob)
    pos = SkyCoord(*mean_xyz, representation_type=CartesianRepresentation)
    pos.representation_type = UnitSphericalRepresentation
    return pos


def posterior_max(prob, nest=False):
    npix = len(prob)
    nside = ah.npix_to_nside(npix)
    i = np.argmax(prob)
    return SkyCoord(
        *hp.pix2ang(nside, i, nest=nest, lonlat=True), unit=u.deg)
