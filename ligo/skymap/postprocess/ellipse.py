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

from astropy.wcs import WCS
import healpy as hp
import numpy as np

from .. import moc

__all__ = ('find_ellipse',)


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
    c = np.cov(xy, aweights=prob, rowvar=False)

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
