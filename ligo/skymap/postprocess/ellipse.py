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

from astropy.wcs import WCS
import healpy as hp
import numpy as np

from .. import moc
from ..extern.quantile import quantile

__all__ = ('find_ellipse',)


def find_ellipse(prob, cl=90, projection='ARC', nest=False):
    """For a HEALPix map, find an ellipse that contains a given probability.

    The orientation is defined as the angle of the semimajor axis
    counterclockwise from west on the plane of the sky. If you think of the
    semimajor distance as the width of the ellipse, then the orientation is the
    clockwise rotation relative to the image x-axis. Equivalently, the
    orientation is the position angle of the semi-minor axis.

    These conventions match the definitions used in DS9 region files [1]_ and
    Aladin drawing commands [2]_.

    Parameters
    ----------
    prob : np.ndarray, astropy.table.Table
        The HEALPix probability map, either as a full rank explicit array
        or as a multi-order map.
    cl : float
        The desired credible level (default: 90).
    projection : str, optional
        The WCS projection (default: 'ARC', or zenithal equidistant).
        For a list of possible values, see the Astropy documentation [3]_.
    nest : bool
        HEALPix pixel ordering (default: False, or ring ordering).

    Returns
    -------
    ra : float
        The ellipse center right ascension in degrees.
    dec : float
        The ellipse center right ascension in degrees.
    a : float
        The lenth of the semimajor axis in degrees.
    b : float
        The length of the semiminor axis in degrees.
    pa : float
        The orientation of the ellipse axis on the plane of the sky in degrees.
    area : float
        The area of the ellipse in square degrees.

    Notes
    -----

    The center of the ellipse is the median a posteriori sky position. The
    length and orientation of the semi-major and semi-minor axes are measured
    as follows:

    1. The sky map is transformed to a WCS projection that may be specified by
       the caller. The default projection is ``ARC`` (zenithal equidistant), in
       which radial distances are proportional to the physical angular
       separation from the center point.
    2. A 1-sigma ellipse is estimated by calculating the covariance matrix in
       the projected image plane using three rounds of sigma clipping to reject
       distant outlier points.
    3. The 1-sigma ellipse is inflated until it encloses an integrated
       probability of ``cl`` (default: 90%).

    The function returns a tuple of the right ascension, declination,
    semi-major distance, semi-minor distance, and orientation angle, all in
    degrees.

    References
    ----------

    .. [1] http://ds9.si.edu/doc/ref/region.html
    .. [2] http://aladin.u-strasbg.fr/java/AladinScriptManual.gml#draw
    .. [3] http://docs.astropy.org/en/stable/wcs/index.html#supported-projections

    Examples
    --------

    **Example 1**

    First, we need some imports.

    >>> from astropy.io import fits
    >>> from astropy.utils.data import download_file
    >>> from astropy.wcs import WCS
    >>> import healpy as hp
    >>> from reproject import reproject_from_healpix
    >>> import subprocess

    Next, we download the BAYESTAR sky map for GW170817 from the
    LIGO Document Control Center.

    >>> url = 'https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz'  # doctest: +SKIP
    >>> filename = download_file(url, cache=True, show_progress=False)  # doctest: +SKIP
    >>> _, healpix_hdu = fits.open(filename)  # doctest: +SKIP
    >>> prob = hp.read_map(healpix_hdu, verbose=False)  # doctest: +SKIP

    Then, we calculate ellipse and write it to a DS9 region file.

    >>> ra, dec, a, b, pa, area = find_ellipse(prob)  # doctest: +SKIP
    >>> print(*np.around([ra, dec, a, b, pa, area], 5))  # doctest: +SKIP
    195.03732 -19.29358 8.66545 1.1793 63.61698 32.07665
    >>> s = 'fk5;ellipse({},{},{},{},{})'.format(ra, dec, a, b, pa)  # doctest: +SKIP
    >>> open('ds9.reg', 'w').write(s)  # doctest: +SKIP

    Then, we reproject a small patch of the HEALPix map, and save it to a file.

    >>> wcs = WCS()  # doctest: +SKIP
    >>> wcs.wcs.ctype = ['RA---ARC', 'DEC--ARC']  # doctest: +SKIP
    >>> wcs.wcs.crval = [ra, dec]  # doctest: +SKIP
    >>> wcs.wcs.crpix = [128, 128]  # doctest: +SKIP
    >>> wcs.wcs.cdelt = [-0.1, 0.1]  # doctest: +SKIP
    >>> img, _ = reproject_from_healpix(healpix_hdu, wcs, [256, 256])  # doctest: +SKIP
    >>> img_hdu = fits.ImageHDU(img, wcs.to_header())  # doctest: +SKIP
    >>> img_hdu.writeto('skymap.fits')  # doctest: +SKIP

    Now open the image and region file in DS9. You should find that the ellipse
    encloses the probability hot spot. You can load the sky map and region file
    from the command line:

    .. code-block:: sh

        $ ds9 skymap.fits -region ds9.reg

    Or you can do this manually:

        1. Open DS9.
        2. Open the sky map: select "File->Open..." and choose ``skymap.fits``
           from the dialog box.
        3. Open the region file: select "Regions->Load Regions..." and choose
           ``ds9.reg`` from the dialog box.

    Now open the image and region file in Aladin.

        1. Open Aladin.
        2. Open the sky map: select "File->Load Local File..." and choose
           ``skymap.fits`` from the dialog box.
        3. Open the sky map: select "File->Load Local File..." and choose
           ``ds9.reg`` from the dialog box.

    You can also compare the original HEALPix file with the ellipse in Aladin:

        1. Open Aladin.
        2. Open the HEALPix file by pasting the URL from the top of this
           example in the Command field at the top of the window and hitting
           return, or by selecting "File->Load Direct URL...", pasting the URL,
           and clicking "Submit."
        3. Open the sky map: select "File->Load Local File..." and choose
           ``ds9.reg`` from the dialog box.

    **Example 2**

    This example shows that we get approximately the same answer for GW171087
    if we read it in as a multi-order map.

    >>> from ..io import read_sky_map  # doctest: +SKIP
    >>> skymap_moc = read_sky_map(healpix_hdu, moc=True)  # doctest: +SKIP
    >>> ellipse = find_ellipse(skymap_moc)  # doctest: +SKIP
    >>> print(*np.around(ellipse, 5))  # doctest: +SKIP
    195.03709 -19.27589 8.67611 1.18167 63.60454 32.08015

    **Example 3**

    I'm not showing the `ra` or `pa` output from the examples below because
    the right ascension is arbitary when dec=90° and the position angle is
    arbitrary when a=b; their arbitrary values may vary depending on your math
    library. Also, I add 0.0 to the outputs because on some platforms you tend
    to get values of dec or pa that get rounded to -0.0, which is within
    numerical precision but would break the doctests (see
    https://stackoverflow.com/questions/11010683).

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
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, area)
    89.90863 0.87034 0.87034 2.37888

    >>> prob = make_uniform_in_sin_theta(10)
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, area)
    89.90828 9.02485 9.02484 255.11972

    >>> prob = make_uniform_in_sin_theta(120)
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, area)
    90.0 107.9745 107.9745 26988.70467

    **Example 4**

    These are approximately Gaussian distributions.

    >>> from scipy import stats
    >>> def make_gaussian(mean, cov, nside=512):
    ...     npix = hp.nside2npix(nside)
    ...     xyz = np.transpose(hp.pix2vec(nside, np.arange(npix)))
    ...     dist = stats.multivariate_normal(mean, cov)
    ...     prob = dist.pdf(xyz)
    ...     return prob / prob.sum()
    ...

    This one is centered at RA=45°, Dec=0° and has a standard deviation of ~1°.

    >>> prob = make_gaussian(
    ...     [1/np.sqrt(2), 1/np.sqrt(2), 0],
    ...     np.square(np.deg2rad(1)))
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(ra, dec, a, b, area)
    45.0 0.0 2.14241 2.14208 14.4677

    This one is centered at RA=45°, Dec=0°, and is elongated in the north-south
    direction.

    >>> prob = make_gaussian(
    ...     [1/np.sqrt(2), 1/np.sqrt(2), 0],
    ...     np.diag(np.square(np.deg2rad([1, 1, 10]))))
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(ra, dec, a, b, pa, area)
    45.0 0.0 13.58769 2.08298 90.0 88.57797

    This one is centered at RA=0°, Dec=0°, and is elongated in the east-west
    direction.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     np.diag(np.square(np.deg2rad([1, 10, 1]))))
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, pa, area)
    0.0 13.58392 2.08238 0.0 88.54623

    This one is centered at RA=0°, Dec=0°, and has its long axis tilted about
    10° to the west of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 0.1, -0.15],
    ...      [0, -0.15, 1]])
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, pa, area)
    0.0 64.77133 33.50754 80.78231 6372.34466

    This one is centered at RA=0°, Dec=0°, and has its long axis tilted about
    10° to the east of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 0.1, 0.15],
    ...      [0, 0.15, 1]])
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, pa, area)
    0.0 64.77133 33.50754 99.21769 6372.34466

    This one is centered at RA=0°, Dec=0°, and has its long axis tilted about
    80° to the east of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 1, 0.15],
    ...      [0, 0.15, 0.1]])
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, pa, area)
    0.0 64.77564 33.50986 170.78252 6372.42573

    This one is centered at RA=0°, Dec=0°, and has its long axis tilted about
    80° to the west of north.

    >>> prob = make_gaussian(
    ...     [1, 0, 0],
    ...     [[0.1, 0, 0],
    ...      [0, 1, -0.15],
    ...      [0, -0.15, 0.1]])
    ...
    >>> ra, dec, a, b, pa, area = np.around(find_ellipse(prob), 5) + 0
    >>> print(dec, a, b, pa, area)
    0.0 64.77564 33.50986 9.21748 6372.42573
    """  # noqa: E501
    try:
        prob['UNIQ']
    except (IndexError, KeyError, ValueError):
        npix = len(prob)
        nside = hp.npix2nside(npix)
        ipix = range(npix)
        area = hp.nside2pixarea(nside, degrees=True)
    else:
        order, ipix = moc.uniq2nest(prob['UNIQ'])
        nside = 1 << order.astype(int)
        ipix = ipix.astype(int)
        area = hp.nside2pixarea(nside)
        prob = prob['PROBDENSITY'] * area
        area *= np.square(180 / np.pi)
        nest = True

    # Find median a posteriori sky position.
    xyz0 = [quantile(x, 0.5, weights=prob)
            for x in hp.pix2vec(nside, ipix, nest=nest)]
    (ra,), (dec,) = hp.vec2ang(np.asarray(xyz0), lonlat=True)

    # Construct WCS with the specified projection
    # and centered on mean direction.
    w = WCS()
    w.wcs.crval = [ra, dec]
    w.wcs.ctype = ['RA---' + projection, 'DEC--' + projection]

    # Transform HEALPix to zenithal equidistant coordinates.
    xy = w.wcs_world2pix(
        np.transpose(
            hp.pix2ang(
                nside, ipix, nest=nest, lonlat=True)), 1)

    # Keep only values that were inside the projection.
    keep = np.logical_and.reduce(np.isfinite(xy), axis=1)
    xy = xy[keep]
    prob = prob[keep]
    if not np.isscalar(area):
        area = area[keep]

    # Find covariance matrix, performing three rounds of sigma-clipping
    # to reject outliers.
    keep = np.ones(len(xy), dtype=bool)
    for _ in range(3):
        c = np.cov(xy[keep], aweights=prob[keep], rowvar=False)
        nsigmas = np.sqrt(np.sum(xy.T * np.linalg.solve(c, xy.T), axis=0))
        keep &= (nsigmas < 3)

    # Find the number of sigma that enclose the cl% credible level.
    i = np.argsort(nsigmas)
    nsigmas = nsigmas[i]
    cls = np.cumsum(prob[i])
    if np.isscalar(area):
        careas = np.arange(1, len(i) + 1) * area
    else:
        careas = np.cumsum(area[i])
    nsigma = np.interp(1e-2 * cl, cls, nsigmas)
    area = np.interp(1e-2 * cl, cls, careas)

    # If the credible level is not within the projection,
    # then stop here and return all nans.
    if 1e-2 * cl > cls[-1]:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    # Find the eigendecomposition of the covariance matrix.
    w, v = np.linalg.eigh(c)

    # Find the semi-minor and semi-major axes.
    b, a = nsigma * np.sqrt(w)

    # Find the position angle.
    pa = np.rad2deg(np.arctan2(*v[0]))

    # An ellipse is symmetric under rotations of 180°.
    # Return the smallest possible positive position angle.
    pa %= 180

    # Done!
    return ra, dec, a, b, pa, area
