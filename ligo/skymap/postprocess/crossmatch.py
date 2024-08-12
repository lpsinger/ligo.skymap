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
"""Catalog cross matching for HEALPix sky maps."""
from collections import namedtuple

import astropy_healpix as ah
from astropy.coordinates import ICRS, SkyCoord, SphericalRepresentation
from astropy import units as u
import healpy as hp
import numpy as np

from .. import distance
from .. import moc

from .cosmology import dVC_dVL_for_DL

__all__ = ('crossmatch', 'CrossmatchResult')


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

    WARNING: The input array is clobbered in the process.
    """
    npix = len(m)
    nside = ah.npix_to_nside(npix)
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
    unit sphere.
    """
    cos_angle_distance = (
        np.cos(phi1 - phi0) * np.sin(theta0) * np.sin(theta1) +
        np.cos(theta0) * np.cos(theta1))
    return np.clip(cos_angle_distance, -1, 1)


def angle_distance(theta0, phi0, theta1, phi1):
    """Angular separation in radians between two points on the unit sphere."""
    return np.arccos(cos_angle_distance(theta0, phi0, theta1, phi1))


# Class to hold return value of find_injection method
CrossmatchResult = namedtuple(
    'CrossmatchResult',
    'searched_area searched_prob offset searched_modes contour_areas '
    'area_probs contour_modes searched_prob_dist contour_dists '
    'searched_vol searched_prob_vol contour_vols probdensity probdensity_vol')
"""Cross match result as returned by
:func:`~ligo.skymap.postprocess.crossmatch.crossmatch`.

Notes
-----
 - All probabilities returned are between 0 and 1.
 - All angles returned are in degrees.
 - All areas returned are in square degrees.
 - All distances are luminosity distances in units of Mpc.
 - All volumes are in units of Mpc³. If :func:`.crossmatch` was run with
   ``cosmology=False``, then all volumes are Euclidean volumes in luminosity
   distance. If :func:`.crossmatch` was run with ``cosmology=True``, then all
   volumes are comoving volumes.

"""
_same_length_as_coordinates = ''' \
Same length as the `coordinates` argument passed to \
:func:`~ligo.skymap.postprocess.crossmatch.crossmatch`.'''
_same_length_as_contours = ''' \
of the probabilities specified by the `contour` argument passed to \
:func:`~ligo.skymap.postprocess.crossmatch.crossmatch`.'''
_same_length_as_areas = ''' \
of the areas specified by the `areas` argument passed to
:func:`~ligo.skymap.postprocess.crossmatch.crossmatch`.'''
CrossmatchResult.searched_area.__doc__ = '''\
Area within the 2D credible region containing each target \
position.''' + _same_length_as_coordinates
CrossmatchResult.searched_prob.__doc__ = '''\
Probability within the 2D credible region containing each target \
position.''' + _same_length_as_coordinates
CrossmatchResult.offset.__doc__ = '''\
Angles on the sky between the target positions and the maximum a posteriori \
position.''' + _same_length_as_coordinates
CrossmatchResult.searched_modes.__doc__ = '''\
Number of disconnected regions within the 2D credible regions \
containing each target position.''' + _same_length_as_coordinates
CrossmatchResult.contour_areas.__doc__ = '''\
Area within the 2D credible regions''' + _same_length_as_contours
CrossmatchResult.area_probs.__doc__ = '''\
Probability within the 2D credible regions''' + _same_length_as_areas
CrossmatchResult.contour_modes.__doc__ = '''\
Number of disconnected regions within the 2D credible \
regions''' + _same_length_as_contours
CrossmatchResult.searched_prob_dist.__doc__ = '''\
Cumulative CDF of distance, marginalized over sky position, at the distance \
of each of the targets.''' + _same_length_as_coordinates
CrossmatchResult.contour_dists.__doc__ = '''\
Distance credible interval, marginalized over sky \
position,''' + _same_length_as_coordinates
CrossmatchResult.searched_vol.__doc__ = '''\
Volume within the 3D credible region containing each target \
position.''' + _same_length_as_coordinates
CrossmatchResult.searched_prob_vol.__doc__ = '''\
Probability within the 3D credible region containing each target \
position.''' + _same_length_as_coordinates
CrossmatchResult.contour_vols.__doc__ = '''\
Volume within the 3D credible regions''' + _same_length_as_contours
CrossmatchResult.probdensity.__doc__ = '''\
2D probability density per steradian at the positions of each of the \
targets.''' + _same_length_as_coordinates
CrossmatchResult.probdensity_vol.__doc__ = '''\
3D probability density per cubic megaparsec at the positions of each of the \
targets.''' + _same_length_as_coordinates


def crossmatch(sky_map, coordinates=None,
               contours=(), areas=(), modes=False, cosmology=False):
    """Cross match a sky map with a catalog of points.

    Given a sky map and the true right ascension and declination (in radians),
    find the smallest area in deg^2 that would have to be searched to find the
    source, the smallest posterior mass, and the angular offset in degrees from
    the true location to the maximum (mode) of the posterior. Optionally, also
    compute the areas of and numbers of modes within the smallest contours
    containing a given total probability.

    Parameters
    ----------
    sky_map : :class:`astropy.table.Table`
        A multiresolution sky map, as returned by
        :func:`ligo.skymap.io.fits.read_sky_map` called with the keyword
        argument ``moc=True``.

    coordinates : :class:`astropy.coordinates.SkyCoord`, optional
        The catalog of target positions to match against.

    contours : :class:`tuple`, optional
        Credible levels between 0 and 1. If this argument is present, then
        calculate the areas and volumes of the 2D and 3D credible regions that
        contain these probabilities. For example, for ``contours=(0.5, 0.9)``,
        then areas and volumes of the 50% and 90% credible regions.

    areas : :class:`tuple`, optional
        Credible areas in square degrees. If this argument is present, then
        calculate the probability contained in the 2D credible levels that have
        these areas. For example, for ``areas=(20, 100)``, then compute the
        probability within the smallest credible levels of 20 deg² and 100
        deg², respectively.

    modes : :class:`bool`, optional
        If True, then enable calculation of the number of distinct modes or
        islands of probability. Note that this option may be computationally
        expensive.

    cosmology : :class:`bool`, optional
        If True, then search space by descending probability density per unit
        comoving volume. If False, then search space by descending probability
        per luminosity distance cubed.

    Returns
    -------
    result : :class:`~ligo.skymap.postprocess.crossmatch.CrossmatchResult`

    Notes
    -----
    This function is also be used for injection finding; see
    :doc:`/tool/ligo_skymap_stats`.

    Examples
    --------
    First, some imports:

    >>> from astroquery.vizier import VizierClass
    >>> from astropy.coordinates import SkyCoord
    >>> from ligo.skymap.io import read_sky_map
    >>> from ligo.skymap.postprocess import crossmatch

    Next, retrieve the GLADE catalog using Astroquery and get the coordinates
    of all its entries:

    >>> vizier = VizierClass(
    ...     row_limit=-1,
    ...     columns=['recno', 'GWGC', '_RAJ2000', '_DEJ2000', 'Dist'])
    >>> cat, = vizier.get_catalogs('VII/281/glade2')
    >>> cat.sort('recno')  # sort catalog so that doctest output is stable
    >>> del cat['recno']
    >>> coordinates = SkyCoord(cat['_RAJ2000'], cat['_DEJ2000'], cat['Dist'])

    Load the multiresolution sky map for S190814bv:

    >>> url = 'https://gracedb.ligo.org/api/superevents/S190814bv/files/bayestar.multiorder.fits'
    >>> skymap = read_sky_map(url, moc=True)

    Perform the cross match:

    >>> result = crossmatch(skymap, coordinates)

    Using the cross match results, we can list the galaxies within the 90%
    credible volume:

    >>> print(cat[result.searched_prob_vol < 0.9])
         _RAJ2000          _DEJ2000        GWGC            Dist
           deg               deg                           Mpc
    ----------------- ----------------- ---------- --------------------
      9.3396700000000 -19.9342460000000    NGC0171    57.56212553960000
     20.2009090000000 -31.1146050000000        ---   137.16022925600001
      8.9144680000000 -20.1252980000000 ESO540-003    49.07809291930000
     10.6762720000000 -21.7740820000000        ---   276.46938505499998
     13.5855170000000 -23.5523850000000        ---   138.44550704800000
     20.6362970000000 -29.9825150000000        ---   160.23313164900000
                  ...               ...        ...                  ...
     10.6939000000000 -25.6778300000000        ---   323.59399999999999
     15.4935000000000 -26.0305000000000        ---   304.78899999999999
     15.2794000000000 -27.0411000000000        ---   320.62700000000001
     14.8324000000000 -27.0460000000000        ---   320.62700000000001
     14.5341000000000 -26.0949000000000        ---   307.61000000000001
     23.1281000000000 -31.1109200000000        ---   320.62700000000001
    Length = 1479 rows

    """  # noqa: E501
    # Astropy coordinates that are constructed without distance have
    # a distance field that is unity (dimensionless).
    if coordinates is None:
        true_ra = true_dec = true_dist = None
    else:
        # Ensure that coordinates are in proper frame and representation
        coordinates = SkyCoord(coordinates,
                               representation_type=SphericalRepresentation,
                               frame=ICRS)
        true_ra = coordinates.ra.rad
        true_dec = coordinates.dec.rad
        if np.any(coordinates.distance != 1):
            true_dist = coordinates.distance.to_value(u.Mpc)
        else:
            true_dist = None

    contours = np.asarray(contours)

    # Sort the pixels by descending posterior probability.
    sky_map = np.flipud(np.sort(sky_map, order='PROBDENSITY'))

    # Find the pixel that contains the injection.
    order, ipix = moc.uniq2nest(sky_map['UNIQ'])
    max_order = np.max(order)
    max_nside = ah.level_to_nside(max_order)
    max_ipix = ipix << np.int64(2 * (max_order - order))
    if true_ra is not None:
        true_theta = 0.5 * np.pi - true_dec
        true_phi = true_ra
        true_pix = hp.ang2pix(max_nside, true_theta, true_phi, nest=True)
        i = np.argsort(max_ipix)
        true_idx = i[np.digitize(true_pix, max_ipix[i]) - 1]

    # Find the angular offset between the mode and true locations.
    mode_theta, mode_phi = hp.pix2ang(
        ah.level_to_nside(order[0]), ipix[0], nest=True)
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

    if true_ra is None:
        searched_area = searched_prob = probdensity = np.nan
    else:
        # Find the smallest area that would have to be searched to find
        # the true location.
        searched_area = area[true_idx]

        # Find the smallest posterior mass that would have to be searched to
        # find the true location.
        searched_prob = prob[true_idx]

        # Find the probability density.
        probdensity = sky_map['PROBDENSITY'][true_idx]

    # Find the contours of the given credible levels.
    contour_idxs = np.digitize(contours, prob) - 1

    # For each of the given confidence levels, compute the area of the
    # smallest region containing that probability.
    contour_areas = np.interp(
        contours, prob, area, left=0, right=4*180**2/np.pi).tolist()

    # For each listed area, find the probability contained within the
    # smallest credible region of that area.
    area_probs = np.interp(areas, area, prob, left=0, right=1).tolist()

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
        dP_dA = sky_map['PROBDENSITY']
        mu = sky_map['DISTMU']
        sigma = sky_map['DISTSIGMA']
        norm = sky_map['DISTNORM']

        # Set up distance grid.
        n_r = 1000
        distmean, _ = distance.parameters_to_marginal_moments(dP, mu, sigma)
        max_r = 6 * distmean
        if true_dist is not None and np.size(true_dist) != 0 \
                and np.max(true_dist) > max_r:
            max_r = np.max(true_dist)
        d_r = max_r / n_r

        # Calculate searched_prob_dist and contour_dists.
        r = d_r * np.arange(1, n_r)
        P_r = distance.marginal_cdf(r, dP, mu, sigma, norm)
        if true_dist is None:
            searched_prob_dist = np.nan
        else:
            searched_prob_dist = np.interp(true_dist, r, P_r, left=0, right=1)
        if len(contours) == 0:
            contour_dists = []
        else:
            lo, hi = np.interp(
                np.vstack((
                    0.5 * (1 - contours),
                    0.5 * (1 + contours)
                )), P_r, r, left=0, right=np.inf)
            contour_dists = (hi - lo).tolist()

        # Calculate volume of each voxel, defined as the region within the
        # HEALPix pixel and contained within the two centric spherical shells
        # with radii (r - d_r / 2) and (r + d_r / 2).
        dV = (np.square(r) + np.square(d_r) / 12) * d_r * dA.reshape(-1, 1)

        # Calculate probability within each voxel.
        dP = np.exp(
            -0.5 * np.square(
                (r.reshape(1, -1) - mu.reshape(-1, 1)) / sigma.reshape(-1, 1)
            )
        ) * (dP_dA * norm / (sigma * np.sqrt(2 * np.pi))).reshape(-1, 1) * dV
        dP[np.isnan(dP)] = 0  # Suppress invalid values

        # Calculate probability density per unit volume.

        if cosmology:
            dV *= dVC_dVL_for_DL(r)
        dP_dV = dP / dV
        i = np.flipud(np.argsort(dP_dV.ravel()))

        P_flat = np.cumsum(dP.ravel()[i])
        V_flat = np.cumsum(dV.ravel()[i])

        contour_vols = np.interp(
            contours, P_flat, V_flat, left=0, right=np.inf).tolist()
        P = np.empty_like(P_flat)
        V = np.empty_like(V_flat)
        P[i] = P_flat
        V[i] = V_flat
        P = P.reshape(dP.shape)
        V = V.reshape(dV.shape)
        if true_dist is None:
            searched_vol = searched_prob_vol = probdensity_vol = np.nan
        else:
            i_radec = true_idx
            i_dist = np.digitize(true_dist, r) - 1
            probdensity_vol = dP_dV[i_radec, i_dist]
            searched_prob_vol = P[i_radec, i_dist]
            searched_vol = V[i_radec, i_dist]
    else:
        searched_vol = searched_prob_vol = searched_prob_dist \
            = probdensity_vol = np.nan
        contour_dists = [np.nan] * len(contours)
        contour_vols = [np.nan] * len(contours)

    # Done.
    return CrossmatchResult(
        searched_area, searched_prob, offset, searched_modes, contour_areas,
        area_probs, contour_modes, searched_prob_dist, contour_dists,
        searched_vol, searched_prob_vol, contour_vols, probdensity,
        probdensity_vol)
