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
"""Tools for transforming sky maps to detector frames. These functions were
written for Fig. 2 of the GW150914 localization and follow-up paper [1]_.

References
----------

.. [1] LSC/Virgo et al., 2016. "Localization and Broadband Follow-up of the
       Gravitational-wave Transient GW150914." ApJL 826, L13
       <https://doi.org/10.3847/2041-8205/826/1/L13>.
"""

from astropy import constants
from astropy.utils import deprecated
import healpy as hp
import numpy as np

__all__ = ('get_detector_pair_axis', 'rotate_map_to_axis', 'polar_profile')


def _get_detector_location(ifo):
    if isinstance(ifo, str):
        try:
            import lalsimulation
        except ImportError:
            raise RuntimeError('Looking up detectors by name '
                               'requires the lalsimulation package.')
        ifo = lalsimulation.DetectorPrefixToLALDetector(ifo)
    try:
        ifo = ifo.location
    except AttributeError:
        pass
    return ifo


@deprecated('0.0.6')
def get_detector_pair_axis(ifo1, ifo2, gmst):
    """Find the sky position where the line between two detectors pierces the
    celestial sphere.

    Parameters
    ----------
    ifo1, ifo2 : str or `~lal.Detector` or `numpy.ndarray`
        The detector positions. Can be specifed by name (e.g. `'H1'`), or by
        `lal.Detector` object (e.g., as returned by
        `lalsimulation.DetectorPrefixToLALDetector('H1')` or by the geocentric
        Cartesian position in meters.

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

    Examples
    --------

    >>> ret = get_detector_pair_axis('H1', 'L1', 41758.384193753656)
    >>> print(*np.around(ret, 5))
    5.94546 -0.47622 0.01001
    """
    ifo1 = _get_detector_location(ifo1)
    ifo2 = _get_detector_location(ifo2)

    n = ifo2 - ifo1
    light_travel_time = np.sqrt(np.sum(np.square(n))) / constants.c.value
    (theta,), (phi,) = hp.vec2ang(n)
    pole_ra = (gmst + phi) % (2 * np.pi)
    pole_dec = 0.5 * np.pi - theta
    return pole_ra, pole_dec, light_travel_time


@deprecated('0.0.6')
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


@deprecated('0.0.6')
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
