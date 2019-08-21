# -*- coding: utf-8 -*-
#
# Copyright (C) 2013-2019  Leo Singer
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
Injection finding for HEALPix sky maps.

This module is deprecated. It is provided forward backwards compatibility. Use
:mod:`ligo.skymap.postprocess.crossmatch` instead.
"""

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.utils import deprecated

from . import crossmatch, CrossmatchResult


__all__ = ('find_injection_moc', 'FoundInjection')


@deprecated(
    '0.1.10',
    alternative='ligo.skymap.postprocess.crossmatch.CrossmatchResult')
class FoundInjection(CrossmatchResult):
    pass


@deprecated(
    '0.1.10', alternative='ligo.skymap.postprocess.crossmatch.crossmatch')
def find_injection_moc(sky_map, true_ra=None, true_dec=None, true_dist=None,
                       contours=(), areas=(), modes=False, cosmology=False):
    """
    Given a sky map and the true right ascension and declination (in radians),
    find the smallest area in deg^2 that would have to be searched to find the
    source, the smallest posterior mass, and the angular offset in degrees from
    the true location to the maximum (mode) of the posterior. Optionally, also
    compute the areas of and numbers of modes within the smallest contours
    containing a given total probability.
    """
    if true_ra is None and true_dec is None:
        coordinates = None
    elif true_ra is None or true_dec is None:
        raise ValueError('Both true_ra and true_dec must be provided or None')
    elif true_dist is None:
        coordinates = SkyCoord(true_ra * u.rad, true_dec * u.rad)
    else:
        coordinates = SkyCoord(true_ra * u.rad, true_dec * u.rad,
                               true_dist * u.Mpc)
    return crossmatch(sky_map, coordinates, contours=contours, areas=areas,
                      modes=modes, cosmology=cosmology)
