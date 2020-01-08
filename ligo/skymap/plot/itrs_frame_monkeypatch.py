#
# Copyright (C) 2020  Leo Singer
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
"""Workaround for https://github.com/astropy/astropy/pull/9609."""
from astropy.coordinates import ITRS, SphericalRepresentation
from astropy.wcs.utils import _wcs_to_celestial_frame_builtin
from astropy.wcs.utils import WCS_FRAME_MAPPINGS


def wcs_to_celestial_frame(*args, **kwargs):
    frame = _wcs_to_celestial_frame_builtin(*args, **kwargs)
    if isinstance(frame, ITRS):
        frame = ITRS(obstime=frame.obstime,
                     representation_type=SphericalRepresentation)
    return frame


def install():
    WCS_FRAME_MAPPINGS[0] = [wcs_to_celestial_frame]
