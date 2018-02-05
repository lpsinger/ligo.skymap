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
Angle utilities
"""
import numpy as np

__all__ = ('reference_angle', 'reference_angle_deg',
           'wrapped_angle', 'wrapped_angle_deg')


def reference_angle(a):
    """Convert an angle to a reference angle between -pi and pi."""
    a = np.mod(a, 2 * np.pi)
    return np.where(a <= np.pi, a, a - 2 * np.pi)


def reference_angle_deg(a):
    """Convert an angle to a reference angle between -180 and 180 degrees."""
    a = np.mod(a, 360)
    return np.where(a <= 180, a, a - 360)


def wrapped_angle(a):
    """Convert an angle to a reference angle between 0 and 2*pi."""
    return np.mod(a, 2 * np.pi)


def wrapped_angle_deg(a):
    """Convert an angle to a reference angle between 0 and 2*pi."""
    return np.mod(a, 360)
