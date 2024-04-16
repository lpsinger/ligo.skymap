#
# Copyright (C) 2014-2018  Leo Singer
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
"""RGB data for the "Cylon red" color map.

A print- and screen-friendly color map designed specifically for plotting
LSC/Virgo/KAGRA sky maps. The color map is constructed in CIE Lab space,
following a linear ramp in lightness (the `l` coordinate) and a cubic spline in
color components (the `a` and `b` coordinates).

This particular color map was selected from 20 random realizations of this
construction."""

if __name__ == '__main__':  # pragma: no cover
    from colormath.color_conversions import convert_color
    from colormath.color_objects import LabColor, sRGBColor
    from scipy.interpolate import interp1d
    import numpy as np

    def lab_to_rgb(*args):
        """Convert Lab color to sRGB, with components clipped to (0, 1)."""
        Lab = LabColor(*args)
        sRGB = convert_color(Lab, sRGBColor)
        return np.clip(sRGB.get_value_tuple(), 0, 1)

    L_samples = np.linspace(100, 0, 5)

    a_samples = (
        33.34664938,
        98.09940562,
        84.48361516,
        76.62970841,
        21.43276891)

    b_samples = (
        62.73345997,
        2.09003022,
        37.28252236,
        76.22507582,
        16.24862535)

    L = np.linspace(100, 0, 255)
    a = interp1d(L_samples, a_samples[::-1], 'cubic')(L)
    b = interp1d(L_samples, b_samples[::-1], 'cubic')(L)

    for line in __doc__.splitlines():
        print('#', line)
    for L, a, b in zip(L, a, b):
        print(*lab_to_rgb(L, a, b), sep=',')
