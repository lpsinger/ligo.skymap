# -*- coding: UTF-8 -*-
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
Sub-sample interpolation for matched filter time series.
"""

import numpy as np
from scipy import optimize

__all__ = ('interpolate_max',)


#
# Lanczos interpolation
#


def lanczos(t, a):
    """The Lanczos kernel."""
    return np.sinc(t) * np.sinc(t / a)


def lanczos_interpolant(t, y):
    """An interpolant constructed by convolution of the Lanczos kernel with
    a set of discrete samples at unit intervals."""
    a = len(y) // 2
    return sum(lanczos(t - i + a, a) * yi for i, yi in enumerate(y))


def lanczos_interpolant_utility_func(t, y):
    """Utility function for Lanczos interpolation."""
    return -abs2(lanczos_interpolant(t, y))


def interpolate_max_lanczos(imax, y, window_length):
    """Find the time and maximum absolute value of a time series by Lanczos
    interpolation."""
    yi = y[imax-window_length:imax+window_length+1]
    tmax = optimize.fminbound(
        lanczos_interpolant_utility_func, -1., 1., (yi,), xtol=1e-5)
    tmax = np.asscalar(tmax)
    ymax = np.asscalar(lanczos_interpolant(tmax, yi))
    return imax + tmax, ymax


#
# Catmull-Rom spline interpolation
#


def poly_catmull_rom(y):
    return np.poly1d([
        -0.5 * y[0] + 1.5 * y[1] - 1.5 * y[2] + 0.5 * y[3],
        y[0] - 2.5 * y[1] + 2 * y[2] - 0.5 * y[3],
        -0.5 * y[0] + 0.5 * y[2],
        y[1]
    ])


def interpolate_max_catmull_rom_even(y):

    # Construct Catmull-Rom interpolating polynomials for
    # real and imaginary parts
    poly_re = poly_catmull_rom(y.real)
    poly_im = poly_catmull_rom(y.imag)

    # Find the roots of d(|y|^2)/dt as approximated
    roots = (poly_re * poly_re.deriv() + poly_im * poly_im.deriv()).r

    # Find which of the two matched interior points has a greater magnitude
    t_max = 0.
    y_max = y[1]
    y_max_abs2 = abs2(y_max)

    new_t_max = 1.
    new_y_max = y[2]
    new_y_max_abs2 = abs2(new_y_max)

    if new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs2 = new_y_max_abs2

    # Find any real root in (0, 1) that has a magnitude greater than the
    # greatest endpoint
    for root in roots:
        if np.isreal(root) and 0 < root < 1:
            new_t_max = root
            new_y_max = poly_re(new_t_max) + poly_im(new_t_max) * 1j
            new_y_max_abs2 = abs2(new_y_max)
            if new_y_max_abs2 > y_max_abs2:
                t_max = new_t_max
                y_max = new_y_max
                y_max_abs2 = new_y_max_abs2

    # Done
    return t_max, y_max


def interpolate_max_catmull_rom(imax, y, window_length):
    t_max, y_max = interpolate_max_catmull_rom_even(y[imax - 2:imax + 2])
    y_max_abs2 = abs2(y_max)
    t_max = t_max - 1

    new_t_max, new_y_max = interpolate_max_catmull_rom_even(
        y[imax - 1:imax + 3])
    new_y_max_abs2 = abs2(new_y_max)

    if new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs2 = new_y_max_abs2

    return imax + t_max, y_max


#
# Quadratic fit
#


def interpolate_max_quadratic_fit(imax, y, window_length):
    """Quadratic fit to absolute value of y. Note that this one does not alter
    the value at the maximum."""

    poly = np.polyfit(
        np.arange(-window_length, window_length + 1.),
        np.abs(y[imax - window_length:imax + window_length + 1]),
        2)

    # Find which of the two matched interior points has a greater magnitude
    t_max = -1.
    y_max = y[imax - 1]
    y_max_abs2 = abs2(y_max)

    new_t_max = 1.
    new_y_max = y[imax + 1]
    new_y_max_abs2 = abs2(new_y_max)

    if new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs2 = new_y_max_abs2

    # Determine if the global extremum of the polynomial is a
    # local maximum in (-1, 1)
    A, B, C = poly
    new_t_max = -0.5 * B / A
    new_y_max_abs2 = np.square(np.polyval(poly, new_t_max))
    if -1 < new_t_max < 1 and new_y_max_abs2 > y_max_abs2:
        t_max = new_t_max

    return imax + t_max, y[imax]


#
# Nearest neighbor interpolation
#


def interpolate_max_nearest_neighbor(imax, y, window_length):
    """Trivial, nearest-neighbor interpolation"""
    return imax, y[imax]


#
# Set default interpolation scheme
#


_interpolants = {
    'catmull-rom': interpolate_max_catmull_rom,
    'lanczos': interpolate_max_lanczos,
    'nearest-neighbor': interpolate_max_nearest_neighbor,
    'quadratic-fit': interpolate_max_quadratic_fit}


def interpolate_max(imax, y, window_length, method='catmull-rom'):
    """Perform sub-sample interpolation to find the phase and amplitude
    at the maximum of the absolute value of a complex series.

    Parameters
    ----------
    imax : int
        The index of the maximum sample in the series.
    y : `numpy.ndarray`
        The complex series.
    window_length : int
        The window of the interpolation function. The interpolation will
        consider a sliding window of `2 * window_length + 1` samples centered
        on `imax`.
    method : {'catmull-rom', 'lanczos', 'nearest-neighbor', 'quadratic-fit'}
        The interpolation method:
        * `catmull-rom`: Catmull-Rom cubic splines
        * `lanczos`: Lanczos filter interpolation
        * `nearest-neighbor`: Nearest neighbor (e.g., no interpolation)
        * `quadratic-fit`: Fit the absolute value of the SNR to a quadratic
          function.

    Returns
    -------
    imax_interp : float
        The interpolated index of the maximum sample, which should be between
        `imax - 0.5` and `imax + 0.5`.
    ymax_interp : complex
        The interpolated value at the maximum.
    """
    return _interpolants[method](imax, y, window_length)
