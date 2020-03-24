#
# Copyright (C) 2013-2020  Leo Singer
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
"""Sub-sample interpolation for matched filter time series.

Example
-------
.. plot::
   :context: reset
   :include-source:
   :align: center

    from ligo.skymap.bayestar.interpolation import interpolate_max
    from matplotlib import pyplot as plt
    import numpy as np

    z = np.asarray([ 9.135017 -2.8185585j,  9.995214 -1.1222992j,
                    10.682851 +0.8188147j, 10.645139 +3.0268786j,
                     9.713133 +5.5589147j,  7.9043484+7.9039335j,
                     5.511646 +9.333084j ,  2.905198 +9.715742j ,
                     0.5302934+9.544538j ])

    amp = np.abs(z)
    arg = np.rad2deg(np.unwrap(np.angle(z)))
    arg -= (np.median(arg) // 360) * 360
    imax = np.argmax(amp)
    window = 4

    fig, (ax_amp, ax_arg) = plt.subplots(2, 1, figsize=(5, 6), sharex=True)
    ax_arg.set_xlabel('Sample index')
    ax_amp.set_ylabel('Amplitude')
    ax_arg.set_ylabel('Phase')
    args, kwargs = ('.-',), dict(color='lightgray', label='data')
    ax_amp.plot(amp, *args, **kwargs)
    ax_arg.plot(arg, *args, **kwargs)
    for method in ['lanczos', 'catmull-rom',
                   'quadratic-fit', 'nearest-neighbor']:
        i, y = interpolate_max(imax, z, window, method)
        amp = np.abs(y)
        arg = np.rad2deg(np.angle(y))
        args, kwargs = ('o',), dict(mfc='none', label=method)
        ax_amp.plot(i, amp, *args, **kwargs)
        ax_arg.plot(i, arg, *args, **kwargs)
    ax_arg.legend()
    fig.tight_layout()

"""
import numpy as np
from scipy import optimize

from .filter import abs2, exp_i, unwrap

__all__ = ('interpolate_max',)


#
# Lanczos interpolation
#


def lanczos(t, a):
    """The Lanczos kernel."""
    return np.where(np.abs(t) < a, np.sinc(t) * np.sinc(t / a), 0)


def lanczos_interpolant(t, y):
    """An interpolant constructed by convolution of the Lanczos kernel with
    a set of discrete samples at unit intervals.
    """
    a = len(y) // 2
    return sum(lanczos(t - i + a, a) * yi for i, yi in enumerate(y))


def lanczos_interpolant_utility_func(t, y):
    """Utility function for Lanczos interpolation."""
    return -abs2(lanczos_interpolant(t, y))


def interpolate_max_lanczos(imax, y, window_length):
    """Find the time and maximum absolute value of a time series by Lanczos
    interpolation.
    """
    yi = y[(imax - window_length):(imax + window_length + 1)]
    tmax = optimize.fminbound(
        lanczos_interpolant_utility_func, -1., 1., (yi,), xtol=1e-5)
    tmax = tmax.item()
    ymax = lanczos_interpolant(tmax, yi).item()
    return imax + tmax, ymax


#
# Catmull-Rom spline interpolation, real and imaginary parts
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
            new_t_max = np.real(root)
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
# Catmull-Rom spline interpolation, amplitude and phase
#


def interpolate_max_catmull_rom_amp_phase_even(y):

    # Construct Catmull-Rom interpolating polynomials for
    # real and imaginary parts
    poly_abs = poly_catmull_rom(np.abs(y))
    poly_arg = poly_catmull_rom(unwrap(np.angle(y)))

    # Find the roots of d(|y|)/dt as approximated
    roots = poly_abs.r

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
            new_t_max = np.real(root)
            new_y_max = poly_abs(new_t_max) * exp_i(poly_arg(new_t_max))
            new_y_max_abs2 = abs2(new_y_max)
            if new_y_max_abs2 > y_max_abs2:
                t_max = new_t_max
                y_max = new_y_max
                y_max_abs2 = new_y_max_abs2

    # Done
    return t_max, y_max


def interpolate_max_catmull_rom_amp_phase(imax, y, window_length):
    t_max, y_max = interpolate_max_catmull_rom_amp_phase_even(
        y[imax - 2:imax + 2])
    y_max_abs2 = abs2(y_max)
    t_max = t_max - 1

    new_t_max, new_y_max = interpolate_max_catmull_rom_amp_phase_even(
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
    the value at the maximum.
    """
    t = np.arange(-window_length, window_length + 1.)
    y = y[imax - window_length:imax + window_length + 1]
    y_abs = np.abs(y)
    a, b, c = np.polyfit(t, y_abs, 2)

    # Find which of the two matched interior points has a greater magnitude
    t_max = -1.
    y_max = y[window_length - 1]
    y_max_abs = y_abs[window_length - 1]

    new_t_max = 1.
    new_y_max = y[window_length + 1]
    new_y_max_abs = y_abs[window_length + 1]

    if new_y_max_abs > y_max_abs:
        t_max = new_t_max
        y_max = new_y_max
        y_max_abs = new_y_max_abs

    # Determine if the global extremum of the polynomial is a
    # local maximum in (-1, 1)
    new_t_max = -0.5 * b / a
    new_y_max_abs = c - 0.25 * np.square(b) / a
    if -1 < new_t_max < 1 and new_y_max_abs > y_max_abs:
        t_max = new_t_max
        y_max_abs = new_y_max_abs
        y_phase = np.interp(t_max, t, np.unwrap(np.angle(y)))
        y_max = y_max_abs * exp_i(y_phase)

    return imax + t_max, y_max


#
# Nearest neighbor interpolation
#


def interpolate_max_nearest_neighbor(imax, y, window_length):
    """Trivial, nearest-neighbor interpolation."""
    return imax, y[imax]


#
# Set default interpolation scheme
#


_interpolants = {
    'catmull-rom-amp-phase': interpolate_max_catmull_rom_amp_phase,
    'catmull-rom': interpolate_max_catmull_rom,
    'lanczos': interpolate_max_lanczos,
    'nearest-neighbor': interpolate_max_nearest_neighbor,
    'quadratic-fit': interpolate_max_quadratic_fit}


def interpolate_max(imax, y, window_length, method='catmull-rom-amp-phase'):
    """Perform sub-sample interpolation to find the phase and amplitude
    at the maximum of the absolute value of a complex series.

    Parameters
    ----------
    imax : int
        The index of the maximum sample in the series.
    y : `numpy.ndarray`
        The complex series.
    window_length : int
        The window of the interpolation function for the `lanczos` and
        `quadratic-fit` methods. The interpolation will consider a sliding
        window of `2 * window_length + 1` samples centered on `imax`.
    method : {'catmull-rom-amp-phase', 'catmull-rom', 'lanczos', 'nearest-neighbor', 'quadratic-fit'}
        The interpolation method. Can be any of the following:

        * `catmull-rom-amp-phase`: Catmull-Rom cubic splines on amplitude and phase
          The `window_length` parameter is ignored (understood to be 2).
        * `catmull-rom`: Catmull-Rom cubic splines on real and imaginary parts
          The `window_length` parameter is ignored (understood to be 2).
        * `lanczos`: Lanczos filter interpolation
        * `nearest-neighbor`: Nearest neighbor (e.g., no interpolation).
          The `window_length` parameter is ignored (understood to be 0).
        * `quadratic-fit`: Fit the absolute value of the SNR to a quadratic
          function and the phase to a linear function.

    Returns
    -------
    imax_interp : float
        The interpolated index of the maximum sample, which should be between
        `imax - 0.5` and `imax + 0.5`.
    ymax_interp : complex
        The interpolated value at the maximum.

    """  # noqa: E501
    return _interpolants[method](imax, np.asarray(y), window_length)
