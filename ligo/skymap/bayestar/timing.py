# -*- coding: utf-8
#
# Copyright (C) 2013-2015  Leo Singer
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
from __future__ import division
"""
Functions for predicting timing accuracy of matched filters.
"""


import logging
import numpy as np
import lal
import lalsimulation
from scipy import interpolate
from scipy import linalg


log = logging.getLogger('BAYESTAR')


_noise_psd_funcs = {}


class vectorize_swig_psd_func(object):
    """Create a vectorized Numpy function from a SWIG-wrapped PSD function.
    SWIG does not provide enough information for Numpy to determine the number
    of input arguments, so we can't just use np.vectorize."""

    def __init__(self, str):
        self.__func = getattr(lalsimulation, str + 'Ptr')
        self.__npyfunc = np.frompyfunc(getattr(lalsimulation, str), 1, 1)

    def __call__(self, f):
        fa = np.asarray(f)
        df = np.diff(fa)
        if fa.ndim == 1 and df.size > 1 and np.all(df[0] == df[1:]):
            fa = np.concatenate((fa, [fa[-1] + df[0]]))
            ret = lal.CreateREAL8FrequencySeries(
                None, 0, fa[0], df[0], lal.DimensionlessUnit, fa.size)
            lalsimulation.SimNoisePSD(ret, 0, self.__func)
            ret = ret.data.data[:-1]
        else:
            ret = self.__npyfunc(f)
        if not np.isscalar(ret):
            ret = ret.astype(float)
        return ret


for _ifos, _func in (
    (("H1", "H2", "L1", "I1"), 'SimNoisePSDaLIGOZeroDetHighPower'),
    (("V1",), 'SimNoisePSDAdvVirgo'),
    (("K1"), 'SimNoisePSDKAGRA')
):
    _func = vectorize_swig_psd_func(_func)
    for _ifo in _ifos:
        _noise_psd_funcs[_ifo] = _func


def get_noise_psd_func(ifo):
    """Find a function that describes the given interferometer's noise PSD."""
    return _noise_psd_funcs[ifo]


class InterpolatedPSD(interpolate.interp1d):
    """Create a (linear in log-log) interpolating function for a discretely
    sampled power spectrum S(f)."""

    def __init__(self, f, S, f_high_truncate=1.0):
        assert f_high_truncate <= 1.0
        f = np.asarray(f)
        S = np.asarray(S)

        # Exclude DC if present
        if f[0] == 0:
            f = f[1:]
            S = S[1:]
        # FIXME: This is a hack to fix an issue with the detection pipeline's
        # PSD conditioning. Remove this when the issue is fixed upstream.
        if f_high_truncate < 1.0:
            log.warn(
                'Truncating PSD at %g of maximum frequency to suppress '
                'rolloff artifacts. This option may be removed in the future.',
                f_high_truncate)
            keep = (f <= f_high_truncate * max(f))
            f = f[keep]
            S = S[keep]
        super(InterpolatedPSD, self).__init__(
            np.log(f), np.log(S),
            kind='linear', bounds_error=False, fill_value=np.inf)
        self._f_min = min(f)
        self._f_max = max(f)

    def __call__(self, f):
        f_min = np.min(f)
        f_max = np.max(f)
        if f_min < self._f_min:
            log.warn('Assuming PSD is infinite at %g Hz because PSD is only '
                     'sampled down to %g Hz', f_min, self._f_min)
        if f_max > self._f_max:
            log.warn('Assuming PSD is infinite at %g Hz because PSD is only '
                     'sampled up to %g Hz', f_max, self._f_max)
        return np.where(
            (f >= self._f_min) & (f <= self._f_max),
            np.exp(super(InterpolatedPSD, self).__call__(np.log(f))), np.inf)


class SignalModel(object):
    """Class to speed up computation of signal/noise-weighted integrals and
    Barankin and Cramér-Rao lower bounds on time and phase estimation.


    Note that the autocorrelation series and the moments are related,
    as shown below.

    Create signal model:
    >>> from . import filter
    >>> sngl = lambda: None
    >>> H = filter.sngl_inspiral_psd(
    ...     'TaylorF2threePointFivePN', mass1=1.4, mass2=1.4)
    >>> S = get_noise_psd_func('H1')
    >>> W = filter.signal_psd_series(H, S)
    >>> sm = SignalModel(W)

    Compute one-sided autocorrelation function:
    >>> out_duration = 0.1
    >>> a, sample_rate = filter.autocorrelation(W, out_duration)

    Restore negative time lags using symmetry:
    >>> a = np.concatenate((a[:0:-1].conj(), a))

    Compute the first 2 frequency moments by taking derivatives of the
    autocorrelation sequence using centered finite differences.
    The nth frequency moment should be given by (-1j)^n a^(n)(t).
    >>> acor_moments = []
    >>> for i in range(2):
    ...     acor_moments.append(a[len(a) // 2])
    ...     a = -0.5j * sample_rate * (a[2:] - a[:-2])
    >>> assert np.all(np.isreal(acor_moments))
    >>> acor_moments = np.real(acor_moments)

    Compute the first 2 frequency moments using this class.
    >>> quad_moments = [sm.get_sn_moment(i) for i in range(2)]

    Compare them.
    >>> for i, (am, qm) in enumerate(zip(acor_moments, quad_moments)):
    ...     assert np.allclose(am, qm, rtol=0.05)
    """

    def __init__(self, h):
        """Create a TaylorF2 signal model with the given masses, PSD function
        S(f), PN amplitude order, and low-frequency cutoff."""

        # Find indices of first and last nonzero samples.
        nonzero = np.flatnonzero(h.data.data)
        first_nonzero = nonzero[0]
        last_nonzero = nonzero[-1]

        # Frequency sample points
        self.dw = 2 * np.pi * h.deltaF
        f = h.f0 + h.deltaF * np.arange(first_nonzero, last_nonzero + 1)
        self.w = 2 * np.pi * f

        # Throw away leading and trailing zeros.
        h = h.data.data[first_nonzero:last_nonzero + 1]

        self.denom_integrand = 4 / (2 * np.pi) * h
        self.den = np.trapz(self.denom_integrand, dx=self.dw)

    def get_horizon_distance(self, snr_thresh=1):
        return np.sqrt(self.den) / snr_thresh

    def get_sn_average(self, func):
        """Get the average of a function of angular frequency, weighted by the
        signal to noise per unit angular frequency."""
        num = np.trapz(func(self.w) * self.denom_integrand, dx=self.dw)
        return num / self.den

    def get_sn_moment(self, power):
        """Get the average of angular frequency to the given power, weighted by
        the signal to noise per unit frequency."""
        return self.get_sn_average(lambda w: w**power)

    def get_crb(self, snr):
        """Get the Cramér-Rao bound, or inverse Fisher information matrix,
        describing the phase and time estimation covariance."""
        w1 = self.get_sn_moment(1)
        w2 = self.get_sn_moment(2)
        I = np.asarray(((1, -w1), (-w1, w2)))
        return linalg.inv(I) / np.square(snr)

    # FIXME: np.vectorize doesn't work on unbound instance methods. The
    # excluded keyword, added in Numpy 1.7, could be used here to exclude the
    # zeroth argument, self.
    def __get_crb_toa_uncert(self, snr):
        return np.sqrt(self.get_crb(snr)[1, 1])

    def get_crb_toa_uncert(self, snr):
        return np.frompyfunc(self.__get_crb_toa_uncert, 1, 1)(snr)
