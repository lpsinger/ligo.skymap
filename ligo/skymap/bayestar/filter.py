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
Utility functions for BAYESTAR that are related to matched filtering.
"""

import logging
import math

import lal
import lalsimulation
import numpy as np
from scipy import interpolate
from scipy import fftpack as fft
from scipy import linalg

log = logging.getLogger('BAYESTAR')


def ceil_pow_2(n):
    """Return the least integer power of 2 that is greater than or equal to n.

    >>> ceil_pow_2(128.0)
    128.0
    >>> ceil_pow_2(0.125)
    0.125
    >>> ceil_pow_2(129.0)
    256.0
    >>> ceil_pow_2(0.126)
    0.25
    >>> ceil_pow_2(1.0)
    1.0
    """
    # frexp splits floats into mantissa and exponent, ldexp does the opposite.
    # For positive numbers, mantissa is in [0.5, 1.).
    mantissa, exponent = math.frexp(n)
    return math.ldexp(
        1 if mantissa >= 0 else float('nan'),
        exponent - 1 if mantissa == 0.5 else exponent
    )


def abscissa(series):
    """Produce the independent variable for a lal TimeSeries or
    FrequencySeries."""
    try:
        delta = series.deltaT
        x0 = float(series.epoch)
    except AttributeError:
        delta = series.deltaF
        x0 = series.f0
    return x0 + delta * np.arange(len(series.data.data))


def exp_i(phi):
    return np.cos(phi) + np.sin(phi) * 1j


def truncated_ifft(y, nsamples_out=None):
    r"""Truncated inverse FFT.

    See http://www.fftw.org/pruned.html for a discussion of related algorithms.

    Perform inverse FFT to obtain truncated autocorrelation time series.
    This makes use of a folded DFT for a speedup of::

        log(nsamples)/log(nsamples_out)

    over directly computing the inverse FFT and truncating. Here is how it
    works. Say we have a frequency-domain signal X[k], for 0 ≤ k ≤ N - 1. We
    want to compute its DFT x[n], for 0 ≤ n ≤ M, where N is divisible by M:
    N = cM, for some integer c. The DFT is::

               N - 1
               ______
               \           2 π i k n
        x[n] =  \     exp[-----------] Y[k]
               /               N
              /------
               k = 0

               c - 1   M - 1
               ______  ______
               \       \           2 π i n (m c + j)
             =  \       \     exp[------------------] Y[m c + j]
               /       /                 c M
              /------ /------
               j = 0   m = 0

               c - 1                     M - 1
               ______                    ______
               \           2 π i n j     \           2 π i n m
             =  \     exp[-----------]    \     exp[-----------] Y[m c + j]
               /               N         /               M
              /------                   /------
               j = 0                     m = 0

    So: we split the frequency series into c deinterlaced sub-signals, each of
    length M, compute the DFT of each sub-signal, and add them back together
    with complex weights.

    Parameters
    ----------
    y : `numpy.ndarray`
        Complex input vector.
    nsamples_out : int, optional
        Length of output vector. By default, same as length of input vector.

    Returns
    -------
    x : `numpy.ndarray`
        The first nsamples_out samples of the IFFT of x, zero-padded if

    Examples
    --------

    First generate the IFFT of a random signal:

    >>> nsamples_out = 1024
    >>> y = np.random.randn(nsamples_out) + np.random.randn(nsamples_out) * 1j
    >>> x = fft.ifft(y)

    Now check that the truncated IFFT agrees:

    >>> np.allclose(x, truncated_ifft(y), rtol=1e-15)
    True
    >>> np.allclose(x, truncated_ifft(y, 1024), rtol=1e-15)
    True
    >>> np.allclose(x[:128], truncated_ifft(y, 128), rtol=1e-15)
    True
    >>> np.allclose(x[:1], truncated_ifft(y, 1), rtol=1e-15)
    True
    >>> np.allclose(x[:32], truncated_ifft(y, 32), rtol=1e-15)
    True
    >>> np.allclose(x[:63], truncated_ifft(y, 63), rtol=1e-15)
    True
    >>> np.allclose(x[:25], truncated_ifft(y, 25), rtol=1e-15)
    True
    >>> truncated_ifft(y, 1025)
    Traceback (most recent call last):
      ...
    ValueError: Input is too short: you gave me an input of length 1024, but you asked for an IFFT of length 1025.
    """  # noqa: E501
    nsamples = len(y)
    if nsamples_out is None:
        nsamples_out = nsamples
    elif nsamples_out > nsamples:
        raise ValueError(
            'Input is too short: you gave me an input of length {0}, '
            'but you asked for an IFFT of length {1}.'.format(
                nsamples, nsamples_out))
    elif nsamples & (nsamples - 1):
        raise NotImplementedError(
            'I am too lazy to implement for nsamples that is '
            'not a power of 2.')

    # Find number of FFTs.
    # FIXME: only works if nsamples is a power of 2.
    # Would be better to find the smallest divisor of nsamples that is
    # greater than or equal to nsamples_out.
    nsamples_batch = int(ceil_pow_2(nsamples_out))
    c = nsamples // nsamples_batch

    # FIXME: Implement for real-to-complex FFTs as well.
    twiddle = exp_i(2 * np.pi * np.arange(nsamples_batch) / nsamples)

    x = fft.ifft(y.reshape(nsamples_batch, c).T)

    result = x[-1]
    for row in x[-2::-1]:
        result *= twiddle  # FIXME: check stability of this recurrence relation
        result += row

    # Now need to truncate remaining samples.
    if nsamples_out < nsamples_batch:
        result = result[:nsamples_out]

    return result / c


def get_approximant_and_orders_from_string(s):
    """Determine the approximant, amplitude order, and phase order for a string
    of the form "TaylorT4threePointFivePN". In this example, the waveform is
    "TaylorT4" and the phase order is 7 (twice 3.5). If the input contains the
    substring "restricted" or "Restricted", then the amplitude order is taken
    to be 0. Otherwise, the amplitude order is the same as the phase order."""
    # SWIG-wrapped functions apparently do not understand Unicode, but
    # often the input argument will come from a Unicode XML file.
    s = str(s)
    approximant = lalsimulation.GetApproximantFromString(s)
    try:
        phase_order = lalsimulation.GetOrderFromString(s)
    except RuntimeError:
        phase_order = -1
    if 'restricted' in s or 'Restricted' in s:
        amplitude_order = 0
    else:
        amplitude_order = phase_order
    return approximant, amplitude_order, phase_order


def get_f_lso(mass1, mass2):
    """Calculate the GW frequency during the last stable orbit of a compact
    binary."""
    return 1 / (6 ** 1.5 * np.pi * (mass1 + mass2) * lal.MTSUN_SI)


def sngl_inspiral_psd(waveform, mass1, mass2, f_min=10, f_max=2048, f_ref=0,
                      **kwargs):
    # FIXME: uberbank mass criterion. Should find a way to get this from
    # pipeline output metadata.
    if waveform == 'o1-uberbank':
        log.warning('Template is unspecified; '
                    'using ER8/O1 uberbank criterion')
        if mass1 + mass2 < 4:
            waveform = 'TaylorF2threePointFivePN'
        else:
            waveform = 'SEOBNRv2_ROM_DoubleSpin'
    elif waveform == 'o2-uberbank':
        log.warning('Template is unspecified; '
                    'using ER10/O2 uberbank criterion')
        if mass1 + mass2 < 4:
            waveform = 'TaylorF2threePointFivePN'
        else:
            waveform = 'SEOBNRv4_ROM'
    approx, ampo, phaseo = get_approximant_and_orders_from_string(waveform)
    log.info('Selected template: %s', waveform)

    # Generate conditioned template.
    params = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(params, phaseo)
    lalsimulation.SimInspiralWaveformParamsInsertPNAmplitudeOrder(params, ampo)
    hplus, hcross = lalsimulation.SimInspiralFD(
        m1=float(mass1) * lal.MSUN_SI, m2=float(mass2) * lal.MSUN_SI,
        S1x=float(kwargs.get('spin1x') or 0),
        S1y=float(kwargs.get('spin1y') or 0),
        S1z=float(kwargs.get('spin1z') or 0),
        S2x=float(kwargs.get('spin2x') or 0),
        S2y=float(kwargs.get('spin2y') or 0),
        S2z=float(kwargs.get('spin2z') or 0),
        distance=1e6 * lal.PC_SI, inclination=0, phiRef=0,
        longAscNodes=0, eccentricity=0, meanPerAno=0,
        deltaF=0, f_min=f_min, f_max=f_max, f_ref=f_ref,
        LALparams=params, approximant=approx)

    # Force `plus' and `cross' waveform to be in quadrature.
    h = 0.5 * (hplus.data.data + 1j * hcross.data.data)

    # For inspiral-only waveforms, nullify frequencies beyond ISCO.
    # FIXME: the waveform generation functions pick the end frequency
    # automatically. Shouldn't SimInspiralFD?
    inspiral_only_waveforms = (
        lalsimulation.TaylorF2,
        lalsimulation.SpinTaylorF2,
        lalsimulation.TaylorF2RedSpin,
        lalsimulation.TaylorF2RedSpinTidal,
        lalsimulation.SpinTaylorT4Fourier)
    if approx in inspiral_only_waveforms:
        h[abscissa(hplus) >= get_f_lso(mass1, mass2)] = 0

    # Drop Nyquist frequency.
    if len(h) % 2:
        h = h[:-1]

    # Create output frequency series.
    psd = lal.CreateREAL8FrequencySeries(
        'signal PSD', 0, hplus.f0, hcross.deltaF, hplus.sampleUnits**2, len(h))
    psd.data.data = abs2(h)

    # Done!
    return psd


def signal_psd_series(H, S):
    n = H.data.data.size
    f = H.f0 + np.arange(1, n) * H.deltaF
    ret = lal.CreateREAL8FrequencySeries(
        'signal PSD / noise PSD', 0, H.f0, H.deltaF, lal.DimensionlessUnit, n)
    ret.data.data[0] = 0
    ret.data.data[1:] = H.data.data[1:] / S(f)
    return ret


def autocorrelation(H, out_duration, normalize=True):
    """
    Calculate the complex autocorrelation sequence a(t), for t >= 0, of an
    inspiral signal.

    Parameters
    ----------
    H : lal.REAL8FrequencySeries
        Signal PSD series.
    S : callable
        Noise power spectral density function.

    Returns
    -------
    acor : `numpy.ndarray`
        The complex-valued autocorrelation sequence.
    sample_rate : float
        The sample rate.
    """

    # Compute duration of template, rounded up to a power of 2.
    H_len = H.data.data.size
    nsamples = 2 * H_len
    sample_rate = nsamples * H.deltaF

    # Compute autopower spectral density.
    power = np.empty(nsamples, H.data.data.dtype)
    power[:H_len] = H.data.data
    power[H_len:] = 0

    # Determine length of output FFT.
    nsamples_out = int(np.ceil(out_duration * sample_rate))

    acor = truncated_ifft(power, nsamples_out)
    if normalize:
        acor /= np.abs(acor[0])

    # If we have done this right, then the zeroth sample represents lag 0
    if np.all(np.isreal(H.data.data)):
        assert np.argmax(np.abs(acor)) == 0
        assert np.isreal(acor[0])

    # Done!
    return acor, float(sample_rate)


def abs2(y):
    """Return the absolute value squared, :math:`|z|^2` ,for a complex number
    :math:`z`, without performing a square root."""
    return np.square(y.real) + np.square(y.imag)


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


class InterpolatedPSD(interpolate.interp1d):
    """Create a (linear in log-log) interpolating function for a discretely
    sampled power spectrum S(f)."""

    def __init__(self, f, S, f_high_truncate=1.0, fill_value=np.inf):
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
            log.warning(
                'Truncating PSD at %g of maximum frequency to suppress '
                'rolloff artifacts. This option may be removed in the future.',
                f_high_truncate)
            keep = (f <= f_high_truncate * max(f))
            f = f[keep]
            S = S[keep]
        super(InterpolatedPSD, self).__init__(
            np.log(f), np.log(S),
            kind='linear', bounds_error=False, fill_value=np.log(fill_value))
        self._f_min = min(f)
        self._f_max = max(f)

    def __call__(self, f):
        f_min = np.min(f)
        f_max = np.max(f)
        if f_min < self._f_min:
            log.warning('Assuming PSD is infinite at %g Hz because PSD is '
                        'only sampled down to %g Hz', f_min, self._f_min)
        if f_max > self._f_max:
            log.warning('Assuming PSD is infinite at %g Hz because PSD is '
                        'only sampled up to %g Hz', f_max, self._f_max)
        return np.where(
            (f >= self._f_min) & (f <= self._f_max),
            np.exp(super(InterpolatedPSD, self).__call__(np.log(f))),
            np.exp(self.fill_value))


class SignalModel(object):
    """Class to speed up computation of signal/noise-weighted integrals and
    Barankin and Cramér-Rao lower bounds on time and phase estimation.


    Note that the autocorrelation series and the moments are related,
    as shown below.

    Examples
    --------

    Create signal model:

    >>> from . import filter
    >>> sngl = lambda: None
    >>> H = filter.sngl_inspiral_psd(
    ...     'TaylorF2threePointFivePN', mass1=1.4, mass2=1.4)
    >>> S = vectorize_swig_psd_func('SimNoisePSDaLIGOZeroDetHighPower')
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
        fisher = np.asarray(((1, -w1), (-w1, w2)))
        return linalg.inv(fisher) / np.square(snr)

    # FIXME: np.vectorize doesn't work on unbound instance methods. The
    # excluded keyword, added in Numpy 1.7, could be used here to exclude the
    # zeroth argument, self.
    def __get_crb_toa_uncert(self, snr):
        return np.sqrt(self.get_crb(snr)[1, 1])

    def get_crb_toa_uncert(self, snr):
        return np.frompyfunc(self.__get_crb_toa_uncert, 1, 1)(snr)
