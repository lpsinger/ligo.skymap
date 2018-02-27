# -*- coding: UTF-8 -*-
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
"""
Functions related to matched filtering.
"""

# General imports
import logging
import numpy as np
import math
from scipy import optimize

# LAL imports
import lal
import lalsimulation

# My own imports
from .decorator import memoized


log = logging.getLogger('BAYESTAR')


# Memoize FFT plans
CreateReverseCOMPLEX16FFTPlan = memoized(lal.CreateReverseCOMPLEX16FFTPlan)


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
    """Truncated inverse FFT.

    See http://www.fftw.org/pruned.html for a discussion of related algorithms.

    Perform inverse FFT to obtain truncated autocorrelation time series.
    This makes use of a folded DFT for a speedup of

        log(nsamples)/log(nsamples_out)

    over directly computing the inverse FFT and truncating. Here is how it
    works. Say we have a frequency-domain signal X[k], for 0 ≤ k ≤ N - 1. We
    want to compute its DFT x[n], for 0 ≤ n ≤ M, where N is divisible by M:
    N = cM, for some integer c. The DFT is:

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


    First generate the IFFT of a random signal:
    >>> nsamples_out = 1024
    >>> y = np.random.randn(nsamples_out) + np.random.randn(nsamples_out) * 1j
    >>> plan = CreateReverseCOMPLEX16FFTPlan(nsamples_out, 0)
    >>> freq = lal.CreateCOMPLEX16Vector(nsamples_out)
    >>> freq.data = y
    >>> time = lal.CreateCOMPLEX16Vector(nsamples_out)
    >>> _ = lal.COMPLEX16VectorFFT(time, freq, plan)
    >>> x = time.data

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
    """
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
    plan = CreateReverseCOMPLEX16FFTPlan(nsamples_batch, 0)
    input_workspace = lal.CreateCOMPLEX16Vector(nsamples_batch)
    output_workspace = lal.CreateCOMPLEX16Vector(nsamples_batch)
    twiddle = exp_i(2 * np.pi * np.arange(nsamples_batch) / nsamples)

    j = c - 1
    input_workspace.data = y[j::c]
    lal.COMPLEX16VectorFFT(output_workspace, input_workspace, plan)
    x = output_workspace.data.copy()  # Make sure this is a deep copy

    for j in range(c-2, -1, -1):
        input_workspace.data = y[j::c]
        lal.COMPLEX16VectorFFT(output_workspace, input_workspace, plan)
        x *= twiddle  # FIXME: check stability of this recurrence relation.
        x += output_workspace.data

    # Now need to truncate remaining samples.
    if nsamples_out < nsamples_batch:
        x = x[:nsamples_out]

    return x


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


def sngl_inspiral_psd(waveform, mass1, mass2, f_min=10, f_max=2048, f_ref=0, **kwargs):
    # FIXME: uberbank mass criterion. Should find a way to get this from
    # pipeline output metadata.
    if waveform == 'o1-uberbank':
        log.warn('Template is unspecified; using ER8/O1 uberbank criterion')
        if mass1 + mass2 < 4:
            waveform = 'TaylorF2threePointFivePN'
        else:
            waveform = 'SEOBNRv2_ROM_DoubleSpin'
    elif waveform == 'o2-uberbank':
        log.warn('Template is unspecified; using ER10/O2 uberbank criterion')
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
        distance=1e6*lal.PC_SI, inclination=0, phiRef=0,
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


def autocorrelation(H, out_duration):
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
    acor /= np.abs(acor[0])

    # If we have done this right, then the zeroth sample represents lag 0
    assert np.argmax(np.abs(acor)) == 0
    assert np.isreal(acor[0])

    # Done!
    return acor, float(sample_rate)


def abs2(y):
    """Return the |z|^2 for a complex number z."""
    return np.square(y.real) + np.square(y.imag)