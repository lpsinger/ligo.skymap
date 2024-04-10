#
# Copyright (C) 2013-2024  Leo Singer
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
"""
Rapid sky localization with BAYESTAR [1]_.

References
----------
.. [1] Singer & Price, 2016. "Rapid Bayesian position reconstruction for
       gravitational-wave transients." PRD, 93, 024013.
       :doi:`10.1103/PhysRevD.93.024013`

"""
import inspect
import logging
import os
import sys
from textwrap import wrap

from astropy.table import Column, Table
from astropy import units as u
import lal
import lalsimulation
import numpy as np

from .. import distance
from . import filter
from ..io.hdf5 import write_samples
from ..io.fits import metadata_for_version_module
from ..io.events.base import Event
from . import filter  # noqa
from ..kde import Clustered2Plus1DSkyKDE
from .. import moc
from .. import healpix_tree
from .. import version
from .. import core
from ..core import (antenna_factor, signal_amplitude_model,
                    log_posterior_toa_phoa_snr as _log_posterior_toa_phoa_snr)
from ..util.numpy import require_contiguous_aligned
from ..util.stopwatch import Stopwatch
from .ez_emcee import ez_emcee

__all__ = ('derasterize', 'localize', 'rasterize', 'antenna_factor',
           'signal_amplitude_model')

log = logging.getLogger('BAYESTAR')

_RESCALE_LOGLIKELIHOOD = 0.83

signal_amplitude_model = require_contiguous_aligned(signal_amplitude_model)
_log_posterior_toa_phoa_snr = require_contiguous_aligned(
    _log_posterior_toa_phoa_snr)


# Wrap so that ufunc parameter names are known
def log_posterior_toa_phoa_snr(
        ra, sin_dec, distance, u, twopsi, t, min_distance, max_distance,
        prior_distance_power, cosmology, gmst, sample_rate, epochs, snrs,
        responses, locations, horizons,
        rescale_loglikelihood=_RESCALE_LOGLIKELIHOOD):
    return _log_posterior_toa_phoa_snr(
        ra, sin_dec, distance, u, twopsi, t, min_distance, max_distance,
        prior_distance_power, cosmology, gmst, sample_rate, epochs, snrs,
        responses, locations, horizons, rescale_loglikelihood)


# Wrap so that fixed parameter values are pulled from keyword arguments
def log_post(params, *args, **kwargs):
    # Names of parameters
    keys = ('ra', 'sin_dec', 'distance', 'u', 'twopsi', 't')

    params = list(params.T)
    kwargs = dict(kwargs)

    return log_posterior_toa_phoa_snr(
        *(kwargs.pop(key) if key in kwargs else params.pop(0) for key in keys),
        *args, **kwargs)


def localize_emcee(args, xmin, xmax, chain_dump=None):
    # Gather posterior samples
    chain = ez_emcee(log_post, xmin, xmax, args=args, vectorize=True)

    # Transform back from sin_dec to dec and cos_inclination to inclination
    chain[:, 1] = np.arcsin(chain[:, 1])
    chain[:, 3] = np.arccos(chain[:, 3])

    # Optionally save posterior sample chain to file.
    if chain_dump:
        _, ndim = chain.shape
        names = 'ra dec distance inclination twopsi time'.split()[:ndim]
        write_samples(Table(rows=chain, names=names, copy=False), chain_dump,
                      path='/bayestar/posterior_samples', overwrite=True)

    # Pass a random subset of 1000 points to the KDE, to save time.
    pts = np.random.permutation(chain)[:1000, :3]
    ckde = Clustered2Plus1DSkyKDE(pts)
    return ckde.as_healpix()


def condition(
        event, waveform='o2-uberbank', f_low=30.0,
        enable_snr_series=True, f_high_truncate=0.95):

    if len(event.singles) == 0:
        raise ValueError('Cannot localize an event with zero detectors.')

    singles = event.singles
    if not enable_snr_series:
        singles = [single for single in singles if single.snr is not None]

    ifos = [single.detector for single in singles]

    # Extract SNRs from table.
    snrs = np.ma.asarray([
        np.ma.masked if single.snr is None else single.snr
        for single in singles])

    # Look up physical parameters for detector.
    detectors = [lalsimulation.DetectorPrefixToLALDetector(str(ifo))
                 for ifo in ifos]
    responses = np.asarray([det.response for det in detectors])
    locations = np.asarray([det.location for det in detectors]) / lal.C_SI

    # Power spectra for each detector.
    psds = [single.psd for single in singles]
    psds = [filter.InterpolatedPSD(filter.abscissa(psd), psd.data.data,
                                   f_high_truncate=f_high_truncate)
            for psd in psds]

    log.debug('calculating templates')
    H = filter.sngl_inspiral_psd(waveform, f_min=f_low, **event.template_args)

    log.debug('calculating noise PSDs')
    HS = [filter.signal_psd_series(H, S) for S in psds]

    # Signal models for each detector.
    log.debug('calculating Fisher matrix elements')
    signal_models = [filter.SignalModel(_) for _ in HS]

    # Get SNR=1 horizon distances for each detector.
    horizons = np.asarray([signal_model.get_horizon_distance()
                           for signal_model in signal_models])

    weights = np.ma.asarray([
        1 / np.square(signal_model.get_crb_toa_uncert(snr))
        for signal_model, snr in zip(signal_models, snrs)])

    # Center detector array.
    locations -= (np.sum(locations * weights.reshape(-1, 1), axis=0) /
                  np.sum(weights))

    if enable_snr_series:
        snr_series = [single.snr_series for single in singles]
        if all(s is None for s in snr_series):
            snr_series = None
    else:
        snr_series = None

    # Maximum barycentered arrival time error:
    # |distance from array barycenter to furthest detector| / c + 5 ms.
    # For LHO+LLO, this is 15.0 ms.
    # For an arbitrary terrestrial detector network, the maximum is 26.3 ms.
    max_abs_t = np.max(
        np.sqrt(np.sum(np.square(locations), axis=1))) + 0.005

    # Calculate a minimum arrival time prior for low-bandwidth signals.
    # This is important for early-warning events, for which the light travel
    # time between detectors may be only a tiny fraction of the template
    # autocorrelation time.
    #
    # The period of the autocorrelation function is 2 pi times the timing
    # uncertainty at SNR=1 (see Eq. (24) of [1]_). We want a quarter period:
    # there is a factor of a half because we want the first zero crossing,
    # and another factor of a half because the SNR time series will be cropped
    # to a two-sided interval.
    max_acor_t = max(0.5 * np.pi * signal_model.get_crb_toa_uncert(1)
                     for signal_model in signal_models)
    max_abs_t = max(max_abs_t, max_acor_t)

    if snr_series is None:
        log.warning("No SNR time series found, so we are creating a "
                    "zero-noise SNR time series from the whitened template's "
                    "autocorrelation sequence. The sky localization "
                    "uncertainty may be underestimated.")

        acors, sample_rates = zip(
            *[filter.autocorrelation(_, max_abs_t) for _ in HS])
        sample_rate = sample_rates[0]
        deltaT = 1 / sample_rate
        nsamples = len(acors[0])
        assert all(sample_rate == _ for _ in sample_rates)
        assert all(nsamples == len(_) for _ in acors)
        nsamples = nsamples * 2 - 1

        snr_series = []
        for acor, single in zip(acors, singles):
            series = lal.CreateCOMPLEX8TimeSeries(
                'fake SNR', 0, 0, deltaT, lal.StrainUnit, nsamples)
            series.epoch = single.time - 0.5 * (nsamples - 1) * deltaT
            acor = np.concatenate((np.conj(acor[:0:-1]), acor))
            series.data.data = single.snr * filter.exp_i(single.phase) * acor
            snr_series.append(series)

    # Ensure that all of the SNR time series have the same sample rate.
    # FIXME: for now, the Python wrapper expects all of the SNR time series to
    # also be the same length.
    deltaT = snr_series[0].deltaT
    sample_rate = 1 / deltaT
    if any(deltaT != series.deltaT for series in snr_series):
        raise ValueError('BAYESTAR does not yet support SNR time series with '
                         'mixed sample rates')

    # Ensure that all of the SNR time series have odd lengths.
    if any(len(series.data.data) % 2 == 0 for series in snr_series):
        raise ValueError('SNR time series must have odd lengths')

    # Trim time series to the desired length.
    max_abs_n = int(np.ceil(max_abs_t * sample_rate))
    desired_length = 2 * max_abs_n - 1
    for i, series in enumerate(snr_series):
        length = len(series.data.data)
        if length > desired_length:
            snr_series[i] = lal.CutCOMPLEX8TimeSeries(
                series, length // 2 + 1 - max_abs_n, desired_length)

    # FIXME: for now, the Python wrapper expects all of the SNR time sries to
    # also be the same length.
    nsamples = len(snr_series[0].data.data)
    if any(nsamples != len(series.data.data) for series in snr_series):
        raise ValueError('BAYESTAR does not yet support SNR time series of '
                         'mixed lengths')

    # Perform sanity checks that the middle sample of the SNR time series match
    # the sngl_inspiral records to the nearest sample (plus the smallest
    # representable LIGOTimeGPS difference of 1 nanosecond).
    for ifo, single, series in zip(ifos, singles, snr_series):
        shift = np.abs(0.5 * (nsamples - 1) * series.deltaT +
                       float(series.epoch - single.time))
        if shift >= deltaT + 1e-8:
            raise ValueError('BAYESTAR expects the SNR time series to be '
                             'centered on the single-detector trigger times, '
                             'but {} was off by {} s'.format(ifo, shift))

    # Extract the TOAs in GPS nanoseconds from the SNR time series, assuming
    # that the trigger happened in the middle.
    toas_ns = [series.epoch.ns() + 1e9 * 0.5 * (len(series.data.data) - 1) *
               series.deltaT for series in snr_series]

    # Collect all of the SNR series in one array.
    snr_series = np.vstack([series.data.data for series in snr_series])

    # Center times of arrival and compute GMST at mean arrival time.
    # Pre-center in integer nanoseconds to preserve precision of
    # initial datatype.
    epoch = sum(toas_ns) // len(toas_ns)
    toas = 1e-9 * (np.asarray(toas_ns) - epoch)
    mean_toa = np.average(toas, weights=weights)
    toas -= mean_toa
    epoch += int(np.round(1e9 * mean_toa))
    epoch = lal.LIGOTimeGPS(0, int(epoch))

    # Translate SNR time series back to time of first sample.
    toas -= 0.5 * (nsamples - 1) * deltaT

    # Convert complex SNRS to amplitude and phase
    snrs_abs = np.abs(snr_series)
    snrs_arg = filter.unwrap(np.angle(snr_series))
    snrs = np.stack((snrs_abs, snrs_arg), axis=-1)

    return epoch, sample_rate, toas, snrs, responses, locations, horizons


def condition_prior(horizons, min_distance=None, max_distance=None,
                    prior_distance_power=None, cosmology=False):
    if cosmology:
        log.warning('Enabling cosmological prior. This feature is UNREVIEWED.')

    # If minimum distance is not specified, then default to 0 Mpc.
    if min_distance is None:
        min_distance = 0

    # If maximum distance is not specified, then default to the SNR=4
    # horizon distance of the most sensitive detector.
    if max_distance is None:
        max_distance = max(horizons) / 4

    # If prior_distance_power is not specified, then default to 2
    # (p(r) ~ r^2, uniform in volume).
    if prior_distance_power is None:
        prior_distance_power = 2

    # Raise an exception if 0 Mpc is the minimum effective distance and the
    # prior is of the form r**k for k<0
    if min_distance == 0 and prior_distance_power < 0:
        raise ValueError(('Prior is a power law r^k with k={}, '
                          'undefined at min_distance=0').format(
                              prior_distance_power))

    return min_distance, max_distance, prior_distance_power, cosmology


def localize(
        event, waveform='o2-uberbank', f_low=30.0,
        min_distance=None, max_distance=None, prior_distance_power=None,
        cosmology=False, mcmc=False, chain_dump=None,
        enable_snr_series=True, f_high_truncate=0.95,
        rescale_loglikelihood=_RESCALE_LOGLIKELIHOOD):
    """Localize a compact binary signal using the BAYESTAR algorithm.

    Parameters
    ----------
    event : `ligo.skymap.io.events.Event`
        The event candidate.
    waveform : str, optional
        The name of the waveform approximant.
    f_low : float, optional
        The low frequency cutoff.
    min_distance, max_distance : float, optional
        The limits of integration over luminosity distance, in Mpc
        (default: determine automatically from detector sensitivity).
    prior_distance_power : int, optional
        The power of distance that appears in the prior
        (default: 2, uniform in volume).
    cosmology: bool, optional
        Set to enable a uniform in comoving volume prior (default: false).
    mcmc : bool, optional
        Set to use MCMC sampling rather than more accurate Gaussian quadrature.
    chain_dump : str, optional
        Save posterior samples to this filename if `mcmc` is set.
    enable_snr_series : bool, optional
        Set to False to disable SNR time series.
    f_high_truncate : float, optional
        Truncate the noise power spectral densities at this factor times the
        highest sampled frequency to suppress artifacts caused by incorrect
        PSD conditioning by some matched filter pipelines.

    Returns
    -------
    skymap : `astropy.table.Table`
        A 3D sky map in multi-order HEALPix format.

    """
    # Hide event parameters, but show all other arguments
    def formatvalue(value):
        if isinstance(value, Event):
            return '=...'
        else:
            return '=' + repr(value)

    frame = inspect.currentframe()
    argstr = inspect.formatargvalues(*inspect.getargvalues(frame),
                                     formatvalue=formatvalue)

    stopwatch = Stopwatch()
    stopwatch.start()

    epoch, sample_rate, toas, snrs, responses, locations, horizons = \
        condition(event, waveform=waveform, f_low=f_low,
                  enable_snr_series=enable_snr_series,
                  f_high_truncate=f_high_truncate)

    min_distance, max_distance, prior_distance_power, cosmology = \
        condition_prior(horizons, min_distance, max_distance,
                        prior_distance_power, cosmology)

    gmst = lal.GreenwichMeanSiderealTime(epoch)

    # Time and run sky localization.
    log.debug('starting computationally-intensive section')
    args = (min_distance, max_distance, prior_distance_power, cosmology, gmst,
            sample_rate, toas, snrs, responses, locations, horizons,
            rescale_loglikelihood)
    if mcmc:
        max_abs_t = 2 * snrs.data.shape[1] / sample_rate
        skymap = localize_emcee(
            args=args,
            xmin=[0, -1, min_distance, -1, 0, 0],
            xmax=[2 * np.pi, 1, max_distance, 1, 2 * np.pi, 2 * max_abs_t],
            chain_dump=chain_dump)
    else:
        skymap, log_bci, log_bsn = core.toa_phoa_snr(*args)
        skymap = Table(skymap, copy=False)
        skymap.meta['log_bci'] = log_bci
        skymap.meta['log_bsn'] = log_bsn

    # Convert distance moments to parameters
    try:
        distmean = skymap.columns.pop('DISTMEAN')
        diststd = skymap.columns.pop('DISTSTD')
    except KeyError:
        distmean, diststd, _ = distance.parameters_to_moments(
            skymap['DISTMU'], skymap['DISTSIGMA'])
    else:
        skymap['DISTMU'], skymap['DISTSIGMA'], skymap['DISTNORM'] = \
            distance.moments_to_parameters(distmean, diststd)

    # Add marginal distance moments
    good = np.isfinite(distmean) & np.isfinite(diststd)
    prob = (moc.uniq2pixarea(skymap['UNIQ']) * skymap['PROBDENSITY'])[good]
    distmean = distmean[good]
    diststd = diststd[good]
    rbar = (prob * distmean).sum()
    r2bar = (prob * (np.square(diststd) + np.square(distmean))).sum()
    skymap.meta['distmean'] = rbar
    skymap.meta['diststd'] = np.sqrt(r2bar - np.square(rbar))

    stopwatch.stop()
    end_time = lal.GPSTimeNow()
    log.info('finished computationally-intensive section in %s', stopwatch)

    # Fill in metadata and return.
    program, _ = os.path.splitext(os.path.basename(sys.argv[0]))
    skymap.meta.update(metadata_for_version_module(version))
    skymap.meta['creator'] = 'BAYESTAR'
    skymap.meta['origin'] = 'LIGO/Virgo/KAGRA'
    skymap.meta['gps_time'] = float(epoch)
    skymap.meta['runtime'] = stopwatch.real
    skymap.meta['instruments'] = {single.detector for single in event.singles}
    skymap.meta['gps_creation_time'] = end_time
    skymap.meta['history'] = [
        '',
        'Generated by calling the following Python function:',
        *wrap('{}.{}{}'.format(__name__, frame.f_code.co_name, argstr), 72),
        '',
        'This was the command line that started the program:',
        *wrap(' '.join([program] + sys.argv[1:]), 72)]

    return skymap


def rasterize(skymap, order=None):
    orig_order, _ = moc.uniq2nest(skymap['UNIQ'].max())

    # Determine whether we need to do nontrivial downsampling.
    downsampling = (order is not None and 0 <= order < orig_order
                    and 'DISTMU' in skymap.dtype.fields.keys())

    # If we are downsampling, then convert from distance parameters to
    # distance moments times probability density so that the averaging
    # that is automatically done by moc.rasterize() correctly marginalizes
    # the moments.
    if downsampling:
        skymap = Table(skymap, copy=True, meta=skymap.meta)

        probdensity = skymap['PROBDENSITY']
        distmu = skymap.columns.pop('DISTMU')
        distsigma = skymap.columns.pop('DISTSIGMA')

        bad = ~(np.isfinite(distmu) & np.isfinite(distsigma))
        distmean, diststd, _ = distance.parameters_to_moments(
            distmu, distsigma)
        distmean[bad] = np.nan
        diststd[bad] = np.nan
        skymap['DISTMEAN'] = probdensity * distmean
        skymap['DISTVAR'] = probdensity * (
            np.square(diststd) + np.square(distmean))

    skymap = Table(moc.rasterize(skymap, order=order),
                   meta=skymap.meta, copy=False)

    # If we are downsampling, then convert back to distance parameters.
    if downsampling:
        distmean = skymap.columns.pop('DISTMEAN') / skymap['PROBDENSITY']
        diststd = np.sqrt(
            skymap.columns.pop('DISTVAR') / skymap['PROBDENSITY']
            - np.square(distmean))
        skymap['DISTMU'], skymap['DISTSIGMA'], skymap['DISTNORM'] = \
            distance.moments_to_parameters(
                distmean, diststd)

    skymap.rename_column('PROBDENSITY', 'PROB')
    skymap['PROB'] *= 4 * np.pi / len(skymap)
    skymap['PROB'].unit = u.pixel ** -1
    return skymap


def derasterize(skymap):
    skymap.rename_column('PROB', 'PROBDENSITY')
    skymap['PROBDENSITY'] *= len(skymap) / (4 * np.pi)
    skymap['PROBDENSITY'].unit = u.steradian ** -1
    nside, _, ipix, _, _, value = zip(
        *healpix_tree.reconstruct_nested(skymap))
    nside = np.asarray(nside)
    ipix = np.asarray(ipix)
    value = np.stack(value)
    uniq = (4 * np.square(nside) + ipix)
    old_units = [column.unit for column in skymap.columns.values()]
    skymap = Table(value, meta=skymap.meta, copy=False)
    for old_unit, column in zip(old_units, skymap.columns.values()):
        column.unit = old_unit
    skymap.add_column(Column(uniq, name='UNIQ'), 0)
    skymap.sort('UNIQ')
    return skymap


def test():
    """Run BAYESTAR C unit tests.

    Examples
    --------
    >>> test()
    0

    """
    return int(core.test())


del require_contiguous_aligned
