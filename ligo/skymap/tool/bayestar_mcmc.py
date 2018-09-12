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
Markov-Chain Monte Carlo sky localization.
"""

from . import (
    ArgumentParser, FileType, waveform_parser,
    prior_parser, random_parser, mkpath)


def parser():
    parser = ArgumentParser(parents=[waveform_parser, prior_parser,
                                     random_parser])
    parser.add_argument(
        'input', metavar='INPUT.{hdf,xml,xml.gz,sqlite}', default='-',
        nargs='+', type=FileType('rb'),
        help='Input LIGO-LW XML file, SQLite file, or PyCBC HDF5 files. '
             'For PyCBC, you must supply the coincidence file '
             '(e.g. "H1L1-HDFINJFIND.hdf" or "H1L1-STATMAP.hdf"), '
             'the template bank file (e.g. H1L1-BANK2HDF.hdf), '
             'the single-detector merged PSD files '
             '(e.g. "H1-MERGE_PSDS.hdf" and "L1-MERGE_PSDS.hdf"), '
             'and the single-detector merged trigger files '
             '(e.g. "H1-HDF_TRIGGER_MERGE.hdf" and '
             '"L1-HDF_TRIGGER_MERGE.hdf"), '
             'in any order.')
    parser.add_argument(
        '--pycbc-sample', default='foreground',
        help='(PyCBC only) sample population')
    parser.add_argument(
        '--coinc-event-id', type=int, nargs='*',
        help='run on only these specified events')
    parser.add_argument(
        '--output', '-o', default='.',
        help='output directory')
    parser.add_argument(
        '--condor-submit', action='store_true',
        help='submit to Condor instead of running locally')

    group = parser.add_argument_group(
        'fixed parameter options',
        'Options to hold certain parameters constant')
    group.add_argument('--ra', type=float, metavar='DEG',
                       help='Right ascension')
    group.add_argument('--dec', type=float, metavar='DEG',
                       help='Declination')
    group.add_argument('--distance', type=float, metavar='Mpc',
                       help='Luminosity distance')

    return parser


def identity(x):
    return x


def main(args=None):
    opts = parser().parse_args(args)

    import logging
    log = logging.getLogger('BAYESTAR')

    # BAYESTAR imports.
    from ..io import events, hdf5
    from ..bayestar import condition, condition_prior, ez_emcee, log_post

    # Other imports.
    from astropy.table import Table
    import numpy as np
    import os
    from collections import OrderedDict
    import subprocess
    import sys

    # Squelch annoying and uniformative LAL log messages.
    import lal
    lal.ClobberDebugLevel(lal.LALNDEBUG)

    # Read coinc file.
    log.info(
        '%s:reading input files', ','.join(file.name for file in opts.input))
    event_source = events.open(*opts.input, sample=opts.pycbc_sample)

    mkpath(opts.output)

    if opts.condor_submit:
        if opts.seed is not None:
            raise NotImplementedError(
                '--seed does not yet work with --condor-submit')
        if opts.coinc_event_id:
            raise ValueError(
                'must not set --coinc-event-id with --condor-submit')
        with subprocess.Popen(['condor_submit'],
                              # FIXME: use text=True instead in Python >= 3.7
                              encoding=sys.stdin.encoding,
                              stdin=subprocess.PIPE) as proc:
            f = proc.stdin
            print('''
                  accounting_group = ligo.dev.o3.cbc.pe.bayestar
                  on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
                  on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
                  on_exit_hold_reason = (ExitBySignal == True \
                    ? strcat("The job exited with signal ", ExitSignal) \
                    : strcat("The job exited with code ", ExitCode))
                  request_memory = 1000 MB
                  universe = vanilla
                  getenv = true
                  executable = /usr/bin/env
                  JobBatchName = BAYESTAR
                  environment = "OMP_NUM_THREADS=1"
                  ''', file=f)
            print('error =', os.path.join(opts.output, '$(cid).err'), file=f)
            print('log =', os.path.join(opts.output, '$(cid).log'), file=f)
            print('arguments = "',
                  *(arg for arg in sys.argv if arg != '--condor-submit'),
                  '--coinc-event-id $(cid)"', file=f)
            print('queue cid in', *event_source, file=f)
        sys.exit(proc.returncode)

    if opts.coinc_event_id:
        event_source = OrderedDict(
            (key, event_source[key]) for key in opts.coinc_event_id)

    # Loop over all sngl_inspiral <-> sngl_inspiral coincs.
    for int_coinc_event_id, event in event_source.items():
        coinc_event_id = 'coinc_event:coinc_event_id:{}'.format(
            int_coinc_event_id)

        log.info('%s:preparing', coinc_event_id)

        epoch, sample_rate, epochs, snrs, responses, locations, horizons = \
            condition(event, waveform=opts.waveform, f_low=opts.f_low,
                      enable_snr_series=opts.enable_snr_series,
                      f_high_truncate=opts.f_high_truncate)

        min_distance, max_distance, prior_distance_power, cosmology = \
            condition_prior(horizons, opts.min_distance, opts.max_distance,
                            opts.prior_distance_power, opts.cosmology)

        gmst = lal.GreenwichMeanSiderealTime(epoch)

        max_abs_t = 2 * snrs.data.shape[1] / sample_rate
        xmin = [0, -1, min_distance, -1, 0, 0]
        xmax = [2 * np.pi, 1, max_distance, 1, 2 * np.pi, 2 * max_abs_t]
        names = 'ra dec distance inclination twopsi time'.split()
        transformed_names = 'ra sin_dec distance u twopsi time'.split()
        forward_transforms = [identity, np.sin, identity,
                              np.cos, identity, identity]
        reverse_transforms = [identity, np.arcsin, identity,
                              np.arccos, identity, identity]
        kwargs = dict(min_distance=min_distance, max_distance=max_distance,
                      prior_distance_power=prior_distance_power,
                      cosmology=cosmology, gmst=gmst, sample_rate=sample_rate,
                      epochs=epochs, snrs=snrs, responses=responses,
                      locations=locations, horizons=horizons)

        # Fix parameters
        for i, key in reversed(list(enumerate(['ra', 'dec', 'distance']))):
            value = getattr(opts, key)
            if value is None:
                continue

            if key in ['ra', 'dec']:
                # FIXME: figure out a more elegant way to address different
                # units in command line arguments and posterior samples
                value = np.deg2rad(value)

            kwargs[transformed_names[i]] = forward_transforms[i](value)
            del (xmin[i], xmax[i], names[i], transformed_names[i],
                 forward_transforms[i], reverse_transforms[i])

        log.info('%s:sampling', coinc_event_id)

        # Run MCMC
        chain = ez_emcee(log_post, xmin, xmax, kwargs=kwargs, vectorize=True)

        # Transform back from sin_dec to dec and cos_inclination to inclination
        for i, func in enumerate(reverse_transforms):
            chain[:, i] = func(chain[:, i])

        # Create Astropy table
        chain = Table(rows=chain, names=names)

        log.info('%s:saving posterior samples', coinc_event_id)

        hdf5.write_samples(
            chain,
            os.path.join(opts.output, '{}.hdf5'.format(int_coinc_event_id)),
            path='/bayestar/posterior_samples', overwrite=True)
