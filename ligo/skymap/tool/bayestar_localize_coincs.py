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
Produce GW sky maps for all coincidences in search pipeline output database
in LIGO-LW XML, LIGO-LW SQLite, or PyCBC HDF5 format.

The distance prior is controlled by the ``--prior-distance-power`` argument.
If you set ``--prior-distance-power`` to k, then the distance prior is
proportional to r^k. The default is 2, uniform in volume.

If the ``--min-distance`` argument is omitted, it defaults to zero. If the
``--max-distance argument`` is omitted, it defaults to the SNR=4 horizon
distance of the most sensitive detector.

A FITS file is created for each sky map, having a filename of the form
``X.fits`` where X is the integer LIGO-LW row ID of the coinc. The ``OBJECT``
card in the FITS header is also set to the integer row ID.
"""

from . import (
    ArgumentParser, FileType, get_waveform_parser, get_posterior_parser,
    get_mcmc_parser, get_random_parser)


ROW_ID_COMMENT = [
    '',
    'The integer value in the OBJECT card in this FITS header is a row ID',
    'that refers to a coinc_event table row in the input LIGO-LW document.',
    '']


def parser():
    parser = ArgumentParser(
        parents=[get_waveform_parser(), get_posterior_parser(),
                 get_mcmc_parser(), get_random_parser()])
    parser.add_argument(
        '-d', '--disable-detector', metavar='X1', type=str, nargs='+',
        help='disable certain detectors')
    parser.add_argument(
        '--keep-going', '-k', default=False, action='store_true',
        help='Keep processing events if a sky map fails to converge')
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
    return parser


def main(args=None):
    with parser().parse_args(args) as opts:
        import logging
        log = logging.getLogger('BAYESTAR')

        # BAYESTAR imports.
        from .. import omp
        from ..io import fits, events
        from ..bayestar import localize

        # Other imports.
        import os
        import subprocess
        import sys

        log.info('Using %d OpenMP thread(s)', omp.num_threads)

        # Read coinc file.
        log.info(
            '%s:reading input files', ','.join(
                file.name for file in opts.input))
        event_source = events.open(*opts.input, sample=opts.pycbc_sample)

        if opts.disable_detector:
            event_source = events.detector_disabled.open(
                event_source, opts.disable_detector)

        os.makedirs(opts.output, exist_ok=True)

        if opts.condor_submit:
            if opts.seed is not None:
                raise NotImplementedError(
                    '--seed does not yet work with --condor-submit')
            if opts.coinc_event_id:
                raise ValueError(
                    'must not set --coinc-event-id with --condor-submit')
            with subprocess.Popen(['condor_submit'],
                                  text=True, stdin=subprocess.PIPE) as proc:
                f = proc.stdin
                print('''
                    accounting_group = ligo.dev.o4.cbc.pe.bayestar
                    on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)
                    on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
                    on_exit_hold_reason = (ExitBySignal == True \
                        ? strcat("The job exited with signal ", ExitSignal) \
                        : strcat("The job exited with code ", ExitCode))
                    request_memory = 2000 MB
                    request_disk = 100 MB
                    universe = vanilla
                    getenv = true
                    executable = /usr/bin/env
                    JobBatchName = BAYESTAR
                    environment = "OMP_NUM_THREADS=1"
                    ''', file=f)
                print(
                    'error =', os.path.join(opts.output, '$(cid).err'), file=f)
                print(
                    'arguments = "',
                    *(arg for arg in sys.argv if arg != '--condor-submit'),
                    '--coinc-event-id $(cid)"', file=f)
                print('queue cid in', *event_source, file=f)
            sys.exit(proc.returncode)

        if opts.coinc_event_id:
            event_source = {
                key: event_source[key] for key in opts.coinc_event_id}

        count_sky_maps_failed = 0

        # Loop over all sngl_inspiral <-> sngl_inspiral coincs.
        for coinc_event_id, event in event_source.items():
            # Loop over sky localization methods
            log.info('%d:computing sky map', coinc_event_id)
            if opts.chain_dump:
                chain_dump = f'{coinc_event_id}.hdf5'
            else:
                chain_dump = None
            try:
                sky_map = localize(
                    event, opts.waveform, opts.f_low, opts.min_distance,
                    opts.max_distance, opts.prior_distance_power,
                    opts.cosmology, mcmc=opts.mcmc, chain_dump=chain_dump,
                    enable_snr_series=opts.enable_snr_series,
                    f_high_truncate=opts.f_high_truncate,
                    rescale_loglikelihood=opts.rescale_loglikelihood)
                sky_map.meta['objid'] = coinc_event_id
                sky_map.meta['comment'] = ROW_ID_COMMENT

            except (ArithmeticError, ValueError):
                log.exception('%d:sky localization failed', coinc_event_id)
                count_sky_maps_failed += 1
                if not opts.keep_going:
                    raise
            else:
                log.info('%d:saving sky map', coinc_event_id)
                filename = f'{coinc_event_id}.fits'
                fits.write_sky_map(
                    os.path.join(opts.output, filename), sky_map, nest=True)

        if count_sky_maps_failed > 0:
            raise RuntimeError("{0} sky map{1} did not converge".format(
                count_sky_maps_failed,
                's' if count_sky_maps_failed > 1 else ''))
