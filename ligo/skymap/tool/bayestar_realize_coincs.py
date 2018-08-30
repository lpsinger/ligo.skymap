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
Synthesize triggers for simulated sources by realizing Gaussian measurement
errors in SNR and time of arrival. The input file (or stdin if the input file
is omitted) should be an optionally gzip-compressed LIGO-LW XML file of the
form produced by lalapps_inspinj. The output file (or stdout if omitted) will
be an optionally gzip-compressed LIGO-LW XML file containing single-detector
triggers and coincidences.

The root-mean square measurement error depends on the SNR of the signal, so
there is a choice for how to generate perturbed time and phase measurements:

 - `zero-noise`: no measurement error at all
 - `gaussian-noise`: measurement error for a matched filter in Gaussian noise
"""

from argparse import FileType
import functools
import os

import lal
from ..util import lal as lal_pickle_support  # noqa: F401
import numpy as np

from . import ArgumentParser, EnableAction, random_parser, register_to_xmldoc


def parser():
    # Determine list of known detectors for command line arguments.
    available_ifos = sorted(det.frDetector.prefix
                            for det in lal.CachedDetectors)

    parser = ArgumentParser(parents=[random_parser])
    parser.add_argument(
        'input', metavar='IN.xml[.gz]', type=FileType('rb'),
        default='-', help='Name of input file')
    parser.add_argument(
        '-o', '--output', metavar='OUT.xml[.gz]', type=FileType('wb'),
        default='-', help='Name of output file')
    parser.add_argument(
        '-j', '--jobs', type=int, default=1, const=None, nargs='?',
        help='Number of threads')
    parser.add_argument(
        '--detector', metavar='|'.join(available_ifos), nargs='+',
        help='Detectors to use', choices=available_ifos, required=True)
    parser.add_argument(
        '--waveform',
        help='Waveform to use for injections (overrides values in '
        'sim_inspiral table)')
    parser.add_argument(
        '--snr-threshold', type=float, default=4.,
        help='Single-detector SNR threshold')
    parser.add_argument(
        '--net-snr-threshold', type=float, default=12.,
        help='Network SNR threshold')
    parser.add_argument(
        '--keep-subthreshold', action='store_true',
        help='Keep sub-threshold triggers that do not contribute to the '
        'network SNR')
    parser.add_argument(
        '--min-triggers', type=int, default=2,
        help='Emit coincidences only when at least this many triggers '
        'are found')
    parser.add_argument(
        '--min-distance', type=float, default=0.0,
        help='Skip events with distance less than this value')
    parser.add_argument(
        '--max-distance', type=float, default=float('inf'),
        help='Skip events with distance greater than this value')
    parser.add_argument(
        '--measurement-error', default='zero-noise',
        choices=('zero-noise', 'gaussian-noise'),
        help='How to compute the measurement error')
    parser.add_argument(
        '--enable-snr-series', action=EnableAction,
        help='Enable output of SNR time series')
    parser.add_argument(
        '--reference-psd', metavar='PSD.xml[.gz]', type=FileType('rb'),
        required=True, help='Name of PSD file')
    parser.add_argument(
        '--f-low', type=float,
        help='Override low frequency cutoff found in sim_inspiral table')
    parser.add_argument(
        '--duty-cycle', type=float, default=1.0,
        help='Single-detector duty cycle')
    return parser


def simulate_snr(ra, dec, psi, inc, distance, epoch, gmst, H, S,
                 response, location, measurement_error='zero-noise'):
    from scipy.interpolate import interp1d

    from ..bayestar import filter
    from ..bayestar.interpolation import interpolate_max

    duration = 0.1

    # Calculate whitened template autocorrelation sequence.
    HS = filter.signal_psd_series(H, S)
    n = len(HS.data.data)
    acor, sample_rate = filter.autocorrelation(HS, duration)

    # Calculate time, amplitude, and phase.
    u = np.cos(inc)
    u2 = np.square(u)
    signal_model = filter.SignalModel(HS)
    horizon = signal_model.get_horizon_distance()
    Fplus, Fcross = lal.ComputeDetAMResponse(response, ra, dec, psi, gmst)
    toa = lal.TimeDelayFromEarthCenter(location, ra, dec, epoch)
    z = (0.5 * (1 + u2) * Fplus + 1j * u * Fcross) * horizon / distance

    # Calculate complex autocorrelation sequence.
    snr_series = z * np.concatenate((acor[:0:-1].conj(), acor))

    # If requested, add noise.
    if measurement_error == 'gaussian-noise':
        sigmasq = 4 * np.sum(HS.deltaF * np.abs(HS.data.data))
        amp = 4 * n * HS.deltaF**0.5 * np.sqrt(HS.data.data / sigmasq)
        N = lal.CreateCOMPLEX16FrequencySeries(
            '', HS.epoch, HS.f0, HS.deltaF, HS.sampleUnits, n)
        N.data.data = amp * (
            np.random.randn(n) + 1j * np.random.randn(n))
        noise_term, sample_rate_2 = filter.autocorrelation(
            N, 2 * duration - 1 / sample_rate, normalize=False)
        assert sample_rate == sample_rate_2
        snr_series += noise_term

    # Shift SNR series to the nearest sample.
    int_samples, frac_samples = divmod(
        (1e-9 * epoch.gpsNanoSeconds + toa) * sample_rate, 1)
    if frac_samples > 0.5:
        int_samples += 1
        frac_samples -= 1
    epoch = lal.LIGOTimeGPS(epoch.gpsSeconds, 0)
    n = len(acor) - 1
    mprime = np.arange(-n, n + 1)
    m = mprime + frac_samples
    re, im = (
        interp1d(m, x, kind='cubic', bounds_error=False, fill_value=0)(mprime)
        for x in (snr_series.real, snr_series.imag))
    snr_series = re + 1j * im

    # Find the trigger values.
    i_nearest = np.argmax(np.abs(snr_series[n-n//2:n+n//2+1])) + n-n//2
    i_interp, z_interp = interpolate_max(i_nearest, snr_series,
                                         n // 2, method='lanczos')
    toa = epoch + (int_samples + i_interp - n) / sample_rate
    snr = np.abs(z_interp)
    phase = np.angle(z_interp)

    # Shift and truncate the SNR time series.
    epoch += (int_samples + i_nearest - n - n // 2) / sample_rate
    snr_series = snr_series[(i_nearest - n // 2):(i_nearest + n // 2 + 1)]
    tseries = lal.CreateCOMPLEX8TimeSeries(
        'snr', epoch, 0, 1 / sample_rate,
        lal.DimensionlessUnit, len(snr_series))
    tseries.data.data = snr_series
    return horizon, snr, phase, toa, tseries


def simulate(seed_sim_inspiral, psds, responses, locations, measurement_error):
    from ..bayestar import filter

    seed, sim_inspiral = seed_sim_inspiral
    np.random.seed(seed)

    # Unpack some values from the row in the table.
    DL = sim_inspiral.distance
    ra = sim_inspiral.longitude
    dec = sim_inspiral.latitude
    inc = sim_inspiral.inclination
    # phi = sim_inspiral.coa_phase  # arbitrary?
    psi = sim_inspiral.polarization
    epoch = sim_inspiral.time_geocent
    gmst = lal.GreenwichMeanSiderealTime(epoch)

    # Signal models for each detector.
    H = filter.sngl_inspiral_psd(
        sim_inspiral.waveform,
        mass1=sim_inspiral.mass1,
        mass2=sim_inspiral.mass2,
        spin1x=sim_inspiral.spin1x,
        spin1y=sim_inspiral.spin1y,
        spin1z=sim_inspiral.spin1z,
        spin2x=sim_inspiral.spin2x,
        spin2y=sim_inspiral.spin2y,
        spin2z=sim_inspiral.spin2z,
        f_min=sim_inspiral.f_lower)

    return [
        simulate_snr(
            ra, dec, psi, inc, DL, epoch, gmst, H, S, response, location,
            measurement_error=measurement_error)
        for S, response, location in zip(psds, responses, locations)]


def main(args=None):
    p = parser()
    opts = p.parse_args(args)

    # LIGO-LW XML imports.
    from glue.ligolw import ligolw
    from glue.ligolw import param as ligolw_param
    from glue.ligolw.utils import process as ligolw_process
    from glue.ligolw.utils import search_summary as ligolw_search_summary
    from glue.ligolw import utils as ligolw_utils
    from glue.ligolw import lsctables

    # glue, LAL and pylal imports.
    from glue import segments
    import glue.lal
    import lal.series
    import lalsimulation
    from lalinspiral.thinca import InspiralCoincDef
    from tqdm import tqdm

    # FIXME: disable progress bar monitor thread.
    #
    # I was getting error messages that look like this:
    #
    # Traceback (most recent call last):
    #   File "/tqdm/_tqdm.py", line 885, in __del__
    #     self.close()
    #   File "/tqdm/_tqdm.py", line 1090, in close
    #     self._decr_instances(self)
    #   File "/tqdm/_tqdm.py", line 454, in _decr_instances
    #     cls.monitor.exit()
    #   File "/tqdm/_monitor.py", line 52, in exit
    #     self.join()
    #   File "/usr/lib64/python3.6/threading.py", line 1053, in join
    #     raise RuntimeError("cannot join current thread")
    # RuntimeError: cannot join current thread
    #
    # I don't know what causes this... maybe a race condition in tqdm's cleanup
    # code. Anyway, this should disable the tqdm monitor thread entirely.
    tqdm.monitor_interval = 0

    # BAYESTAR imports.
    from ..io.events.ligolw import ContentHandler
    from ..bayestar import filter

    # Open output file.
    out_xmldoc = ligolw.Document()
    out_xmldoc.appendChild(ligolw.LIGO_LW())

    # Write process metadata to output file.
    process = register_to_xmldoc(
        out_xmldoc, p, opts, ifos=opts.detector,
        comment="Simulated coincidences")

    # Add search summary to output file.
    all_time = segments.segment(
        [glue.lal.LIGOTimeGPS(0), glue.lal.LIGOTimeGPS(2e9)])
    search_summary_table = lsctables.New(lsctables.SearchSummaryTable)
    out_xmldoc.childNodes[0].appendChild(search_summary_table)
    ligolw_search_summary.append_search_summary(
        out_xmldoc, process, inseg=all_time, outseg=all_time)

    # Read PSDs.
    xmldoc, _ = ligolw_utils.load_fileobj(
        opts.reference_psd, contenthandler=lal.series.PSDContentHandler)
    psds = lal.series.read_psd_xmldoc(xmldoc, root_name=None)
    psds = {
        key: filter.InterpolatedPSD(filter.abscissa(psd), psd.data.data)
        for key, psd in psds.items() if psd is not None}
    psds = [psds[ifo] for ifo in opts.detector]

    # Read injection file.
    xmldoc, _ = ligolw_utils.load_fileobj(
        opts.input, contenthandler=ContentHandler)

    # Extract simulation table from injection file.
    sim_inspiral_table = lsctables.SimInspiralTable.get_table(xmldoc)

    # Create a SnglInspiral table and initialize its row ID counter.
    sngl_inspiral_table = lsctables.New(lsctables.SnglInspiralTable)
    out_xmldoc.childNodes[0].appendChild(sngl_inspiral_table)
    sngl_inspiral_table.set_next_id(lsctables.SnglInspiralID(0))

    # Create a time slide entry.  Needed for coinc_event rows.
    time_slide_table = lsctables.New(lsctables.TimeSlideTable)
    out_xmldoc.childNodes[0].appendChild(time_slide_table)
    time_slide_id = time_slide_table.get_time_slide_id(
        {ifo: 0 for ifo in opts.detector}, create_new=process)

    # Create a CoincDef table and record a CoincDef row for
    # sngl_inspiral <-> sngl_inspiral coincidences.
    coinc_def_table = lsctables.New(lsctables.CoincDefTable)
    out_xmldoc.childNodes[0].appendChild(coinc_def_table)
    coinc_def = InspiralCoincDef
    coinc_def_id = coinc_def_table.get_next_id()
    coinc_def.coinc_def_id = coinc_def_id
    coinc_def_table.append(coinc_def)

    # Create a CoincMap table.
    coinc_map_table = lsctables.New(lsctables.CoincMapTable)
    out_xmldoc.childNodes[0].appendChild(coinc_map_table)

    # Create a CoincEvent table.
    coinc_table = lsctables.New(lsctables.CoincTable)
    out_xmldoc.childNodes[0].appendChild(coinc_table)

    # Create a CoincInspiral table.
    coinc_inspiral_table = lsctables.New(lsctables.CoincInspiralTable)
    out_xmldoc.childNodes[0].appendChild(coinc_inspiral_table)

    # Precompute values that are common to all simulations.
    detectors = [lalsimulation.DetectorPrefixToLALDetector(ifo)
                 for ifo in opts.detector]
    responses = [det.response for det in detectors]
    locations = [det.location for det in detectors]

    # Fix up sim_inspiral table with values from command line options.
    sim_inspiral_table[:] = [
        row for row in sim_inspiral_table
        if opts.min_distance <= row.distance <= opts.max_distance]
    if opts.f_low is not None:
        for row in sim_inspiral_table:
            row.f_lower = opts.f_low
    if opts.waveform is not None:
        for row in sim_inspiral_table:
            row.waveform = opts.waveform
    # FIXME: Set transverse spin components to 0.
    for row in sim_inspiral_table:
        row.spin1x = row.spin1y = row.spin2x = row.spin2y = 0

    if opts.jobs == 1:
        pool_map = map
    else:
        from .. import omp
        from multiprocessing import Pool
        omp.num_threads = 1  # disable OpenMP parallelism
        pool_map = Pool(opts.jobs).imap

    func = functools.partial(simulate, psds=psds,
                             responses=responses, locations=locations,
                             measurement_error=opts.measurement_error)

    # Make sure that each thread gets a different random number state.
    # We start by drawing a random integer s in the main thread, and
    # then the i'th subprocess will seed itself with the integer i + s.
    #
    # The seed must be an unsigned 32-bit integer, so if there are n
    # threads, then s must be drawn from the interval [0, 2**32 - n).
    #
    # Note that *we* are thread 0, so there are a total of
    # n=1+len(sim_inspiral_table) threads.
    seed = np.random.randint(0, 2 ** 32 - len(sim_inspiral_table) - 1)
    np.random.seed(seed)

    for sim_inspiral, simulation in zip(
            sim_inspiral_table,
            tqdm(pool_map(func,
                          zip(seed + 1 + np.arange(len(sim_inspiral_table)),
                              sim_inspiral_table)),
                 total=len(sim_inspiral_table))):

        sngl_inspirals = []
        used_snr_series = []
        net_snr = 0.0
        count_triggers = 0

        # Loop over individual detectors and create SnglInspiral entries.
        for ifo, (horizon, abs_snr, arg_snr, toa, series) \
                in zip(opts.detector, simulation):

            if np.random.uniform() > opts.duty_cycle:
                continue
            elif abs_snr >= opts.snr_threshold:
                # If SNR < threshold, then the injection is not found. Skip it.
                count_triggers += 1
                net_snr += np.square(abs_snr)
            elif not opts.keep_subthreshold:
                continue

            # Create SnglInspiral entry.
            sngl_inspiral = lsctables.SnglInspiral()
            for validcolumn in lsctables.SnglInspiralTable.validcolumns.keys():
                setattr(sngl_inspiral, validcolumn, None)
            sngl_inspiral.process_id = process.process_id
            sngl_inspiral.ifo = ifo
            sngl_inspiral.mass1 = sim_inspiral.mass1
            sngl_inspiral.mass2 = sim_inspiral.mass2
            sngl_inspiral.spin1x = sim_inspiral.spin1x
            sngl_inspiral.spin1y = sim_inspiral.spin1y
            sngl_inspiral.spin1z = sim_inspiral.spin1z
            sngl_inspiral.spin2x = sim_inspiral.spin2x
            sngl_inspiral.spin2y = sim_inspiral.spin2y
            sngl_inspiral.spin2z = sim_inspiral.spin2z
            sngl_inspiral.end = toa
            sngl_inspiral.snr = abs_snr
            sngl_inspiral.coa_phase = arg_snr
            sngl_inspiral.eff_distance = horizon / sngl_inspiral.snr
            sngl_inspirals.append(sngl_inspiral)
            used_snr_series.append(series)

        net_snr = np.sqrt(net_snr)

        # If too few triggers were found, then skip this event.
        if count_triggers < opts.min_triggers:
            continue

        # If network SNR < threshold, then the injection is not found. Skip it.
        if net_snr < opts.net_snr_threshold:
            continue

        # Add Coinc table entry.
        coinc = lsctables.Coinc()
        coinc.coinc_event_id = coinc_table.get_next_id()
        coinc.process_id = process.process_id
        coinc.coinc_def_id = coinc_def_id
        coinc.time_slide_id = time_slide_id
        coinc.set_instruments(opts.detector)
        coinc.nevents = len(opts.detector)
        coinc.likelihood = None
        coinc_table.append(coinc)

        # Add CoincInspiral table entry.
        coinc_inspiral = lsctables.CoincInspiral()
        coinc_inspiral.coinc_event_id = coinc.coinc_event_id
        coinc_inspiral.ifos = lsctables.instrumentsproperty.set(
            sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals)
        # FIXME: should only be detected sngls
        coinc_inspiral.end = lal.LIGOTimeGPS(
            sum(sngl_inspiral.end.ns() for sngl_inspiral in sngl_inspirals) //
            len(sngl_inspirals) * 1e-9)
        coinc_inspiral.mass = sim_inspiral.mass1 + sim_inspiral.mass2
        coinc_inspiral.mchirp = sim_inspiral.mchirp
        coinc_inspiral.combined_far = 0.0  # Not provided
        coinc_inspiral.false_alarm_rate = 0.0  # Not provided
        coinc_inspiral.minimum_duration = None  # Not provided
        coinc_inspiral.snr = net_snr
        coinc_inspiral_table.append(coinc_inspiral)

        # Record all sngl_inspiral records and associate them with coincs.
        for sngl_inspiral, series in zip(sngl_inspirals, used_snr_series):
            # Give this sngl_inspiral record an id and add it to the table.
            sngl_inspiral.event_id = sngl_inspiral_table.get_next_id()
            sngl_inspiral_table.append(sngl_inspiral)

            if opts.enable_snr_series:
                elem = lal.series.build_COMPLEX8TimeSeries(series)
                elem.appendChild(
                    ligolw_param.Param.from_pyvalue(
                        u'event_id', sngl_inspiral.event_id))
                out_xmldoc.childNodes[0].appendChild(elem)

            # Add CoincMap entry.
            coinc_map = lsctables.CoincMap()
            coinc_map.coinc_event_id = coinc.coinc_event_id
            coinc_map.table_name = sngl_inspiral_table.tableName
            coinc_map.event_id = sngl_inspiral.event_id
            coinc_map_table.append(coinc_map)

    # Record process end time.
    ligolw_process.set_process_end_time(process)

    # Write output file.
    with ligolw_utils.SignalsTrap():
        ligolw_utils.write_fileobj(
            out_xmldoc, opts.output,
            gz=(os.path.splitext(opts.output.name)[-1] == ".gz"))
