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
"""Synthesize triggers for simulated sources by realizing Gaussian measurement
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

import copy
import functools

import lal
import numpy as np

from . import (
    ArgumentParser, EnableAction, FileType, get_random_parser,
    register_to_xmldoc, write_fileobj)


def parser():
    # Determine list of known detectors for command line arguments.
    available_ifos = sorted(det.frDetector.prefix
                            for det in lal.CachedDetectors)

    parser = ArgumentParser(parents=[get_random_parser()])
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
        '--enable-snr-series', action=EnableAction, default=True,
        help='Enable output of SNR time series')
    parser.add_argument(
        '--reference-psd', metavar='PSD.xml[.gz]', type=FileType('rb'),
        required=True, help='Name of PSD file')
    parser.add_argument(
        '--f-low', type=float,
        help='Override low frequency cutoff found in sim_inspiral table')
    parser.add_argument(
        '--f-high', type=float,
        help='Set high frequency cutoff to simulate early warning')
    parser.add_argument(
        '--duty-cycle', type=float, default=1.0,
        help='Single-detector duty cycle')
    parser.add_argument(
        '-P', '--preserve-ids', action='store_true',
        help='Preserve original simulation_ids')
    return parser


def simulate_snr(ra, dec, psi, inc, distance, epoch, gmst, H, S,
                 response, location, measurement_error='zero-noise',
                 duration=0.1):
    from scipy.interpolate import interp1d

    from ..bayestar import filter
    from ..bayestar.interpolation import interpolate_max

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
            N, len(snr_series) / sample_rate, normalize=False)
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


def simulate(seed, sim_inspiral, psds, responses, locations, measurement_error,
             f_low=None, f_high=None, waveform=None):
    from ..bayestar import filter

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

    f_low = f_low or sim_inspiral.f_lower
    waveform = waveform or sim_inspiral.waveform

    # Signal models for each detector.
    H = filter.sngl_inspiral_psd(
        waveform,
        mass1=sim_inspiral.mass1,
        mass2=sim_inspiral.mass2,
        spin1x=sim_inspiral.spin1x,
        spin1y=sim_inspiral.spin1y,
        spin1z=sim_inspiral.spin1z,
        spin2x=sim_inspiral.spin2x,
        spin2y=sim_inspiral.spin2y,
        spin2z=sim_inspiral.spin2z,
        f_min=f_low,
        f_final=f_high)

    return [
        simulate_snr(
            ra, dec, psi, inc, DL, epoch, gmst, H, S, response, location,
            measurement_error=measurement_error)
        for S, response, location in zip(psds, responses, locations)]


def main(args=None):
    p = parser()
    with p.parse_args(args) as opts:
        # LIGO-LW XML imports.
        from ligo.lw import ligolw
        from ligo.lw.param import Param
        from ligo.lw.utils.search_summary import append_search_summary
        from ligo.lw import utils as ligolw_utils
        from ligo.lw.lsctables import (
            New, CoincDefTable, CoincID, CoincInspiralTable, CoincMapTable,
            CoincTable, ProcessParamsTable, ProcessTable, SimInspiralTable,
            SnglInspiralTable, TimeSlideTable)

        # glue, LAL and pylal imports.
        from ligo import segments
        import lal
        import lal.series
        import lalsimulation
        from lalinspiral.inspinjfind import InspiralSCExactCoincDef
        from lalinspiral.thinca import InspiralCoincDef
        from tqdm import tqdm

        # BAYESTAR imports.
        from ..io.events.ligolw import ContentHandler
        from ..bayestar import filter
        from ..util.progress import progress_map

        # Read PSDs.
        xmldoc = ligolw_utils.load_fileobj(
            opts.reference_psd, contenthandler=lal.series.PSDContentHandler)
        psds = lal.series.read_psd_xmldoc(xmldoc, root_name=None)
        psds = {
            key: filter.InterpolatedPSD(filter.abscissa(psd), psd.data.data)
            for key, psd in psds.items() if psd is not None}
        psds = [psds[ifo] for ifo in opts.detector]

        # Extract simulation table from injection file.
        inj_xmldoc = ligolw_utils.load_fileobj(
            opts.input, contenthandler=ContentHandler)
        orig_sim_inspiral_table = SimInspiralTable.get_table(inj_xmldoc)

        # Prune injections that are outside distance limits.
        orig_sim_inspiral_table[:] = [
            row for row in orig_sim_inspiral_table
            if opts.min_distance <= row.distance <= opts.max_distance]

        # Open output file.
        xmldoc = ligolw.Document()
        xmlroot = xmldoc.appendChild(ligolw.LIGO_LW())

        # Create tables. Process and ProcessParams tables are copied from the
        # injection file.
        coinc_def_table = xmlroot.appendChild(New(CoincDefTable))
        coinc_inspiral_table = xmlroot.appendChild(New(CoincInspiralTable))
        coinc_map_table = xmlroot.appendChild(New(CoincMapTable))
        coinc_table = xmlroot.appendChild(New(CoincTable))
        xmlroot.appendChild(ProcessParamsTable.get_table(inj_xmldoc))
        xmlroot.appendChild(ProcessTable.get_table(inj_xmldoc))
        sim_inspiral_table = xmlroot.appendChild(New(SimInspiralTable))
        sngl_inspiral_table = xmlroot.appendChild(New(SnglInspiralTable))
        time_slide_table = xmlroot.appendChild(New(TimeSlideTable))

        # Write process metadata to output file.
        process = register_to_xmldoc(
            xmldoc, p, opts, instruments=opts.detector,
            comment="Simulated coincidences")

        # Add search summary to output file.
        all_time = segments.segment([lal.LIGOTimeGPS(0), lal.LIGOTimeGPS(2e9)])
        append_search_summary(xmldoc, process, inseg=all_time, outseg=all_time)

        # Create a time slide entry.  Needed for coinc_event rows.
        time_slide_id = time_slide_table.get_time_slide_id(
            {ifo: 0 for ifo in opts.detector}, create_new=process)

        # Populate CoincDef table.
        inspiral_coinc_def = copy.copy(InspiralCoincDef)
        inspiral_coinc_def.coinc_def_id = coinc_def_table.get_next_id()
        coinc_def_table.append(inspiral_coinc_def)
        found_coinc_def = copy.copy(InspiralSCExactCoincDef)
        found_coinc_def.coinc_def_id = coinc_def_table.get_next_id()
        coinc_def_table.append(found_coinc_def)

        # Precompute values that are common to all simulations.
        detectors = [lalsimulation.DetectorPrefixToLALDetector(ifo)
                     for ifo in opts.detector]
        responses = [det.response for det in detectors]
        locations = [det.location for det in detectors]

        func = functools.partial(simulate, psds=psds,
                                 responses=responses, locations=locations,
                                 measurement_error=opts.measurement_error,
                                 f_low=opts.f_low, f_high=opts.f_high,
                                 waveform=opts.waveform)

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

        with tqdm(desc='accepted') as progress:
            for sim_inspiral, simulation in zip(
                    orig_sim_inspiral_table,
                    progress_map(
                        func,
                        np.arange(len(orig_sim_inspiral_table)) + seed + 1,
                        orig_sim_inspiral_table, jobs=opts.jobs)):

                sngl_inspirals = []
                used_snr_series = []
                net_snr = 0.0
                count_triggers = 0

                # Loop over individual detectors and create SnglInspiral
                # entries.
                for ifo, (horizon, abs_snr, arg_snr, toa, series) \
                        in zip(opts.detector, simulation):

                    if np.random.uniform() > opts.duty_cycle:
                        continue
                    elif abs_snr >= opts.snr_threshold:
                        # If SNR < threshold, then the injection is not found.
                        # Skip it.
                        count_triggers += 1
                        net_snr += np.square(abs_snr)
                    elif not opts.keep_subthreshold:
                        continue

                    # Create SnglInspiral entry.
                    used_snr_series.append(series)
                    sngl_inspirals.append(
                        sngl_inspiral_table.RowType(**dict(
                            dict.fromkeys(
                                sngl_inspiral_table.validcolumns, None),
                            process_id=process.process_id,
                            ifo=ifo,
                            mass1=sim_inspiral.mass1,
                            mass2=sim_inspiral.mass2,
                            spin1x=sim_inspiral.spin1x,
                            spin1y=sim_inspiral.spin1y,
                            spin1z=sim_inspiral.spin1z,
                            spin2x=sim_inspiral.spin2x,
                            spin2y=sim_inspiral.spin2y,
                            spin2z=sim_inspiral.spin2z,
                            end=toa,
                            snr=abs_snr,
                            coa_phase=arg_snr,
                            f_final=opts.f_high,
                            eff_distance=horizon / abs_snr)))

                net_snr = np.sqrt(net_snr)

                # If too few triggers were found, then skip this event.
                if count_triggers < opts.min_triggers:
                    continue

                # If network SNR < threshold, then the injection is not found.
                # Skip it.
                if net_snr < opts.net_snr_threshold:
                    continue

                # Add Coinc table entry.
                coinc = coinc_table.appendRow(
                    coinc_event_id=coinc_table.get_next_id(),
                    process_id=process.process_id,
                    coinc_def_id=inspiral_coinc_def.coinc_def_id,
                    time_slide_id=time_slide_id,
                    insts=opts.detector,
                    nevents=len(opts.detector),
                    likelihood=None)

                # Add CoincInspiral table entry.
                coinc_inspiral_table.appendRow(
                    coinc_event_id=coinc.coinc_event_id,
                    instruments=[
                        sngl_inspiral.ifo for sngl_inspiral in sngl_inspirals],
                    end=lal.LIGOTimeGPS(1e-9 * np.mean([
                        sngl_inspiral.end.ns()
                        for sngl_inspiral in sngl_inspirals
                        if sngl_inspiral.end is not None])),
                    mass=sim_inspiral.mass1 + sim_inspiral.mass2,
                    mchirp=sim_inspiral.mchirp,
                    combined_far=0.0,  # Not provided
                    false_alarm_rate=0.0,  # Not provided
                    minimum_duration=None,  # Not provided
                    snr=net_snr)

                # Record all sngl_inspiral records and associate them with
                # coincs.
                for sngl_inspiral, series in zip(
                        sngl_inspirals, used_snr_series):
                    # Give this sngl_inspiral record an id and add it to the
                    # table.
                    sngl_inspiral.event_id = sngl_inspiral_table.get_next_id()
                    sngl_inspiral_table.append(sngl_inspiral)

                    if opts.enable_snr_series:
                        elem = lal.series.build_COMPLEX8TimeSeries(series)
                        elem.appendChild(
                            Param.from_pyvalue(
                                'event_id', sngl_inspiral.event_id))
                        xmlroot.appendChild(elem)

                    # Add CoincMap entry.
                    coinc_map_table.appendRow(
                        coinc_event_id=coinc.coinc_event_id,
                        table_name=sngl_inspiral_table.tableName,
                        event_id=sngl_inspiral.event_id)

                # Record injection
                if not opts.preserve_ids:
                    sim_inspiral.simulation_id \
                        = sim_inspiral_table.get_next_id()
                sim_inspiral_table.append(sim_inspiral)

                progress.update()

        # Record coincidence associating injections with events.
        for i, sim_inspiral in enumerate(sim_inspiral_table):
            coinc = coinc_table.appendRow(
                coinc_event_id=coinc_table.get_next_id(),
                process_id=process.process_id,
                coinc_def_id=found_coinc_def.coinc_def_id,
                time_slide_id=time_slide_id,
                instruments=None,
                nevents=None,
                likelihood=None)
            coinc_map_table.appendRow(
                coinc_event_id=coinc.coinc_event_id,
                table_name=sim_inspiral_table.tableName,
                event_id=sim_inspiral.simulation_id)
            coinc_map_table.appendRow(
                coinc_event_id=coinc.coinc_event_id,
                table_name=coinc_table.tableName,
                event_id=CoincID(i))

        # Record process end time.
        process.set_end_time_now()

        # Write output file.
        write_fileobj(xmldoc, opts.output)
