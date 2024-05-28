#
# Copyright (C) 2019-2024  Leo Singer
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
"""Rough-cut injection tool.

The idea is to efficiently sample events, uniformly in "sensitive volume"
(differential comoving volume divided by 1 + z), and from a distribution of
masses and spins, such that later detection cuts will not reject an excessive
number of events.

This occurs in two steps. First, we divide the intrinsic parameter space into a
very coarse 10x10x10x10 grid and calculate the maximum horizon distance in each
grid cell. Second, we directly sample injections jointly from the mass and spin
distribution and a uniform and isotropic spatial distribution with a redshift
cutoff that is piecewise constant in the masses and spins.
"""

from functools import partial

from astropy import cosmology
from astropy.cosmology._utils import vectorize_redshift_method
from astropy import units
from astropy.units import dimensionless_unscaled
import numpy as np
from scipy.integrate import quad, fixed_quad
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.ndimage import maximum_filter

from ..util import progress_map
from ..bayestar.filter import sngl_inspiral_psd
from . import (
    ArgumentParser, FileType, get_random_parser, register_to_xmldoc,
    write_fileobj)


def get_decisive_snr(snrs, min_triggers):
    """Return the SNR for the trigger that decides if an event is detectable.

    Parameters
    ----------
    snrs : list
        List of SNRs (floats).
    min_triggers : int
        Minimum number of triggers to form a coincidence.

    Returns
    -------
    decisive_snr : float

    """
    return sorted(snrs)[-min_triggers]


def lo_hi_nonzero(x):
    nonzero = np.flatnonzero(x)
    return nonzero[0], nonzero[-1]


class GWCosmo:
    """Evaluate GW distance figures of merit for a given cosmology.

    Parameters
    ----------
    cosmo : :class:`astropy.cosmology.FLRW`
        The cosmological model.

    """

    def __init__(self, cosmology):
        self.cosmo = cosmology

    def z_at_snr(self, psds, waveform, f_low, snr_threshold, min_triggers,
                 mass1, mass2, spin1z, spin2z):
        """
        Get redshift at which a waveform attains a given SNR.

        Parameters
        ----------
        psds : list
            List of :class:`lal.REAL8FrequencySeries` objects.
        waveform : str
            Waveform approximant name.
        f_low : float
            Low-frequency cutoff for template.
        snr_threshold : float
            Minimum single-detector SNR.
        min_triggers : int
            Minimum number of triggers to form a coincidence.
        params : list
            List of waveform parameters: mass1, mass2, spin1z, spin2z.

        Returns
        -------
        comoving_distance : float
            Comoving distance in Mpc.

        """
        # Construct waveform
        series = sngl_inspiral_psd(waveform, f_low=f_low,
                                   mass1=mass1, mass2=mass2,
                                   spin1z=spin1z, spin2z=spin2z)
        i_lo, i_hi = lo_hi_nonzero(series.data.data)
        log_f = np.log(series.f0 + series.deltaF * np.arange(i_lo, i_hi + 1))
        log_f_lo = log_f[0]
        log_f_hi = log_f[-1]
        num = interp1d(
            log_f, np.log(series.data.data[i_lo:i_hi + 1]),
            fill_value=-np.inf, bounds_error=False, assume_sorted=True)

        denoms = []
        for series in psds:
            i_lo, i_hi = lo_hi_nonzero(
                np.isfinite(series.data.data) & (series.data.data != 0))
            log_f = np.log(
                series.f0 + series.deltaF * np.arange(i_lo, i_hi + 1))
            denom = interp1d(
                log_f, log_f - np.log(series.data.data[i_lo:i_hi + 1]),
                fill_value=-np.inf, bounds_error=False, assume_sorted=True)
            denoms.append(denom)

        def snr_at_z(z):
            logzp1 = np.log(z + 1)
            integrand = lambda log_f: [
                np.exp(num(log_f + logzp1) + denom(log_f)) for denom in denoms]
            integrals, _ = fixed_quad(
                integrand, log_f_lo, log_f_hi - logzp1, n=1024)
            snr = get_decisive_snr(np.sqrt(4 * integrals), min_triggers)
            with np.errstate(divide='ignore'):
                snr /= self.cosmo.angular_diameter_distance(z).to_value(
                    units.Mpc)
            return snr

        def root_func(z):
            return snr_at_z(z) - snr_threshold

        return root_scalar(root_func, bracket=(0, 1e3)).root

    def get_max_z(self, psds, waveform, f_low, snr_threshold, min_triggers,
                  mass1, mass2, spin1z, spin2z, jobs=1):
        # Calculate the maximum distance on the grid.
        params = [mass1, mass2, spin1z, spin2z]
        shape = np.broadcast_shapes(*(param.shape for param in params))
        result = list(progress_map(
            partial(self.z_at_snr, psds, waveform, f_low,
                    snr_threshold, min_triggers),
            *(param.ravel() for param in params),
            jobs=jobs))
        result = np.reshape(result, shape)

        assert np.all(result >= 0), 'some redshifts are negative'
        assert np.all(np.isfinite(result)), 'some redshifts are not finite'
        return result

    @vectorize_redshift_method
    def _sensitive_volume_integral(self, z):
        dh3_sr = self.cosmo.hubble_distance**3 / units.sr

        def integrand(z):
            result = self.cosmo.differential_comoving_volume(z)
            result /= (1 + z) * dh3_sr
            return result.to_value(dimensionless_unscaled)

        result, _ = quad(integrand, 0, z)
        return result

    def sensitive_volume(self, z):
        """Sensitive volume :math:`V(z)` out to redshift :math:`z`.

        Given a population of events that occur at a constant rate density
        :math:`R` per unit comoving volume per unit proper time, the number of
        observed events out to a redshift :math:`N(z)` over an observation time
        :math:`T` is :math:`N(z) = R T V(z)`.
        """
        dh3 = self.cosmo.hubble_distance**3
        return 4 * np.pi * dh3 * self._sensitive_volume_integral(z)

    def sensitive_distance(self, z):
        r"""Sensitive distance as a function of redshift :math:`z`.

        The sensitive distance is the distance :math:`d_s(z)` defined such that
        :math:`V(z) = 4/3\pi {d_s(z)}^3`, where :math:`V(z)` is the sensitive
        volume.
        """
        dh = self.cosmo.hubble_distance
        return dh * np.cbrt(3 * self._sensitive_volume_integral(z))


def cell_max(values):
    r"""
    Find pairwise max of consecutive elements across all axes of an array.

    Parameters
    ----------
    values : :class:`numpy.ndarray`
        An input array of :math:`n` dimensions,
        :math:`(m_0, m_1, \dots, m_{n-1})`.

    Returns
    -------
    maxima : :class:`numpy.ndarray`
        An input array of :math:`n` dimensions, each with a length 1 less than
        the input array,
        :math:`(m_0 - 1, m_1 - 1, \dots, m_{n-1} - 1)`.

    """
    maxima = maximum_filter(values, size=2, mode='constant')
    indices = (slice(1, None),) * np.ndim(values)
    return maxima[indices]


def assert_not_reached():  # pragma: no cover
    raise AssertionError('This line should not be reached.')


def parser():
    parser = ArgumentParser(parents=[get_random_parser()])
    parser.add_argument(
        '--cosmology', choices=cosmology.available,
        default='Planck15', help='Cosmological model')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--distribution', help='Use a preset distribution', choices=(
            'bns_astro', 'bns_broad', 'nsbh_astro', 'nsbh_broad',
            'bbh_astro', 'bbh_broad'))
    group.add_argument(
        '--distribution-samples',
        help='Load samples of the intrinsic mass and spin distribution from '
             'any file that can be read as an Astropy table. The table '
             'columns should be mass1, mass2, spin1z, and spin2z.')

    parser.add_argument(
        '--reference-psd', type=FileType('rb'), metavar='PSD.xml[.gz]',
        required=True, help='PSD file')
    parser.add_argument(
        '--f-low', type=float, default=25.0,
        help='Low frequency cutoff in Hz')
    parser.add_argument(
        '--snr-threshold', type=float, default=4.,
        help='Single-detector SNR threshold')
    parser.add_argument(
        '--min-triggers', type=int, default=2,
        help='Emit coincidences only when at least this many triggers '
        'are found')
    parser.add_argument(
        '--min-snr', type=float,
        help='Minimum decisive SNR of injections given the reference PSDs. '
             'Deprecated; use the synonymous --snr-threshold option instead.')
    parser.add_argument(
        '--max-distance', type=float, metavar='Mpc',
        help='Maximum luminosity distance for injections')
    parser.add_argument(
        '--waveform', default='o2-uberbank',
        help='Waveform approximant')
    parser.add_argument(
        '--nsamples', type=int, default=100000,
        help='Output this many injections')
    parser.add_argument(
        '-o', '--output', type=FileType('wb'), default='-',
        metavar='INJ.xml[.gz]', help='Output file, optionally gzip-compressed')
    parser.add_argument(
        '-j', '--jobs', type=int, default=1, const=None, nargs='?',
        help='Number of threads')
    return parser


def main(args=None):
    p = parser()
    with p.parse_args(args) as args:
        import warnings

        from astropy.table import Table
        from ligo.lw import lsctables
        from ligo.lw import utils as ligolw_utils
        from ligo.lw import ligolw
        import lal.series
        from scipy import stats

        if args.min_snr is not None:
            warnings.warn(
                'The --min-snr threshold option is deprecated. '
                'Please use the synonymous --snr-threshold option instead.',
                UserWarning)
            args.snr_threshold = args.min_snr

        xmldoc = ligolw.Document()
        xmlroot = xmldoc.appendChild(ligolw.LIGO_LW())
        process = register_to_xmldoc(xmldoc, p, args)

        # Read PSDs
        psds = list(
            lal.series.read_psd_xmldoc(
                ligolw_utils.load_fileobj(
                    args.reference_psd,
                    contenthandler=lal.series.PSDContentHandler)).values())

        if len(psds) < args.min_triggers:
            parser.error(
                f'The number of PSDs ({len(psds)}) must be greater than or '
                f'equal to the value of --min-triggers ({args.min_triggers}).')

        gwcosmo = GWCosmo(getattr(cosmology, args.cosmology))

        if args.distribution:
            ns_mass_min = 1.0
            ns_mass_max = 2.0
            bh_mass_min = 5.0
            bh_mass_max = 50.0

            ns_astro_spin_min = -0.05
            ns_astro_spin_max = +0.05
            ns_astro_mass_dist = stats.norm(1.33, 0.09)
            ns_astro_spin_dist = stats.uniform(
                ns_astro_spin_min, ns_astro_spin_max - ns_astro_spin_min)

            ns_broad_spin_min = -0.4
            ns_broad_spin_max = +0.4
            ns_broad_mass_dist = stats.uniform(
                ns_mass_min, ns_mass_max - ns_mass_min)
            ns_broad_spin_dist = stats.uniform(
                ns_broad_spin_min, ns_broad_spin_max - ns_broad_spin_min)

            bh_astro_spin_min = -0.99
            bh_astro_spin_max = +0.99
            bh_astro_mass_dist = stats.pareto(b=1.3)
            bh_astro_spin_dist = stats.uniform(
                bh_astro_spin_min, bh_astro_spin_max - bh_astro_spin_min)

            bh_broad_spin_min = -0.99
            bh_broad_spin_max = +0.99
            bh_broad_mass_dist = stats.reciprocal(bh_mass_min, bh_mass_max)
            bh_broad_spin_dist = stats.uniform(
                bh_broad_spin_min, bh_broad_spin_max - bh_broad_spin_min)

            if args.distribution.startswith('bns_'):
                m1_min = m2_min = ns_mass_min
                m1_max = m2_max = ns_mass_max
                if args.distribution.endswith('_astro'):
                    x1_min = x2_min = ns_astro_spin_min
                    x1_max = x2_max = ns_astro_spin_max
                    m1_dist = m2_dist = ns_astro_mass_dist
                    x1_dist = x2_dist = ns_astro_spin_dist
                elif args.distribution.endswith('_broad'):
                    x1_min = x2_min = ns_broad_spin_min
                    x1_max = x2_max = ns_broad_spin_max
                    m1_dist = m2_dist = ns_broad_mass_dist
                    x1_dist = x2_dist = ns_broad_spin_dist
                else:  # pragma: no cover
                    assert_not_reached()
            elif args.distribution.startswith('nsbh_'):
                m1_min = bh_mass_min
                m1_max = bh_mass_max
                m2_min = ns_mass_min
                m2_max = ns_mass_max
                if args.distribution.endswith('_astro'):
                    x1_min = bh_astro_spin_min
                    x1_max = bh_astro_spin_max
                    x2_min = ns_astro_spin_min
                    x2_max = ns_astro_spin_max
                    m1_dist = bh_astro_mass_dist
                    m2_dist = ns_astro_mass_dist
                    x1_dist = bh_astro_spin_dist
                    x2_dist = ns_astro_spin_dist
                elif args.distribution.endswith('_broad'):
                    x1_min = bh_broad_spin_min
                    x1_max = bh_broad_spin_max
                    x2_min = ns_broad_spin_min
                    x2_max = ns_broad_spin_max
                    m1_dist = bh_broad_mass_dist
                    m2_dist = ns_broad_mass_dist
                    x1_dist = bh_broad_spin_dist
                    x2_dist = ns_broad_spin_dist
                else:  # pragma: no cover
                    assert_not_reached()
            elif args.distribution.startswith('bbh_'):
                m1_min = m2_min = bh_mass_min
                m1_max = m2_max = bh_mass_max
                if args.distribution.endswith('_astro'):
                    x1_min = x2_min = bh_astro_spin_min
                    x1_max = x2_max = bh_astro_spin_max
                    m1_dist = m2_dist = bh_astro_mass_dist
                    x1_dist = x2_dist = bh_astro_spin_dist
                elif args.distribution.endswith('_broad'):
                    x1_min = x2_min = bh_broad_spin_min
                    x1_max = x2_max = bh_broad_spin_max
                    m1_dist = m2_dist = bh_broad_mass_dist
                    x1_dist = x2_dist = bh_broad_spin_dist
                else:  # pragma: no cover
                    assert_not_reached()
            else:  # pragma: no cover
                assert_not_reached()

            dists = (m1_dist, m2_dist, x1_dist, x2_dist)

            # Construct mass1, mass2, spin1z, spin2z grid.
            m1 = np.geomspace(m1_min, m1_max, 10)
            m2 = np.geomspace(m2_min, m2_max, 10)
            x1 = np.linspace(x1_min, x1_max, 10)
            x2 = np.linspace(x2_min, x2_max, 10)
            params = m1, m2, x1, x2

            # Calculate the maximum distance on the grid.
            max_z = gwcosmo.get_max_z(
                psds, args.waveform, args.f_low,
                args.snr_threshold, args.min_triggers,
                *np.meshgrid(m1, m2, x1, x2, indexing='ij'), jobs=args.jobs)
            if args.max_distance is not None:
                new_max_z = cosmology.z_at_value(
                    gwcosmo.cosmo.luminosity_distance,
                    args.max_distance * units.Mpc)
                max_z[max_z > new_max_z] = new_max_z
            max_distance = gwcosmo.sensitive_distance(
                max_z).to_value(units.Mpc)

            # Find piecewise constant approximate upper bound on distance.
            max_distance = cell_max(max_distance)

            # Calculate V * T in each grid cell
            cdfs = [dist.cdf(param) for param, dist in zip(params, dists)]
            cdf_los = [cdf[:-1] for cdf in cdfs]
            cdfs = [np.diff(cdf) for cdf in cdfs]
            probs = np.prod(np.meshgrid(*cdfs, indexing='ij'), axis=0)
            probs /= probs.sum()
            probs *= 4/3*np.pi*max_distance**3
            volume = probs.sum()
            probs /= volume
            probs = probs.ravel()

            # Draw random grid cells
            dist = stats.rv_discrete(values=(np.arange(len(probs)), probs))
            indices = np.unravel_index(
                dist.rvs(size=args.nsamples), max_distance.shape)

            # Draw random intrinsic params from each cell
            cols = {}
            cols['mass1'], cols['mass2'], cols['spin1z'], cols['spin2z'] = [
                dist.ppf(
                    stats.uniform(cdf_lo[i], cdf[i]).rvs(size=args.nsamples))
                for i, dist, cdf_lo, cdf in zip(indices, dists, cdf_los, cdfs)]
        elif args.distribution_samples:
            # Load distribution samples.
            samples = Table.read(args.distribution_samples)

            # Calculate the maximum sensitive distance for each sample.
            max_z = gwcosmo.get_max_z(
                psds, args.waveform, args.f_low,
                args.snr_threshold, args.min_triggers,
                samples['mass1'], samples['mass2'],
                samples['spin1z'], samples['spin2z'], jobs=args.jobs)
            if args.max_distance is not None:
                new_max_z = cosmology.z_at_value(
                    gwcosmo.cosmo.luminosity_distance,
                    args.max_distance * units.Mpc)
                max_z[max_z > new_max_z] = new_max_z
            max_distance = gwcosmo.sensitive_distance(
                max_z).to_value(units.Mpc)

            # Calculate V * T for each sample.
            probs = 1 / len(max_distance)
            probs *= 4/3*np.pi*max_distance**3
            volume = probs.sum()
            probs /= volume

            # Draw weighted samples for the simulated events.
            dist = stats.rv_discrete(values=(np.arange(len(probs)), probs))
            # Note that we do this in small batches because
            # stats.rv_discrete.rvs has quadratic memory usage, number of
            # values times number of samples, which might cause us to run out
            # of memory if we did it all at once.
            n_batches = max(args.nsamples * len(probs) // 1_000_000_000, 1)
            batch_sizes = [len(subarray) for subarray in
                           np.array_split(np.empty(args.nsamples), n_batches)]
            indices = np.concatenate([dist.rvs(size=batch_size)
                                      for batch_size in batch_sizes])

            cols = {key: samples[key][indices]
                    for key in ['mass1', 'mass2', 'spin1z', 'spin2z']}
        else:
            assert_not_reached()

        volumetric_rate = (
            args.nsamples / volume * units.year**-1 * units.Mpc**-3)

        # Swap binary components as needed to ensure that mass1 >= mass2.
        # Note that the .copy() is important.
        # See https://github.com/numpy/numpy/issues/14428
        swap = cols['mass1'] < cols['mass2']
        cols['mass1'][swap], cols['mass2'][swap] = \
            cols['mass2'][swap].copy(), cols['mass1'][swap].copy()
        cols['spin1z'][swap], cols['spin2z'][swap] = \
            cols['spin2z'][swap].copy(), cols['spin1z'][swap].copy()

        # Draw random extrinsic parameters
        cols['distance'] = stats.powerlaw(
            a=3, scale=max_distance[indices]).rvs(size=args.nsamples)
        cols['longitude'] = stats.uniform(0, 2 * np.pi).rvs(
            size=args.nsamples)
        cols['latitude'] = np.arcsin(stats.uniform(-1, 2).rvs(
            size=args.nsamples))
        cols['inclination'] = np.arccos(stats.uniform(-1, 2).rvs(
            size=args.nsamples))
        cols['polarization'] = stats.uniform(0, 2 * np.pi).rvs(
            size=args.nsamples)
        cols['coa_phase'] = stats.uniform(-np.pi, 2 * np.pi).rvs(
            size=args.nsamples)
        cols['time_geocent'] = stats.uniform(
            1e9, units.year.to(units.second)).rvs(size=args.nsamples)

        # Convert from sensitive distance to redshift and comoving distance.
        # FIXME: Replace this brute-force lookup table with a solver.
        z = np.linspace(0, max_z.max(), 10000)
        ds = gwcosmo.sensitive_distance(z).to_value(units.Mpc)
        dc = gwcosmo.cosmo.comoving_distance(z).to_value(units.Mpc)
        z_for_ds = interp1d(ds, z, kind='cubic', assume_sorted=True)
        dc_for_ds = interp1d(ds, dc, kind='cubic', assume_sorted=True)
        zp1 = 1 + z_for_ds(cols['distance'])
        cols['distance'] = dc_for_ds(cols['distance'])

        # Apply redshift factor to convert from comoving distance and source
        # frame masses to luminosity distance and observer frame masses.
        for key in ['distance', 'mass1', 'mass2']:
            cols[key] *= zp1

        # Populate sim_inspiral table
        sims = xmlroot.appendChild(lsctables.New(lsctables.SimInspiralTable))
        for row in zip(*cols.values()):
            sims.appendRow(
                **dict(
                    dict.fromkeys(sims.validcolumns, None),
                    process_id=process.process_id,
                    simulation_id=sims.get_next_id(),
                    waveform=args.waveform,
                    f_lower=args.f_low,
                    **dict(zip(cols.keys(), row))))

        # Record process end time.
        process.comment = str(volumetric_rate)
        process.set_end_time_now()

        # Write output file.
        write_fileobj(xmldoc, args.output)
