#
# Copyright (C) 2019-2020  Leo Singer
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

The idea is to efficiently sample events, uniformly in comoving volume, and
from a distribution of masses and spins, such that later detection cuts will
not reject an excessive number of events. We divide the intrinsic parameter
space into a very coarse grid and we calculate the maximum horizon distance in
each grid cell.
"""

from functools import partial

from astropy import cosmology
from astropy.cosmology.core import vectorize_if_needed
from astropy import units
from astropy.units import dimensionless_unscaled
import lal
import numpy as np
from scipy.integrate import quad, fixed_quad
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar
from scipy.ndimage import maximum_filter

from ..util import progress_map
from ..bayestar.filter import sngl_inspiral_psd
from . import (
    ArgumentParser, FileType, random_parser, register_to_xmldoc, write_fileobj)

lal.ClobberDebugLevel(lal.LALNDEBUG)


def get_decisive_snr(snrs):
    """Return the SNR for the trigger that decides if an event is detectable.

    If there are two or more detectors, then the decisive SNR is the SNR of the
    second loudest detector (since a coincidence of two or more events is
    required). If there is only one detector, then the decisive SNR is just the
    SNR of that detector. If there are no detectors, then 0 is returned.

    Parameters
    ----------
    snrs : list
        List of SNRs (floats).

    Returns
    -------
    decisive_snr : float

    """
    if len(snrs) > 1:
        return sorted(snrs)[-2]
    elif len(snrs) == 1:
        return snrs[0]
    else:
        return 0.0


def lo_hi_nonzero(x):
    nonzero = np.flatnonzero(x)
    return nonzero[0], nonzero[-1]


def z_at_snr(cosmo, psds, waveform, f_low, snr, mass1, mass2, spin1z, spin2z):
    """
    Get redshift at which a waveform attains a given SNR.

    Parameters
    ----------
    cosmo : :class:`astropy.cosmology.FLRW`
        The cosmological model.
    psds : list
        List of :class:`lal.REAL8FrequencySeries` objects.
    waveform : str
        Waveform approximant name.
    f_low : float
        Low-frequency cutoff for template.
    snr : float
        Target SNR.
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
        log_f = np.log(series.f0 + series.deltaF * np.arange(i_lo, i_hi + 1))
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
        snr = get_decisive_snr(np.sqrt(4 * integrals))
        with np.errstate(divide='ignore'):
            snr /= cosmo.angular_diameter_distance(z).to_value(units.Mpc)
        return snr

    return root_scalar(lambda z: snr_at_z(z) - snr, bracket=(0, 1e3)).root


def get_max_z(cosmo, psds, waveform, f_low, snr, mass1, mass2, spin1z, spin2z,
              jobs=1):
    # Calculate the maximum distance on the grid.
    params = [mass1, mass2, spin1z, spin2z]
    result = progress_map(
        partial(z_at_snr, cosmo, psds, waveform, f_low, snr),
        *(param.ravel() for param in np.meshgrid(*params, indexing='ij')),
        jobs=jobs)
    result = np.reshape(result, tuple(len(param) for param in params))

    assert np.all(result >= 0), 'some redshifts are negative'
    assert np.all(np.isfinite(result)), 'some redshifts are not finite'
    return result


def _sensitive_volume_integral(cosmo, z):
    dh3_sr = cosmo.hubble_distance**3 / units.sr

    def integrand(z):
        result = cosmo.differential_comoving_volume(z)
        result /= (1 + z) * dh3_sr
        return result.to_value(dimensionless_unscaled)

    def integral(z):
        result, _ = quad(integrand, 0, z)
        return result

    return vectorize_if_needed(integral, z)


def sensitive_volume(cosmo, z):
    """Sensitive volume :math:`V(z)` out to redshift :math:`z`.

    Given a population of events that occur at a constant rate density
    :math:`R` per unit comoving volume per unit proper time, the number of
    observed events out to a redshift :math:`N(z)` over an observation time
    :math:`T` is :math:`N(z) = R T V(z)`.
    """
    dh3 = cosmo.hubble_distance**3
    return 4 * np.pi * dh3 * _sensitive_volume_integral(cosmo, z)


def sensitive_distance(cosmo, z):
    r"""Sensitive distance as a function of redshift :math:`z`.

    The sensitive distance is the distance :math:`d_s(z)` defined such that
    :math:`V(z) = 4/3\pi {d_s(z)}^3`, where :math:`V(z)` is the sensitive
    volume.
    """
    dh = cosmo.hubble_distance
    return dh * np.cbrt(3 * _sensitive_volume_integral(cosmo, z))


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
    parser = ArgumentParser(parents=[random_parser])
    parser.add_argument(
        '--cosmology', choices=cosmology.parameters.available,
        default='Planck15', help='Cosmological model')
    parser.add_argument(
        '--distribution', required=True, choices=(
            'bns_astro', 'bns_broad', 'nsbh_astro', 'nsbh_broad',
            'bbh_astro', 'bbh_broad'))
    parser.add_argument('--reference-psd', type=FileType('rb'), required=True)
    parser.add_argument('--f-low', type=float, default=25.0)
    parser.add_argument('--min-snr', type=float, default=4)
    parser.add_argument('--waveform', default='o2-uberbank')
    parser.add_argument('--nsamples', type=int, default=100000)
    parser.add_argument('-o', '--output', type=FileType('wb'), default='-')
    parser.add_argument(
        '-j', '--jobs', type=int, default=1, const=None, nargs='?',
        help='Number of threads')
    return parser


def main(args=None):
    from glue.ligolw import lsctables
    from glue.ligolw.utils import process as ligolw_process
    from glue.ligolw import utils as ligolw_utils
    from glue.ligolw import ligolw
    import lal.series
    from scipy import stats

    p = parser()
    args = p.parse_args(args)

    xmldoc = ligolw.Document()
    xmlroot = xmldoc.appendChild(ligolw.LIGO_LW())
    process = register_to_xmldoc(xmldoc, p, args)

    cosmo = cosmology.default_cosmology.get_cosmology_from_string(
        args.cosmology)

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
    ns_broad_mass_dist = stats.uniform(ns_mass_min, ns_mass_max - ns_mass_min)
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

    # Read PSDs
    psds = list(
        lal.series.read_psd_xmldoc(
            ligolw_utils.load_fileobj(
                args.reference_psd,
                contenthandler=lal.series.PSDContentHandler)[0]).values())

    # Construct mass1, mass2, spin1z, spin2z grid.
    m1 = np.geomspace(m1_min, m1_max, 10)
    m2 = np.geomspace(m2_min, m2_max, 10)
    x1 = np.linspace(x1_min, x1_max, 10)
    x2 = np.linspace(x2_min, x2_max, 10)
    params = m1, m2, x1, x2

    # Calculate the maximum distance on the grid.
    max_z = get_max_z(
        cosmo, psds, args.waveform, args.f_low, args.min_snr, m1, m2, x1, x2,
        jobs=args.jobs)
    max_distance = sensitive_distance(cosmo, max_z).to_value(units.Mpc)

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

    volumetric_rate = args.nsamples / volume * units.year**-1 * units.Mpc**-3

    # Draw random grid cells
    dist = stats.rv_discrete(values=(np.arange(len(probs)), probs))
    indices = np.unravel_index(
        dist.rvs(size=args.nsamples), max_distance.shape)

    # Draw random intrinsic params from each cell
    cols = {}
    cols['mass1'], cols['mass2'], cols['spin1z'], cols['spin2z'] = [
        dist.ppf(stats.uniform(cdf_lo[i], cdf[i]).rvs(size=args.nsamples))
        for i, dist, cdf_lo, cdf in zip(indices, dists, cdf_los, cdfs)]

    # Draw random extrinsic parameters
    cols['distance'] = stats.powerlaw(a=3, scale=max_distance[indices]).rvs(
        size=args.nsamples)
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
    cols['time_geocent'] = stats.uniform(1e9, units.year.to(units.second)).rvs(
        size=args.nsamples)

    # Convert from sensitive distance to redshift and comoving distance.
    # FIXME: Replace this brute-force lookup table with a solver.
    z = np.linspace(0, max_z.max(), 10000)
    ds = sensitive_distance(cosmo, z).to_value(units.Mpc)
    dc = cosmo.comoving_distance(z).to_value(units.Mpc)
    z_for_ds = interp1d(ds, z, kind='cubic', assume_sorted=True)
    dc_for_ds = interp1d(ds, dc, kind='cubic', assume_sorted=True)
    zp1 = 1 + z_for_ds(cols['distance'])
    cols['distance'] = dc_for_ds(cols['distance'])

    # Apply redshift factor to convert from comoving distance and source frame
    # masses to luminosity distance and observer frame masses.
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
    ligolw_process.set_process_end_time(process)

    # Write output file.
    write_fileobj(xmldoc, args.output)
