#
# Copyright (C) 2019  Leo Singer
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
"""Rough-cut injection tool.

The idea is to efficiently sample events, uniformly in comoving volume, and
from a distribution of masses and spins, such that later detection cuts will
not reject an excessive number of events. We divide the intrinsic parameter
space into a very coarse grid and we calculate the maximum horizon distance in
each grid cell."""

import functools

from astropy import cosmology
from astropy import units
import lal
import numpy as np

from ..bayestar.filter import (
    InterpolatedPSD, abscissa, signal_psd_series, sngl_inspiral_psd)
from . import (
    ArgumentParser, FileType, random_parser, register_to_xmldoc)

lal.ClobberDebugLevel(lal.LALNDEBUG)


def get_snr_at_z(cosmo, psds, H, z):
    f = abscissa(H)
    interp = InterpolatedPSD(f, H.data.data, fill_value=0)
    Hinterp = lal.CreateREAL8FrequencySeries(
        H.name, H.epoch, H.f0, H.deltaF, H.sampleUnits, len(H.data.data))
    Hinterp.data.data[:] = interp(f * (1 + z))
    HSs = [signal_psd_series(Hinterp, psd) for psd in psds]
    snrs = sorted(np.sqrt(4 * np.trapz(HS.data.data, dx=HS.deltaF))
                  for HS in HSs if np.any(HS.data.data != 0))
    if len(snrs) > 1:
        result = snrs[-2]
    elif len(snrs) == 1:
        result = snrs[-1]
    else:
        result = 0.0
    return result / cosmo.angular_diameter_distance(z).value


def get_max_comoving_distance(cosmo, psds, waveform, f_low, min_snr, params):
    mass1, mass2, spin1z, spin2z = params
    H = sngl_inspiral_psd(waveform, f_low=f_low,
                          mass1=mass1, mass2=mass2,
                          spin1z=spin1z, spin2z=spin2z)
    nonzero = np.flatnonzero(H.data.data)
    if len(nonzero) == 0:
        return 0.0
    H = lal.CutREAL8FrequencySeries(
        H, int(nonzero[0]), int(nonzero[-1] - nonzero[0] + 1))
    func = functools.partial(get_snr_at_z, cosmo, psds, H)
    z = cosmology.z_at_value(func, min_snr, zmax=100)
    return cosmo.comoving_distance(z).value


def z_for_comoving_distance(cosmo, distance):
    if distance == 0.0:
        return 0.0
    return cosmology.z_at_value(
        cosmo.comoving_distance, distance * units.Mpc, zmax=100)


def parser():
    parser = ArgumentParser(parents=[random_parser])
    parser.add_argument(
        '--cosmology', choices=cosmology.parameters.available,
        default='WMAP9', help='Cosmological model')
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
    return parser


def main(args=None):
    import itertools
    import os

    from astropy.utils.console import ProgressBar
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
        else:
            raise AssertionError('This line hould not be reached.')
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
        else:
            raise AssertionError('This line hould not be reached.')
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
        else:
            raise AssertionError('This line hould not be reached.')
    else:
        raise AssertionError('This line hould not be reached.')

    dists = (m1_dist, m2_dist, x1_dist, x2_dist)

    # Read PSDs
    psds = list(
        lal.series.read_psd_xmldoc(
            ligolw_utils.load_fileobj(
                args.reference_psd,
                contenthandler=lal.series.PSDContentHandler)[0]).values())
    psds = tuple(InterpolatedPSD(abscissa(psd), psd.data.data) for psd in psds)

    # Construct mass1, mass2, spin1z, spin2z grid.
    m1 = np.geomspace(m1_min, m1_max, 5)
    m2 = np.geomspace(m2_min, m2_max, 5)
    x1 = np.linspace(x1_min, x1_max, 5)
    x2 = np.linspace(x2_min, x2_max, 5)
    params = m1, m2, x1, x2

    # Calculate the maximum distance on the grid.
    shape = tuple(len(param) for param in params)
    max_distance = np.reshape(
        ProgressBar.map(
            functools.partial(
                get_max_comoving_distance, cosmo, psds,
                args.waveform, args.f_low, args.min_snr),
            np.column_stack([param.ravel() for param
                             in np.meshgrid(*params, indexing='ij')]),
            multiprocess=True),
        shape)

    # Make sure that we filled in all entries
    assert np.all(max_distance >= 0)

    # Find piecewise constant approximate upper bound on distance:
    # Calculate approximate gradient at each grid point, then approximate
    # function with a plane at that point, and find the maximum of that plane
    # in a square patch around that point
    max_distance_grad = np.asarray(np.gradient(max_distance, *params))
    param_edges = [
        np.concatenate(((p[0],), 0.5 * (p[1:] + p[:-1]), (p[-1],)))
        for p in params]
    param_los = [param_edge[:-1] for param_edge in param_edges]
    param_his = [param_edge[1:] for param_edge in param_edges]
    lo_hi_deltas = [((param_lo, param_hi) - param)
                    for param_lo, param_hi, param
                    in zip(param_los, param_his, params)]
    corner_deltas = np.asarray([np.meshgrid(*delta, indexing='ij')
                                for delta in itertools.product(*lo_hi_deltas)])
    max_distance += (corner_deltas * max_distance_grad).sum(1).max(0)

    # Truncate maximum distance at the particle horizon.
    max_distance = np.minimum(
        max_distance, cosmo.comoving_distance(np.inf).value)

    # Calculate V * T in each grid cell
    cdf_los = [dist.cdf(param_lo) for param_lo, dist in zip(param_los, dists)]
    cdf_his = [dist.cdf(param_hi) for param_hi, dist in zip(param_his, dists)]
    cdfs = [cdf_hi - cdf_lo for cdf_lo, cdf_hi in zip(cdf_los, cdf_his)]
    probs = np.prod(np.meshgrid(*cdfs, indexing='ij'), axis=0)
    probs /= probs.sum()
    probs *= 4/3*np.pi*max_distance**3
    volume = probs.sum()
    probs /= volume
    probs = probs.ravel()

    volumetric_rate = args.nsamples / volume * units.year**-1 * units.Mpc**-3

    # Draw random grid cells
    dist = stats.rv_discrete(values=(np.arange(len(probs)), probs))
    indices = np.unravel_index(dist.rvs(size=args.nsamples), shape)

    # Draw random intrinsic params from each cell
    values = [
        dist.ppf(stats.uniform(cdf_lo[i], cdf[i]).rvs(size=args.nsamples))
        for i, dist, cdf_lo, cdf in zip(indices, dists, cdf_los, cdfs)]

    # Draw random extrinsic parameters for each cell
    dist = stats.powerlaw(a=3, scale=max_distance[indices])
    values.append(dist.rvs(size=args.nsamples))
    dist = stats.uniform(0, 2 * np.pi)
    values.append(dist.rvs(size=args.nsamples))
    dist = stats.uniform(-1, 2)
    values.append(np.arcsin(dist.rvs(size=args.nsamples)))
    dist = stats.uniform(-1, 2)
    values.append(np.arccos(dist.rvs(size=args.nsamples)))
    dist = stats.uniform(0, 2 * np.pi)
    values.append(dist.rvs(size=args.nsamples))
    dist = stats.uniform(-np.pi, 2 * np.pi)
    values.append(dist.rvs(size=args.nsamples))
    dist = stats.uniform(1e9, units.year.to(units.second))
    values.append(np.sort(dist.rvs(size=args.nsamples)))

    # Populate sim_inspiral table
    sims = xmlroot.appendChild(lsctables.New(lsctables.SimInspiralTable))
    keys = ('mass1', 'mass2', 'spin1z', 'spin2z',
            'distance', 'longitude', 'latitude',
            'inclination', 'polarization', 'coa_phase', 'time_geocent')
    for row in zip(*values):
        sims.appendRow(
            **dict(
                dict.fromkeys(sims.validcolumns, None),
                process_id=process.process_id,
                simulation_id=sims.get_next_id(),
                waveform=args.waveform,
                f_lower=args.f_low,
                **dict(zip(keys, row))))

    # Apply redshift factor
    colnames = ['distance', 'mass1', 'mass2']
    columns = [sims.getColumnByName(colname) for colname in colnames]
    zp1 = 1 + np.asarray(ProgressBar.map(
        functools.partial(z_for_comoving_distance, cosmo),
        np.asarray(columns[0]),
        multiprocess=True))
    for column in columns:
        column[:] = np.asarray(column) * zp1

    # Record process end time.
    process.comment = str(volumetric_rate)
    ligolw_process.set_process_end_time(process)

    # Write output file.
    with ligolw_utils.SignalsTrap():
        ligolw_utils.write_fileobj(
            xmldoc, args.output,
            gz=(os.path.splitext(args.output.name)[-1] == ".gz"))
