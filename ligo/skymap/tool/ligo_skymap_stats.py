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
Calculate summary statistics for a batch of sky maps.

The filenames of the sky maps may be provided as positional command line
arguments, and may also be provided as globs (such as ``*.fits.gz``). If
supplied with the optional ``--database`` argument, then also match sky maps
with injections from an inspinjfind-style sqlite database.

All angular separations are in degrees, all areas are in square degrees, and
all volumes are in cubic megaparsecs. The output is written as tab-separated
values with the following columns:

+-----------------------------------------------------------------------------+
| **From the search pipeline database**                                       |
+------------------------+----------------------------------------------------+
| ``coinc_event_id``     | event identifier                                   |
+------------------------+----------------------------------------------------+
| ``simulation_id``      | injection identifier                               |
+------------------------+----------------------------------------------------+
| ``far``                | false alarm rate                                   |
+------------------------+----------------------------------------------------+
| ``snr``                | signal to noise ratio                              |
+------------------------+----------------------------------------------------+
| **Injection finding**                                                       |
+------------------------+----------------------------------------------------+
| ``searched_area``      | area of 2D credible region containing the true sky |
|                        | location                                           |
+------------------------+----------------------------------------------------+
| ``searched_prob``      | probability in that 2D credible region             |
+------------------------+----------------------------------------------------+
| ``searched_prob_dist`` | marginal distance CDF at the true distance         |
+------------------------+----------------------------------------------------+
| ``searched_vol``       | volume of 3D credible region containing the true   |
|                        | position                                           |
+------------------------+----------------------------------------------------+
| ``searched_prob_vol``  | probability contained in that volume               |
+------------------------+----------------------------------------------------+
| ``offset``             | angular separation between the maximum             |
|                        | *a posteriori* position and the true sky position  |
+------------------------+----------------------------------------------------+
| **Additional metadata from the sky maps**                                   |
+------------------------+----------------------------------------------------+
| ``runtime``            | wall clock run time to generate sky map            |
+------------------------+----------------------------------------------------+
| ``distmean``           | mean *a posteriori* distance                       |
+------------------------+----------------------------------------------------+
| ``diststd``            | *a posteriori* standard deviation of distance      |
+------------------------+----------------------------------------------------+
| ``log_bci``            | natural log Bayes factor, coherent vs. incoherent  |
+------------------------+----------------------------------------------------+
| ``log_bsn``            | natural log Bayes factor, signal vs. noise         |
+------------------------+----------------------------------------------------+
| **Credible levels** (if ``--area`` or ``--contour`` options present)        |
+------------------------+----------------------------------------------------+
| ``area(P)``            | area of the *P* percent 2D credible region         |
+------------------------+----------------------------------------------------+
| ``prob(A)``            | probability contained within the 2D credible level |
|                        | of area *A*                                        |
+------------------------+----------------------------------------------------+
| ``dist(P)``            | distance for a cumulative marginal probability of  |
|                        | *P* percent                                        |
+------------------------+----------------------------------------------------+
| ``vol(P)``             | volume of the *P* percent 3D credible region       |
+------------------------+----------------------------------------------------+
| **Modes** (if ``--modes`` option is present)                                |
+------------------------+----------------------------------------------------+
| ``searched_modes``     | number of simply connected figures in the 2D       |
|                        | credible region containing the true sky location   |
+------------------------+----------------------------------------------------+
| ``modes(P)``           | number of simply connected figures in the *P*      |
|                        | percent 2D credible region                         |
+------------------------+----------------------------------------------------+

"""

import sqlite3

import numpy as np

from . import ArgumentParser, FileType, SQLiteType
from ..io import fits
from ..postprocess import find_injection_moc
from ..util import sqlite


def parser():
    parser = ArgumentParser()
    parser.add_argument(
        '-o', '--output', metavar='OUT.dat', type=FileType('w'), default='-',
        help='Name of output file')
    parser.add_argument(
        '-j', '--jobs', type=int, default=1, const=None, nargs='?',
        help='Number of threads')
    parser.add_argument(
        '-p', '--contour', default=[], nargs='+', type=float,
        metavar='PERCENT',
        help='Report the area of the smallest contour and the number of modes '
        'containing this much probability.')
    parser.add_argument(
        '-a', '--area', default=[], nargs='+', type=float, metavar='DEG2',
        help='Report the largest probability contained within any region '
        'of this area in square degrees. Can be repeated multiple times.')
    parser.add_argument(
        '--modes', action='store_true',
        help='Compute number of disjoint modes')
    parser.add_argument(
        '-d', '--database', type=SQLiteType('r'), metavar='DB.sqlite',
        help='Input SQLite database from search pipeline')
    parser.add_argument(
        'fitsfilenames', metavar='GLOB.fits[.gz]', nargs='+', action='glob',
        help='Input FITS filenames and/or globs')
    parser.add_argument(
        '--cosmology', action='store_true',
        help='Report volume localizations as comoving volumes.')
    return parser


def startup(dbfilename, opts_contour, opts_modes, opts_area, opts_cosmology):
    global db, contours, modes, areas, cosmology
    if dbfilename is None:
        db = None
    else:
        db = sqlite3.connect(dbfilename)
    contours = opts_contour
    modes = opts_modes
    areas = opts_area
    cosmology = opts_cosmology


def process(fitsfilename):
    sky_map = fits.read_sky_map(fitsfilename, moc=True)

    coinc_event_id = sky_map.meta.get('objid')
    try:
        runtime = sky_map.meta['runtime']
    except KeyError:
        runtime = float('nan')

    contour_pvalues = 0.01 * np.asarray(contours)

    if db is None:
        simulation_id = true_ra = true_dec = true_dist = far = snr = None
    else:
        row = db.execute(
            """
            SELECT DISTINCT sim.simulation_id AS simulation_id,
            sim.longitude AS ra, sim.latitude AS dec, sim.distance AS distance,
            ci.combined_far AS far, ci.snr AS snr
            FROM coinc_event_map AS cem1 INNER JOIN coinc_event_map AS cem2
            ON (cem1.coinc_event_id = cem2.coinc_event_id)
            INNER JOIN sim_inspiral AS sim
            ON (cem1.event_id = sim.simulation_id)
            INNER JOIN coinc_inspiral AS ci
            ON (cem2.event_id = ci.coinc_event_id)
            WHERE cem1.table_name = 'sim_inspiral'
            AND cem2.table_name = 'coinc_event' AND cem2.event_id = ?
            """, (coinc_event_id,)).fetchone()
        if row is None:
            return None
        simulation_id, true_ra, true_dec, true_dist, far, snr = row

    (
        searched_area, searched_prob, offset, searched_modes, contour_areas,
        area_probs, contour_modes, searched_prob_dist, contour_dists,
        searched_vol, searched_prob_vol, contour_vols
    ) = find_injection_moc(
        sky_map, true_ra, true_dec, true_dist,
        contours=contour_pvalues, areas=areas, modes=modes, cosmology=cosmology
    )

    if snr is None:
        snr = np.nan
    if far is None:
        far = np.nan
    distmean = sky_map.meta.get('distmean', np.nan)
    diststd = sky_map.meta.get('diststd', np.nan)
    log_bci = sky_map.meta.get('log_bci', np.nan)
    log_bsn = sky_map.meta.get('log_bsn', np.nan)

    ret = [coinc_event_id]
    if db is not None:
        ret += [
            simulation_id, far, snr, searched_area,
            searched_prob, searched_prob_dist, searched_vol, searched_prob_vol,
            offset]
    ret += [runtime, distmean, diststd, log_bci, log_bsn]
    ret += contour_areas + area_probs + contour_dists + contour_vols
    if modes:
        if db is not None:
            ret += [searched_modes]
        ret += contour_modes
    return ret


def main(args=None):
    opts = parser().parse_args(args)

    from tqdm import tqdm
    from .. import omp

    if opts.database is None:
        dbfilename = None
    else:
        dbfilename = sqlite.get_filename(opts.database)
    args = (dbfilename, opts.contour, opts.modes, opts.area, opts.cosmology)
    if opts.jobs == 1:
        pool_map = map
    else:
        omp.num_threads = 1  # disable OpenMP parallelism

        from multiprocessing import Pool
        pool_map = Pool(opts.jobs, startup, args).imap
    startup(*args)

    colnames = ['coinc_event_id']
    if opts.database is not None:
        colnames += ['simulation_id', 'far', 'snr', 'searched_area',
                     'searched_prob', 'searched_prob_dist', 'searched_vol',
                     'searched_prob_vol', 'offset']
    colnames += ['runtime', 'distmean', 'diststd', 'log_bci', 'log_bsn']
    colnames += ['area({0:g})'.format(_) for _ in contours]
    colnames += ['prob({0:g})'.format(_) for _ in areas]
    colnames += ['dist({0:g})'.format(_) for _ in contours]
    colnames += ['vol({0:g})'.format(_) for _ in contours]
    if modes:
        if opts.database is not None:
            colnames += ['searched_modes']
        colnames += ["modes({0:g})".format(p) for p in contours]
    print(*colnames, sep="\t", file=opts.output)

    with tqdm(total=len(opts.fitsfilenames)) as progress:
        for record in pool_map(process, opts.fitsfilenames):
            if record is not None:
                print(*record, sep="\t", file=opts.output)
            progress.update()
