#
# Copyright (C) 2011-2018  Will M. Farr <will.farr@ligo.org>
#                          Leo P. Singer <leo.singer@ligo.org>
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
Generate a FITS sky map file from posterior samples using clustering and
kernel density estimation.

The output consist of two files:

*  ``skypost.obj``, a :mod:`pickle` representation of the kernel density
   estimator
*  ``skymap.fits.gz``, a 3D localization in HEALPix/FITS format
"""

from argparse import SUPPRESS

from . import ArgumentParser, DirType, EnableAction, FileType, random_parser


def parser():
    # Command line interface.
    parser = ArgumentParser(parents=[random_parser])
    parser.add_argument('samples', type=FileType('rb'), metavar='SAMPLES.hdf5',
                        help='posterior samples file')
    # Only present for backward compatibility with --samples syntax
    parser.add_argument('--samples', action='store_false', dest='_ignored',
                        help=SUPPRESS)
    parser.add_argument('--outdir', '-o', default='.',
                        type=DirType(create=True), help='output directory')
    parser.add_argument('--fitsoutname', default='skymap.fits',
                        metavar='SKYMAP.fits[.gz]',
                        help='filename for the FITS file')
    parser.add_argument('--loadpost', type=FileType('rb'),
                        metavar='SKYPOST.obj',
                        help='filename for pickled posterior state')
    parser.add_argument('--maxpts', type=int,
                        help='maximum number of posterior points to use')
    parser.add_argument('--trials', type=int, default=5,
                        help='number of trials at each clustering number')
    parser.add_argument('--enable-distance-map', action=EnableAction,
                        help='generate HEALPix map of distance estimates')
    parser.add_argument('--enable-multiresolution', action=EnableAction,
                        default=True,
                        help='generate a multiresolution HEALPix map')
    parser.add_argument('-j', '--jobs', action='store_true',
                        help='Use multiple threads')
    parser.add_argument('--instruments', metavar='H1|L1|V1|...', nargs='+',
                        help='instruments to store in FITS header')
    parser.add_argument('--objid', help='event ID to store in FITS header')
    return parser


def main(args=None):
    _parser = parser()
    args = _parser.parse_args(args)

    # Late imports
    from .. import io
    from ..bayestar import rasterize
    from .. import version
    from astropy.table import Table
    from astropy.time import Time
    import numpy as np
    import os
    import sys
    import pickle
    from ..kde import Clustered2Plus1DSkyKDE, Clustered2DSkyKDE
    import logging

    log = logging.getLogger()

    log.info('reading samples')
    try:
        data = io.read_samples(args.samples.name)
    except IOError:
        # FIXME: remove this code path once we support only HDF5
        data = Table.read(args.samples, format='ascii')

    if args.maxpts is not None:
        # Shuffle the data and take a random subsample.
        # Note that if arg.maxpts > len(data), then this
        # just selects the entire array.
        log.info('taking random subsample of chain')
        data = data[np.random.choice(len(data), args.maxpts, replace=False)]
    try:
        dist = data['dist']
    except KeyError:
        try:
            dist = data['distance']
        except KeyError:
            dist = None

    if args.loadpost is None:
        if dist is None:
            if args.enable_distance_map:
                _parser.error("The posterior samples file '{0}' does not have "
                              "a distance column named 'dist' or 'distance'. "
                              "Cannot generate distance map.".format(
                                  args.samples.name))
            pts = np.column_stack((data['ra'], data['dec']))
        else:
            pts = np.column_stack((data['ra'], data['dec'], dist))
        if args.enable_distance_map:
            cls = Clustered2Plus1DSkyKDE
        else:
            cls = Clustered2DSkyKDE
        skypost = cls(pts, trials=args.trials, multiprocess=args.jobs)

        log.info('pickling')
        with open(os.path.join(args.outdir, 'skypost.obj'), 'wb') as out:
            pickle.dump(skypost, out)
    else:
        skypost = pickle.load(args.loadpost)
        skypost.multiprocess = args.j

    log.info('making skymap')
    hpmap = skypost.as_healpix()
    if not args.enable_multiresolution:
        hpmap = rasterize(hpmap)
    hpmap.meta.update(io.fits.metadata_for_version_module(version))
    hpmap.meta['creator'] = _parser.prog
    hpmap.meta['origin'] = 'LIGO/Virgo'
    hpmap.meta['gps_creation_time'] = Time.now().gps
    hpmap.meta['history'] = [
        '', 'Generated by running the following script:',
        ' '.join([_parser.prog] + sys.argv[1:])]
    if args.objid is not None:
        hpmap.meta['objid'] = args.objid
    if args.instruments:
        hpmap.meta['instruments'] = args.instruments
    if args.enable_distance_map:
        hpmap.meta['distmean'] = np.mean(dist)
        hpmap.meta['diststd'] = np.std(dist)

    keys = ['time', 'time_mean', 'time_maxl']
    for key in keys:
        try:
            time = data[key]
        except KeyError:
            continue
        else:
            hpmap.meta['gps_time'] = time.mean()
            break
    else:
        log.warning(
            'Cannot determine the event time from any of the columns %r', keys)

    io.write_sky_map(os.path.join(args.outdir, args.fitsoutname),
                     hpmap, nest=True)
