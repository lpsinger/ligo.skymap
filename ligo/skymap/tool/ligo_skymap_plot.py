# -*- coding: utf-8 -*-
#
# Copyright (C) 2011-2018  Leo Singer
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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Plot an all-sky map on a Mollweide projection.
By default, plot in celestial coordinates (RA, Dec).

To plot in geographic coordinates (longitude, latitude) with major
coastlines overlaid, provide the ``--geo`` flag.

Public-domain cartographic data is courtesy of `Natural Earth
<http://www.naturalearthdata.com>`_ and processed with `MapShaper
<http://www.mapshaper.org>`_.
"""

from . import ArgumentParser, FileType, SQLiteType, figure_parser


def parser():
    parser = ArgumentParser(parents=[figure_parser])
    parser.add_argument(
        '--annotate', default=False, action='store_true',
        help='annotate plot with information about the event')
    parser.add_argument(
        '--contour', metavar='PERCENT', type=float, nargs='+',
        help='plot contour enclosing this percentage of'
        ' probability mass [may be specified multiple times, default: none]')
    parser.add_argument(
        '--colorbar', default=False, action='store_true',
        help='Show colorbar')
    parser.add_argument(
        '--radec', nargs=2, metavar='deg', type=float, action='append',
        default=[], help='right ascension (deg) and declination (deg) to mark')
    parser.add_argument(
        '--inj-database', metavar='FILE.sqlite', type=SQLiteType('r'),
        help='read injection positions from database')
    parser.add_argument(
        '--geo', action='store_true',
        help='Plot in geographic coordinates, (lat, lon) instead of (RA, Dec)')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    return parser


def main(args=None):
    opts = parser().parse_args(args)

    # Late imports

    import os
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    import healpy as hp
    from ..io import fits
    from .. import plot
    from .. import postprocess
    from astropy.coordinates import SkyCoord
    from astropy.time import Time

    skymap, metadata = fits.read_sky_map(opts.input.name, nest=None)
    nside = hp.npix2nside(len(skymap))

    # Convert sky map from probability to probability per square degree.
    deg2perpix = hp.nside2pixarea(nside, degrees=True)
    probperdeg2 = skymap / deg2perpix

    if opts.geo:
        obstime = Time(metadata['gps_time'], format='gps').utc.isot
        ax = plt.axes(projection='geo degrees mollweide', obstime=obstime)
    else:
        ax = plt.axes(projection='astro hours mollweide')
    ax.grid()

    # Plot sky map.
    vmax = probperdeg2.max()
    img = ax.imshow_hpx(
        (probperdeg2, 'ICRS'), nested=metadata['nest'], vmin=0., vmax=vmax)

    # Add colorbar.
    if opts.colorbar:
        cb = plot.colorbar(img)
        cb.set_label(r'prob. per deg$^2$')

    # Add contours.
    if opts.contour:
        cls = 100 * postprocess.find_greedy_credible_levels(skymap)
        cs = ax.contour_hpx(
            (cls, 'ICRS'), nested=metadata['nest'],
            colors='k', linewidths=0.5, levels=opts.contour)
        fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
        plt.clabel(cs, fmt=fmt, fontsize=6, inline=True)

    # Add continents.
    if opts.geo:
        geojson_filename = os.path.join(
            os.path.dirname(plot.__file__), 'ne_simplified_coastline.json')
        with open(geojson_filename, 'r') as geojson_file:
            geoms = json.load(geojson_file)['geometries']
        verts = [coord for geom in geoms
                 for coord in zip(*geom['coordinates'])]
        plt.plot(*verts, color='0.5', linewidth=0.5,
                 transform=ax.get_transform('world'))

    radecs = opts.radec
    if opts.inj_database:
        query = '''SELECT DISTINCT longitude, latitude FROM sim_inspiral AS si
                   INNER JOIN coinc_event_map AS cm1
                   ON (si.simulation_id = cm1.event_id)
                   INNER JOIN coinc_event_map AS cm2
                   ON (cm1.coinc_event_id = cm2.coinc_event_id)
                   WHERE cm2.event_id = ?'''
        (ra, dec), = opts.inj_database.execute(
            query, (metadata['objid'],)).fetchall()
        radecs.append(np.rad2deg([ra, dec]).tolist())

    # Add markers (e.g., for injections or external triggers).
    for ra, dec in radecs:
        ax.plot_coord(
            SkyCoord(ra, dec, unit='deg'), '*',
            markerfacecolor='white', markeredgecolor='black', markersize=10)

    # Add a white outline to all text to make it stand out from the background.
    plot.outline_text(ax)

    if opts.annotate:
        text = []
        try:
            objid = metadata['objid']
        except KeyError:
            pass
        else:
            text.append('event ID: {}'.format(objid))
        if opts.contour:
            pp = np.round(opts.contour).astype(int)
            ii = np.round(np.searchsorted(np.sort(cls), opts.contour) *
                          deg2perpix).astype(int)
            for i, p in zip(ii, pp):
                # FIXME: use Unicode symbol instead of TeX '$^2$'
                # because of broken fonts on Scientific Linux 7.
                text.append(
                    u'{:d}% area: {:d} deg²'.format(p, i, grouping=True))
        ax.text(1, 1, '\n'.join(text), transform=ax.transAxes, ha='right')

    # Show or save output.
    opts.output()
