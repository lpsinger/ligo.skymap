#
# Copyright (C) 2011-2025  Leo Singer
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

import numpy as np

from . import ArgumentParser, FileType, SQLiteType
from .matplotlib import get_figure_parser


def parser():
    parser = ArgumentParser(parents=[get_figure_parser()])
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
        help='Use a terrestrial reference frame with coordinates (lon, lat) '
             'instead of the celestial frame with coordinates (RA, dec) '
             'and draw continents on the map')
    parser.add_argument(
        '--projection', type=str, default='mollweide',
        choices=['mollweide', 'aitoff', 'globe', 'zoom'],
        help='Projection style')
    parser.add_argument(
        '--projection-center', metavar='CENTER',
        help='Specify the center for globe and zoom projections, e.g. 14h 10d')
    parser.add_argument(
        '--zoom-radius', metavar='RADIUS',
        help='Specify the radius for zoom projections, e.g. 4deg')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    return parser


def _format_area(area):
    """Format area values for plot annotations.

    Values are printed to at least three significant figures. Values with three
    or more figures to the left of the decimal point are printed as integers.
    A comma separator is added for the thousands place.

    >>> _format_area(0.345678)
    '0.346'
    >>> _format_area(3.45678)
    '3.46'
    >>> _format_area(34.5678)
    '34.6'
    >>> _format_area(345.678)
    '346'
    >>> _format_area(3456.78)
    '3,457'
    >>> _format_area(34567.8)
    '34,568'
    """
    if area <= 100:
        return np.format_float_positional(
            area, precision=3, fractional=False, trim='-')
    else:
        return f'{np.round(area).astype(int):,d}'


def main(args=None):
    with parser().parse_args(args) as opts:
        # Late imports
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import rcParams
        from ..io import fits
        from .. import moc
        from .. import plot
        from .. import postprocess
        from astropy.coordinates import SkyCoord
        from astropy.table import Table
        from astropy.time import Time
        from astropy import units as u

        sr_to_deg2 = u.sr.to(u.deg**2)

        skymap = fits.read_sky_map(opts.input.name, moc=True)

        axes_args = {}
        if opts.geo:
            axes_args['projection'] = 'geo'
            obstime = Time(skymap.meta['gps_time'], format='gps').utc.isot
            axes_args['obstime'] = obstime
        else:
            axes_args['projection'] = 'astro'
        axes_args['projection'] += ' ' + opts.projection
        if opts.projection_center is not None:
            axes_args['center'] = SkyCoord(opts.projection_center)
        if opts.zoom_radius is not None:
            axes_args['radius'] = opts.zoom_radius
        ax = plt.axes(**axes_args)
        ax.grid()

        # Add contours.
        if opts.contour:
            dA = moc.uniq2pixarea(skymap['UNIQ'])
            dP = skymap['PROBDENSITY'] * dA
            cls = 100 * postprocess.find_greedy_credible_levels(
                dP, skymap['PROBDENSITY'])

            table = Table({'UNIQ': skymap['UNIQ'], 'CLS': cls})
            cs = ax.contour_hpx(
                (table, 'ICRS'), colors='k', linewidths=0.5,
                levels=opts.contour, order='nearest-neighbor')
            fmt = r'%g\%%' if rcParams['text.usetex'] else '%g%%'
            plt.clabel(cs, fmt=fmt, fontsize=6, inline=True)

        # Add continents.
        if opts.geo:
            plt.plot(*plot.coastlines(), color='0.5', linewidth=0.5,
                     transform=ax.get_transform('world'))

        radecs = opts.radec
        if opts.inj_database:
            query = '''SELECT DISTINCT longitude, latitude
                       FROM sim_inspiral AS si
                       INNER JOIN coinc_event_map AS cm1
                       ON (si.simulation_id = cm1.event_id)
                       INNER JOIN coinc_event_map AS cm2
                       ON (cm1.coinc_event_id = cm2.coinc_event_id)
                       WHERE cm2.event_id = ?
                       AND cm1.table_name = 'sim_inspiral'
                       AND cm2.table_name = 'coinc_event'
                    '''
            (ra, dec), = opts.inj_database.execute(
                query, (skymap.meta['objid'],)).fetchall()
            radecs.append(np.rad2deg([ra, dec]).tolist())

        # Add markers (e.g., for injections or external triggers).
        for ra, dec in radecs:
            ax.plot_coord(
                SkyCoord(ra, dec, unit=u.deg), '*',
                markerfacecolor='white', markeredgecolor='black',
                markersize=10)

        # Add a white outline to all text to make it stand out from the
        # background.
        plot.outline_text(ax)

        if opts.annotate:
            text = []
            try:
                objid = skymap.meta['objid']
            except KeyError:
                pass
            else:
                text.append('event ID: {}'.format(objid))
            if opts.contour:
                i = np.flipud(np.argsort(skymap['PROBDENSITY']))
                areas = postprocess.interp_greedy_credible_levels(
                    opts.contour, cls[i], np.cumsum(dA[i]), right=4*np.pi)
                pp = np.round(opts.contour).astype(int)
                ii = areas * sr_to_deg2
                for i, p in zip(ii, pp):
                    # FIXME: use Unicode symbol instead of TeX '$^2$'
                    # because of broken fonts on Scientific Linux 7.
                    text.append(
                        '{:d}% area: {} deg²'.format(p, _format_area(i)))
            ax.text(1, 1, '\n'.join(text), transform=ax.transAxes, ha='right')

        # Plot sky map.
        skymap['PROBDENSITY'] *= 1 / sr_to_deg2
        img = ax.imshow_hpx((skymap, 'ICRS'), vmin=0, order='nearest-neighbor')

        # Add colorbar.
        if opts.colorbar:
            cb = plot.colorbar(img)
            cb.set_label(r'prob. per deg$^2$')

        # Show or save output.
        opts.output()
