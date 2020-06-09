#
# Copyright (C) 2011-2020  Leo Singer
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
        '--marker', metavar='matplotlib_marker', default='*', type=str,
        help='symbol used when marking --radec (has to be one of the matplotlib marker strings)')
    parser.add_argument(
        '--marker-color', metavar='color_name', default='white', type=str,
        help='color for the --radec marker')
    parser.add_argument(
        '--marker-ecolor', metavar='color_name', default='black', type=str,
        help='edge color for the --radec marker')
    parser.add_argument(
        '--marker-size', metavar='size_int', default=10, type=int,
        help='size of the --radec marker symbol')
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
    
    # separate display group for PTA related options
    pta_group = parser.add_argument_group('PTA options', 
                        'Options for plotting the pulsars of a Pulsar Timing Array')
    pta_group.add_argument(
        '--pta-file', metavar='pulsars.csv', type=str,
        help='path to PTA pulsars file with columns ra, dec (in degrees), '
        'and noise rms (in seconds). The rms column can be left out if '
        'specifying the option --pta-fix-marker-size.')
    pta_group.add_argument(
        '--pta-size', metavar='int', type=int,
        help='Optional when using --pta-file, only plot the first INT pulsar')
    pta_group.add_argument(
        '--pta-marker', metavar='matplotlib_marker', default='*', type=str,
        help='symbol used when marking pulsars from --pta-file')
    pta_group.add_argument(
        '--pta-color', metavar='color_name', default='white', type=str,
        help='color for the --pta-file markers')
    pta_group.add_argument(
        '--pta-ecolor', metavar='color_name', default='black', type=str,
        help='edge color for the --pta-file markers')
    pta_group.add_argument(
        '--pta-fix-marker-size', metavar='size_int', default=False,
        help="Optional when using --pta-file. By default, the pulsar noise rms "
        "area read from the PTA file and used to scale the pulsar markers' sizes. "
        "Using this option fixes the markers' sizes instead.")
    
    return parser


def main(args=None):
    opts = parser().parse_args(args)

    # Late imports

    import os
    import json
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    from ..io import fits
    from .. import plot
    from .. import postprocess
    import astropy_healpix as ah
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astropy import units as u

    skymap, metadata = fits.read_sky_map(opts.input.name, nest=None)
    nside = ah.npix_to_nside(len(skymap))

    # Convert sky map from probability to probability per square degree.
    deg2perpix = ah.nside_to_pixel_area(nside).to_value(u.deg**2)
    probperdeg2 = skymap / deg2perpix

    axes_args = {}
    if opts.geo:
        axes_args['projection'] = 'geo'
        obstime = Time(metadata['gps_time'], format='gps').utc.isot
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
                   WHERE cm2.event_id = ?
                   AND cm1.table_name = 'sim_inspiral'
                   AND cm2.table_name = 'coinc_event_id'
                   '''
        (ra, dec), = opts.inj_database.execute(
            query, (metadata['objid'],)).fetchall()
        radecs.append(np.rad2deg([ra, dec]).tolist())

    # Add markers (e.g., for injections or external triggers).
    for ra, dec in radecs:
        ax.plot_coord(
            SkyCoord(ra, dec, unit='deg'), 
            opts.marker,
            markerfacecolor=opts.marker_color, 
            markeredgecolor=opts.marker_ecolor,
            markersize=opts.marker_size)
        
    # If given, add additional markers for the pulsars from a PTA
    if opts.pta_file:
        try:
            with open(opts.pta_file, 'r') as f:
                pta_data = np.loadtxt(f, comments='#', delimiter=',')
            
        except IOError:
            print('Could not find PTA file, skipping')
        except ValueError:
            print('Could not read points from PTA file, skipping')
            
        
        if opts.pta_size:
            pta_data = pta_data[:opts.pta_size]
        
        # fixed pta marker sizes
        if opts.pta_fix_marker_size:
            pulsar_sizes = np.full(len(pta_data), opts.pta_fix_marker_size)
        # scale pta marker size with noise rms: bigger markers for lower noise
        else:
            pulsar_sizes = 10 * (pta_data[:, 2] / 1.e-7)**(-0.4)
            
        for i, pulsar in enumerate(pta_data):
            ra, dec = pulsar[:2]
            ax.plot_coord(
                SkyCoord(ra, dec, unit='deg'),
                opts.pta_marker,
                markerfacecolor=opts.pta_color,
                markeredgecolor=opts.pta_ecolor,
                markersize=pulsar_sizes[i])

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
                text.append('{:d}% area: {:,d} deg²'.format(p, i))
        ax.text(1, 1, '\n'.join(text), transform=ax.transAxes, ha='right')

    # Show or save output.
    opts.output()
