# -*- coding: utf-8 -*-
#
# Copyright (C) 2018  Leo Singer
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
Make an airmass chart for a LIGO/Virgo probability sky map.
"""

from argparse import FileType

import numpy as np

from . import ArgumentParser, figure_parser


def parser():
    from astropy.coordinates import EarthLocation
    parser = ArgumentParser(parents=[figure_parser])
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    parser.add_argument(
        '--site', choices=EarthLocation.get_site_names(),
        help='Observatory site', required=True)
    parser.add_argument(
        '--time', help='UTC time')
    return parser


def condition_secz(x):
    """Condition secz airmass formula: values <= 0 are below the horizon,
    which we map to infinite airmass."""
    return np.where(x <= 0, np.inf, x)


def main(args=None):
    opts = parser().parse_args(args)

    # Late imports
    import json

    from astroplan import Observer
    from astroplan.plots import plot_airmass
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from astropy.utils.data import get_pkg_data_fileobj
    from matplotlib import dates
    from matplotlib.patches import Patch
    from matplotlib import pyplot as plt
    from tqdm import tqdm
    import pytz

    from ..io import fits
    from .. import moc
    from .. import plot  # noqa
    from ..extern.quantile import quantile

    # FIXME: Astropy's site registry provides time zones,
    # but the Astropy API does not provide access to them.
    with get_pkg_data_fileobj('coordinates/sites.json',
                              package='astropy') as f:
        sites = json.load(f)
    timezones = {key: value.get('timezone') for key, value in sites.items()}
    for site in sites.values():
        timezones.update(dict.fromkeys(site['aliases'], site.get('timezone')))

    m = fits.read_sky_map(opts.input.name, moc=True)

    # Make an empty airmass chart.
    # FIXME: have to add a dummy target until
    # https://github.com/astropy/astroplan/pull/349
    # is in a release of astroplan
    observer = Observer.at_site(opts.site)
    t0 = Time(opts.time) if opts.time is not None else Time.now()
    t0 = observer.midnight(t0)
    ax = plot_airmass([SkyCoord(0, 0, unit='rad')], observer, t0,
                      brightness_shading=True, altitude_yaxis=True)

    # Remove the fake source and determine times that were used for the plot.
    del ax.lines[:]
    times = Time(np.linspace(*ax.get_xlim()), format='plot_date')

    theta, phi = moc.uniq2ang(m['UNIQ'])
    coords = SkyCoord(phi, 0.5 * np.pi - theta, unit='rad')
    prob = moc.uniq2pixarea(m['UNIQ']) * m['PROBDENSITY']

    levels = np.arange(0.9, 0.0, -0.1)
    nlevels = len(levels)
    quantiles = np.concatenate((0.5 - 0.5 * levels, 0.5 + 0.5 * levels))

    airmass = np.column_stack([
        quantile(
            condition_secz(coords.transform_to(observer.altaz(t)).secz),
            quantiles,
            weights=prob)
        for t in tqdm(times)])

    cmap = plt.get_cmap()
    for level, lo, hi in zip(levels, airmass[:nlevels], airmass[nlevels:]):
        ax.fill_between(times.plot_date, lo, hi, color=cmap(level), zorder=2)

    ax.legend(
        [Patch(facecolor=cmap(level)) for level in levels],
        ['{}%'.format(int(100 * level)) for level in levels])
    # ax.set_title('{} from {}'.format(m.meta['objid'], observer.name))

    # Add local time axis
    if opts.site in timezones:
        timezone = timezones[opts.site]
        tzinfo = pytz.timezone(timezone)
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ax.get_xticks())
        ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M', tz=tzinfo))
        plt.setp(ax2.get_xticklabels(), rotation=-30, ha='right')
        ax2.set_xlabel("Time from {} [{}]".format(
            min(times).to_datetime(tzinfo).date(),
            timezone))

    # Show or save output.
    opts.output()
