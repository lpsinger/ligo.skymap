#
# Copyright (C) 2018-2023  Leo Singer
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
"""Make an airmass chart for a LIGO/Virgo/KAGRA probability sky map."""

import numpy as np

from . import ArgumentParser, FileType, HelpChoicesAction
from .matplotlib import get_figure_parser


def parser():
    from astropy.coordinates import EarthLocation
    site_names = EarthLocation.get_site_names()
    parser = ArgumentParser(parents=[get_figure_parser()])
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Print airmass table to stdout')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    parser.add_argument(
        '--time', help='UTC time')
    parser.add_argument(
        '--site', metavar='SITE', choices=site_names, help='Observatory site')
    parser.add_argument(
        '--help-site', action=HelpChoicesAction, choices=site_names)
    parser.add_argument(
        '--site-longitude', metavar='DEG', type=float,
        help='Observatory longitude on the WGS84 ellipsoid. '
        'Mutually exclusive with --site.')
    parser.add_argument(
        '--site-latitude', metavar='DEG', type=float,
        help='Observatory latitude on the WGS84 ellipsoid. '
        'Mutually exclusive with --site.')
    parser.add_argument(
        '--site-height', metavar='METERS', type=float,
        help='Observatory height from the WGS84 ellipsoid. '
        'Mutually exclusive with --site.')
    parser.add_argument(
        '--site-timezone',
        help='Observatory time zone, e.g. "US/Pacific". '
        'Mutually exclusive with --site.')
    return parser


def condition_secz(x):
    """Condition secz airmass formula: values <= 0 are below the horizon,
    which we map to infinite airmass.
    """
    return np.where(x <= 0, np.inf, x)


def clip_verylarge(x, max=1e300):
    return np.clip(x, -max, max)


def main(args=None):
    p = parser()
    opts = p.parse_args(args)

    # Late imports
    import operator
    import sys

    from astroplan import Observer
    from astroplan.plots import plot_airmass
    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.table import Table
    from astropy.time import Time
    from astropy import units as u
    from matplotlib import dates
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize
    from matplotlib.patches import Patch
    from matplotlib import pyplot as plt
    from tqdm import tqdm
    import pytz

    from ..io import fits
    from .. import moc
    from .. import plot  # noqa
    from ..extern.numpy.quantile import percentile

    if opts.site is None:
        if opts.site_longitude is None or opts.site_latitude is None:
            p.error('must specify either --site or both '
                    '--site-longitude and --site-latitude')
        location = EarthLocation(
            lon=opts.site_longitude * u.deg,
            lat=opts.site_latitude * u.deg,
            height=(opts.site_height or 0) * u.m)
        if opts.site_timezone is not None:
            location.info.meta = {'timezone': opts.site_timezone}
        observer = Observer(location)
    else:
        if not ((opts.site_longitude is None) and
                (opts.site_latitude is None) and
                (opts.site_height is None) and
                (opts.site_timezone is None)):
            p.error('argument --site not allowed with arguments '
                    '--site-longitude, --site-latitude, '
                    '--site-height, or --site-timezone')
        observer = Observer.at_site(opts.site)

    m = fits.read_sky_map(opts.input.name, moc=True)

    # Make an empty airmass chart.
    t0 = Time(opts.time) if opts.time is not None else Time.now()
    t0 = observer.midnight(t0)
    ax = plot_airmass([], observer, t0, altitude_yaxis=True)

    # Remove the fake source and determine times that were used for the plot.
    for artist in ax.lines:
        artist.remove()
    times = Time(np.linspace(*ax.get_xlim()), format='plot_date')

    theta, phi = moc.uniq2ang(m['UNIQ'])
    coords = SkyCoord(phi, 0.5 * np.pi - theta, unit='rad')
    prob = moc.uniq2pixarea(m['UNIQ']) * m['PROBDENSITY']

    levels = np.arange(90, 0, -10)
    nlevels = len(levels)
    percentiles = np.concatenate((50 - 0.5 * levels, 50 + 0.5 * levels))

    airmass = np.column_stack([
        percentile(
            condition_secz(coords.transform_to(observer.altaz(t)).secz),
            percentiles,
            weights=prob)
        for t in tqdm(times)])

    cmap = ScalarMappable(Normalize(0, 100), plt.get_cmap())
    for level, lo, hi in zip(levels, airmass[:nlevels], airmass[nlevels:]):
        ax.fill_between(
            times.plot_date,
            clip_verylarge(lo),  # Clip infinities to large but finite values
            clip_verylarge(hi),  # because fill_between cannot handle inf
            color=cmap.to_rgba(level), zorder=2)

    ax.legend(
        [Patch(facecolor=cmap.to_rgba(level)) for level in levels],
        ['{}%'.format(level) for level in levels])
    # ax.set_title('{} from {}'.format(m.meta['objid'], observer.name))

    # Adapted from astroplan
    start = times[0]
    twilights = [
        (times[0].datetime, 0.0),
        (observer.sun_set_time(
            Time(start), which='next').datetime, 0.0),
        (observer.twilight_evening_civil(
            Time(start), which='next').datetime, 0.1),
        (observer.twilight_evening_nautical(
            Time(start), which='next').datetime, 0.2),
        (observer.twilight_evening_astronomical(
            Time(start), which='next').datetime, 0.3),
        (observer.twilight_morning_astronomical(
            Time(start), which='next').datetime, 0.4),
        (observer.twilight_morning_nautical(
            Time(start), which='next').datetime, 0.3),
        (observer.twilight_morning_civil(
            Time(start), which='next').datetime, 0.2),
        (observer.sun_rise_time(
            Time(start), which='next').datetime, 0.1),
        (times[-1].datetime, 0.0),
    ]

    twilights.sort(key=operator.itemgetter(0))
    for i, twi in enumerate(twilights[1:], 1):
        if twi[1] != 0:
            ax.axvspan(twilights[i - 1][0], twilights[i][0],
                       ymin=0, ymax=1, color='grey', alpha=twi[1], linewidth=0)
        if twi[1] != 0.4:
            ax.axvspan(twilights[i - 1][0], twilights[i][0],
                       ymin=0, ymax=1, color='white', alpha=0.8 - 2 * twi[1],
                       zorder=3, linewidth=0)

    # Add local time axis
    timezone = (observer.location.info.meta or {}).get('timezone')
    if timezone:
        tzinfo = pytz.timezone(timezone)
        ax2 = ax.twiny()
        ax2.set_xlim(ax.get_xlim())
        ax2.set_xticks(ax.get_xticks())
        ax2.xaxis.set_major_formatter(dates.DateFormatter('%H:%M', tz=tzinfo))
        plt.setp(ax2.get_xticklabels(), rotation=-30, ha='right')
        ax2.set_xlabel("Time from {} [{}]".format(
            min(times).to_datetime(tzinfo).date(),
            timezone))

    if opts.verbose:
        # Write airmass table to stdout.
        times.format = 'isot'
        table = Table(masked=True)
        table['time'] = times
        table['sun_alt'] = np.ma.masked_greater_equal(
            observer.sun_altaz(times).alt, 0)
        table['sun_alt'].format = lambda x: '{}'.format(int(np.round(x)))
        for p, data in sorted(zip(percentiles, airmass)):
            table[str(p)] = np.ma.masked_invalid(data)
            table[str(p)].format = lambda x: '{:.01f}'.format(np.around(x, 1))
        table.write(sys.stdout, format='ascii.fixed_width')

    # Show or save output.
    opts.output()
