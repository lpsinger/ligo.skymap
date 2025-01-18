#
# Copyright (C) 2019-2025  Leo Singer
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
"""Make an observability chart for a LIGO/Virgo/KAGRA probability sky map."""

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
        '--max-airmass', default=2.5, type=float, help='Maximum airmass')
    parser.add_argument(
        '--twilight', default='astronomical',
        choices=('astronomical', 'nautical', 'civil'),
        help='Twilight definition: astronomical (-18 degrees), '
        'nautical (-12 degrees), or civil (-6 degrees)')
    parser.add_argument(
        '--site', nargs='*', default=[], metavar='SITE', choices=site_names,
        help='Observatory site')
    parser.add_argument(
        '--help-site', action=HelpChoicesAction, choices=site_names)
    parser.add_argument(
        '--site-name', nargs='*', default=[], help='Observatory name.')
    parser.add_argument(
        '--site-longitude', nargs='*', default=[], metavar='DEG', type=float,
        help='Observatory longitude on the WGS84 ellipsoid.')
    parser.add_argument(
        '--site-latitude', nargs='*', default=[], metavar='DEG', type=float,
        help='Observatory latitude on the WGS84 ellipsoid.')
    parser.add_argument(
        '--site-height', nargs='*', default=[], metavar='METERS', type=float,
        help='Observatory height from the WGS84 ellipsoid.')
    return parser


def condition_secz(x):
    """Condition secz airmass formula: values <= 0 are below the horizon,
    which we map to infinite airmass.
    """
    return np.where(x <= 0, np.inf, x)


def main(args=None):
    p = parser()
    with p.parse_args(args) as opts:
        # Late imports
        from astroplan import (
            AirmassConstraint, AtNightConstraint, Observer,
            is_event_observable)
        from astropy.coordinates import EarthLocation, SkyCoord
        from astropy.time import Time
        from astropy import units as u
        from matplotlib import dates
        from matplotlib import pyplot as plt
        from tqdm import tqdm

        from ..io import fits
        from .. import moc
        from .. import plot  # noqa

        names = ('name', 'longitude', 'latitude', 'height')
        length0, *lengths = (
            len(getattr(opts, 'site_{}'.format(name))) for name in names)
        if not all(length0 == length for length in lengths):
            p.error(
                'these options require equal numbers of arguments: {}'.format(
                    ', '.join('--site-{}'.format(name) for name in names)))

        observers = [Observer.at_site(site) for site in opts.site]
        for name, lon, lat, height in zip(
                opts.site_name, opts.site_longitude, opts.site_latitude,
                opts.site_height):
            location = EarthLocation(
                lon=lon * u.deg,
                lat=lat * u.deg,
                height=(height or 0) * u.m)
            observers.append(Observer(location, name=name))
        observers = list(reversed(observers))

        m = fits.read_sky_map(opts.input.name, moc=True)

        t0 = Time(opts.time) if opts.time is not None else Time.now()
        times = t0 + np.linspace(0, 1) * u.day

        theta, phi = moc.uniq2ang(m['UNIQ'])
        coords = SkyCoord(phi, 0.5 * np.pi - theta, unit=u.rad)
        prob = np.asarray(moc.uniq2pixarea(m['UNIQ']) * m['PROBDENSITY'])

        constraints = [
            getattr(AtNightConstraint, 'twilight_{}'.format(opts.twilight))(),
            AirmassConstraint(opts.max_airmass)]

        fig = plt.figure()
        width, height = fig.get_size_inches()
        fig.set_size_inches(width, (len(observers) + 1) / 16 * width)
        ax = plt.axes()
        locator = dates.AutoDateLocator()
        formatter = dates.DateFormatter('%H:%M')
        ax.set_xlim([times[0].plot_date, times[-1].plot_date])
        ax.xaxis.set_major_formatter(formatter)
        ax.xaxis.set_major_locator(locator)
        ax.set_xlabel("Time from {0} [UTC]".format(min(times).datetime.date()))
        plt.setp(ax.get_xticklabels(), rotation=30, ha='right')
        ax.set_yticks(np.arange(len(observers)))
        ax.set_yticklabels([observer.name for observer in observers])
        ax.yaxis.set_tick_params(left=False)
        ax.grid(axis='x')
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)

        for i, observer in enumerate(tqdm(observers)):
            observable = 100 * np.dot(prob, is_event_observable(
                constraints, observer, coords, times))
            ax.contourf(
                times.plot_date, [i - 0.4, i + 0.4],
                np.tile(observable, (2, 1)), levels=np.arange(10, 110, 10),
                cmap=plt.get_cmap().reversed())

        plt.tight_layout()

        # Show or save output.
        opts.output()
