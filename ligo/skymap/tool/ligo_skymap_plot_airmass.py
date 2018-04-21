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


def main(args=None):
    opts = parser().parse_args(args)

    # Late imports
    from astroplan import Observer
    from astroplan.plots import plot_airmass
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from matplotlib.patches import Patch
    from matplotlib import pyplot as plt
    import numpy as np

    from ..io import fits
    from .. import moc
    from .. import plot  # noqa
    from ..util import quantile

    m = fits.read_sky_map(opts.input.name, moc=True)

    # Make an empty airmass chart.
    # FIXME: have to add a dummy target until
    # https://github.com/astropy/astroplan/pull/349
    # is in a release of astroplan
    observer = Observer.at_site(opts.site)
    ax = plot_airmass([SkyCoord(0, 0, unit='rad')], observer,
                      Time(opts.time) if opts.time is not None else Time.now(),
                      brightness_shading=True)
    del ax.lines[:]

    # Remove the fake source and determine times that were used for the plot.
    t = Time(np.linspace(*ax.get_xlim()), format='plot_date')

    # Calculate airmass at every point in time and for every HEALPix pixel.
    theta, phi = moc.uniq2ang(m['UNIQ'])
    coords = SkyCoord(phi, 0.5 * np.pi - theta, unit='rad')
    prob = moc.uniq2pixarea(m['UNIQ']) * m['PROBDENSITY']
    airmass = coords.transform_to(observer.altaz(t[:, np.newaxis])).secz
    airmass[airmass <= 0] = np.inf

    # Plot probability-weighted airmass contours.
    cmap = plt.get_cmap()
    alphas = np.arange(0.45, 0.0, -0.05)
    for alpha in alphas:
        lo, hi = quantile(airmass, [0.5 - alpha, 0.5 + alpha],
                          weights=prob, axis=1)
        ax.fill_between(t.plot_date, lo, hi, color=cmap(2 * alpha), zorder=2)

    ax.legend(
        [Patch(facecolor=cmap(2 * alpha)) for alpha in alphas],
        ['{}%'.format(int(200 * alpha)) for alpha in alphas])
    ax.set_title('{} from {}'.format(m.meta['objid'], observer.name))

    plt.subplots_adjust(bottom=0.2)

    # Show or save output.
    opts.output()
