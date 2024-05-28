#
# Copyright (C) 2017-2024  Leo Singer
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
"""Create a P-P plot to compare a posterior sample chain with a sky map."""

import re

import numpy as np

from . import ArgumentParser, FileType
from .matplotlib import get_figure_parser


def fmt(x, sigfigs, force_scientific=False):
    """Round and format a number in scientific notation."""
    places = sigfigs - int(np.floor(np.log10(x)))
    x_rounded = np.around(x, places)
    if places <= 0 and not force_scientific:
        return '{:d}'.format(int(x_rounded))
    else:
        s = ('{:.' + str(sigfigs) + 'e}').format(x_rounded)
        return re.sub(r'^(.*)e\+?(-?)0*(\d+)$', r'$\1 \\times 10^{\2\3}$', s)


def parser():
    parser = ArgumentParser(parents=[get_figure_parser()])
    parser.add_argument(
        'skymap', metavar='SKYMAP.fits[.gz]', type=FileType('rb'),
        help='FITS sky map file')
    parser.add_argument(
        'samples', metavar='SAMPLES.hdf5', type=FileType('rb'),
        help='HDF5 posterior sample chain file')
    parser.add_argument(
        '-m', '--max-points', type=int,
        help='Maximum number of posterior sample points '
        '[default: all of them]')
    parser.add_argument(
        '-p', '--contour', default=[], nargs='+', type=float,
        metavar='PERCENT',
        help='Report the areas and volumes within the smallest contours '
        'containing this much probability.')
    return parser


def main(args=None):
    with parser().parse_args(args) as args:
        # Late imports.
        from astropy.table import Table
        from matplotlib import pyplot as plt
        from scipy.interpolate import interp1d
        from .. import io
        from .. import plot as _  # noqa: F401
        from ..postprocess import find_injection_moc

        # Read input.
        skymap = io.read_sky_map(args.skymap.name, moc=True)
        chain = io.read_samples(args.samples.name)

        # If required, downselect to a smaller number of posterior samples.
        if args.max_points is not None:
            chain = Table(np.random.permutation(chain)[:args.max_points],
                          copy=False)

        # Calculate P-P plot.
        contours = np.asarray(args.contour)
        result = find_injection_moc(skymap,
                                    chain['ra'], chain['dec'], chain['dist'],
                                    contours=1e-2 * contours)

        # Make Matplotlib figure.
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111, projection='pp_plot')
        ax.add_diagonal()
        ax.add_series(result.searched_prob, label='R.A., Dec.')
        searched_area_func = interp1d(np.linspace(0, 1, len(chain)),
                                      np.sort(result.searched_area),
                                      bounds_error=False)
        if 'DISTMU' in skymap.colnames:
            ax.add_series(result.searched_prob_dist, label='Distance')
            ax.add_series(result.searched_prob_vol, label='Volume')
            searched_vol_func = interp1d(np.linspace(0, 1, len(chain)),
                                         np.sort(result.searched_vol),
                                         bounds_error=False)
        for p, area, vol in zip(
                args.contour, result.contour_areas, result.contour_vols):
            text = '{:g}%\n{} deg$^2$'.format(p, fmt(area, 2))
            if 'DISTMU' in skymap.colnames:
                text += '\n{} Mpc$^3$'.format(
                    fmt(vol, 2, force_scientific=True))
            ax.annotate(
                text, (1e-2 * p, 1e-2 * p), (0, -150),
                xycoords='data', textcoords='offset points',
                horizontalalignment='right', backgroundcolor='white',
                arrowprops=dict(connectionstyle='bar,angle=0,fraction=0',
                                arrowstyle='-|>', linewidth=2, color='black'))
            area = searched_area_func(1e-2 * p)
            text = '{:g}%\n{} deg$^2$'.format(p, fmt(area, 2))
            if 'DISTMU' in skymap.colnames:
                vol = searched_vol_func(1e-2 * p)
                text += '\n{} Mpc$^3$'.format(
                    fmt(vol, 2, force_scientific=True))
            ax.annotate(
                text, (1e-2 * p, 1e-2 * p), (-75, 0),
                xycoords='data', textcoords='offset points',
                horizontalalignment='right', verticalalignment='center',
                backgroundcolor='white',
                arrowprops=dict(connectionstyle='bar,angle=0,fraction=0',
                                arrowstyle='-|>', linewidth=2, color='black'))
        ax.set_xlabel('searched probability')
        ax.set_ylabel('cumulative fraction of posterior samples')
        ax.set_title(args.skymap.name)
        ax.legend()
        ax.grid()

        # Show or save output.
        args.output()
