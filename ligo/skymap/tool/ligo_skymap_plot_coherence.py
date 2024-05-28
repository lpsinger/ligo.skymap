#
# Copyright (C) 2011-2024  Leo Singer
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
Show a sky map's Bayes factor for coherence vs. incoherence as a bullet chart.
"""

from . import ArgumentParser, FileType
from .matplotlib import get_figure_parser


def parser():
    parser = ArgumentParser(parents=[get_figure_parser()])
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    parser.set_defaults(colormap='RdYlBu')
    return parser


def main(args=None):
    with parser().parse_args(args) as opts:
        # Late imports
        from astropy.io import fits
        import numpy as np
        from ..plot import plot_bayes_factor

        header = fits.getheader(opts.input, 1)
        logb = header['LOGBCI']
        objid = header.get('OBJECT')

        title = 'Coherence'
        if objid:
            title += f' of {objid}'
        logb_string = np.format_float_positional(logb, 1, trim='0', sign=True)
        title += fr' $[\ln\,B = {logb_string}]$'

        plot_bayes_factor(logb, title=title, palette=opts.colormap)

        # Show or save output.
        opts.output()
