#
# Copyright (C) 2013-2024  Giuseppe Greco, Leo Singer, and CDS team.
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
# along with this program. If not, see <http://www.gnu.org/licenses/>.
"""
Create a contour for a credible level of an all-sky probability map. The input
is a HEALPix FITS probability map. The output is a `Multi-Order Coverage (MOC)
<http://ivoa.net/documents/MOC/>`_ FITS file.
"""

from . import ArgumentParser, FileType


def parser():
    parser = ArgumentParser()

    parser.add_argument(
        '-o', '--output', metavar='FILE.fits', required=True,
        help='output file')
    parser.add_argument(
        '-c', '--contour', metavar='PERCENT', type=float, required=True,
        help='MOC region enclosing this percentage of probability \
              [range is 0-100]')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input multi-order or flatten \
                                      HEALPix FITS file')

    return parser


def main(args=None):
    p = parser()
    with parser().parse_args(args) as opts:
        import astropy_healpix as ah
        import astropy.units as u

        try:
            from mocpy import MOC
        except ImportError:
            p.error('This command-line tool requires mocpy >= 0.8.2. '
                    'Please install it by running "pip install mocpy".')

        from ..io import read_sky_map

        # Read multi-order sky map
        skymap = read_sky_map(opts.input.name, moc=True)

        uniq = skymap['UNIQ']
        probdensity = skymap['PROBDENSITY']

        level, ipix = ah.uniq_to_level_ipix(uniq)
        area = ah.nside_to_pixel_area(
            ah.level_to_nside(level)).to_value(u.steradian)

        prob = probdensity * area

        # Create MOC
        contour_decimal = opts.contour / 100
        moc = MOC.from_valued_healpix_cells(
            uniq, prob, max_depth=level.max(),
            cumul_from=0.0, cumul_to=contour_decimal)

        # Write MOC
        moc.write(opts.output, format='fits', overwrite=True)
