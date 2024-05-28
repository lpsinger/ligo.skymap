#
# Copyright (C) 2018-2024  Leo Singer
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
"""Convert a HEALPix FITS file from multi-resolution UNIQ indexing to the more
common IMPLICIT indexing.
"""

from . import ArgumentParser, FileType


def parser():
    parser = ArgumentParser()
    parser.add_argument('input', metavar='INPUT.fits[.gz]',
                        type=FileType('rb'), help='Input FITS file')
    parser.add_argument('output', metavar='OUTPUT.fits',
                        type=FileType('wb'), help='Output FITS file')
    parser.add_argument('--nside', type=int, help='Output HEALPix resolution')
    return parser


def main(args=None):
    with parser().parse_args(args) as args:
        import logging
        import warnings
        import astropy_healpix as ah
        from astropy.io import fits
        from ..io import read_sky_map, write_sky_map
        from ..bayestar import rasterize

        log = logging.getLogger()

        if args.nside is None:
            order = None
        else:
            order = ah.nside_to_level(args.nside)

        log.info('reading FITS file %s', args.input.name)
        hdus = fits.open(args.input)
        ordering = hdus[1].header['ORDERING']
        expected_ordering = 'NUNIQ'
        if ordering != expected_ordering:
            msg = 'Expected the FITS file {} to have ordering {}, but it is {}'
            warnings.warn(
                msg.format(args.input.name, expected_ordering, ordering))
        log.debug('converting original FITS file to Astropy table')
        table = read_sky_map(hdus, moc=True)
        log.debug('flattening HEALPix tree')
        table = rasterize(table, order=order)
        log.info('writing FITS file %s', args.output.name)
        write_sky_map(args.output.name, table, nest=True)
        log.debug('done')
