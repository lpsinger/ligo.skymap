#
# Copyright (C) 2019-2024  Leo Singer
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
List most likely constellations for a localization.

Just for fun, and for public outreach purposes.
"""

from . import ArgumentParser, FileType


def parser():
    parser = ArgumentParser()
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    parser.add_argument(
        '-o', '--output', metavar='OUT.dat', type=FileType('w'), default='-',
        help='Name of output file')
    return parser


def main(args=None):
    with parser().parse_args(args) as opts:
        # Late imports
        from ..io import fits
        import astropy_healpix as ah
        from astropy.coordinates import SkyCoord
        from astropy.table import Table
        from astropy import units as u
        import healpy as hp
        import numpy as np

        prob, meta = fits.read_sky_map(opts.input.name, nest=None)
        npix = len(prob)
        nside = ah.npix_to_nside(npix)
        ipix = np.arange(npix)
        ra, dec = hp.pix2ang(nside, ipix, lonlat=True, nest=meta['nest'])
        coord = SkyCoord(ra * u.deg, dec * u.deg)
        table = Table(
            {'prob': prob, 'constellation': coord.get_constellation()},
            copy=False)
        table = table.group_by('constellation').groups.aggregate(np.sum)
        table.sort('prob')
        table.reverse()
        table.write(opts.output, format='ascii.tab')
