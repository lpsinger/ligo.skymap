#
# Copyright (C) 2023-2024  Leo Singer
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import astropy_healpix as ah
import numpy as np
from reproject.utils import parse_output_projection
from reproject.healpix.utils import parse_coord_system

from ..healpix_tree import HEALPIX_MACHINE_ORDER, HEALPIX_MACHINE_NSIDE
from ..moc import uniq2nest


def reproject_from_healpix_moc(
    input_data, output_projection, shape_out=None
):
    """
    Reproject multiorder HEALPIX data to a standard projection.
    Adapted from :meth:`reproject.reproject_from_healpix`.

    Parameters
    ----------
    input_data : tuple
        A tuple consisting of the following:
        -   A multi-order HEALPix dataset stored as an Astropy table whose
            first column is called UNIQ and contains the NUNIQ pixel index.
            Every point on the unit sphere must be contained in exactly one
            pixel in the dataset.
        -   An instance of `~astropy.coordinates.BaseCoordinateFrame` or a
            string alias for a coordinate frame.
    output_projection : `~astropy.wcs.WCS` or `~astropy.io.fits.Header`
        The output projection, which can be either a `~astropy.wcs.WCS`
        or a `~astropy.io.fits.Header` instance.
    shape_out : tuple, optional
        If ``output_projection`` is a `~astropy.wcs.WCS` instance, the
        shape of the output data should be specified separately.

    Returns
    -------
    array_new : `~numpy.ndarray`
        The reprojected array.
    footprint : `~numpy.ndarray`
        Footprint of the input array in the output array. Values of 0 indicate
        no coverage or valid values in the input image, while values of 1
        indicate valid values.
    """
    array_in, coord_system_in = input_data
    coord_system_in = parse_coord_system(coord_system_in)
    wcs_out, shape_out = parse_output_projection(
        output_projection, shape_out=shape_out)

    # Look up lon, lat of pixels in reference system and convert celestial
    # coordinates
    yinds, xinds = np.indices(shape_out)
    world_in = wcs_out.pixel_to_world(xinds, yinds).transform_to(
        coord_system_in)
    world_in_cart = world_in.represent_as("cartesian").xyz.value

    hpx_in = ah.xyz_to_healpix(
        *world_in_cart, nside=HEALPIX_MACHINE_NSIDE, order='nested')
    order, hpx_data = uniq2nest(array_in['UNIQ'])
    hpx_data <<= 2 * (HEALPIX_MACHINE_ORDER - order)
    sorter = np.argsort(hpx_data)
    i = np.searchsorted(hpx_data, hpx_in, 'right', sorter=sorter) - 1
    data = array_in.columns[1][sorter][i]

    footprint = (hpx_in != -1).astype(float)

    return data, footprint
