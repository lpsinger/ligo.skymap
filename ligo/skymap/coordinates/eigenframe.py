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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""
Astropy coordinate frame for eigendecomposition of a cloud of points or a 3D
sky map.
"""

from astropy.coordinates import (
    BaseCoordinateFrame, CartesianRepresentation, DynamicMatrixTransform,
    frame_transform_graph, ICRS, SphericalRepresentation)
from astropy.coordinates import CartesianRepresentationAttribute
from astropy.units import dimensionless_unscaled
import numpy as np

from ..distance import principal_axes

__all__ = ('EigenFrame',)


class EigenFrame(BaseCoordinateFrame):
    """A coordinate frame that has its axes aligned with the principal
    components of a cloud of points.
    """

    e_x = CartesianRepresentationAttribute(
        default=CartesianRepresentation(1, 0, 0, unit=dimensionless_unscaled),
        unit=dimensionless_unscaled)
    e_y = CartesianRepresentationAttribute(
        default=CartesianRepresentation(0, 1, 0, unit=dimensionless_unscaled),
        unit=dimensionless_unscaled)
    e_z = CartesianRepresentationAttribute(
        default=CartesianRepresentation(0, 0, 1, unit=dimensionless_unscaled),
        unit=dimensionless_unscaled)

    default_representation = SphericalRepresentation

    @classmethod
    def for_coords(cls, coords):
        """Create a coordinate frame that has its axes aligned with the
        principal components of a cloud of points.

        Parameters
        ----------
        coords : `astropy.coordinates.SkyCoord`
            A cloud of points

        Returns
        -------
        frame : `EigenFrame`
            A new coordinate frame

        """
        v = coords.icrs.cartesian.xyz.value
        _, r = np.linalg.eigh(np.dot(v, v.T))
        r = r[:, ::-1]  # Order by descending eigenvalue
        e_x, e_y, e_z = CartesianRepresentation(r, unit=dimensionless_unscaled)
        return cls(e_x=e_x, e_y=e_y, e_z=e_z)

    @classmethod
    def for_skymap(cls, prob, distmu, distsigma, nest=False):
        """Create a coordinate frame that has its axes aligned with the
        principal components of a 3D sky map.

        Parameters
        ----------
        prob : `numpy.ndarray`
            Marginal probability (pix^-2)
        distmu : `numpy.ndarray`
            Distance location parameter (Mpc)
        distsigma : `numpy.ndarray`
            Distance scale parameter (Mpc)
        distnorm : `numpy.ndarray`
            Distance normalization factor (Mpc^-2)
        nest : bool, default=False
            Indicates whether the input sky map is in nested rather than
            ring-indexed HEALPix coordinates (default: ring).

        Returns
        -------
        frame : `EigenFrame`
            A new coordinate frame

        """
        r = principal_axes(prob, distmu, distsigma, nest=nest)
        r = r[:, ::-1]  # Order by descending eigenvalue
        e_x, e_y, e_z = CartesianRepresentation(r, unit=dimensionless_unscaled)
        return cls(e_x=e_x, e_y=e_y, e_z=e_z)


@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, EigenFrame)
def icrs_to_eigenframe(from_coo, to_frame):
    return np.vstack((to_frame.e_x.xyz.value,
                      to_frame.e_y.xyz.value,
                      to_frame.e_z.xyz.value))


@frame_transform_graph.transform(DynamicMatrixTransform, EigenFrame, ICRS)
def eigenframe_to_icrs(from_coo, to_frame):
    return np.column_stack((from_coo.e_x.xyz.value,
                            from_coo.e_y.xyz.value,
                            from_coo.e_z.xyz.value))
