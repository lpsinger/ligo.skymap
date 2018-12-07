#
# Copyright (C) 2018  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""Astropy coordinate frames to visualize triangulation rings from pairs of
gravitational-wave detectors. These are useful for generating plots similar to
Fig. 2 of the GW150914 localization and follow-up paper [1]_.

Example
-------

.. plot::
   :context: reset
   :include-source:
   :align: center

    from astropy.coordinates import EarthLocation
    from astropy.time import Time
    from ligo.skymap.coordinates import DetectorFrame
    from ligo.skymap.io import read_sky_map
    import ligo.skymap.plot
    from matplotlib import pyplot as plt

    # Download GW150914 localization
    url = 'https://dcc.ligo.org/public/0122/P1500227/012/bayestar_gstlal_C01.fits.gz'
    m, meta = ligo.skymap.io.read_sky_map(url)

    # Plot sky map on an orthographic projection
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(
        111, projection='astro globe', center='130d -70d')
    ax.imshow_hpx(m, cmap='cylon')

    # Hide the original ('RA', 'Dec') ticks
    for coord in ax.coords:
        coord.set_ticks_visible(False)
        coord.set_ticklabel_visible(False)

    # Construct Hanford-Livingston detector frame at the time of the event
    frame = DetectorFrame(site_1=EarthLocation.of_site('H1'),
                          site_2=EarthLocation.of_site('L1'),
                          obstime=Time(meta['gps_time'], format='gps'))

    # Draw grid for detector frame
    ax.get_coords_overlay(frame).grid()

References
----------

.. [1] LSC/Virgo et al., 2016. "Localization and Broadband Follow-up of the
       Gravitational-wave Transient GW150914." ApJL 826, L13
       <https://doi.org/10.3847/2041-8205/826/1/L13>.
"""  # noqa: E501

from astropy.coordinates import (
    CartesianRepresentation, DynamicMatrixTransform, EarthLocation,
    EarthLocationAttribute, frame_transform_graph, ITRS,
    SphericalRepresentation)
from astropy.coordinates.matrix_utilities import matrix_transpose
from astropy import units as u
import lal
import numpy as np

__all__ = ('DetectorFrame',)

# Add gravitational-wave detectors to site registry
registry = EarthLocation._get_site_registry()
for detector in lal.CachedDetectors:
    names = [detector.frDetector.name, detector.frDetector.prefix]
    location = EarthLocation(*detector.location, unit=u.m)
    registry.add_site(names, location)
    del names, detector
del lal, registry


class DetectorFrame(ITRS):
    """A coordinate frames to visualize triangulation rings from pairs of
    gravitational-wave detectors."""

    site_1 = EarthLocationAttribute()
    site_2 = EarthLocationAttribute()

    default_representation = SphericalRepresentation


@frame_transform_graph.transform(DynamicMatrixTransform, ITRS, DetectorFrame)
def itrs_to_detectorframe(from_coo, to_frame):
    e_z = CartesianRepresentation(u.Quantity(to_frame.site_1.geocentric) -
                                  u.Quantity(to_frame.site_2.geocentric))
    e_z /= e_z.norm()
    e_x = CartesianRepresentation(0, 0, 1).cross(e_z)
    e_x /= e_x.norm()
    e_y = e_z.cross(e_x)

    return np.row_stack((e_x.xyz.value,
                         e_y.xyz.value,
                         e_z.xyz.value))


@frame_transform_graph.transform(DynamicMatrixTransform, DetectorFrame, ITRS)
def detectorframe_to_itrs(from_coo, to_frame):
    return matrix_transpose(itrs_to_detectorframe(to_frame, from_coo))
