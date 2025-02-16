#
# Copyright (C) 2012-2024  Leo Singer
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
Axes subclasses for astronomical mapmaking.

This module adds several :class:`astropy.visualization.wcsaxes.WCSAxes`
subclasses to the Matplotlib projection registry. The projections have names of
the form :samp:`{astro_or_geo_or_galactic} [{lon_units}] {projection}`.

:samp:`{astro_or_geo_or_galactic}` may be ``astro``, ``geo``, or ``galactic``.
It controls the reference frame, either celestial (ICRS), terrestrial (ITRS),
or galactic.

:samp:`{lon_units}` may be ``hours`` or ``degrees``. It controls the units of
the longitude axis. If omitted, ``astro`` implies ``hours`` and ``geo`` implies
degrees.

:samp:`{projection}` may be any of the following:

* ``aitoff`` for the Aitoff all-sky projection

* ``mollweide`` for the Mollweide all-sky projection

* ``globe`` for an orthographic projection, like the three-dimensional view of
  the Earth from a distant satellite

* ``zoom`` for a gnomonic projection suitable for visualizing small zoomed-in
  patches

All projections support the ``center`` argument, while some support additional
arguments. The ``globe`` projections also support the ``rotate`` argument, and
the ``zoom`` projections also supports the ``radius`` and ``rotate`` arguments.

Examples
--------
.. plot::
   :context: reset
   :include-source:
   :align: center

    import ligo.skymap.plot
    from matplotlib import pyplot as plt
    ax = plt.axes(projection='astro hours mollweide')
    ax.grid()

.. plot::
   :context: reset
   :include-source:
   :align: center

    import ligo.skymap.plot
    from matplotlib import pyplot as plt
    ax = plt.axes(projection='geo aitoff')
    ax.grid()

.. plot::
   :context: reset
   :include-source:
   :align: center

    import ligo.skymap.plot
    from matplotlib import pyplot as plt
    ax = plt.axes(projection='astro zoom',
                  center='5h -32d', radius='5 deg', rotate='20 deg')
    ax.grid()

.. plot::
   :context: reset
   :include-source:
   :align: center

    import ligo.skymap.plot
    from matplotlib import pyplot as plt
    ax = plt.axes(projection='geo globe', center='-50d +23d')
    ax.grid()

Insets
------
You can use insets to link zoom-in views between axes. There are two supported
styles of insets: rectangular and circular (loupe). The example below shows
both kinds of insets.

.. plot::
   :context: reset
   :include-source:
   :align: center

    import ligo.skymap.plot
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(9, 4), dpi=100)

    ax_globe = plt.axes(
        [0.1, 0.1, 0.8, 0.8],
        projection='astro degrees globe',
        center='120d +23d')

    ax_zoom_rect = plt.axes(
        [0.0, 0.2, 0.4, 0.4],
        projection='astro degrees zoom',
        center='150d +30d',
        radius='9 deg')

    ax_zoom_circle = plt.axes(
        [0.55, 0.1, 0.6, 0.6],
        projection='astro degrees zoom',
        center='120d +10d',
        radius='5 deg')

    ax_globe.mark_inset_axes(ax_zoom_rect)
    ax_globe.connect_inset_axes(ax_zoom_rect, 'upper left')
    ax_globe.connect_inset_axes(ax_zoom_rect, 'lower right')

    ax_globe.mark_inset_circle(ax_zoom_circle, '120d +10d', '4 deg')
    ax_globe.connect_inset_circle(ax_zoom_circle, '120d +10d', '4 deg')

    ax_globe.grid()
    ax_zoom_rect.grid()
    ax_zoom_circle.grid()

    for ax in [ax_globe, ax_zoom_rect, ax_zoom_circle]:
        ax.set_facecolor('none')
        for key in ['ra', 'dec']:
            ax.coords[key].set_auto_axislabel(False)

Complete Example
----------------
The following example demonstrates most of the features of this module.

.. plot::
   :context: reset
   :include-source:
   :align: center

    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy import units as u
    import ligo.skymap.plot
    from matplotlib import pyplot as plt

    url = 'https://dcc.ligo.org/public/0146/G1701985/001/bayestar_no_virgo.fits.gz'
    center = SkyCoord.from_name('NGC 4993')

    fig = plt.figure(figsize=(4, 4), dpi=100)

    ax = plt.axes(
        [0.05, 0.05, 0.9, 0.9],
        projection='astro globe',
        center=center)

    ax_inset = plt.axes(
        [0.59, 0.3, 0.4, 0.4],
        projection='astro zoom',
        center=center,
        radius=10*u.deg)

    for key in ['ra', 'dec']:
        ax_inset.coords[key].set_ticklabel_visible(False)
        ax_inset.coords[key].set_ticks_visible(False)
    ax.grid()
    ax.mark_inset_axes(ax_inset)
    ax.connect_inset_axes(ax_inset, 'upper left')
    ax.connect_inset_axes(ax_inset, 'lower left')
    ax_inset.scalebar((0.1, 0.1), 5 * u.deg).label()
    ax_inset.compass(0.9, 0.1, 0.2)

    ax.imshow_hpx(url, cmap='cylon')
    ax_inset.imshow_hpx(url, cmap='cylon')
    ax_inset.plot(
        center.ra.deg, center.dec.deg,
        transform=ax_inset.get_transform('world'),
        marker=ligo.skymap.plot.reticle(),
        markersize=30,
        markeredgewidth=3)

"""  # noqa: E501
from itertools import product
from warnings import warn

from astropy.coordinates import SkyCoord, UnitSphericalRepresentation
from astropy.io.fits import Header
from astropy.time import Time
from astropy.visualization.wcsaxes import SphericalCircle, WCSAxes
from astropy.visualization.wcsaxes.formatter_locator import (
    AngleFormatterLocator)
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from astropy import units as u
from matplotlib import rcParams
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, PathPatch
from matplotlib.path import Path
from matplotlib.projections import projection_registry
import numpy as np
from reproject import reproject_from_healpix
from scipy.optimize import minimize_scalar

from .angle import reference_angle_deg, wrapped_angle_deg
from .reproject_from_healpix_moc import reproject_from_healpix_moc

__all__ = ['AutoScaledWCSAxes', 'ScaleBar']


class WCSInsetPatch(PathPatch):
    """Subclass of `matplotlib.patches.PathPatch` for marking the outline of
    one `astropy.visualization.wcsaxes.WCSAxes` inside another.
    """

    def __init__(self, ax, *args, **kwargs):
        self._ax = ax
        super().__init__(
            None, *args, fill=False,
            edgecolor=ax.coords.frame.get_color(),
            linewidth=ax.coords.frame.get_linewidth(),
            **kwargs)

    def get_path(self):
        frame = self._ax.coords.frame
        return frame.patch.get_path().interpolated(50).transformed(
            frame.transform)


class WCSInsetConnectionPatch(ConnectionPatch):
    """Patch to connect an inset WCS axes inside another WCS axes.

    Notes
    -----
    FIXME: This class assumes that the projection of the circle in figure-inch
    coordinates *is* a circle. It will have noticeable artifacts if the
    projection is very distorted."""

    _corners_map = {1: 3, 2: 1, 3: 0, 4: 2}

    def __init__(self, ax, ax_inset, loc, **kwargs):
        try:
            loc = AnchoredOffsetbox.codes[loc]
        except KeyError:
            loc = int(loc)
        corners = ax_inset.viewLim.corners()
        transform = (ax_inset.coords.frame.transform +
                     ax.coords.frame.transform.inverted())
        xy_inset = corners[self._corners_map[loc]]
        xy = transform.transform_point(xy_inset)
        super().__init__(
            xy, xy_inset, 'data', 'data', axesA=ax, axesB=ax_inset,
            color=ax_inset.coords.frame.get_color(),
            linewidth=ax_inset.coords.frame.get_linewidth(),
            **kwargs)


class WCSCircleInsetConnectionPatch(PathPatch):
    """Patch to connect a circular inset WCS axes inside another WCS axes."""

    def __init__(self, ax1, ax2, coord, radius, sign, *args, **kwargs):
        self._axs = (ax1, ax2)
        self._coord = coord.icrs
        self._radius = radius
        self._sign = sign
        super().__init__(None, *args, **kwargs, clip_on=False, transform=None)

    def get_path(self):
        # Calculate the position and radius of the inset in figure-inch
        # coordinates.
        offset = self._coord.directional_offset_by(0 * u.deg, self._radius)
        transforms = [ax.get_transform('world') for ax in self._axs]
        centers = np.asarray([
            tx.transform_point((self._coord.ra.deg, self._coord.dec.deg))
            for tx in transforms])
        offsets = np.asarray([
            tx.transform_point((offset.ra.deg, offset.dec.deg))
            for tx in transforms])

        # Plot outer tangents.
        r0, r1 = np.sqrt(np.sum(np.square(centers - offsets), axis=-1))
        dx, dy = np.diff(centers, axis=0).ravel()
        gamma = -np.arctan(dy / dx)
        beta = np.arcsin((r1 - r0) / np.sqrt(np.square(dx) + np.square(dy)))
        alpha = gamma - self._sign * beta
        p0 = centers[0] + self._sign * np.asarray([
            r0 * np.sin(alpha), r0 * np.cos(alpha)])
        p1 = centers[1] + self._sign * np.asarray([
            r1 * np.sin(alpha), r1 * np.cos(alpha)])
        return Path(np.vstack((p0, p1)), np.asarray([
            Path.MOVETO, Path.LINETO]))


class AutoScaledWCSAxes(WCSAxes):
    """Axes base class. The pixel scale is adjusted to the DPI of the image,
    and there are a variety of convenience methods.
    """

    name = 'astro wcs'

    def __init__(self, *args, header, obstime=None, **kwargs):
        super().__init__(*args, aspect=1, **kwargs)
        h = Header(header, copy=True)
        naxis1 = h['NAXIS1']
        naxis2 = h['NAXIS2']
        scale = min(self.bbox.width / naxis1, self.bbox.height / naxis2)
        h['NAXIS1'] = int(np.ceil(naxis1 * scale))
        h['NAXIS2'] = int(np.ceil(naxis2 * scale))
        scale1 = h['NAXIS1'] / naxis1
        scale2 = h['NAXIS2'] / naxis2
        h['CRPIX1'] = (h['CRPIX1'] - 1) * (h['NAXIS1'] - 1) / (naxis1 - 1) + 1
        h['CRPIX2'] = (h['CRPIX2'] - 1) * (h['NAXIS2'] - 1) / (naxis2 - 1) + 1
        h['CDELT1'] /= scale1
        h['CDELT2'] /= scale2
        if obstime is not None:
            h['MJD-OBS'] = Time(obstime).utc.mjd
            h['DATE-OBS'] = Time(obstime).utc.isot
        self.reset_wcs(WCS(h))
        self.set_xlim(-0.5, h['NAXIS1'] - 0.5)
        self.set_ylim(-0.5, h['NAXIS2'] - 0.5)
        self._header = h

    @property
    def header(self):
        return self._header

    def mark_inset_axes(self, ax, *args, **kwargs):
        """Outline the footprint of another WCSAxes inside this one.

        Parameters
        ----------
        ax : `astropy.visualization.wcsaxes.WCSAxes`
            The other axes.

        Other parameters
        ----------------
        args :
            Extra arguments for `matplotlib.patches.PathPatch`
        kwargs :
            Extra keyword arguments for `matplotlib.patches.PathPatch`

        Returns
        -------
        patch : `matplotlib.patches.PathPatch`

        """
        return self.add_patch(WCSInsetPatch(
            ax, *args, transform=self.get_transform('world'), **kwargs))

    def mark_inset_circle(self, ax, center, radius, *args, **kwargs):
        """Outline a circle in this and another Axes to create a loupe.

        Parameters
        ----------
        ax : `astropy.visualization.wcsaxes.WCSAxes`
            The other axes.
        coord : `astropy.coordinates.SkyCoord`
            The center of the circle.
        radius : `astropy.units.Quantity`
            The radius of the circle in units that are compatible with degrees.

        Other parameters
        ----------------
        args :
            Extra arguments for `matplotlib.patches.PathPatch`
        kwargs :
            Extra keyword arguments for `matplotlib.patches.PathPatch`

        Returns
        -------
        patch1 : `matplotlib.patches.PathPatch`
            The outline of the circle in these Axes.
        patch2 : `matplotlib.patches.PathPatch`
            The outline of the circle in the other Axes.
        """
        center = SkyCoord(
            center, representation_type=UnitSphericalRepresentation).icrs
        radius = u.Quantity(radius)
        args = ((center.ra, center.dec), radius, *args)
        kwargs = {'facecolor': 'none',
                  'edgecolor': rcParams['axes.edgecolor'],
                  'linewidth': rcParams['axes.linewidth'],
                  **kwargs}
        for ax in (self, ax):
            ax.add_patch(SphericalCircle(*args, **kwargs,
                                         transform=ax.get_transform('world')))

    def connect_inset_axes(self, ax, loc, *args, **kwargs):
        """Connect a corner of another WCSAxes to the matching point inside
        this one.

        Parameters
        ----------
        ax : `astropy.visualization.wcsaxes.WCSAxes`
            The other axes.
        loc : int, str
            Which corner to connect. For valid values, see
            `matplotlib.offsetbox.AnchoredOffsetbox`.

        Other parameters
        ----------------
        args :
            Extra arguments for `matplotlib.patches.ConnectionPatch`
        kwargs :
            Extra keyword arguments for `matplotlib.patches.ConnectionPatch`

        Returns
        -------
        patch : `matplotlib.patches.ConnectionPatch`

        """
        return self.add_patch(WCSInsetConnectionPatch(
            self, ax, loc, *args, **kwargs))

    def connect_inset_circle(self, ax, center, radius, *args, **kwargs):
        """Connect a circle in this and another Axes to create a loupe.

        Parameters
        ----------
        ax : `astropy.visualization.wcsaxes.WCSAxes`
            The other axes.
        coord : `astropy.coordinates.SkyCoord`
            The center of the circle.
        radius : `astropy.units.Quantity`
            The radius of the circle in units that are compatible with degrees.

        Other parameters
        ----------------
        args :
            Extra arguments for `matplotlib.patches.PathPatch`
        kwargs :
            Extra keyword arguments for `matplotlib.patches.PathPatch`

        Returns
        -------
        patch1, patch2 : `matplotlib.patches.ConnectionPatch`
            The two connecting patches.
        """
        center = SkyCoord(
            center, representation_type=UnitSphericalRepresentation).icrs
        radius = u.Quantity(radius)
        kwargs = {'color': rcParams['axes.edgecolor'],
                  'linewidth': rcParams['axes.linewidth'],
                  **kwargs}
        for sign in (-1, 1):
            self.add_patch(WCSCircleInsetConnectionPatch(
                self, ax, center, radius, sign, *args, **kwargs))

    def compass(self, x, y, size):
        """Add a compass to indicate the north and east directions.

        Parameters
        ----------
        x, y : float
            Position of compass vertex in axes coordinates.
        size : float
            Size of compass in axes coordinates.

        """
        xy = x, y
        scale = self.wcs.pixel_scale_matrix
        scale /= np.sqrt(np.abs(np.linalg.det(scale)))
        return [self.annotate(label, xy, xy + size * n,
                              self.transAxes, self.transAxes,
                              ha='center', va='center',
                              arrowprops=dict(arrowstyle='<-',
                                              shrinkA=0.0, shrinkB=0.0))
                for n, label, ha, va in zip(scale, 'EN',
                                            ['right', 'center'],
                                            ['center', 'bottom'])]

    def scalebar(self, *args, **kwargs):
        """Add scale bar.

        Parameters
        ----------
        xy : tuple
            The axes coordinates of the scale bar.
        length : `astropy.units.Quantity`
            The length of the scale bar in angle-compatible units.

        Other parameters
        ----------------
        args :
            Extra arguments for `matplotlib.patches.FancyArrowPatch`
        kwargs :
            Extra keyword arguments for `matplotlib.patches.FancyArrowPatch`

        Returns
        -------
        patch : `matplotlib.patches.FancyArrowPatch`

        """
        return self.add_patch(ScaleBar(self, *args, **kwargs))

    def _reproject_hpx(self, data, hdu_in=None, order='bilinear',
                       nested=False, field=0, smooth=None):
        if isinstance(data, np.ndarray):
            data = (data, self.header['RADESYS'])

        # It's normal for reproject_from_healpix to produce some Numpy invalid
        # value warnings for points that land outside the projection.
        with np.errstate(invalid='ignore'):
            try:
                # Check if the input is a multiorder sky map
                data[0]['UNIQ']
            except (IndexError, KeyError, TypeError):
                img, mask = reproject_from_healpix(
                    data, self.header, hdu_in=hdu_in, order=order,
                    nested=nested, field=field)
            else:
                if order != 'nearest-neighbor':
                    warn('You requested bilinear interpolation of a '
                         'multi-order sky map, but only nearest-neighbor '
                         'interpolation is currently spported', UserWarning)
                img, mask = reproject_from_healpix_moc(data, self.header)
        img = np.ma.array(img, mask=~mask.astype(bool))

        if smooth is not None:
            # Infrequently used imports
            from astropy.convolution import convolve_fft, Gaussian2DKernel

            pixsize = np.mean(np.abs(self.wcs.wcs.cdelt)) * u.deg
            smooth = (smooth / pixsize).to(u.dimensionless_unscaled).value
            kernel = Gaussian2DKernel(smooth)
            # Ignore divide by zero warnings for pixels that have no valid
            # neighbors.
            with np.errstate(invalid='ignore'):
                img = convolve_fft(img, kernel, fill_value=np.nan)

        return img

    def contour_hpx(self, data, hdu_in=None, order='bilinear', nested=False,
                    field=0, smooth=None, **kwargs):
        """Add contour levels for a HEALPix data set.

        Parameters
        ----------
        data : `numpy.ndarray` or str or `~astropy.io.fits.TableHDU` or `~astropy.io.fits.BinTableHDU` or tuple
            The HEALPix data set. If this is a `numpy.ndarray`, then it is
            interpreted as the HEALPix array in the same coordinate system as
            the axes. Otherwise, the input data can be any type that is
            understood by `reproject.reproject_from_healpix`.
        smooth : `astropy.units.Quantity`, optional
            An optional smoothing length in angle-compatible units.

        Other parameters
        ----------------
        hdu_in, order, nested, field, smooth :
            Extra arguments for `reproject.reproject_from_healpix`
        kwargs :
            Extra keyword arguments for `matplotlib.axes.Axes.contour`

        Returns
        -------
        countours : `matplotlib.contour.QuadContourSet`

        """  # noqa: E501
        img = self._reproject_hpx(data, hdu_in=hdu_in, order=order,
                                  nested=nested, field=field, smooth=smooth)
        return self.contour(img, **kwargs)

    def contourf_hpx(self, data, hdu_in=None, order='bilinear', nested=False,
                     field=0, smooth=None, **kwargs):
        """Add filled contour levels for a HEALPix data set.

        Parameters
        ----------
        data : `numpy.ndarray` or str or `~astropy.io.fits.TableHDU` or `~astropy.io.fits.BinTableHDU` or tuple
            The HEALPix data set. If this is a `numpy.ndarray`, then it is
            interpreted as the HEALPix array in the same coordinate system as
            the axes. Otherwise, the input data can be any type that is
            understood by `reproject.reproject_from_healpix`.
        smooth : `astropy.units.Quantity`, optional
            An optional smoothing length in angle-compatible units.

        Other parameters
        ----------------
        hdu_in, order, nested, field, smooth :
            Extra arguments for `reproject.reproject_from_healpix`
        kwargs :
            Extra keyword arguments for `matplotlib.axes.Axes.contour`

        Returns
        -------
        contours : `matplotlib.contour.QuadContourSet`

        """  # noqa: E501
        img = self._reproject_hpx(data, hdu_in=hdu_in, order=order,
                                  nested=nested, field=field, smooth=smooth)
        return self.contourf(img, **kwargs)

    def imshow_hpx(self, data, hdu_in=None, order='bilinear', nested=False,
                   field=0, smooth=None, **kwargs):
        """Add an image for a HEALPix data set.

        Parameters
        ----------
        data : `numpy.ndarray` or str or `~astropy.io.fits.TableHDU` or `~astropy.io.fits.BinTableHDU` or tuple
            The HEALPix data set. If this is a `numpy.ndarray`, then it is
            interpreted as the HEALPix array in the same coordinate system as
            the axes. Otherwise, the input data can be any type that is
            understood by `reproject.reproject_from_healpix`.
        smooth : `astropy.units.Quantity`, optional
            An optional smoothing length in angle-compatible units.

        Other parameters
        ----------------
        hdu_in, order, nested, field, smooth :
            Extra arguments for `reproject.reproject_from_healpix`
        kwargs :
            Extra keyword arguments for `matplotlib.axes.Axes.contour`

        Returns
        -------
        image : `matplotlib.image.AxesImage`

        """  # noqa: E501
        img = self._reproject_hpx(data, hdu_in=hdu_in, order=order,
                                  nested=nested, field=field, smooth=smooth)
        return self.imshow(img, **kwargs)


class ScaleBar(FancyArrowPatch):

    def _func(self, dx, x, y):
        p1, p2 = self._transAxesToWorld.transform([[x, y], [x + dx, y]])
        p1 = SkyCoord(*p1, unit=u.deg)
        p2 = SkyCoord(*p2, unit=u.deg)
        return np.square((p1.separation(p2) - self._length).value)

    def __init__(self, ax, xy, length, *args, **kwargs):
        x, y = xy
        self._ax = ax
        self._length = u.Quantity(length)
        self._transAxesToWorld = (
            (ax.transAxes - ax.transData) + ax.coords.frame.transform)
        dx = minimize_scalar(
            self._func, args=xy, bounds=[0, 1 - x], method='bounded').x
        custom_kwargs = kwargs
        kwargs = dict(
            capstyle='round',
            color='black',
            linewidth=rcParams['lines.linewidth'],
        )
        kwargs.update(custom_kwargs)
        super().__init__(
            xy, (x + dx, y),
            *args,
            arrowstyle='-',
            shrinkA=0.0,
            shrinkB=0.0,
            transform=ax.transAxes,
            **kwargs)

    def label(self, **kwargs):
        (x0, y), (x1, _) = self._posA_posB
        s = ' {0.value:g}{0.unit:unicode}'.format(self._length)
        return self._ax.text(
            0.5 * (x0 + x1), y, s,
            ha='center', va='bottom', transform=self._ax.transAxes, **kwargs)


class Astro:
    _crval1 = 180
    _xcoord = 'RA--'
    _ycoord = 'DEC-'
    _radesys = 'ICRS'


class GeoAngleFormatterLocator(AngleFormatterLocator):

    def formatter(self, values, spacing):
        return super().formatter(
            reference_angle_deg(values.to(u.deg).value) * u.deg, spacing)


class Geo:
    _crval1 = 0
    _radesys = 'ITRS'
    _xcoord = 'TLON'
    _ycoord = 'TLAT'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.invert_xaxis()
        fl = self.coords[0]._formatter_locator
        self.coords[0]._formatter_locator = GeoAngleFormatterLocator(
            values=fl.values,
            number=fl.number,
            spacing=fl.spacing,
            format=fl.format,
            format_unit=fl.format_unit)


class GalacticAngleFormatterLocator(AngleFormatterLocator):

    def formatter(self, values, spacing):
        return super().formatter(
            wrapped_angle_deg(values.to(u.deg).value) * u.deg, spacing)


class Galactic:
    _crval1 = 0
    _radesys = 'GALACTIC'
    _xcoord = 'GLON'
    _ycoord = 'GLAT'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        fl = self.coords[0]._formatter_locator
        self.coords[0]._formatter_locator = GalacticAngleFormatterLocator(
            values=fl.values,
            number=fl.number,
            spacing=fl.spacing,
            format=fl.format,
            format_unit=fl.format_unit)


class Degrees:
    """WCS axes with longitude axis in degrees."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.coords[0].set_format_unit(u.degree)


class Hours:
    """WCS axes with longitude axis in hour angle."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.coords[0].set_format_unit(u.hourangle)


class Globe(AutoScaledWCSAxes):

    def __init__(self, *args, center='0d 0d', rotate=None, **kwargs):
        center = SkyCoord(
            center, representation_type=UnitSphericalRepresentation).icrs
        header = {
            'NAXIS': 2,
            'NAXIS1': 180,
            'NAXIS2': 180,
            'CRPIX1': 90.5,
            'CRPIX2': 90.5,
            'CRVAL1': center.ra.deg,
            'CRVAL2': center.dec.deg,
            'CDELT1': -2 / np.pi,
            'CDELT2': 2 / np.pi,
            'CTYPE1': self._xcoord + '-SIN',
            'CTYPE2': self._ycoord + '-SIN',
            'RADESYS': self._radesys}
        if rotate is not None:
            header['LONPOLE'] = u.Quantity(rotate).to_value(u.deg)
        super().__init__(
            *args, frame_class=EllipticalFrame, header=header, **kwargs)


class Zoom(AutoScaledWCSAxes):

    def __init__(self, *args, center='0d 0d', radius='1 deg', rotate=None,
                 **kwargs):
        center = SkyCoord(
            center, representation_type=UnitSphericalRepresentation).icrs
        radius = u.Quantity(radius).to(u.deg).value
        header = {
            'NAXIS': 2,
            'NAXIS1': 512,
            'NAXIS2': 512,
            'CRPIX1': 256.5,
            'CRPIX2': 256.5,
            'CRVAL1': center.ra.deg,
            'CRVAL2': center.dec.deg,
            'CDELT1': -radius / 256,
            'CDELT2': radius / 256,
            'CTYPE1': self._xcoord + '-TAN',
            'CTYPE2': self._ycoord + '-TAN',
            'RADESYS': self._radesys}
        if rotate is not None:
            header['LONPOLE'] = u.Quantity(rotate).to_value(u.deg)
        super().__init__(*args, header=header, **kwargs)


class AllSkyAxes(AutoScaledWCSAxes):
    """Base class for a multi-purpose all-sky projection."""

    def __init__(self, *args, center=None, **kwargs):
        if center is None:
            center = f"{self._crval1}d 0d"
        center = SkyCoord(
            center, representation_type=UnitSphericalRepresentation).icrs
        header = {
            'NAXIS': 2,
            'NAXIS1': 360,
            'NAXIS2': 180,
            'CRPIX1': 180.5,
            'CRPIX2': 90.5,
            'CRVAL1': center.ra.deg,
            'CRVAL2': center.dec.deg,
            'CDELT1': -2 * np.sqrt(2) / np.pi,
            'CDELT2': 2 * np.sqrt(2) / np.pi,
            'CTYPE1': self._xcoord + '-' + self._wcsprj,
            'CTYPE2': self._ycoord + '-' + self._wcsprj,
            'RADESYS': self._radesys}
        super().__init__(
            *args, frame_class=EllipticalFrame, header=header, **kwargs)
        self.coords[0].set_ticks(spacing=45 * u.deg)
        self.coords[1].set_ticks(spacing=30 * u.deg)
        self.coords[0].set_ticklabel(exclude_overlapping=True)
        self.coords[1].set_ticklabel(exclude_overlapping=True)


class Aitoff(AllSkyAxes):
    _wcsprj = 'AIT'


class Mollweide(AllSkyAxes):
    _wcsprj = 'MOL'


moddict = globals()

#
# Create subclasses and register all projections:
# '{astro|geo|galactic} {hours|degrees} {aitoff|globe|mollweide|zoom}'
#
bases1 = (Astro, Geo, Galactic)
bases2 = (Hours, Degrees)
bases3 = (Aitoff, Globe, Mollweide, Zoom)
for bases in product(bases1, bases2, bases3):
    class_name = ''.join(cls.__name__ for cls in bases) + 'Axes'
    projection = ' '.join(cls.__name__.lower() for cls in bases)
    new_class = type(class_name, bases, {'name': projection})
    projection_registry.register(new_class)
    moddict[class_name] = new_class
    __all__.append(class_name)

#
# Create some synonyms:
# 'astro' will be short for 'astro hours',
# 'geo' will be short for 'geo degrees'
#

bases2 = (Hours, Degrees, Degrees)
for base1, base2 in zip(bases1, bases2):
    for base3 in (Aitoff, Globe, Mollweide, Zoom):
        bases = (base1, base2, base3)
        orig_class_name = ''.join(cls.__name__ for cls in bases) + 'Axes'
        orig_class = moddict[orig_class_name]
        class_name = ''.join(cls.__name__ for cls in (base1, base3)) + 'Axes'
        projection = ' '.join(cls.__name__.lower() for cls in (base1, base3))
        new_class = type(class_name, (orig_class,), {'name': projection})
        projection_registry.register(new_class)
        moddict[class_name] = new_class
        __all__.append(class_name)

del class_name, moddict, projection, projection_registry, new_class
__all__ = tuple(__all__)
