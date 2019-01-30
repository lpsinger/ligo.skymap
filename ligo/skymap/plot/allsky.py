#
# Copyright (C) 2012-2018  Leo Singer
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
Axes subclasses for all-sky maps. This adds several
`astropy.visualization.wcsaxes.WCSAxes` subclasses to the
Matplotlib projection registry. The projections are:

*  ``astro degrees aitoff``
*  ``astro degrees mollweide``
*  ``astro hours aitoff``
*  ``astro hours mollweide``
*  ``geo degrees aitoff``
*  ``geo hours aitoff``
*  ``geo degrees mollweide``
*  ``geo hours mollweide``
*  ``astro globe`` with option `center`
*  ``astro zoom`` with options `center` and `radius`

Example
-------

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
from astropy.coordinates import SkyCoord
from astropy.io.fits import Header
from astropy.time import Time
from astropy.visualization.wcsaxes import WCSAxes
from astropy.visualization.wcsaxes.formatter_locator import (
    AngleFormatterLocator)
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.wcs import WCS
from astropy import units as u
from matplotlib import rcParams
from matplotlib.offsetbox import AnchoredOffsetbox
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, PathPatch
from matplotlib.projections import projection_registry
import numpy as np
from reproject import reproject_from_healpix
from scipy.ndimage import gaussian_filter
import scipy.optimize
from .angle import reference_angle_deg

__all__ = (
    'AstroDegreesAitoffAllSkyAxes',
    'AstroDegreesMollweideAllSkyAxes',
    'AstroHoursAitoffAllSkyAxes',
    'AstroHoursMollweideAllSkyAxes',
    'AutoScaledWCSAxes',
    'GeoDegreesAitoffAllSkyAxes',
    'GeoDegreesMollweideAllSkyAxes',
    'GeoHoursAitoffAllSkyAxes',
    'GeoHoursMollweideAllSkyAxes',
    'GlobeAxes',
    'ScaleBar',
    'ZoomSkyAxes')


class WCSInsetPatch(PathPatch):
    """Subclass of `matplotlib.patches.PathPatch` for marking the outline of
    one `astropy.visualization.wcsaxes.WCSAxes` inside another."""

    def __init__(self, ax, *args, **kwargs):
        self._ax = ax
        super(WCSInsetPatch, self).__init__(
            None, *args, fill=False,
            edgecolor=ax.coords.frame.get_color(),
            linewidth=ax.coords.frame.get_linewidth(),
            **kwargs)

    def get_path(self):
        frame = self._ax.coords.frame
        return frame.patch.get_path().interpolated(50).transformed(
            frame.transform)


class WCSInsetConnectionPatch(ConnectionPatch):
    """Patch to connect an inset WCS axes inside another WCS axes."""

    _corners_map = {1: 3, 2: 1, 3: 0, 4: 2}

    def __init__(self, ax, ax_inset, loc, *args, **kwargs):
        try:
            loc = AnchoredOffsetbox.codes[loc]
        except KeyError:
            loc = int(loc)
        corners = ax_inset.viewLim.corners()
        transform = (ax_inset.coords.frame.transform +
                     ax.coords.frame.transform.inverted())
        xy_inset = corners[self._corners_map[loc]]
        xy = transform.transform_point(xy_inset)
        super(WCSInsetConnectionPatch, self).__init__(
            xy, xy_inset, 'data', 'data', ax, ax_inset, *args,
            color=ax_inset.coords.frame.get_color(),
            linewidth=ax_inset.coords.frame.get_linewidth(),
            **kwargs)


class AutoScaledWCSAxes(WCSAxes):
    """Axes base class. The pixel scale is adjusted to the DPI of the image,
    and there are a variety of convenience methods."""

    def __init__(self, header, *args, **kwargs):
        super(AutoScaledWCSAxes, self).__init__(*args, aspect=1, **kwargs)
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

    def connect_inset_axes(self, ax, loc, *args, **kwargs):
        """Convenience function to connect a corner of another WCSAxes to the
        matching point inside this one.

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
        img, mask = reproject_from_healpix(data, self.header, hdu_in=hdu_in,
                                           order=order, nested=nested,
                                           field=field)
        img = np.ma.array(img, mask=~mask.astype(bool))
        if smooth is not None:
            pixsize = np.mean(np.abs(self.wcs.wcs.cdelt)) * u.deg
            smooth = (smooth / pixsize).to(u.dimensionless_unscaled).value
            img = gaussian_filter(img, smooth)
        return img

    def contour_hpx(self, data, hdu_in=None, order='bilinear', nested=False,
                    field=0, smooth=None, **kwargs):
        """Add contour levels for a HEALPix data set.

        Parameters
        ----------
        data : `numpy.ndarray` or str or `~astropy.io.fits.TableHDU` or `~astropy.io.fits.BinTableHDU` or tuple
            The HEALPix data set. If this is a `numpy.ndarray`, then it is
            interpreted as the HEALPi array in the same coordinate system as
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
            interpreted as the HEALPi array in the same coordinate system as
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
            interpreted as the HEALPi array in the same coordinate system as
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


class GlobeAxes(AutoScaledWCSAxes):

    name = 'astro globe'

    def __init__(self, *args, **kwargs):
        kwargs = dict(kwargs)
        center = SkyCoord(kwargs.pop('center', '0d 0d')).icrs
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
            'CTYPE1': 'RA---SIN',
            'CTYPE2': 'DEC--SIN',
            'RADESYS': 'ICRS'}
        super(GlobeAxes, self).__init__(
            header, *args, frame_class=EllipticalFrame, **kwargs)


class ZoomSkyAxes(AutoScaledWCSAxes):

    name = 'astro zoom'

    def __init__(self, *args, **kwargs):
        kwargs = dict(kwargs)
        center = SkyCoord(kwargs.pop('center', '0d 0d')).icrs
        radius = u.Quantity(kwargs.pop('radius', '1 deg')).to(u.deg).value
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
            'CTYPE1': 'RA---TAN',
            'CTYPE2': 'DEC--TAN',
            'RADESYS': 'ICRS'}
        super(ZoomSkyAxes, self).__init__(header, *args, **kwargs)


class AllSkyAxes(AutoScaledWCSAxes):
    """Base class for a multi-purpose all-sky projection"""

    def __init__(self, *args, **kwargs):
        kwargs = dict(kwargs)
        obstime = kwargs.pop('obstime', None)
        header = {
            'NAXIS': 2,
            'NAXIS1': 360,
            'NAXIS2': 180,
            'CRPIX1': 180.5,
            'CRPIX2': 90.5,
            'CRVAL1': self._crval1,
            'CRVAL2': 0.0,
            'CDELT1': -2 * np.sqrt(2) / np.pi,
            'CDELT2': 2 * np.sqrt(2) / np.pi,
            'CTYPE1': self._xcoord + '-' + self._wcsprj,
            'CTYPE2': self._ycoord + '-' + self._wcsprj,
            'RADESYS': self._radesys}
        if obstime is not None:
            header['DATE-OBS'] = Time(obstime).utc.isot
        super(AllSkyAxes, self).__init__(
            header, *args, frame_class=EllipticalFrame, **kwargs)
        self.coords[0].set_ticks(spacing=45 * u.deg)
        self.coords[1].set_ticks(spacing=30 * u.deg)
        self.coords[0].set_ticklabel(exclude_overlapping=True)
        self.coords[1].set_ticklabel(exclude_overlapping=True)


class Astro(object):
    _crval1 = 180
    _xcoord = 'RA--'
    _ycoord = 'DEC-'
    _radesys = 'ICRS'


class GeoAngleFormatterLocator(AngleFormatterLocator):

    def formatter(self, values, spacing):
        return super(GeoAngleFormatterLocator, self).formatter(
            reference_angle_deg(values.to(u.deg).value) * u.deg, spacing)


class Geo(WCSAxes):
    _crval1 = 0
    _radesys = 'ITRS'
    _xcoord = 'TLON'
    _ycoord = 'TLAT'

    def __init__(self, *args, **kwargs):
        super(Geo, self).__init__(*args, **kwargs)
        self.invert_xaxis()
        fl = self.coords[0]._formatter_locator
        self.coords[0]._formatter_locator = GeoAngleFormatterLocator(
            values=fl.values,
            number=fl.number,
            spacing=fl.spacing,
            format=fl.format)


class Degrees(WCSAxes):
    """WCS axes with longitude axis in degrees"""

    def __init__(self, *args, **kwargs):
        super(Degrees, self).__init__(*args, **kwargs)
        self.coords[0].set_major_formatter('d')


class Hours(WCSAxes):
    """WCS axes with longitude axis in hour angle"""

    def __init__(self, *args, **kwargs):
        super(Hours, self).__init__(*args, **kwargs)
        self.coords[0].set_major_formatter('hh')


class Aitoff(object):
    _wcsprj = 'AIT'


class Mollweide(object):
    _wcsprj = 'MOL'


class AstroDegreesAitoffAllSkyAxes(Astro, Degrees, Aitoff, AllSkyAxes):
    name = 'astro degrees aitoff'


class AstroDegreesMollweideAllSkyAxes(Astro, Degrees, Mollweide, AllSkyAxes):
    name = 'astro degrees mollweide'


class AstroHoursAitoffAllSkyAxes(Astro, Hours, Aitoff, AllSkyAxes):
    name = 'astro hours aitoff'


class AstroHoursMollweideAllSkyAxes(Astro, Hours, Mollweide, AllSkyAxes):
    name = 'astro hours mollweide'


class GeoDegreesAitoffAllSkyAxes(Geo, Degrees, Aitoff, AllSkyAxes):
    name = 'geo degrees aitoff'


class GeoHoursAitoffAllSkyAxes(Geo, Hours, Aitoff, AllSkyAxes):
    name = 'geo hours aitoff'


class GeoDegreesMollweideAllSkyAxes(Geo, Degrees, Mollweide, AllSkyAxes):
    name = 'geo degrees mollweide'


class GeoHoursMollweideAllSkyAxes(Geo, Hours, Mollweide, AllSkyAxes):
    name = 'geo hours mollweide'


projection_registry.register(
    AstroDegreesAitoffAllSkyAxes,
    AstroDegreesMollweideAllSkyAxes,
    AstroHoursAitoffAllSkyAxes,
    AstroHoursMollweideAllSkyAxes,
    GeoDegreesAitoffAllSkyAxes,
    GeoHoursAitoffAllSkyAxes,
    GeoDegreesMollweideAllSkyAxes,
    GeoHoursMollweideAllSkyAxes,
    GlobeAxes,
    ZoomSkyAxes)


class ScaleBar(FancyArrowPatch):

    def _func(self, dx, x, y):
        p1, p2 = self._transAxesToWorld.transform([[x, y], [x + dx, y]])
        p1 = SkyCoord(*p1, unit=u.deg)
        p2 = SkyCoord(*p2, unit=u.deg)
        return np.square((p1.separation(p2) - self._length).value)

    def __init__(self, ax, xy, length, *args, **kwargs):
        x, y = xy
        self._ax = ax
        self._length = length
        self._transAxesToWorld = (
            (ax.transAxes - ax.transData) + ax.coords.frame.transform)
        dx = scipy.optimize.minimize_scalar(
            self._func, args=xy, bounds=[0, 1 - x], method='bounded').x
        super(ScaleBar, self).__init__(
            xy, (x + dx, y),
            *args,
            arrowstyle='-',
            capstyle='round',
            color='black',
            linewidth=rcParams['lines.linewidth'],
            shrinkA=0.0,
            shrinkB=0.0,
            transform=ax.transAxes,
            **kwargs)

    def label(self):
        (x0, y), (x1, _) = self._posA_posB
        s = u' {0.value:g}{0.unit:unicode}'.format(self._length)
        return self._ax.text(
            0.5 * (x0 + x1), y, s,
            ha='center', va='bottom', transform=self._ax.transAxes)
