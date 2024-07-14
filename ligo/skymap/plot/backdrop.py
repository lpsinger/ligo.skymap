#
# Copyright (C) 2017-2020  Leo Singer
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
"""Backdrops for astronomical plots."""

from importlib.resources import files
import json
import warnings

from astropy.io import fits
from astropy.time import Time
from astropy.utils.data import download_file
from astropy.wcs import WCS
from matplotlib.image import imread
import numpy as np
from PIL.Image import DecompressionBombWarning
from reproject import reproject_interp

__all__ = ('bluemarble', 'blackmarble', 'coastlines', 'mellinger',
           'reproject_interp_rgb')


def big_imread(*args, **kwargs):
    """Wrapper for imread() that suppresses warnings when loading very large
    images (usually tiffs). Most of the all-sky images that we use in this
    module are large enough to trigger this warning:

        DecompressionBombWarning: Image size (91125000 pixels) exceeds limit of
        89478485 pixels, could be decompression bomb DOS attack.
    """
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DecompressionBombWarning)
        img = imread(*args, **kwargs)
    return img


def mellinger():
    """Get the Mellinger Milky Way panorama.

    Retrieve, cache, and return the Mellinger Milky Way panorama. See
    http://www.milkywaysky.com.

    Returns
    -------
    `astropy.io.fits.ImageHDU`
        A FITS WCS image in ICRS coordinates.

    Examples
    --------
    .. plot::
       :context: reset
       :include-source:
       :align: center

        from astropy.visualization import (ImageNormalize,
                                           AsymmetricPercentileInterval)
        from astropy.wcs import WCS
        from matplotlib import pyplot as plt
        from ligo.skymap.plot import mellinger
        from reproject import reproject_interp

        ax = plt.axes(projection='astro hours aitoff')
        backdrop = mellinger()
        backdrop_wcs = WCS(backdrop.header).dropaxis(-1)
        interval = AsymmetricPercentileInterval(45, 98)
        norm = ImageNormalize(backdrop.data, interval)
        backdrop_reprojected = np.asarray([
            reproject_interp((layer, backdrop_wcs), ax.header)[0]
            for layer in norm(backdrop.data)])
        backdrop_reprojected = np.rollaxis(backdrop_reprojected, 0, 3)
        ax.imshow(backdrop_reprojected)

    """
    url = 'https://web.archive.org/web/20160317214047/http://galaxy.phy.cmich.edu/~axel/mwpan2/mwpan2_RGB_3600.fits'  # noqa: E501
    hdu, = fits.open(url, cache=True)
    return hdu


def bluemarble(t, resolution='low'):
    """Get the "Blue Marble" image.

    Retrieve, cache, and return the NASA/NOAO/NPP "Blue Marble" image showing
    landforms and oceans.

    See https://visibleearth.nasa.gov/view.php?id=74117.

    Parameters
    ----------
    t : `astropy.time.Time`
        Time to embed in the WCS header.
    resolution : {'low', 'high'}
        Specify which version to use: the "low" resolution version (5400x2700
        pixels, the default) or the "high" resolution version (21600x10800
        pixels).

    Returns
    -------
    `astropy.io.fits.ImageHDU`
        A FITS WCS image in ICRS coordinates.

    Examples
    --------
    .. plot::
       :context: reset
       :include-source:
       :align: center

        from matplotlib import pyplot as plt
        from ligo.skymap.plot import bluemarble, reproject_interp_rgb

        obstime = '2017-08-17 12:41:04'
        ax = plt.axes(projection='geo degrees aitoff', obstime=obstime)
        ax.imshow(reproject_interp_rgb(bluemarble(obstime), ax.header))

    """
    variants = {
        'low': '5400x2700',
        'high': '21600x10800'
    }

    url = ('https://eoimages.gsfc.nasa.gov/images/imagerecords/74000/74117/'
           'world.200408.3x{}.png'.format(variants[resolution]))
    img = big_imread(download_file(url, cache=True))
    height, width, ndim = img.shape
    gmst_deg = Time(t).sidereal_time('mean', 'greenwich').deg
    header = fits.Header(dict(
        NAXIS=3,
        NAXIS1=ndim, NAXIS2=width, NAXIS3=height,
        CRPIX2=width / 2, CRPIX3=height / 2,
        CRVAL2=gmst_deg % 360, CRVAL3=0,
        CDELT2=360 / width,
        CDELT3=-180 / height,
        CTYPE2='RA---CAR',
        CTYPE3='DEC--CAR',
        RADESYSa='ICRS').items())
    return fits.ImageHDU(img[:, :, :], header)


def blackmarble(t, resolution='low'):
    """Get the "Black Marble" image.

    Get the NASA/NOAO/NPP image showing city lights, at the sidereal time given
    by t. See https://visibleearth.nasa.gov/view.php?id=79765.

    Parameters
    ----------
    t : `astropy.time.Time`
        Time to embed in the WCS header.
    resolution : {'low', 'mid', 'high'}
        Specify which version to use: the "low" resolution version (3600x1800
        pixels, the default), the "mid" resolution version (13500x6750 pixels),
        or the "high" resolution version (54000x27000 pixels).

    Returns
    -------
    `astropy.io.fits.ImageHDU`
        A FITS WCS image in ICRS coordinates.

    Examples
    --------
    .. plot::
       :context: reset
       :include-source:
       :align: center

        from matplotlib import pyplot as plt
        from ligo.skymap.plot import blackmarble, reproject_interp_rgb

        obstime = '2017-08-17 12:41:04'
        ax = plt.axes(projection='geo degrees aitoff', obstime=obstime)
        ax.imshow(reproject_interp_rgb(blackmarble(obstime), ax.header))

    """
    variants = {
        'low': '3600x1800',
        'high': '13500x6750',
        'mid': '54000x27000'
    }

    url = ('http://eoimages.gsfc.nasa.gov/images/imagerecords/79000/79765/'
           'dnb_land_ocean_ice.2012.{}_geo.tif'.format(variants[resolution]))
    img = big_imread(download_file(url, cache=True))
    height, width, ndim = img.shape
    gmst_deg = Time(t).sidereal_time('mean', 'greenwich').deg
    header = fits.Header(dict(
        NAXIS=3,
        NAXIS1=ndim, NAXIS2=width, NAXIS3=height,
        CRPIX2=width / 2, CRPIX3=height / 2,
        CRVAL2=gmst_deg % 360, CRVAL3=0,
        CDELT2=360 / width,
        CDELT3=-180 / height,
        CTYPE2='RA---CAR',
        CTYPE3='DEC--CAR',
        RADESYSa='ICRS').items())
    return fits.ImageHDU(img[:, :, :], header)


def reproject_interp_rgb(input_data, *args, **kwargs):
    data = input_data.data
    wcs = WCS(input_data.header).celestial
    return np.moveaxis(np.stack([
        reproject_interp((data[:, :, i], wcs),
                         *args, **kwargs)[0].astype(data.dtype)
        for i in range(3)]), 0, -1)


def coastlines():
    with files(__package__).joinpath(
            'ne_simplified_coastline.json').open() as f:
        geoms = json.load(f)['geometries']
    return [coord for geom in geoms for coord in zip(*geom['coordinates'])]
