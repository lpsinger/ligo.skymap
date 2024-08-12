from astropy.io.fits import Header
from astropy.table import Table
from astropy_healpix import HEALPix
import numpy as np
from reproject import reproject_from_healpix

from ...moc import nest2uniq
from ..reproject_from_healpix_moc import reproject_from_healpix_moc


def test_reproject_from_healpix_moc():
    header = Header({
        'NAXIS': 2,
        'NAXIS1': 180,
        'NAXIS2': 180,
        'CRPIX1': 90.5,
        'CRPIX2': 90.5,
        'CRVAL1': 180,
        'CRVAL2': 0,
        'CDELT1': 1,
        'CDELT2': 1,
        'CTYPE1': 'RA---CAR',
        'CTYPE2': 'DEC--CAR',
        'RADESYS': 'ICRS',
    })
    hpx = HEALPix(nside=32)
    map_flat = np.arange(hpx.npix)

    map_moc = Table({
        'UNIQ': nest2uniq(np.int8(hpx.level), hpx.ring_to_nested(map_flat)),
        'IMG': map_flat,
    })

    expected = reproject_from_healpix(
        (map_flat, 'ICRS'), header, order='nearest-neighbor', nested=False)
    result = reproject_from_healpix_moc(
        (map_moc, 'ICRS'), header, 'nearest-neighbor')

    np.testing.assert_array_equal(result, expected)
