from astropy.coordinates import EarthLocation, ITRS, SkyCoord
from astropy.time import Time
from astropy import units as u
from pytest import approx

from .. import DetectorFrame


def test_detector_frame():
    site_1 = EarthLocation.of_site('H1')
    site_2 = EarthLocation.of_site('L1')
    t = Time(1e9, format='gps')

    detector_frame = DetectorFrame(site_1=site_1, site_2=site_2, obstime=t)
    itrs_frame = ITRS(obstime=t)

    itrs_coord = SkyCoord(*(u.Quantity(site_1.geocentric) -
                            u.Quantity(site_2.geocentric)).value,
                          frame=itrs_frame)
    assert itrs_coord.transform_to(detector_frame).lat.deg == 90

    detector_coord = SkyCoord(lon=0 * u.deg, lat=90 * u.deg,
                              frame=detector_frame)
    converted = detector_coord.transform_to(itrs_frame)
    assert converted.spherical.lon.deg == approx(itrs_coord.spherical.lon.deg)
    assert converted.spherical.lat.deg == approx(itrs_coord.spherical.lat.deg)
