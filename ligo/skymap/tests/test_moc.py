from astropy import table
import healpy as hp
import numpy as np
import pytest

from .. import moc


def input_skymap(order1, d_order, fraction):
    """Construct a test multi-resolution sky map, with values that are
    proportional to the NESTED pixel index.

    To make the test more interesting by mixing together multiple resolutions,
    part of the sky map is refined to a higher order.

    Parameters
    ----------
    order1 : int
        The HEALPix resolution order.
    d_order : int
        The increase in orer for part of the sky map.
    fraction : float
        The fraction of the original pixels to refine.
    """
    order2 = order1 + d_order
    npix1 = hp.nside2npix(hp.order2nside(order1))
    npix2 = hp.nside2npix(hp.order2nside(order2))
    ipix1 = np.arange(npix1)
    ipix2 = np.arange(npix2)

    data1 = table.Table({
        'UNIQ': moc.nest2uniq(order1, ipix1.astype(np.uint64)),
        'VALUE': ipix1.astype(float),
        'VALUE2': np.pi * ipix1.astype(float)
    })

    data2 = table.Table({
        'UNIQ': moc.nest2uniq(order2, ipix2.astype(np.uint64)),
        'VALUE': np.repeat(ipix1, npix2 // npix1).astype(float),
        'VALUE2': np.pi * np.repeat(ipix1, npix2 // npix1).astype(float)
    })

    n = int(npix1 * (1 - fraction))
    return table.vstack((data1[:n], data2[n * npix2 // npix1:]))


@pytest.mark.parametrize('order_in', [6])
@pytest.mark.parametrize('d_order_in', range(3))
@pytest.mark.parametrize('fraction_in', [0, 0.25, 0.5, 1])
@pytest.mark.parametrize('order_out', range(6))
def test_rasterize_downsample(order_in, d_order_in, fraction_in, order_out):
    npix_in = hp.nside2npix(hp.order2nside(order_in))
    npix_out = hp.nside2npix(hp.order2nside(order_out))
    skymap_in = input_skymap(order_in, d_order_in, fraction_in)
    skymap_out = moc.rasterize(skymap_in, order_out)

    assert len(skymap_out) == npix_out
    reps = npix_in // npix_out
    expected = np.mean(np.arange(npix_in).reshape(-1, reps), axis=1)
    np.testing.assert_array_equal(skymap_out['VALUE'], expected)


@pytest.mark.parametrize('order_in', [2])
@pytest.mark.parametrize('d_order_in', range(3))
@pytest.mark.parametrize('fraction_in', [0, 0.25, 0.5, 1])
@pytest.mark.parametrize('order_out', range(3, 9))
def test_rasterize_upsample(order_in, d_order_in, fraction_in, order_out):
    npix_in = hp.nside2npix(hp.order2nside(order_in))
    npix_out = hp.nside2npix(hp.order2nside(order_out))
    skymap_in = input_skymap(order_in, d_order_in, fraction_in)
    skymap_out = moc.rasterize(skymap_in, order_out)

    assert len(skymap_out) == npix_out
    ipix = np.arange(npix_in)
    reps = npix_out // npix_in
    for i in range(reps):
        np.testing.assert_array_equal(skymap_out['VALUE'][i::reps], ipix)


@pytest.mark.parametrize('order', range(3))
def test_rasterize_default(order):
    npix = hp.nside2npix(hp.order2nside(order))
    skymap_in = input_skymap(order, 0, 0)
    skymap_out = moc.rasterize(skymap_in)
    assert len(skymap_out) == npix
