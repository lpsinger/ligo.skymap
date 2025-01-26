import astropy_healpix as ah
from astropy import table
from astropy import units as u
import numpy as np
import pytest

from .. import moc


@pytest.mark.parametrize('order', [0, 1, 2])
@pytest.mark.parametrize('ipix', [-1, -2, -3])
def test_nest2uniq_invalid(order, ipix):
    """Test nest2uniq for invalid values."""
    assert moc.nest2uniq(np.int8(order), ipix) == -1


@pytest.mark.parametrize('uniq', [-1, 0, 2, 3])
def test_uniq2order_invalid(uniq):
    """Test uniq2order for invalid values."""
    assert moc.uniq2order(uniq) == -1


@pytest.mark.parametrize('uniq', [-1, 0, 2, 3])
def test_uniq2pixarea_invalid(uniq):
    """Test uniq2order for invalid values."""
    assert np.isnan(moc.uniq2pixarea(uniq))


@pytest.mark.parametrize('uniq', [-1, 0, 2, 3])
def test_uniq2nest_invalid(uniq):
    """Test uniq2order for invalid values."""
    order, nest = moc.uniq2nest(uniq)
    assert order == -1
    assert nest == -1


@pytest.mark.parametrize('uniq', [-1, 0, 2, 3])
def test_uniq2ang_invalid(uniq):
    """Test uniq2order for invalid values."""
    theta, phi = moc.uniq2ang(uniq)
    assert np.isnan(theta)
    assert np.isnan(phi)


def input_skymap(order1, d_order, fraction,
                 dtype=np.float64, extra_fields=True):
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
    dtype : numpy.dtype
        Numpy dtype for the data fields.
    extra_fields : boolean
        Add extra fields whose values are not tested, in order to check
        that we get correct results regardless of strut packing in memory.
    """
    order2 = order1 + d_order
    npix1 = ah.nside_to_npix(ah.level_to_nside(order1))
    npix2 = ah.nside_to_npix(ah.level_to_nside(order2))
    ipix1 = np.arange(npix1)
    ipix2 = np.arange(npix2)

    fields = {
        'UNIQ': moc.nest2uniq(order1, ipix1),
        'VALUE': ipix1.astype(dtype),
    }
    if extra_fields:
        fields['VALUE2'] = np.pi * ipix1.astype(dtype)
    data1 = table.Table(fields)

    fields = {
        'UNIQ': moc.nest2uniq(order2, ipix2),
        'VALUE': np.repeat(ipix1, npix2 // npix1).astype(dtype),
    }
    if extra_fields:
        fields['VALUE2'] = np.pi * np.repeat(
            ipix1, npix2 // npix1).astype(dtype)
    data2 = table.Table(fields)

    n = int(npix1 * (1 - fraction))
    return table.vstack((data1[:n], data2[n * npix2 // npix1:]))


def test_rasterize_oom():
    """Test that rasterize() will correctly raise a MemoryError if it runs out
    of memory.
    """
    # A pixel at the highest possible 64-bit HEALPix resolution.
    uniq = moc.nest2uniq(np.int8(29), 0)
    data = table.Table({'UNIQ': [uniq], 'VALUE': [0]})
    with pytest.raises(MemoryError):
        moc._rasterize(data)


@pytest.mark.parametrize('order_in', [6])
@pytest.mark.parametrize('d_order_in', range(3))
@pytest.mark.parametrize('fraction_in', [0, 0.25, 0.5, 1])
@pytest.mark.parametrize('order_out', range(6))
@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('extra_fields', [False, True])
def test_rasterize_downsample(order_in, d_order_in, fraction_in, order_out,
                              dtype, extra_fields):
    npix_in = ah.nside_to_npix(ah.level_to_nside(order_in))
    npix_out = ah.nside_to_npix(ah.level_to_nside(order_out))
    skymap_in = input_skymap(order_in, d_order_in, fraction_in, dtype,
                             extra_fields)
    skymap_out = moc.rasterize(skymap_in, order_out)

    assert len(skymap_out) == npix_out
    reps = npix_in // npix_out
    expected = np.mean(
        np.arange(npix_in).reshape(-1, reps), axis=1).astype(dtype)
    np.testing.assert_array_equal(skymap_out['VALUE'], expected)


@pytest.mark.parametrize('order_in', [2])
@pytest.mark.parametrize('d_order_in', range(3))
@pytest.mark.parametrize('fraction_in', [0, 0.25, 0.5, 1])
@pytest.mark.parametrize('order_out', range(3, 9))
@pytest.mark.parametrize('dtype', [np.float32, np.float64])
@pytest.mark.parametrize('extra_fields', [False, True])
def test_rasterize_upsample(order_in, d_order_in, fraction_in, order_out,
                            dtype, extra_fields):
    npix_in = ah.nside_to_npix(ah.level_to_nside(order_in))
    npix_out = ah.nside_to_npix(ah.level_to_nside(order_out))
    skymap_in = input_skymap(order_in, d_order_in, fraction_in, dtype,
                             extra_fields)
    skymap_out = moc.rasterize(skymap_in, order_out)

    assert len(skymap_out) == npix_out
    ipix = np.arange(npix_in, dtype=dtype)
    reps = npix_out // npix_in
    for i in range(reps):
        np.testing.assert_array_equal(skymap_out['VALUE'][i::reps], ipix)


@pytest.mark.parametrize('order', range(3))
def test_rasterize_default(order):
    npix = ah.nside_to_npix(ah.level_to_nside(order))
    skymap_in = input_skymap(order, 0, 0)
    skymap_out = moc.rasterize(skymap_in)
    assert len(skymap_out) == npix


@pytest.mark.parametrize('order', range(1, 4))
@pytest.mark.parametrize('rounds', range(3))
def test_bayestar_adaptive_grid(order, rounds, benchmark):
    def func(pts):
        ras, decs = pts.T
        return ras + decs

    top_nside = ah.level_to_nside(order)

    skymap_out = benchmark(
        moc.bayestar_adaptive_grid, func, top_nside=top_nside,
        rounds=rounds)

    assert len(skymap_out) == ah.nside_to_npix(top_nside) * (1 + .75 * rounds)

    level, ipix = ah.uniq_to_level_ipix(skymap_out['UNIQ'])
    nside = ah.level_to_nside(level)
    area = ah.nside_to_pixel_area(nside).to_value(u.steradian)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order='nested')
    expected = func(np.column_stack((ra.rad, dec.rad)))
    expected /= (expected * area).sum()
    np.testing.assert_array_equal(skymap_out["PROBDENSITY"], expected)

    assert area.sum() == pytest.approx(4 * np.pi)


@pytest.mark.parametrize('order', range(1, 4))
def test_bayestar_adaptive_grid_refinement(order):
    """Check that bayestar_adaptive_grid refines the correct pixels."""
    top_nside = ah.level_to_nside(order)
    top_hpx = ah.HEALPix(top_nside, order='nested')
    high_prob_ipix = np.random.choice(
        top_hpx.npix, top_hpx.npix // 4, replace=False)
    low_prob_ipix = np.setdiff1d(np.arange(top_hpx.npix), high_prob_ipix)

    high_prob_uniq = ah.level_ipix_to_uniq(
        order + 1, (4 * high_prob_ipix + np.arange(4)[:, np.newaxis]).ravel())
    low_prob_uniq = ah.level_ipix_to_uniq(order, (low_prob_ipix))

    skymap_expected = table.Table({
        'UNIQ': np.concatenate((high_prob_uniq, low_prob_uniq)),
        'PROBDENSITY': np.concatenate(
            (np.ones(high_prob_uniq.shape), np.zeros(low_prob_uniq.shape)))
    })
    skymap_expected['PROBDENSITY'] /= np.pi
    skymap_expected.sort('UNIQ')

    def func(pts):
        ra, dec = pts.T
        ipix = top_hpx.lonlat_to_healpix(ra * u.rad, dec * u.rad)
        return np.isin(ipix, high_prob_ipix).astype(float)

    skymap_out = moc.bayestar_adaptive_grid(
        func, top_nside=top_nside, rounds=1)
    skymap_out.sort('UNIQ')

    for key in ['UNIQ', 'PROBDENSITY']:
        np.testing.assert_array_almost_equal(
            skymap_out[key], skymap_expected[key])
