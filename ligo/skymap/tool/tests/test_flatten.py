from astropy import table
import healpy as hp
import numpy as np
import pytest

from . import run_entry_point
from ...distance import moments_to_parameters, parameters_to_marginal_moments
from ... import moc
from ...io.fits import read_sky_map, write_sky_map


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

    # Create a random sky map.
    area = hp.nside2pixarea(hp.order2nside(order1))
    probdensity = np.random.uniform(0, 1, npix1)
    prob = probdensity * area
    normalization = prob.sum()
    prob /= normalization
    probdensity /= normalization
    distmean = np.random.uniform(100, 110, npix1)
    diststd = np.random.uniform(0, 1 / np.sqrt(3) - 0.1, npix1) * distmean
    distmu, distsigma, distnorm = moments_to_parameters(distmean, diststd)
    assert np.all(np.isfinite(distmu))

    data1 = table.Table({
        'UNIQ': moc.nest2uniq(order1, ipix1.astype(np.uint64)),
        'PROBDENSITY': probdensity,
        'DISTMU': distmu,
        'DISTSIGMA': distsigma,
        'DISTNORM': distnorm
    })

    # Add some upsampled pixels.
    data2 = table.Table(np.repeat(data1, npix2 // npix1))
    data2['UNIQ'] = moc.nest2uniq(order2, ipix2.astype(np.uint64))
    n = int(npix1 * (1 - fraction))
    result = table.vstack((data1[:n], data2[n * npix2 // npix1:]))

    # Add marginal distance mean and standard deviation.
    rbar = (prob * distmean).sum()
    r2bar = (prob * (np.square(diststd) + np.square(distmean))).sum()
    result.meta['distmean'] = rbar
    result.meta['diststd'] = np.sqrt(r2bar - np.square(rbar))

    return result


@pytest.mark.parametrize('order_in', [2])
@pytest.mark.parametrize('d_order_in', range(3))
@pytest.mark.parametrize('fraction_in', [0, 0.25, 0.5, 1])
@pytest.mark.parametrize('nside_out', [None, 1, 2, 4, 8, 512])
def test_flatten(tmpdir, order_in, d_order_in, fraction_in, nside_out):
    """Test ligo-skyamp-flatten."""
    input_filename = str(tmpdir / 'bayestar.fits')
    output_filename = str(tmpdir / 'bayestar.fits.gz')

    skymap = input_skymap(order_in, d_order_in, fraction_in)
    write_sky_map(input_filename, skymap, moc=True)
    expected_distmean = skymap.meta['distmean']
    expected_diststd = skymap.meta['diststd']

    args = ['ligo-skymap-flatten', input_filename, output_filename]
    if nside_out is not None:
        args.extend(['--nside', str(nside_out)])
    run_entry_point(*args)

    (prob, distmu, distsigma, distnorm), _ = read_sky_map(
        output_filename, distances=True)
    distmean, diststd = parameters_to_marginal_moments(
        prob, distmu, distsigma)

    if nside_out is not None:
        assert len(prob) == hp.nside2npix(nside_out)

    assert prob.sum() == pytest.approx(1)
    assert distmean == pytest.approx(expected_distmean)
    assert diststd == pytest.approx(expected_diststd)
