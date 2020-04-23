from astropy.coordinates import CartesianRepresentation, SkyCoord
from astropy.table import Table
from astropy import units as u
import numpy as np
import healpy as hp
import pytest
from scipy import stats

from ...distance import cartesian_kde_to_moments, moments_to_parameters
from ...moc import nest2uniq
from ..crossmatch import crossmatch


@pytest.fixture
def cartesian_gaussian():
    """Generate a 3D Cartesian Gaussian with random mean and covariance."""
    rotation = stats.special_ortho_group(3).rvs()
    eigenvalues = stats.gamma(2).rvs(3)
    cov = rotation @ np.diag(eigenvalues) @ rotation.T
    mean = stats.norm.rvs(size=3)
    return stats.multivariate_normal(mean, cov)


def cartesian_gaussian_to_skymap(level, mean, cov):
    """Convert a 3D Cartesian Gaussian to a 3D sky map."""
    # Set up HEALPix grid.
    nside = 2**level
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    coords = np.column_stack(hp.pix2vec(nside, ipix, nest=True))

    # Create sky map using same method as KDE to convert from Cartesian
    # to spherical representation. This is just a special case where there
    # is only a single cluster and a single KDE sample.
    means = mean[np.newaxis, ..., np.newaxis]
    inv_covs = np.linalg.inv(cov)[np.newaxis, ...]
    weights = np.ones(1)
    probdensity, distmean, diststd = np.transpose([
        cartesian_kde_to_moments(n, means, inv_covs, weights) for n in coords])

    # Create 3D, multi-order sky map.
    uniq = nest2uniq(level, ipix)
    distmu, distsigma, distnorm = moments_to_parameters(distmean, diststd)
    return Table({'UNIQ': uniq,
                  'PROBDENSITY': probdensity,
                  'DISTMU': distmu,
                  'DISTSIGMA': distsigma,
                  'DISTNORM': distnorm})


@pytest.fixture(params=range(3))
def contours(request):
    """Generate random values for the ``contours`` argument."""
    return np.random.uniform(size=request.param)


@pytest.mark.parametrize('n_coordinates', [None, *range(3)])
def test_crossmatch_cartesian_gaussian_distribution(
        cartesian_gaussian, contours, n_coordinates):
    """Test on a Cartesian Gaussian distribution."""
    skymap = cartesian_gaussian_to_skymap(
        5, cartesian_gaussian.mean, cartesian_gaussian.cov)

    if n_coordinates is None:
        coordinates = None
    else:
        coordinates_xyz = cartesian_gaussian.rvs(size=n_coordinates)
        coordinates = SkyCoord(*coordinates_xyz.T * u.Mpc,
                               representation_type=CartesianRepresentation)

    result = crossmatch(skymap, contours=contours, coordinates=coordinates)

    # Test credible volumes against the known analytic expression.
    expected = (4 / 3 * np.pi * stats.chi(3).ppf(contours)**3 *
                np.sqrt(np.linalg.det(cartesian_gaussian.cov)))
    np.testing.assert_allclose(result.contour_vols, expected, rtol=1e-3)

    if coordinates is None:
        assert np.isnan(result.probdensity_vol)
    elif np.size(coordinates) == 0:
        assert np.size(result.probdensity_vol) == 0
    else:
        expected = cartesian_gaussian.pdf(coordinates_xyz)
        np.testing.assert_allclose(result.probdensity_vol, expected,
                                   rtol=1e-1)
