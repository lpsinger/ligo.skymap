from astropy.table import Table
import numpy as np
import healpy as hp
from scipy import stats

from ...distance import cartesian_kde_to_moments, moments_to_parameters
from ...moc import nest2uniq
from ..crossmatch import crossmatch


def test_crossmatch():
    # Construct a random mean vector and covariance matrix
    # for a multivariate Gaussian distribution in Cartesian coordinates.
    rotation = stats.special_ortho_group(3).rvs()
    eigenvalues = stats.gamma(2).rvs(3)
    cart_inv_cov = rotation @ np.diag(1 / eigenvalues) @ rotation.T
    cart_mean = stats.norm.rvs(size=3)

    # Set up HEALPix grid.
    level = 5
    nside = 2**level
    npix = hp.nside2npix(nside)
    ipix = np.arange(npix)
    cart_coords = np.column_stack(hp.pix2vec(nside, ipix, nest=True))

    # Create sky map using same method as KDE.
    cart_means = cart_mean[np.newaxis, ..., np.newaxis]
    cart_inv_covs = cart_inv_cov[np.newaxis, ...]
    weights = np.ones(1)
    probdensity, distmean, diststd = np.transpose([
        cartesian_kde_to_moments(n, cart_means, cart_inv_covs, weights)
        for n in cart_coords])

    # Create 3D, multi-order sky map.
    uniq = nest2uniq(level, ipix)
    distmu, distsigma, distnorm = moments_to_parameters(distmean, diststd)
    skymap = Table({'UNIQ': uniq,
                    'PROBDENSITY': probdensity,
                    'DISTMU': distmu,
                    'DISTSIGMA': distsigma,
                    'DISTNORM': distnorm})

    # Run function under test.
    contours = [0.1, 0.5, 0.9]
    result = crossmatch(skymap, contours=contours)

    # Test credible volumes against the known analytic expression.
    expected = (4 / 3 * np.pi * stats.chi(3).ppf(contours)**3 *
                np.sqrt(np.prod(eigenvalues)))
    np.testing.assert_allclose(result.contour_vols, expected, rtol=1e-3)
