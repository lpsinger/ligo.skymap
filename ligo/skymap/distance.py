#
# Copyright (C) 2017-2024  Leo Singer
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
Distance distribution functions from [1]_, [2]_, [3]_.

References
----------
.. [1] Singer, Chen, & Holz, 2016. "Going the Distance: Mapping Host Galaxies
   of LIGO and Virgo Sources in Three Dimensions Using Local Cosmography and
   Targeted Follow-up." ApJL, 829, L15.
   :doi:`10.3847/2041-8205/829/1/L15`

.. [2] Singer, Chen, & Holz, 2016. "Supplement: 'Going the Distance: Mapping
   Host Galaxies of LIGO and Virgo Sources in Three Dimensions Using Local
   Cosmography and Targeted Follow-up' (2016, ApJL, 829, L15)." ApJS, 226, 10.
   :doi:`10.3847/0067-0049/226/1/10`

.. [3] https://asd.gsfc.nasa.gov/Leo.Singer/going-the-distance

"""

from astropy.table import Table
import astropy_healpix as ah
import numpy as np
import healpy as hp
import scipy.special
from .core import (conditional_pdf, conditional_cdf, conditional_ppf,
                   moments_to_parameters, parameters_to_moments,
                   volume_render as _volume_render, marginal_pdf, marginal_cdf,
                   marginal_ppf)
from .util.numpy import add_newdoc_ufunc, require_contiguous_aligned
from .healpix_tree import HEALPIX_MACHINE_ORDER
from . import moc

__all__ = ('conditional_pdf', 'conditional_cdf', 'conditional_ppf',
           'moments_to_parameters', 'parameters_to_moments', 'volume_render',
           'marginal_pdf', 'marginal_cdf', 'marginal_ppf', 'ud_grade',
           'conditional_kde', 'cartesian_kde_to_moments', 'principal_axes',
           'parameters_to_moments')


add_newdoc_ufunc(conditional_pdf, """\
Conditional distance probability density function (ansatz).

Parameters
----------
r : `numpy.ndarray`
    Distance (Mpc)
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Returns
-------
pdf : `numpy.ndarray`
    Conditional probability density according to ansatz.

""")


add_newdoc_ufunc(conditional_cdf, """\
Cumulative conditional distribution of distance (ansatz).

Parameters
----------
r : `numpy.ndarray`
    Distance (Mpc)
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Returns
-------
pdf : `numpy.ndarray`
    Conditional probability density according to ansatz.

Examples
--------
Test against numerical integral of pdf.

>>> import scipy.integrate
>>> distmu = 10.0
>>> distsigma = 5.0
>>> _, _, distnorm = parameters_to_moments(distmu, distsigma)
>>> r = 8.0
>>> expected, _ = scipy.integrate.quad(
...     conditional_pdf, 0, r,
...     (distmu, distsigma, distnorm))
>>> result = conditional_cdf(
...     r, distmu, distsigma, distnorm)
>>> np.testing.assert_almost_equal(result, expected)

For negative distances, it returns 0.

>>> print(conditional_cdf(-1, distmu, distsigma, distnorm))
0.0

For infinite postive distance, it returns 1 (provided that distnorm normalizes
the distribution).


>>> expected = 1.0
>>> result = conditional_cdf(np.inf, distmu, distsigma, distnorm)
>>> np.testing.assert_almost_equal(result, expected)

""")


add_newdoc_ufunc(conditional_ppf, """\
Point percent function (inverse cdf) of distribution of distance (ansatz).

Parameters
----------
p : `numpy.ndarray`
    The cumulative distribution function
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Returns
-------
r : `numpy.ndarray`
    Distance at which the cdf is equal to `p`.

Examples
--------
Test against numerical estimate.

>>> import scipy.optimize
>>> distmu = 10.0
>>> distsigma = 5.0
>>> distnorm = 1.0
>>> p = 0.16  # "one-sigma" lower limit
>>> expected_r16 = scipy.optimize.brentq(
... lambda r: conditional_cdf(r, distmu, distsigma, distnorm) - p, 0.0, 100.0)
>>> r16 = conditional_ppf(p, distmu, distsigma, distnorm)
>>> np.testing.assert_almost_equal(r16, expected_r16)

""")


add_newdoc_ufunc(moments_to_parameters, """\
Convert ansatz moments to parameters.

This function is the inverse of `parameters_to_moments`.

Parameters
----------
distmean : `numpy.ndarray`
    Conditional mean of distance (Mpc)
diststd : `numpy.ndarray`
    Conditional standard deviation of distance (Mpc)

Returns
-------
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

""")


add_newdoc_ufunc(parameters_to_moments, """\
Convert ansatz parameters to moments.

This function is the inverse of `moments_to_parameters`.

Parameters
----------
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)

Returns
-------
distmean : `numpy.ndarray`
    Conditional mean of distance (Mpc)
diststd : `numpy.ndarray`
    Conditional standard deviation of distance (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Examples
--------
For mu=0, sigma=1, the ansatz is a chi distribution with 3 degrees of
freedom, and the moments have simple expressions.

>>> mean, std, norm = parameters_to_moments(0, 1)
>>> expected_mean = 2 * np.sqrt(2 / np.pi)
>>> expected_std = np.sqrt(3 - expected_mean**2)
>>> expected_norm = 2.0
>>> np.testing.assert_allclose(mean, expected_mean)
>>> np.testing.assert_allclose(std, expected_std)
>>> np.testing.assert_allclose(norm, expected_norm)

Check that the moments scale as expected when we vary sigma.

>>> sigma = np.logspace(-8, 8)
>>> mean, std, norm = parameters_to_moments(0, sigma)
>>> np.testing.assert_allclose(mean, expected_mean * sigma)
>>> np.testing.assert_allclose(std, expected_std * sigma)
>>> np.testing.assert_allclose(norm, expected_norm / sigma**2)

Check some more arbitrary values using numerical quadrature:

>>> import scipy.integrate
>>> sigma = 1.0
>>> for mu in np.linspace(-10, 10):
...     mean, std, norm = parameters_to_moments(mu, sigma)
...     moments = np.empty(3)
...     for k in range(3):
...         moments[k], _ = scipy.integrate.quad(
...             lambda r: r**k * conditional_pdf(r, mu, sigma, 1.0),
...             0, np.inf)
...     expected_norm = 1 / moments[0]
...     expected_mean, r2 = moments[1:] * expected_norm
...     expected_std = np.sqrt(r2 - np.square(expected_mean))
...     np.testing.assert_approx_equal(mean, expected_mean, 5)
...     np.testing.assert_approx_equal(std, expected_std, 5)
...     np.testing.assert_approx_equal(norm, expected_norm, 5)

""")


_volume_render = require_contiguous_aligned(_volume_render)


def volume_render(x, y, max_distance, axis0, axis1, R, skymap):
    """Perform volumetric rendering of a 3D multi-resolution sky map.

    Parameters
    ----------
    x : `numpy.ndarray`
        X-coordinate in rendered image
    y : `numpy.ndarray`
        Y-coordinate in rendered image
    max_distance : float
        Limit of integration from `-max_distance` to `+max_distance`
    axis0 : int
        Index of axis to assign to x-coordinate
    axis1 : int
        Index of axis to assign to y-coordinate
    R : `numpy.ndarray`
        Rotation matrix as provided by `principal_axes`
    skymap : `astropy.table.Table`
        Multi-resolution sky map

    Returns
    -------
    image : `numpy.ndarray`
        Rendered image

    Examples
    --------
    Test volume rendering of a normal unit sphere...
    First, set up the 3D sky map.

    >>> import astropy_healpix as ah
    >>> from astropy.table import Table
    >>> import numpy as np
    >>> level = 5
    >>> nside = 2**level
    >>> npix = ah.nside_to_npix(nside)
    >>> ipix = np.arange(npix)
    >>> skymap = Table({
    ...     'UNIQ': ah.level_ipix_to_uniq(level, ipix),
    ...     'PROBDENSITY': np.ones(npix) / (4 * np.pi),
    ...     'DISTMU': np.zeros(npix),
    ...     'DISTSIGMA': np.ones(npix),
    ...     'DISTNORM': np.ones(npix) * 2.0
    ... })

    The conditional distance distribution should be a chi distribution with
    3 degrees of freedom.

    >>> from scipy.stats import norm, chi
    >>> r = np.linspace(0, 10.0)
    >>> actual = conditional_pdf(
    ...     r, skymap['DISTMU'][0], skymap['DISTSIGMA'][0], skymap['DISTNORM'][0])
    >>> expected = chi(3).pdf(r)
    >>> np.testing.assert_almost_equal(actual, expected)

    Next, run the volume renderer.

    >>> dmax = 4.0
    >>> n = 64
    >>> s = np.logspace(-dmax, dmax, n)
    >>> x, y = np.meshgrid(s, s)
    >>> R = np.eye(3)
    >>> P = volume_render(x, y, dmax, 0, 1, R, skymap)

    Next, integrate analytically.

    >>> P_expected = norm.pdf(x) * norm.pdf(y) * (norm.cdf(dmax) - norm.cdf(-dmax))

    Compare the two.

    >>> np.testing.assert_almost_equal(P, P_expected, decimal=4)

    Last, check that we don't have a coordinate singularity at the origin.

    >>> x = np.concatenate(([0], np.logspace(1 - n, 0, n) * dmax))
    >>> y = 0.0
    >>> P = volume_render(x, y, dmax, 0, 1, R, skymap)
    >>> P_expected = norm.pdf(x) * norm.pdf(y) * (norm.cdf(dmax) - norm.cdf(-dmax))
    >>> np.testing.assert_allclose(P, P_expected, rtol=1e-4)

    """  # noqa: E501
    skymap = Table(skymap)
    uniq = skymap.columns.pop('UNIQ')
    nside = 1 << np.int64(moc.uniq2order(uniq.max()))
    order, nest = moc.uniq2nest(uniq)
    skymap['NEST'] = nest << np.int64(2 * (HEALPIX_MACHINE_ORDER - order))
    skymap.sort('NEST')
    return _volume_render(
        x, y, max_distance, axis0, axis1, R, nside,
        skymap['NEST'], skymap['PROBDENSITY'], skymap['DISTMU'],
        skymap['DISTSIGMA'], skymap['DISTNORM'])


add_newdoc_ufunc(marginal_pdf, """\
Calculate all-sky marginal pdf (ansatz).

Parameters
----------
r : `numpy.ndarray`
    Distance (Mpc)
prob : `numpy.ndarray`
    Marginal probability (pix^-2)
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Returns
-------
pdf : `numpy.ndarray`
    Marginal probability density according to ansatz.

Examples
--------

>>> npix = 12
>>> prob, distmu, distsigma, distnorm = np.random.uniform(size=(4, 12))
>>> r = np.linspace(0, 1)
>>> pdf_expected = np.dot(
...     conditional_pdf(r[:, np.newaxis], distmu, distsigma, distnorm), prob)
>>> pdf = marginal_pdf(r, prob, distmu, distsigma, distnorm)
>>> np.testing.assert_allclose(pdf, pdf_expected, rtol=1e-4)

""")
marginal_pdf = require_contiguous_aligned(marginal_pdf)


add_newdoc_ufunc(marginal_cdf, """\
Calculate all-sky marginal cdf (ansatz).

Parameters
----------
r : `numpy.ndarray`
    Distance (Mpc)
prob : `numpy.ndarray`
    Marginal probability (pix^-2)
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Returns
-------
cdf : `numpy.ndarray`
    Marginal cumulative probability according to ansatz.

Examples
--------

>>> npix = 12
>>> prob, distmu, distsigma = np.random.uniform(size=(3, 12))
>>> prob /= prob.sum()
>>> _, _, distnorm = parameters_to_moments(distmu, distsigma)
>>> r = np.linspace(0, 1)
>>> cdf_expected = np.dot(
...     conditional_cdf(r[:, np.newaxis], distmu, distsigma, distnorm), prob)
>>> cdf = marginal_cdf(r, prob, distmu, distsigma, distnorm)
>>> np.testing.assert_allclose(cdf, cdf_expected, rtol=1e-4)

For negative distances, it returns 0.

>>> print(marginal_cdf(-1, prob, distmu, distsigma, distnorm))
0.0

For infinite postive distance, it returns the sum of prob
(provided that distnorm normalizes the distribution).

>>> expected = 1.0
>>> result = marginal_cdf(np.inf, prob, distmu, distsigma, distnorm)
>>> np.testing.assert_almost_equal(result, expected)

""")
marginal_cdf = require_contiguous_aligned(marginal_cdf)


add_newdoc_ufunc(marginal_ppf, """\
Point percent function (inverse cdf) of marginal distribution of distance
(ansatz).

Parameters
----------
p : `numpy.ndarray`
    The cumulative distribution function
prob : `numpy.ndarray`
    Marginal probability (pix^-2)
distmu : `numpy.ndarray`
    Distance location parameter (Mpc)
distsigma : `numpy.ndarray`
    Distance scale parameter (Mpc)
distnorm : `numpy.ndarray`
    Distance normalization factor (Mpc^-2)

Returns
-------
r : `numpy.ndarray`
    Distance at which the cdf is equal to `p`.

Examples
--------

>>> from astropy.utils.misc import NumpyRNGContext
>>> npix = 12
>>> with NumpyRNGContext(0):
...     prob, distmu, distsigma, distnorm = np.random.uniform(size=(4, 12))
>>> r_expected = np.linspace(0.4, 0.7)
>>> cdf = marginal_cdf(r_expected, prob, distmu, distsigma, distnorm)
>>> r = marginal_ppf(cdf, prob, distmu, distsigma, distnorm)
>>> np.testing.assert_allclose(r, r_expected, rtol=1e-4)

""")
marginal_ppf = require_contiguous_aligned(marginal_ppf)


def ud_grade(prob, distmu, distsigma, *args, **kwargs):
    """
    Upsample or downsample a distance-resolved sky map.

    Parameters
    ----------
    prob : `numpy.ndarray`
        Marginal probability (pix^-2)
    distmu : `numpy.ndarray`
        Distance location parameter (Mpc)
    distsigma : `numpy.ndarray`
        Distance scale parameter (Mpc)
    *args, **kwargs :
        Additional arguments to `healpy.ud_grade` (e.g.,
        `nside`, `order_in`, `order_out`).

    Returns
    -------
    prob : `numpy.ndarray`
        Resampled marginal probability (pix^-2)
    distmu : `numpy.ndarray`
        Resampled distance location parameter (Mpc)
    distsigma : `numpy.ndarray`
        Resampled distance scale parameter (Mpc)
    distnorm : `numpy.ndarray`
        Resampled distance normalization factor (Mpc^-2)

    """
    bad = ~(np.isfinite(distmu) & np.isfinite(distsigma))
    distmean, diststd, _ = parameters_to_moments(distmu, distsigma)
    distmean[bad] = 0
    diststd[bad] = 0
    distmean = hp.ud_grade(prob * distmu, *args, power=-2, **kwargs)
    diststd = hp.ud_grade(prob * np.square(diststd), *args, power=-2, **kwargs)
    prob = hp.ud_grade(prob, *args, power=-2, **kwargs)
    distmean /= prob
    diststd = np.sqrt(diststd / prob)
    bad = ~hp.ud_grade(~bad, *args, power=-2, **kwargs)
    distmean[bad] = np.inf
    diststd[bad] = 1
    distmu, distsigma, distnorm = moments_to_parameters(distmean, diststd)
    return prob, distmu, distsigma, distnorm


def _conditional_kde(n, X, Cinv, W):
    Cinv_n = np.dot(Cinv, n)
    cinv = np.dot(n, Cinv_n)
    x = np.dot(Cinv_n, X) / cinv
    w = W * (0.5 / np.pi) * np.sqrt(np.linalg.det(Cinv) / cinv) * np.exp(
        0.5 * (np.square(x) * cinv - (np.dot(Cinv, X) * X).sum(0)))
    return x, cinv, w


def conditional_kde(n, datasets, inverse_covariances, weights):
    return [
        _conditional_kde(n, X, Cinv, W)
        for X, Cinv, W in zip(datasets, inverse_covariances, weights)]


def cartesian_kde_to_moments(n, datasets, inverse_covariances, weights):
    """
    Calculate the marginal probability, conditional mean, and conditional
    standard deviation of a mixture of three-dimensional kernel density
    estimators (KDEs), in a given direction specified by a unit vector.

    Parameters
    ----------
    n : `numpy.ndarray`
        A unit vector; an array of length 3.
    datasets : list of `numpy.ndarray`
        A list 2D Numpy arrays specifying the sample points of the KDEs.
        The first dimension of each array is 3.
    inverse_covariances: list of `numpy.ndarray`
        An array of 3x3 matrices specifying the inverses of the covariance
        matrices of the KDEs. The list has the same length as the datasets
        parameter.
    weights : list
        A list of floating-point weights.

    Returns
    -------
    prob : float
        The marginal probability in direction n, integrated over all distances.
    mean : float
        The conditional mean in direction n.
    std : float
        The conditional standard deviation in direction n.

    Examples
    --------
    >>> # Some imports
    >>> import scipy.stats
    >>> import scipy.integrate
    >>> # Construct random dataset for KDE
    >>> np.random.seed(0)
    >>> nclusters = 5
    >>> ndata = np.random.randint(0, 1000, nclusters)
    >>> covs = [np.random.uniform(0, 1, size=(3, 3)) for _ in range(nclusters)]
    >>> covs = [_ + _.T + 3 * np.eye(3) for _ in covs]
    >>> means = np.random.uniform(-1, 1, size=(nclusters, 3))
    >>> datasets = [np.random.multivariate_normal(m, c, n).T
    ...     for m, c, n in zip(means, covs, ndata)]
    >>> weights = ndata / float(np.sum(ndata))
    >>>
    >>> # Construct set of KDEs
    >>> kdes = [scipy.stats.gaussian_kde(_) for _ in datasets]
    >>>
    >>> # Random unit vector n
    >>> n = np.random.normal(size=3)
    >>> n /= np.sqrt(np.sum(np.square(n)))
    >>>
    >>> # Analytically evaluate conditional mean and std. dev. in direction n
    >>> datasets = [_.dataset for _ in kdes]
    >>> inverse_covariances = [_.inv_cov for _ in kdes]
    >>> result_prob, result_mean, result_std = cartesian_kde_to_moments(
    ...     n, datasets, inverse_covariances, weights)
    >>>
    >>> # Numerically integrate conditional distance moments
    >>> def rkbar(k):
    ...     def integrand(r):
    ...         return r ** k * np.sum([kde(r * n) * weight
    ...             for kde, weight in zip(kdes, weights)])
    ...     integral, err = scipy.integrate.quad(integrand, 0, np.inf)
    ...     return integral
    ...
    >>> r0bar = rkbar(2)
    >>> r1bar = rkbar(3)
    >>> r2bar = rkbar(4)
    >>>
    >>> # Extract conditional mean and std. dev.
    >>> r1bar /= r0bar
    >>> r2bar /= r0bar
    >>> expected_prob = r0bar
    >>> expected_mean = r1bar
    >>> expected_std = np.sqrt(r2bar - np.square(r1bar))
    >>>
    >>> # Check that the two methods give almost the same result
    >>> np.testing.assert_almost_equal(result_prob, expected_prob)
    >>> np.testing.assert_almost_equal(result_mean, expected_mean)
    >>> np.testing.assert_almost_equal(result_std, expected_std)
    >>>
    >>> # Check that KDE is normalized over unit sphere.
    >>> nside = 32
    >>> npix = ah.nside_to_npix(nside)
    >>> prob, _, _ = np.transpose([cartesian_kde_to_moments(
    ...     np.asarray(hp.pix2vec(nside, ipix)),
    ...     datasets, inverse_covariances, weights)
    ...     for ipix in range(npix)])
    >>> result_integral = prob.sum() * hp.nside2pixarea(nside)
    >>> np.testing.assert_almost_equal(result_integral, 1.0, decimal=4)

    """
    # Initialize moments of conditional KDE.
    r0bar = 0
    r1bar = 0
    r2bar = 0

    # Loop over KDEs.
    for X, Cinv, W in zip(datasets, inverse_covariances, weights):
        x, cinv, w = _conditional_kde(n, X, Cinv, W)

        # Accumulate moments of conditional KDE.
        c = 1 / cinv
        x2 = np.square(x)
        a = scipy.special.ndtr(x * np.sqrt(cinv))
        b = np.sqrt(0.5 / np.pi * c) * np.exp(-0.5 * cinv * x2)
        r0bar_ = (x2 + c) * a + x * b
        r1bar_ = x * (x2 + 3 * c) * a + (x2 + 2 * c) * b
        r2bar_ = (x2 * x2 + 6 * x2 * c + 3 * c * c) * a + x * (x2 + 5 * c) * b
        r0bar += np.mean(w * r0bar_)
        r1bar += np.mean(w * r1bar_)
        r2bar += np.mean(w * r2bar_)

    # Normalize moments.
    with np.errstate(invalid='ignore'):
        r1bar /= r0bar
        r2bar /= r0bar
    var = r2bar - np.square(r1bar)

    # Handle invalid values.
    if var >= 0:
        mean = r1bar
        std = np.sqrt(var)
    else:
        mean = np.inf
        std = 1.0
    prob = r0bar

    # Done!
    return prob, mean, std


def principal_axes(prob, distmu, distsigma, nest=False):
    npix = len(prob)
    nside = ah.npix_to_nside(npix)
    good = np.isfinite(prob) & np.isfinite(distmu) & np.isfinite(distsigma)
    ipix = np.flatnonzero(good)
    distmean, diststd, _ = parameters_to_moments(distmu[good], distsigma[good])
    mass = prob[good] * (np.square(diststd) + np.square(distmean))
    xyz = np.asarray(hp.pix2vec(nside, ipix, nest=nest))
    cov = np.dot(xyz * mass, xyz.T)
    L, V = np.linalg.eigh(cov)
    if np.linalg.det(V) < 0:
        V = -V
    return V


def principal_axes_moc(skymap):
    skymap = skymap[
        np.isfinite(skymap['PROBDENSITY']) &
        np.isfinite(skymap['DISTMU']) &
        np.isfinite(skymap['DISTNORM'])]

    prob = skymap['PROBDENSITY'] * moc.uniq2pixarea(skymap['UNIQ'])
    distmean, diststd, _ = parameters_to_moments(
        skymap['DISTMU'], skymap['DISTSIGMA'])
    mass = prob * (np.square(diststd) + np.square(distmean))

    order, ipix = moc.uniq2nest(skymap['UNIQ'])
    nside = 1 << np.int64(order)
    xyz = np.asarray(hp.pix2vec(nside, ipix, nest=True))

    cov = np.dot(xyz * mass, xyz.T)
    L, V = np.linalg.eigh(cov)
    if np.linalg.det(V) < 0:
        V = -V
    return V


def parameters_to_marginal_moments(prob, distmu, distsigma):
    """Calculate the marginal (integrated all-sky) mean and standard deviation
    of distance from the ansatz parameters.

    Parameters
    ----------
    prob : `numpy.ndarray`
        Marginal probability (pix^-2)
    distmu : `numpy.ndarray`
        Distance location parameter (Mpc)
    distsigma : `numpy.ndarray`
        Distance scale parameter (Mpc)

    Returns
    -------
    distmean : float
        Mean distance (Mpc)
    diststd : float
        Std. deviation of distance (Mpc)

    """
    good = np.isfinite(prob) & np.isfinite(distmu) & np.isfinite(distsigma)
    prob = prob[good]
    distmu = distmu[good]
    distsigma = distsigma[good]
    distmean, diststd, _ = parameters_to_moments(distmu, distsigma)
    rbar = (prob * distmean).sum()
    r2bar = (prob * (np.square(diststd) + np.square(distmean))).sum()
    return rbar, np.sqrt(r2bar - np.square(rbar))


del add_newdoc_ufunc, require_contiguous_aligned
