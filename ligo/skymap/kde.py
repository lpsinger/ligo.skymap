#
# Copyright (C) 2012-2020  Will M. Farr <will.farr@ligo.org>
#                          Leo P. Singer <leo.singer@ligo.org>
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

import copyreg
from functools import partial

from astropy.coordinates import SkyCoord
from astropy.utils.misc import NumpyRNGContext
import healpy as hp
import logging
import numpy as np
from scipy.stats import gaussian_kde

from . import distance
from . import moc
from .coordinates import EigenFrame
from .util import progress_map

log = logging.getLogger()

__all__ = ('BoundedKDE', 'Clustered2DSkyKDE', 'Clustered3DSkyKDE',
           'Clustered2Plus1DSkyKDE')


class BoundedKDE(gaussian_kde):
    """Density estimation using a KDE on bounded domains.

    Bounds can be any combination of low or high (if no bound, set to
    ``float('inf')`` or ``float('-inf')``), and can be periodic or
    non-periodic.  Cannot handle topologies that have
    multi-dimensional periodicities; will only handle topologies that
    are direct products of (arbitrary numbers of) R, [0,1], and S1.

    Parameters
    ----------
    pts : :class:`numpy.ndarray`
        ``(Ndim, Npts)`` shaped array of points (as in :class:`gaussian_kde`).
    low
        Lower bounds; if ``None``, assume no lower bounds.
    high
        Upper bounds; if ``None``, assume no upper bounds.
    periodic
        Boolean array giving periodicity in each dimension; if
        ``None`` assume no dimension is periodic.
    bw_method : optional
        Bandwidth estimation method (see :class:`gaussian_kde`).

    """

    def __init__(self, pts, low=-np.inf, high=np.inf, periodic=False,
                 bw_method=None):

        super().__init__(pts, bw_method=bw_method)
        self._low = np.broadcast_to(
            low, self.d).astype(self.dataset.dtype)
        self._high = np.broadcast_to(
            high, self.d).astype(self.dataset.dtype)
        self._periodic = np.broadcast_to(
            periodic, self.d).astype(bool)

    def evaluate(self, pts):
        """Evaluate the KDE at the given points."""
        pts = np.atleast_2d(pts)
        d, m = pts.shape
        if d != self.d and d == 1 and m == self.d:
            pts = pts.T

        pts_orig = pts
        pts = np.copy(pts_orig)

        den = super().evaluate(pts)

        for i, (low, high, period) in enumerate(zip(self._low, self._high,
                                                    self._periodic)):
            if period:
                p = high - low

                pts[i, :] += p
                den += super().evaluate(pts)

                pts[i, :] -= 2.0 * p
                den += super().evaluate(pts)

                pts[i, :] = pts_orig[i, :]

            else:
                if not np.isneginf(low):
                    pts[i, :] = 2.0 * low - pts[i, :]
                    den += super().evaluate(pts)
                    pts[i, :] = pts_orig[i, :]

                if not np.isposinf(high):
                    pts[i, :] = 2.0 * high - pts[i, :]
                    den += super().evaluate(pts)
                    pts[i, :] = pts_orig[i, :]

        return den

    __call__ = evaluate

    def quantile(self, pt):
        """Quantile of ``pt``, evaluated by a greedy algorithm.

        Parameters
        ----------
        pt
            The point at which the quantile value is to be computed.

        Notes
        -----
        The quantile of ``pt`` is the fraction of points used to construct the
        KDE that have a lower KDE density than ``pt``.

        """
        return np.count_nonzero(self(self.dataset) < self(pt)) / self.n


def km_assign(mus, cov, pts):
    """Implement the assignment step in the k-means algorithm.

    Given a set of centers, ``mus``, a covariance matrix used to produce a
    metric on the space, ``cov``, and a set of points, ``pts`` (shape ``(npts,
    ndim)``), assigns each point to its nearest center, returning an array of
    indices of shape ``(npts,)`` giving the assignments.
    """
    k = mus.shape[0]
    n = pts.shape[0]

    dists = np.zeros((k, n))

    for i, mu in enumerate(mus):
        dx = pts - mu
        try:
            dists[i, :] = np.sum(dx * np.linalg.solve(cov, dx.T).T, axis=1)
        except np.linalg.LinAlgError:
            dists[i, :] = np.nan

    return np.nanargmin(dists, axis=0)


def km_centroids(pts, assign, k):
    """Implement the centroid-update step of the k-means algorithm.

    Given a set of points, ``pts``, of shape ``(npts, ndim)``, and an
    assignment of each point to a region, ``assign``, and the number of means,
    ``k``, returns an array of shape ``(k, ndim)`` giving the centroid of each
    region.
    """
    mus = np.zeros((k, pts.shape[1]))
    for i in range(k):
        sel = assign == i
        if np.sum(sel) > 0:
            mus[i, :] = np.mean(pts[sel, :], axis=0)
        else:
            mus[i, :] = pts[np.random.randint(pts.shape[0]), :]

    return mus


def k_means(pts, k):
    """Perform k-means clustering on the set of points.

    Parameters
    ----------
    pts
        Array of shape ``(npts, ndim)`` giving the points on which k-means is
        to operate.
    k
        Positive integer giving the number of regions.

    Returns
    -------
    centroids
        An ``(k, ndim)`` array giving the centroid of each region.
    assign
        An ``(npts,)`` array of integers between 0 (inclusive) and k
        (exclusive) indicating the assignment of each point to a region.

    """
    assert pts.shape[0] > k, 'must have more points than means'

    cov = np.cov(pts, rowvar=0)

    mus = np.random.permutation(pts)[:k, :]
    assign = km_assign(mus, cov, pts)
    while True:
        old_assign = assign

        mus = km_centroids(pts, assign, k)
        assign = km_assign(mus, cov, pts)

        if np.all(assign == old_assign):
            break

    return mus, assign


def _cluster(cls, pts, trials, i, seed, jobs):
    k = i // trials
    if k == 0:
        raise ValueError('Expected at least one cluster')
    try:
        if k == 1:
            assign = np.zeros(len(pts), dtype=np.intp)
        else:
            with NumpyRNGContext(i + seed):
                _, assign = k_means(pts, k)
        obj = cls(pts, assign=assign)
    except np.linalg.LinAlgError:
        return -np.inf,
    else:
        return obj.bic, k, obj.kdes


class ClusteredKDE:

    def __init__(self, pts, max_k=40, trials=5, assign=None, jobs=1):
        self.jobs = jobs
        if assign is None:
            log.info('clustering ...')
            # Make sure that each thread gets a different random number state.
            # We start by drawing a random integer s in the main thread, and
            # then the i'th subprocess will seed itself with the integer i + s.
            #
            # The seed must be an unsigned 32-bit integer, so if there are n
            # threads, then s must be drawn from the interval [0, 2**32 - n).
            seed = np.random.randint(0, 2**32 - max_k * trials)
            func = partial(_cluster, type(self), pts, trials, seed=seed,
                           jobs=jobs)
            self.bic, self.k, self.kdes = max(
                self._map(func, range(trials, (max_k + 1) * trials)),
                key=lambda items: items[:2])
        else:
            # Build KDEs for each cluster, skipping degenerate clusters
            self.kdes = []
            npts, ndim = pts.shape
            self.k = assign.max() + 1
            for i in range(self.k):
                sel = (assign == i)
                cluster_pts = pts[sel, :]
                # Equivalent to but faster than len(set(pts))
                nuniq = len(np.unique(cluster_pts, axis=0))
                # Skip if there are fewer unique points than dimensions
                if nuniq <= ndim:
                    continue
                try:
                    kde = gaussian_kde(cluster_pts.T)
                except (np.linalg.LinAlgError, ValueError):
                    # If there are fewer unique points than degrees of freedom,
                    # then the KDE will fail because the covariance matrix is
                    # singular. In that case, don't bother adding that cluster.
                    pass
                else:
                    self.kdes.append(kde)

            # Calculate BIC
            # The number of parameters is:
            #
            # * ndim for each centroid location
            #
            # * (ndim+1)*ndim/2 Kernel covariances for each cluster
            #
            # * one weighting factor for the cluster (minus one for the
            #   overall constraint that the weights must sum to one)
            nparams = (self.k * ndim +
                       0.5 * self.k * (ndim + 1) * ndim + self.k - 1)
            with np.errstate(divide='ignore'):
                self.bic = (
                    np.sum(np.log(self.eval_kdes(pts))) -
                    0.5 * nparams * np.log(npts))

    def eval_kdes(self, pts):
        pts = pts.T
        return sum(w * kde(pts) for w, kde in zip(self.weights, self.kdes))

    def __call__(self, pts):
        return self.eval_kdes(pts)

    @property
    def weights(self):
        """Get the cluster weights: the fraction of the points within each
        cluster.
        """
        w = np.asarray([kde.n for kde in self.kdes])
        return w / np.sum(w)

    def _map(self, func, items):
        return progress_map(func, items, jobs=self.jobs)


class SkyKDE(ClusteredKDE):

    @classmethod
    def transform(cls, pts):
        """Override in sub-classes to transform points."""
        raise NotImplementedError

    def __init__(self, pts, max_k=40, trials=5, assign=None, jobs=1):
        if assign is None:
            pts = self.transform(pts)
        super().__init__(
            pts, max_k=max_k, trials=trials, assign=assign, jobs=jobs)

    def __call__(self, pts):
        return super().__call__(self.transform(pts))

    def as_healpix(self, top_nside=16, rounds=8):
        return moc.bayestar_adaptive_grid(self, top_nside=top_nside,
                                          rounds=rounds)


# We have to put in some hooks to make instances of Clustered2DSkyKDE picklable
# because we dynamically create subclasses with different values of the 'frame'
# class variable. This gets even trickier because we need both the class and
# instance objects to be picklable.


class _Clustered2DSkyKDEMeta(type):  # noqa: N802
    """Metaclass to make dynamically created subclasses of Clustered2DSkyKDE
    picklable.
    """


def _Clustered2DSkyKDEMeta_pickle(cls):  # noqa: N802
    """Pickle dynamically created subclasses of Clustered2DSkyKDE."""
    return type, (cls.__name__, cls.__bases__, {'frame': cls.frame})


# Register function to pickle subclasses of Clustered2DSkyKDE.
copyreg.pickle(_Clustered2DSkyKDEMeta, _Clustered2DSkyKDEMeta_pickle)


def _Clustered2DSkyKDE_factory(name, frame):  # noqa: N802
    """Unpickle instances of dynamically created subclasses of
    Clustered2DSkyKDE.

    FIXME: In Python 3, we could make this a class method of Clustered2DSkyKDE.
    Unfortunately, Python 2 is picky about pickling bound class methods.
    """
    new_cls = type(name, (Clustered2DSkyKDE,), {'frame': frame})
    return super(Clustered2DSkyKDE, Clustered2DSkyKDE).__new__(new_cls)


class Clustered2DSkyKDE(SkyKDE, metaclass=_Clustered2DSkyKDEMeta):
    r"""Represents a kernel-density estimate of a sky-position PDF that has
    been decomposed into clusters, using a different kernel for each
    cluster.

    The estimated PDF is

    .. math::

      p\left( \vec{\theta} \right) = \sum_{i = 0}^{k-1} \frac{N_i}{N}
      \sum_{\vec{x} \in C_i} N\left[\vec{x}, \Sigma_i\right]\left( \vec{\theta}
      \right)

    where :math:`C_i` is the set of points belonging to cluster
    :math:`i`, :math:`N_i` is the number of points in this cluster,
    :math:`\Sigma_i` is the optimally-converging KDE covariance
    associated to cluster :math:`i`.

    The number of clusters, :math:`k` is chosen to maximize the `BIC
    <http://en.wikipedia.org/wiki/Bayesian_information_criterion>`_
    for the given set of points being drawn from the clustered KDE.
    The points are assigned to clusters using the k-means algorithm,
    with a decorrelated metric.  The overall clustering behavior is
    similar to the well-known `X-Means
    <http://www.cs.cmu.edu/~dpelleg/download/xmeans.pdf>`_ algorithm.
    """

    frame = None

    @classmethod
    def transform(cls, pts):
        pts = SkyCoord(*pts.T, unit='rad').transform_to(cls.frame).spherical
        return np.column_stack((pts.lon.rad, np.sin(pts.lat.rad)))

    def __new__(cls, pts, *args, **kwargs):
        frame = EigenFrame.for_coords(SkyCoord(*pts.T, unit='rad'))
        name = '{:s}_{:x}'.format(cls.__name__, id(frame))
        new_cls = type(name, (cls,), {'frame': frame})
        return super().__new__(new_cls)

    def __reduce__(self):
        """Pickle instances of dynamically created subclasses of
        Clustered2DSkyKDE.
        """
        factory_args = self.__class__.__name__, self.frame
        return _Clustered2DSkyKDE_factory, factory_args, self.__dict__

    def eval_kdes(self, pts):
        base = super().eval_kdes
        dphis = (0.0, 2 * np.pi, -2 * np.pi)
        phi, z = pts.T
        return sum(base(np.column_stack((phi + dphi, z))) for dphi in dphis)


class Clustered3DSkyKDE(SkyKDE):
    """Like :class:`Clustered2DSkyKDE`, but clusters in 3D
    space.  Can compute volumetric posterior density (per cubic Mpc),
    and also produce Healpix maps of the mean and standard deviation
    of the log-distance.
    """

    @classmethod
    def transform(cls, pts):
        return SkyCoord(*pts.T, unit='rad').cartesian.xyz.value.T

    def __call__(self, pts, distances=False):
        """Given an array of positions in RA, DEC, compute the marginal sky
        posterior and optinally the conditional distance parameters.
        """
        func = partial(distance.cartesian_kde_to_moments,
                       datasets=[_.dataset for _ in self.kdes],
                       inverse_covariances=[_.inv_cov for _ in self.kdes],
                       weights=self.weights)
        probdensity, mean, std = zip(*self._map(func, self.transform(pts)))
        if distances:
            mu, sigma, norm = distance.moments_to_parameters(mean, std)
            return probdensity, mu, sigma, norm
        else:
            return probdensity

    def posterior_spherical(self, pts):
        """Evaluate the posterior probability density in spherical polar
        coordinates, as a function of (ra, dec, distance).
        """
        return super().__call__(pts)

    def as_healpix(self, top_nside=16):
        """Return a HEALPix multi-order map of the posterior density and
        conditional distance distribution parameters.
        """
        m = super().as_healpix(top_nside=top_nside)
        order, ipix = moc.uniq2nest(m['UNIQ'])
        nside = 2 ** order.astype(int)
        theta, phi = hp.pix2ang(nside, ipix, nest=True)
        p = np.column_stack((phi, 0.5 * np.pi - theta))
        print('evaluating distance layers ...')
        _, m['DISTMU'], m['DISTSIGMA'], m['DISTNORM'] = self(p, distances=True)
        return m


class Clustered2Plus1DSkyKDE(Clustered3DSkyKDE):
    """A hybrid sky map estimator that uses a 2D clustered KDE for the marginal
    distribution as a function of (RA, Dec) and a 3D clustered KDE for the
    conditional distance distribution.
    """

    def __init__(self, pts, max_k=40, trials=5, assign=None, jobs=1):
        if assign is None:
            self.twod = Clustered2DSkyKDE(
                pts, max_k=max_k, trials=trials, assign=assign, jobs=jobs)
        super().__init__(
            pts, max_k=max_k, trials=trials, assign=assign, jobs=jobs)

    def __call__(self, pts, distances=False):
        probdensity = self.twod(pts)
        if distances:
            _, distmu, distsigma, distnorm = super().__call__(
                pts, distances=True)
            return probdensity, distmu, distsigma, distnorm
        else:
            return probdensity

    def posterior_spherical(self, pts):
        """Evaluate the posterior probability density in spherical polar
        coordinates, as a function of (ra, dec, distance).
        """
        return self(pts) * super().posterior_spherical(pts) / super().__call__(
            pts)
