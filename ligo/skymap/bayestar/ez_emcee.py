# Copyright (C) 2018  Leo Singer
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
import numpy as np
from tqdm import tqdm

__all__ = ('ez_emcee',)


class LogProbFunctionWithPriorBounds:
    """Wrap a log probability function to enforce parameter bounds."""

    def __init__(self, log_prob_fn, lo, hi):
        self._inner = log_prob_fn
        self._lo = lo
        self._hi = hi

    def __call__(self, y, *args, **kwargs):
        x = transform_backward(y, self._lo, self._hi)
        return self._inner(x, *args, **kwargs) - np.log(jac(x, self._lo, self._hi))


def transform_forward(x, xmin, xmax):
    """Transform from problem coordinates in (lo, hi) sampling coordinates in
    (-inf, inf)."""
    return np.log((x - xmin) / (xmax - x))


def jac(x, xmin, xmax):
    return np.prod(1 / (x - xmin) + 1 / (xmax - x), axis=-1)


def transform_backward(y, xmin, xmax):
    """Transform from sampling coordinates in (-inf, inf) to problem
    coordinates in (lo, hi)."""
    return xmin / (np.exp(y) + 1) + xmax / (np.exp(-y) + 1)


def ez_emcee(log_prob_fn, lo, hi, nindep=200, nwalkers=None, nburnin=500,
             **kwargs):
    """Fire-and-forget MCMC sampling using `emcee.EnsembleSampler`, featuring
    automated convergence monitoring, progress tracking, and thinning.

    The parameters are bounded in the finite interval described by ``lo`` and
    ``hi``. (Currently, the bounds must be finite, but we will add infinite and
    half-infinite paramter bounds in the future.)

    If run in an interactive terminal, live progress is shown including the
    current sample number, the total required number of samples, time elapsed
    and estimated time remaining, acceptance fraction, and autocorrelation
    length.

    Sampling terminates when all chains have accumulated the requested number
    of independent samples.

    Parameters
    ----------
    log_prob_fn : callable
        The log probability function. It should take as its argument the
        parameter vector as an of length ``ndim``, or if it is vectorized, an
        2D array with ``ndim`` columns.
    lo : list, `numpy.ndarray`
        List of lower limits of parameters, of length ``ndim``.
    hi : list, `numpy.ndarray`
        List of upper limits of parameters, of length ``ndim``.
    nindep : int, optional
        Minimum number of independent samples.
    nwalkers : int, optional
        Number of walkers. The default is 4 times the number of dimensions.
    nburnin : int, optional
        Number of samples to discard during burn-in phase.

    Returns
    -------
    chain : `numpy.ndarray`
        The thinned and flattened posterior sample chain,
        with at least ``nindep`` * ``nwalkers`` rows
        and exactly ``ndim`` columns.

    Other parameters
    ----------------
    kwargs :
        Extra keyword arguments for `emcee.EnsembleSampler`.
        *Tip:* Consider setting the `pool` or `vectorized` keyword arguments in
        order to speed up likelihood evaluations.

    Notes
    -----
    The autocorrelation length, which has a complexity of :math:`O(N \log N)`
    in the number of samples, is recalulated at geometrically progressing
    intervals so that its amortized complexity per sample is constant. (In
    simpler terms, as the chains grow longer and the autocorrelation length
    takes longer to compute, we update it less frequently so that it is never
    more expensive than sampling the chain in the first place.)
    """

    # Optional dependencies
    from emcee import EnsembleSampler
    from emcee.autocorr import AutocorrError

    lo = np.asarray(lo)
    hi = np.asarray(hi)
    ndim = len(lo)
    log_prob_fn = LogProbFunctionWithPriorBounds(log_prob_fn, lo, hi)

    if nwalkers is None:
        nwalkers = 4 * ndim

    nsteps = 64

    with tqdm(total=nburnin + nindep * nsteps) as progress:

        sampler = EnsembleSampler(nwalkers, ndim, log_prob_fn, **kwargs)
        pos = np.random.uniform(lo, hi, (nwalkers, ndim))
        pos = transform_forward(pos, lo, hi)

        # Burn in
        progress.set_description('Burning in')
        for pos, log_prob, rstate in sampler.sample(
                pos, iterations=nburnin, store=False):
            progress.update()

        acl_bad = True
        acl = 0
        while acl_bad or sampler.iteration < nindep * acl:

            # Advance the chain
            progress.total = nburnin + max(sampler.iteration + nsteps,
                                           nindep * acl)
            progress.set_description('Sampling')
            for pos, log_prob, rstate in sampler.sample(
                    pos, log_prob, rstate, iterations=nsteps):
                progress.update()

            # Refresh convergence statistics
            progress.set_description('Checking convergence')
            try:
                acl = sampler.get_autocorr_time()
                acl_bad = False
            except AutocorrError as e:
                acl = e.tau
                acl_bad = True
            acl = np.max(acl)
            if np.isfinite(acl):
                acl = int(np.ceil(acl))
            else:
                acl_bad = True
            accept = np.mean(sampler.acceptance_fraction)
            progress.set_postfix(acl=acl, accept=accept)

            # The autocorrelation time calculation has complexity N log N in
            # the number of posterior samples. Only refresh the autocorrelation
            # length estimate on logarithmically spaced samples so that the
            # amortized complexity per sample is constant.
            nsteps *= 2

    chain = sampler.get_chain(thin=acl, flat=True)
    return transform_backward(chain, lo, hi)
