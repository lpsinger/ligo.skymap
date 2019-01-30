# Copied from https://github.com/numpy/numpy/pull/9211
#
# Copyright (c) 2005-2017, NumPy Developers.
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
#     * Redistributions of source code must retain the above copyright
#        notice, this list of conditions and the following disclaimer.
# 
#     * Redistributions in binary form must reproduce the above
#        copyright notice, this list of conditions and the following
#        disclaimer in the documentation and/or other materials provided
#        with the distribution.
# 
#     * Neither the name of the NumPy Developers nor the names of any
#        contributors may be used to endorse or promote products derived
#        from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import warnings

import numpy as np
from numpy.lib.function_base import _ureduce, array_function_dispatch

__all__ = ('percentile', 'quantile')


def _percentile_dispatcher(a, q, axis=None, weights=None, out=None,
                           overwrite_input=None, interpolation=None,
                           keepdims=None):
    return (a, q, out)


@array_function_dispatch(_percentile_dispatcher)
def percentile(a, q, axis=None, weights=None, out=None,
               overwrite_input=False, interpolation='linear', keepdims=False):
    q = np.true_divide(q, 100.0)  # handles the asarray for us too
    if not _quantile_is_valid(q):
        raise ValueError("Percentiles must be in the range [0, 100]")
    return _quantile_unchecked(
        a, q, axis, weights, out, overwrite_input, interpolation, keepdims)


def _quantile_dispatcher(a, q, axis=None, weights=None, out=None,
                         overwrite_input=None, interpolation=None,
                         keepdims=None):
    return (a, q, out)


@array_function_dispatch(_quantile_dispatcher)
def quantile(a, q, axis=None, weights=None, out=None,
             overwrite_input=False, interpolation='linear', keepdims=False):
    q = np.asanyarray(q)
    if not _quantile_is_valid(q):
        raise ValueError("Quantiles must be in the range [0, 1]")
    return _quantile_unchecked(
        a, q, axis, weights, out, overwrite_input, interpolation, keepdims)


def _quantile_unchecked(a, q, axis=None, weights=None, out=None,
                        overwrite_input=False, interpolation='linear',
                        keepdims=False):
    """Assumes that q is in [0, 1], and is an ndarray"""
    wgt = _validate_weights(a, q, axis, weights)

    r, k = _ureduce(a, func=_quantile_ureduce_func, q=q, axis=axis,
                    weights=wgt, out=out,
                    overwrite_input=overwrite_input,
                    interpolation=interpolation)
    if keepdims:
        return r.reshape(q.shape + k)
    else:
        return r


def _quantile_is_valid(q):
    # avoid expensive reductions, relevant for arrays with < O(1000) elements
    if q.ndim == 1 and q.size < 10:
        for i in range(q.size):
            if q[i] < 0.0 or q[i] > 1.0:
                return False
    else:
        # faster than any()
        if np.count_nonzero(q < 0.0) or np.count_nonzero(q > 1.0):
            return False
    return True


def _validate_weights(a, q, axis, weights):
    if weights is None:
        return None

    a = asanyarray(a)
    wgt = np.asanyarray(weights)

    if not np.issubdtype(wgt.dtype, np.number):
        # convert every element to np.float
        # raises ValueError if any element fails to convert
        wgt = wgt.astype(float)

    if np.isnan(wgt).any():
        raise ValueError("No weight can be NaN.")

    if (wgt < 0).any():
        raise ValueError("Negative weight not allowed.")

    if issubclass(a.dtype.type, np.bool_):
        raise TypeError("Boolean weights not supported")
    elif issubclass(a.dtype.type, np.integer):
        result_dtype = np.result_type(a.dtype, wgt.dtype, 'f8')
    else:
        result_dtype = np.result_type(a.dtype, wgt.dtype)

    # determine if the weights can be applied to the array
    if a.shape == wgt.shape or wgt.size == 1:
        broadcastable = True
    elif a.ndim == wgt.ndim:  # same ndim, but some dims can be size = 1
        broadcastable = all([a_dim == w_dim or w_dim == 1
                             for a_dim, w_dim
                             in zip(a.shape, wgt.shape)])
    else:
        broadcastable = False

    if not broadcastable:  # 1-D array needs preprocessing before broadcast
        if axis is None:
            raise TypeError(
                "Axis must be specified when shapes of a and weights "
                "differ and not broadcastable.")
        if wgt.ndim != 1:
            raise TypeError(
                "1D weights expected when shapes of a and weights differ "
                " and not broadcastable.")
        if wgt.shape[0] != a.shape[axis]:
            raise ValueError(
                "Length of weights not compatible with specified axis.")
        # set up wgt to broadcast along axis
        wgt = np.broadcast_to(wgt, (a.ndim-1)*(1,) + wgt.shape)
        wgt = wgt.swapaxes(-1, axis)
    else:  # same shape, or at least broadcastable
        if axis is None:
            # to avoid ap = a.flatten() in _quantile_ureduce_func,
            # we explicitly set axis
            axis = tuple(range(a.ndim))

    scl = wgt.sum(axis=axis, dtype=result_dtype)
    if np.any(scl == 0.0):
        raise ZeroDivisionError(
            "Weights sum to zero, can't be normalized")

    # Obtain a weights array of the same shape as reduced a
    wgt = np.broadcast_to(wgt, a.shape)
    wgt, _ = _ureduce(wgt, func=lambda x, **kwargs: x, axis=axis)

    return wgt


def _quantile_ureduce_func(a, q, axis=None, weights=None, out=None,
                           overwrite_input=False, interpolation='linear',
                           keepdims=False):
    a = asarray(a)

    if q.ndim == 0:
        # Do not allow 0-d arrays because following code fails for scalar
        zerod = True
        q = q[None]
    else:
        zerod = False

    if a.size == 1:  # all quantiles point to the same value
        return np.repeat(a, q.size)

    # prepare a for partioning
    if overwrite_input:
        if axis is None:
            ap = a.ravel()
        else:
            ap = a
    else:
        if axis is None:
            ap = a.flatten()
        else:
            ap = a.copy()

    if axis is None:
        axis = 0

    if weights is None:
        Nx = ap.shape[axis]
        indices = q * (Nx - 1)
    else:
        # axis will be moved to -1 and indices will be floats
        ap, indices, Nx = _indices_for_weighted_quantile(ap, q, axis, weights)

    if interpolation == 'lower':
        indices = np.floor(indices).astype(intp)
    elif interpolation == 'higher':
        indices = np.ceil(indices).astype(intp)
    elif interpolation == 'midpoint':
        indices = 0.5 * (np.floor(indices) + np.ceil(indices))
    elif interpolation == 'nearest':
        indices = np.around(indices).astype(intp)
    elif interpolation == 'linear':
        pass  # keep index as fraction and interpolate
    else:
        raise ValueError(
            "interpolation can only be 'linear', 'lower' 'higher', "
            "'midpoint', or 'nearest'")

    n = np.array(False, dtype=bool)  # check for nan's flag
    inexact = np.issubdtype(a.dtype, np.inexact)  # if array could have nan's

    if indices.dtype == intp:  # take the points along axis

        if inexact:
            indices = concatenate((indices, [-1]))  # to move nan's to end

        ap.partition(indices, axis=axis)
        # ensure axis with q-th is first
        ap = np.moveaxis(ap, axis, 0)
        axis = 0

        if inexact:
            indices = indices[:-1]
            n = np.isnan(ap[-1:, ...])

        if zerod:
            indices = indices[0]

        r = np.take(ap, indices, axis=axis, out=out)

    else:  # weight the points above and below the indices
        indices_below = np.floor(indices).astype(intp)
        indices_above = indices_below + 1
        indices_above[indices_above > Nx - 1] = Nx - 1

        weights_above = indices - indices_below
        weights_below = 1.0 - weights_above

        if weights is None:
            if inexact:
                indices_above = concatenate((indices_above, [-1]))

            weights_shape = [1, ] * ap.ndim
            weights_shape[axis] = len(indices)
            weights_below.shape = weights_shape
            weights_above.shape = weights_shape

            ap.partition(concatenate((indices_below, indices_above)),
                         axis=axis)

            # ensure axis with qth is first
            ap = np.moveaxis(ap, axis, 0)
            weights_below = np.moveaxis(weights_below, axis, 0)
            weights_above = np.moveaxis(weights_above, axis, 0)
            axis = 0

            if inexact:
                indices_above = indices_above[:-1]
                n = np.isnan(ap[-1:, ...])

            x1 = np.take(ap, indices_below, axis=axis) * weights_below
            x2 = np.take(ap, indices_above, axis=axis) * weights_above

        else:  # weights were specificed, that means axis was moved to -1

            x1 = np.take_along_axis(ap, indices_below, axis=-1) * weights_below
            x2 = np.take_along_axis(ap, indices_above, axis=-1) * weights_above

            if x1.ndim > 1:  # move the axis back
                x1 = np.swapaxes(x1, axis, -1)
                x2 = np.swapaxes(x2, axis, -1)

        # ensure axis with q-th is first
        x1 = np.moveaxis(x1, axis, 0)
        x2 = np.moveaxis(x2, axis, 0)

        if zerod:
            x1 = x1.squeeze(0)
            x2 = x2.squeeze(0)

        if out is not None:
            r = add(x1, x2, out=out)
        else:
            r = add(x1, x2)

    if np.any(n):
        warnings.warn("Invalid value encountered in quantile",
                      RuntimeWarning, stacklevel=3)
        if zerod:
            if ap.ndim == 1:
                if out is not None:
                    out[...] = a.dtype.type(np.nan)
                    r = out
                else:
                    r = a.dtype.type(np.nan)
            else:
                r[n.squeeze(axis)] = a.dtype.type(np.nan)
        else:
            if r.ndim == 1:
                r[:] = a.dtype.type(np.nan)
            else:
                r[..., n.squeeze(axis)] = a.dtype.type(np.nan)

    return r


def _indices_for_weighted_quantile(ap, q, axis, weights):
    # first move to axis -1 for np.vectorize
    ap = np.swapaxes(ap.astype('f8'), axis, -1)
    weights = np.swapaxes(weights.astype('f8'), axis, -1)

    # (1) values with weight=0 are made nan and later sorted to end of array
    if (weights == 0).any():
        ap[weights == 0] = np.nan

    def _sort_by_index(vector, vec_indices):
        return vector[vec_indices]
    # this func vectorizes sort along axis
    arraysort = np.vectorize(_sort_by_index, signature='(i),(i)->(i)')

    # (2) compute the sorted data array.  NaNs are sorted to the end of array
    ind_sorted = np.argsort(ap, axis=-1)  # sort values long axis
    ap_sorted = arraysort(ap, ind_sorted)

    # align the weights to the sorted data array
    ws_sorted = arraysort(weights, ind_sorted)

    # along axis -1,
    # (3) normalize weights by the smallest non-zero weight
    # (4) if more than 1 weights, subtract the first one by unit weight
    def _normalize_and_set_boundary(ws):  # ws is a 1-D vector
        w_vec = ws.copy()
        inds = w_vec > 0
        if (w_vec[inds] < 1).any():
            w_vec[inds] = w_vec[inds] / w_vec[inds].min()
        if len(w_vec) > 1:
            w_vec[0] -= 1
        return w_vec

    normalize_and_bound = np.vectorize(_normalize_and_set_boundary,
                                       signature='(i)->(i)')
    ws_sorted = normalize_and_bound(ws_sorted)

    nonzero_w_inds = ws_sorted > 0

    # compute the cumulative weights
    cum_w = ws_sorted.cumsum(axis=-1)
    cum_w_max = cum_w.max(axis=-1, keepdims=True)

    # (5) compute the upper bounds of probability bands
    # upper bound is just the normalized cumulative weight
    prob_band_upper_bounds = cum_w / cum_w_max

    # (6) finding lower bound requires computing the transition band width,
    # from the prior band's right edge to the current band's left edge,
    # which is (current band upper bound - prior band upper bound) / weight
    # first need the prior band upper bound
    prior_band_upper_bounds = np.roll(prob_band_upper_bounds, 1, axis=-1)
    prior_band_upper_bounds[..., 0] = 0  # nothing prior to left-most band

    # the weight denominator could contain 0's (due to input or nans)
    # so the division is only performed where weight > 0 (nonzero_w_inds)
    transition_band_width = ws_sorted  # to initialize
    # transition_band_width inherits the fact that
    # band width == 0 where weight == 0
    transition_band_width[nonzero_w_inds] =\
        ((prob_band_upper_bounds[nonzero_w_inds] -
          prior_band_upper_bounds[nonzero_w_inds]) /
         ws_sorted[nonzero_w_inds])  # ok to overwrite ws_sorted

    # the transition_band_width computed should align with the left band
    # it bounds to, so it's rolled to the left by 1 position
    transition_band_width = np.roll(transition_band_width, -1, axis=-1)

    # (7) the lower bound to the current probability band is equal to
    # prior band's upper bound + prior band's transition band width,
    # so we sum and then roll to the right by 1 position
    prob_band_lower_bounds =\
        np.roll(prob_band_upper_bounds + transition_band_width, 1, axis=-1)
    prob_band_lower_bounds[..., 0] = 0  # boundary condition on the left

    # combine and riffle shuffle
    quantile_bounds = np.concatenate([prob_band_upper_bounds,
                                      prob_band_lower_bounds],
                                     axis=-1)
    quantile_bounds[..., 0::2] = prob_band_lower_bounds
    quantile_bounds[..., 1::2] = prob_band_upper_bounds

    ap = np.concatenate([ap_sorted, ap_sorted], axis=-1)
    ap[..., 0::2] = ap_sorted
    ap[..., 1::2] = ap_sorted

    # (8) interpolate for indices where q sits along axis of ap
    Nx = ap.shape[-1]
    indices_hard = np.arange(Nx)
    vec_interp_func = np.vectorize(np.interp, signature='(n),(m),(m)->(n)')
    indices = vec_interp_func(q, quantile_bounds, indices_hard)

    return ap, indices, Nx
