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

try:
    import numpy._core.numeric as _nx
except ImportError:
    # FIXME: Remove after we require Numpy >=2.0.0.
    import numpy.core.numeric as _nx


def _ureduce(a, func, **kwargs):
    a = np.asanyarray(a)
    axis = kwargs.get('axis', None)
    if axis is not None:
        keepdim = list(a.shape)
        nd = a.ndim
        axis = _nx.normalize_axis_tuple(axis, nd)

        for ax in axis:
            keepdim[ax] = 1

        if len(axis) == 1:
            kwargs['axis'] = axis[0]
        else:
            keep = set(range(nd)) - set(axis)
            nkeep = len(keep)
            # swap axis that should not be reduced to front
            for i, s in enumerate(sorted(keep)):
                a = a.swapaxes(i, s)
            # merge reduced axis
            a = a.reshape(a.shape[:nkeep] + (-1,))
            kwargs['axis'] = -1
        keepdim = tuple(keepdim)
    else:
        keepdim = (1,) * a.ndim

    r = func(a, **kwargs)
    return r, keepdim


def percentile(a, q, axis=None, weights=None, out=None,
               overwrite_input=False, interpolation='linear', keepdims=False):
    q = np.true_divide(q, 100.0)  # handles the asarray for us too
    if not _quantile_is_valid(q):
        raise ValueError("Percentiles must be in the range [0, 100]")
    return _quantile_unchecked(
        a, q, axis, weights, out, overwrite_input, interpolation, keepdims)


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
        wgt = None
    else:
        a = np.asanyarray(a)
        wgt = np.asanyarray(weights)

        if issubclass(a.dtype.type, (np.integer, np.bool_)):
            result_dtype = np.result_type(a.dtype, wgt.dtype, 'f8')
        else:
            result_dtype = np.result_type(a.dtype, wgt.dtype)

        broadcastable = False
        if a.ndim == wgt.ndim and a.shape != wgt.shape:
            broadcastable = all([a_dim == w_dim or w_dim == 1
                                 for a_dim, w_dim
                                 in zip(a.shape, wgt.shape)])

        if a.shape != wgt.shape and not broadcastable:
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
            if not np.issubdtype(wgt.dtype, np.number):
                raise ValueError("All weight entries must be numeric.")

            nan_ws = wgt[np.isnan(wgt)]
            if nan_ws.size > 0:
                raise ValueError("No weight can be NaN.")

            negative_ws = wgt[wgt < 0]
            if negative_ws.size > 0:
                raise ValueError("Negative weight not allowed.")

            # setup wgt to broadcast along axis
            wgt = np.broadcast_to(wgt, (a.ndim-1)*(1,) + wgt.shape)
            wgt = wgt.swapaxes(-1, axis)
        else:  # same shape, or at least broadcastable
            if axis is None:
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
    a = np.asarray(a)

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

    # move axis to -1 for np.vectorize() operations.  Will move back later.
    ap = np.swapaxes(ap, axis, -1)

    if weights is None:
        Nx = ap.shape[-1]
        indices = q * (Nx - 1)
    else:
        # need a copy of weights for later array assignment
        weights = np.swapaxes(weights.astype('f8'), axis, -1)
        weights[weights < 0.] = 0.  # negative weights are treated as 0
        # values with weight=0 are assigned minimum value and later moved left
        abs_min = np.amin(ap)
        ap[weights == 0.] = abs_min - 1.

        def _sort_by_index(vector, vec_indices):
            return vector[vec_indices]
        # this func vectorizes sort along axis
        arraysort = np.vectorize(_sort_by_index, signature='(i),(i)->(i)')

        ind_sorted = np.argsort(ap, axis=-1)  # sort values long axis
        ap_sorted = arraysort(ap, ind_sorted).astype('f8')

        n = np.isnan(ap_sorted[..., -1:])
        if n.ndim > 1:
            n = np.swapaxes(n, axis, -1)

        ws_sorted = arraysort(weights, ind_sorted).astype('f8')
        ws_sorted[np.isnan(ap_sorted)] = 0.  # neglect nans from calculation
        nonzero_w_inds = ws_sorted > 0.

        cum_w = ws_sorted.cumsum(axis=-1)
        cum_w_max = cum_w.max(axis=-1)

        # some manipulation to get lower/upper percentage bounds
        normalized_w_upper = (cum_w.T / cum_w_max.T).T
        prior_cum_w = np.roll(normalized_w_upper, 1, axis=-1)
        prior_cum_w[..., 0] = 0.

        w_slice = ws_sorted  # .copy()
        # in case any input weight is less than 1, we renormalize by min
        if True in (ws_sorted[nonzero_w_inds] < 1.0):
            ws_sorted[nonzero_w_inds] =\
                ws_sorted[nonzero_w_inds] / ws_sorted[nonzero_w_inds].min()

        w_slice[nonzero_w_inds] = ((normalized_w_upper[nonzero_w_inds] -
                                    prior_cum_w[nonzero_w_inds]) /
                                   ws_sorted[nonzero_w_inds])

        w_slice = np.roll(w_slice, -1, axis=-1)
        # now create the lower percentage bound
        normalized_w_lower = np.roll(normalized_w_upper + w_slice, 1, axis=-1)
        normalized_w_lower[..., 0] = 0.0

        # now we subtract by left-most w_slice value
        new_w_upper = (normalized_w_upper.T - w_slice[..., 0].T).T
        new_w_upper[new_w_upper < 0.0] = 0.0
        new_w_lower = (normalized_w_lower.T - w_slice[..., 0].T).T
        new_w_lower[new_w_lower < 0.0] = 0.0
        new_w_lower[..., 0] = 0.0

        # renormalize by right-most bound
        normalized_w_upper = (new_w_upper.T / new_w_upper[..., -1].T).T
        normalized_w_lower = (new_w_lower.T / new_w_upper[..., -1].T).T

        # combine and resort
        cum_w_bands = np.concatenate([normalized_w_upper, normalized_w_lower],
                                     axis=-1)
        inds_resort = np.argsort(cum_w_bands, axis=-1)
        cum_w_bands = arraysort(cum_w_bands, inds_resort)

        ap = np.concatenate([ap_sorted, ap_sorted], axis=-1)
        ap = arraysort(ap, inds_resort)

        # interpolate
        Nx = ap.shape[-1]
        indices_hard = np.arange(Nx)
        vec_interp_func = np.vectorize(np.interp, signature='(n),(m),(m)->(n)')
        indices = vec_interp_func(q, cum_w_bands, indices_hard)

    if interpolation == 'lower':
        indices = np.floor(indices).astype(np.intp)
    elif interpolation == 'higher':
        indices = np.ceil(indices).astype(np.intp)
    elif interpolation == 'midpoint':
        indices = 0.5 * (np.floor(indices) + np.ceil(indices))
    elif interpolation == 'nearest':
        indices = np.around(indices).astype(np.intp)
    elif interpolation == 'linear':
        pass  # keep index as fraction and interpolate
    else:
        raise ValueError(
            "interpolation can only be 'linear', 'lower' 'higher', "
            "'midpoint', or 'nearest'")

    inexact = np.issubdtype(a.dtype, np.inexact)

    if indices.dtype == np.intp:

        if weights is None:
            if inexact:
                indices = np.concatenate((indices, [-1]))  # to move nan's to end
            ap.partition(indices, axis=-1)
            n = np.isnan(ap[..., -1:])
            if inexact:
                indices = indices[:-1]
            if n.ndim > 1:
                n = np.swapaxes(n, axis, -1)

        r = np.take(ap, indices, axis=-1)

        if r.ndim > 1:
            r = np.swapaxes(r, axis, -1)  # move the axis back

        r = np.moveaxis(r, axis, 0)

        if zerod:
            r = r.squeeze(0)

        if out is not None:
            r = np.add(r, 0, out=out)

    else:  # weight the points above and below the indices
        indices_below = np.floor(indices).astype(np.intp)
        indices_above = indices_below + 1
        indices_above[indices_above > Nx - 1] = Nx - 1

        if weights is None:
            if inexact:
                indices_above = np.concatenate((indices_above, [-1]))
            ap.partition(np.concatenate((indices_below, indices_above)), axis=-1)
            n = np.isnan(ap[..., -1:])
            if inexact:
                indices_above = indices_above[:-1]
            if n.ndim > 1:
                n = np.swapaxes(n, axis, -1)

        weights_above = indices - indices_below
        weights_below = 1.0 - weights_above

        def _take1d(vec, inds, wts):
            return np.take(vec, inds) * wts

        vec_take = np.vectorize(_take1d, signature='(n),(m),(m)->(m)')

        x1 = vec_take(ap, indices_below, weights_below)
        x2 = vec_take(ap, indices_above, weights_above)

        if x1.ndim > 1:  # move the axis back
            x1 = np.swapaxes(x1, axis, -1)
            x2 = np.swapaxes(x2, axis, -1)

        x1 = np.moveaxis(x1, axis, 0)
        x2 = np.moveaxis(x2, axis, 0)

        if zerod:
            x1 = x1.squeeze(0)
            x2 = x2.squeeze(0)

        if out is not None:
            r = np.add(x1, x2, out=out)
        else:
            r = np.add(x1, x2)

    if np.any(n):
        warnings.warn("Invalid value encountered in percentile",
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
