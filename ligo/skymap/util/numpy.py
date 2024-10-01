#
# Copyright (C) 2018-2024  Leo Singer
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
import functools
import numpy as np

try:
    from numpy._core.umath import _add_newdoc_ufunc
except ImportError:
    # FIXME: Remove after we require Numpy >=2.0.0.
    from numpy.core.umath import _add_newdoc_ufunc

__all__ = ('add_newdoc_ufunc', 'require_contiguous_aligned')


def add_newdoc_ufunc(func, doc):  # pragma: no cover
    """Set the docstring for a Numpy ufunc.

    The function :func:`numpy._core.umath._add_newdoc_ufunc` can only change a
    ufunc's docstring if it is `NULL`. This workaround avoids an exception when
    the user tries to `reload()` this module.

    Notes
    -----
    :func:`numpy._core.umath._add_newdoc_ufunc` is not part of Numpy's public
    API, but according to upstream developers it is unlikely to go away any
    time soon.

    See https://github.com/numpy/numpy/issues/26233.
    """
    try:
        _add_newdoc_ufunc(func, doc)
    except ValueError as e:
        msg = 'Cannot change docstring of ufunc with non-NULL docstring'
        if e.args[0] == msg:
            pass


def require_contiguous_aligned(func):
    """Wrap a Numpy ufunc to guarantee that all of its inputs are
    C-contiguous arrays.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        n = func.nin
        args = [arg if i >= n or np.isscalar(arg)
                else np.require(arg, requirements={'CONTIGUOUS', 'ALIGNED'})
                for i, arg in enumerate(args)]
        return func(*args, **kwargs)
    return wrapper
