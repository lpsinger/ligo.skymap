#
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
import functools
import numpy as np

__all__ = ('add_newdoc_ufunc', 'require_contiguous')


def add_newdoc_ufunc(func, doc):  # pragma: no cover
    """The function `np.lib.add_newdoc_ufunc` can only change a ufunc's
    docstring if it is `NULL`. This workaround avoids an exception when the
    user tries to `reload()` this module."""
    try:
        np.lib.add_newdoc_ufunc(func, doc)
    except ValueError as e:
        msg = 'Cannot change docstring of ufunc with non-NULL docstring'
        if e.args[0] == msg:
            pass


def require_contiguous(func):
    """Wrap a Numpy ufunc to guarantee that all of its inputs are
    C-contiguous arrays."""
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        n = func.nin
        args = [arg if i >= n or np.isscalar(arg)
                else np.ascontiguousarray(arg)
                for i, arg in enumerate(args)]
        return func(*args, **kwargs)
    return wrapper
