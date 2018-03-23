#
# Copyright (C) 2013-2018  Leo Singer
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
"""
Collection of Python decorators.
"""

from functools import lru_cache, wraps
from astropy.utils.misc import NumpyRNGContext

__all__ = ('memoized', 'with_numpy_random_seed')

memoized = lru_cache(maxsize=None)


def with_numpy_random_seed(func, seed=0):
    """Decorate a function so that it is called with a pre-defined random seed.
    The random seed is restored when the function returns."""

    @wraps(func)
    def wrapped_func(*args, **kwargs):
        with NumpyRNGContext(seed):
            ret = func(*args, **kwargs)
        return ret

    return wrapped_func
