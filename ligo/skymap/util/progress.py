#
# Copyright (C) 2019  Leo Singer
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
"""Tools for progress bars"""

try:
    from billiard import Pool
except ImportError:
    from multiprocessing import Pool
from operator import itemgetter

from tqdm.auto import tqdm

__all__ = ('progress_map',)


class WrappedFunc:

    def __init__(self, func):
        self.func = func

    def __call__(self, i_args):
        i, args = i_args
        return i, self.func(*args)


def progress_map(func, *iterables, jobs=1, **kwargs):
    r"""
    Map a function across iterables of arguments.

    This is comparable to :meth:`astropy.utils.console.ProgressBar.map`, except
    that it is implemented using :mod:`tqdm` and so provides more detailed and
    accurate progress information.
    """
    total = min(len(iterable) for iterable in iterables)
    if jobs == 1:
        return list(tqdm(map(func, *iterables), total=total, **kwargs))
    else:
        with Pool(jobs) as pool:
            return [
                item[1] for item in sorted(
                    tqdm(
                        pool.imap_unordered(
                            WrappedFunc(func),
                            enumerate(zip(*iterables))
                        ),
                        total=total, **kwargs
                    ),
                    key=itemgetter(0)
                )
            ]
