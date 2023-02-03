#
# Copyright (C) 2019-2020  Leo Singer
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
from heapq import heappop, heappush
from operator import length_hint

from tqdm.auto import tqdm

from .. import omp

__all__ = ('progress_map',)


class WrappedFunc:

    def __init__(self, func):
        self.func = func

    def __call__(self, i_args):
        i, args = i_args
        return i, self.func(*args)


def _get_total_estimate(*iterables):
    """Estimate total loop iterations for mapping over multiple iterables."""
    estimates = (length_hint(iterable, -1) for iterable in iterables)
    valid_estimates = (estimate for estimate in estimates if estimate != -1)
    return min(valid_estimates, default=None)


def _results_in_order(completed):
    """Put results back into order and yield them as quickly as they arrive."""
    heap = []
    current = 0
    for i_result in completed:
        i, result = i_result
        if i == current:
            yield result
            current += 1
            while heap and heap[0][0] == current:
                _, result = heappop(heap)
                yield result
                current += 1
        else:
            heappush(heap, i_result)
    assert not heap, 'The heap must be empty'


def _init_process():
    """Disable OpenMP when using multiprocessing."""
    omp.num_threads = 1


def progress_map(func, *iterables, jobs=1, **kwargs):
    """Map a function across iterables of arguments.

    This is comparable to :meth:`astropy.utils.console.ProgressBar.map`, except
    that it is implemented using :mod:`tqdm` and so provides more detailed and
    accurate progress information.
    """
    total = _get_total_estimate(*iterables)
    if jobs == 1:
        yield from tqdm(map(func, *iterables), total=total, **kwargs)
    else:
        with Pool(jobs, _init_process) as pool:
            yield from _results_in_order(
                tqdm(
                    pool.imap_unordered(
                        WrappedFunc(func),
                        enumerate(zip(*iterables))
                    ),
                    total=total, **kwargs
                )
            )
