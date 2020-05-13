#
# Copyright (C) 2020  Leo Singer
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
"""Performance measurement utilities."""
from resource import getrusage, RUSAGE_SELF
from time import perf_counter

import numpy as np


class StopwatchTimes:

    def __init__(self, real=0, user=0, sys=0):
        self.real = real
        self.user = user
        self.sys = sys

    def __iadd__(self, other):
        self.real += other.real
        self.user += other.user
        self.sys += other.sys
        return self

    def __isub__(self, other):
        self.real -= other.real
        self.user -= other.user
        self.sys -= other.sys
        return self

    def __add__(self, other):
        return StopwatchTimes(self.real + other.real,
                              self.user + other.user,
                              self.sys + other.sys)

    def __sub__(self, other):
        return StopwatchTimes(self.real - other.real,
                              self.user - other.user,
                              self.sys - other.sys)

    def __repr__(self):
        return f'{self.__class__.__name__}(real={self.real!r}, user={self.user!r}, sys={self.sys!r})'  # noqa: E501

    def __str__(self):
        real, user, sys = (np.format_float_positional(val, 3, unique=False)
                           for val in (self.real, self.user, self.sys))
        return f'real={real}s, user={user}s, sys={sys}s'

    @classmethod
    def now(cls):
        rusage = getrusage(RUSAGE_SELF)
        return cls(perf_counter(), rusage.ru_utime, rusage.ru_stime)

    def reset(self):
        self.real = self.user = self.sys = 0


class Stopwatch(StopwatchTimes):
    """A code profiling utility that mimics the interface of a stopwatch."""

    def __init__(self):
        super().__init__()
        self._base = None

    def start(self):
        self._base = StopwatchTimes.now()

    def stop(self):
        self += StopwatchTimes.now() - self._base
        self._base = None
        return self

    def lap(self):
        delta = StopwatchTimes.now() - self._base
        self += delta
        return delta
