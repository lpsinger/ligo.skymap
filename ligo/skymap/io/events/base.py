# Copyright (C) 2017  Leo Singer
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
Base classes for reading events from search pipelines.
"""

from abc import ABCMeta, abstractmethod
from collections.abc import Mapping

__all__ = ('EventSource', 'Event', 'SingleEvent')


def _fmt(obj, keys):
    kvs = ', '.join('{}={!r}'.format(key, getattr(obj, key)) for key in keys)
    return '<{}({})>'.format(obj.__class__.__name__, kvs)


class EventSource(Mapping):
    """Abstraction of a source of coincident events.

    This is a mapping from event IDs (which may be any hashable type, but are
    generally integers or strings) to instances of `Event`.
    """

    def __str__(self):
        try:
            length = len(self)
        except (NotImplementedError, TypeError):
            contents = '...'
        else:
            contents = '...{} items...'.format(length)
        return '<{}({{{}}})>'.format(self.__class__.__name__, contents)

    def __repr__(self):
        try:
            len(self)
        except NotImplementedError:
            contents = '...'
        else:
            contents = ', '.join('{}: {!r}'.format(key, value)
                                 for key, value in self.items())
        return '{}({{{}}})'.format(self.__class__.__name__, contents)


class Event(metaclass=ABCMeta):
    """Abstraction of a coincident trigger.

    Attributes
    ----------
    singles : list, tuple
        Sequence of `SingleEvent`
    template_args : dict
        Dictionary of template parameters
    """

    @property
    @abstractmethod
    def singles(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def template_args(self):
        raise NotImplementedError

    __str_keys = ('singles',)

    def __str__(self):
        return _fmt(self, self.__str_keys)

    __repr__ = __str__


class SingleEvent(metaclass=ABCMeta):
    """Abstraction of a single-detector trigger.

    Attributes
    ----------
    detector : str
        Instrument name (e.g. 'H1')
    snr : float
        Signal to noise ratio
    phase : float
        Phase on arrival
    time : float
        GPS time on arrival
    zerolag_time : float
        GPS time on arrival in zero-lag data, without time slides applied
    psd : `REAL8FrequencySeries`
        Power spectral density
    snr_series : `COMPLEX8TimeSeries`
        SNR time series
    """

    @property
    @abstractmethod
    def detector(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def snr(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def phase(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def time(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def zerolag_time(self):
        raise NotImplementedError

    @property
    @abstractmethod
    def psd(self):
        raise NotImplementedError

    @property
    def snr_series(self):
        return None

    __str_keys = ('detector', 'snr', 'phase', 'time')

    def __str__(self):
        keys = self.__str_keys
        if self.time != self.zerolag_time:
            keys += ('zerolag_time',)
        return _fmt(self, keys)

    __repr__ = __str__
