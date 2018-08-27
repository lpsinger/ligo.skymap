# Copyright (C) 2017-2018  Leo Singer
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
from ligo.gracedb import rest

from .base import EventSource
from .ligolw import LigoLWEventSource

__all__ = ('GraceDBEventSource',)


class GraceDBEventSource(EventSource):
    """Read events from GraceDB.

    Parameters
    ----------
    graceids : list
        List of GraceDB ID strings.
    client : `ligo.gracedb.rest.GraceDb`, optional
        Client object
    """

    def __init__(self, graceids, client=None):
        if client is None:
            client = rest.GraceDb()
        self._client = client
        self._graceids = graceids

    def __iter__(self):
        return iter(self._graceids)

    def __getitem__(self, graceid):
        coinc_file = self._client.files(graceid, 'coinc.xml')
        psd_file = self._client.files(graceid, 'psd.xml.gz')
        event, = LigoLWEventSource(
            coinc_file, psd_file=psd_file, coinc_def=None).values()
        return event

    def __len__(self):
        return len(self._graceids)


open = GraceDBEventSource
