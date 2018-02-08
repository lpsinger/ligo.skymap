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
Read events from a list of GraceDB IDs.
"""
from .base import *
from .ligolw import *

__all__ = ('GraceDBEventSource',)


class GraceDBEventSource(EventSource):

    def __init__(self, graceids, client=None):
        if client is None:
            from ligo.gracedb.rest import GraceDb
            client = GraceDb()
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
        try:
            return len(self._graceids)
        except TypeError:
            raise NotImplementedError


open = GraceDBEventSource
