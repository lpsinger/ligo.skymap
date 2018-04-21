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
"""
Read events from GstLal-style SQLite output.
"""
import os
import sqlite3

from glue.ligolw import dbtables

from ...util import sqlite
from .ligolw import LigoLWEventSource

__all__ = ('SQLiteEventSource',)


class SQLiteEventSource(LigoLWEventSource):

    def __init__(self, f, *args, **kwargs):
        if isinstance(f, sqlite3.Connection):
            db = f
            filename = sqlite.get_filename(f)
        else:
            if hasattr(f, 'read'):
                filename = f.name
                f.close()
            else:
                filename = f
            db = sqlite.open(filename, 'r')
        super(SQLiteEventSource, self).__init__(dbtables.get_xml(db),
                                                *args, **kwargs)
        self._fallbackpath = os.path.dirname(filename) if filename else None


open = SQLiteEventSource
