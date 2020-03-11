# Copyright (C) 2017-2020  Leo Singer
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
"""Read events from a GstLal-style SQLite output."""
import os
import sqlite3

from ligo.lw import dbtables

from ...util import sqlite
from .ligolw import LigoLWEventSource

__all__ = ('SQLiteEventSource',)


class SQLiteEventSource(LigoLWEventSource):
    """Read events from LIGO-LW SQLite files.

    Parameters
    ----------
    f : str, file-like object, or `sqlite3.Connection` instance
        The SQLite database.

    Returns
    -------
    `~ligo.skymap.io.events.EventSource`

    """

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
        super().__init__(dbtables.get_xml(db), *args, **kwargs)
        self._fallbackpath = os.path.dirname(filename) if filename else None


open = SQLiteEventSource
