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
"""Tools for reading and writing SQLite databases"""

import os
import sqlite3
import sys
_open = open


def _open_a(string):
    return sqlite3.connect(string)


def _open_r(string):
    return sqlite3.connect('file:{}?mode=ro'.format(string), uri=True)


def _open_w(string):
    with _open(string, 'wb') as f:
        pass
    return sqlite3.connect(string)


_openers = {'a': _open_a, 'r': _open_r, 'w': _open_w}


def open(string, mode):
    """Open an SQLite database with an `open`-style mode flag.

    Parameters
    ----------
    string : str
        Path of the SQLite database file
    mode : {'r', 'w', 'a'}
        Access mode: read only, clobber and overwrite, or modify in place.

    Returns
    -------
    `sqlite3.Connection`

    Raises
    ------
    ValueError
        If the filename is invalid (e.g. ``/dev/stdin``)
    KeyError
        If the requested mode is invalid
    OSError
        If the database could not be opened in the specified mode
    """
    if string in {'-', '/dev/stdin', '/dev/stdout'}:
        raise ValueError('Cannot open stdin/stdout as an SQLite database')
    try:
        opener = _openers[mode]
    except KeyError:
        raise ValueError('Invalid mode "{}". Must be one of "{}".'.format(
            mode, ''.join(_openers.keys())))
    try:
        return opener(string)
    except (OSError, sqlite3.Error) as e:
        raise OSError('Failed to open database {}: {}'.format(string, e))


def get_filename(connection):
    """Get the name of the file associated with an SQLite connection.

    Parameters
    ----------
    connection : `sqlite3.Connection`
        The database connection

    Returns
    -------
    str
        The name of the file that contains the SQLite database

    Raises
    ------
    RuntimeError
        If more than one database is attached to the connection

    Examples
    --------

    >>> from ..util.file import TemporaryDirectory
    >>> with TemporaryDirectory() as d:
    ...     with sqlite3.connect(os.path.join(d, 'test.sqlite')) as db:
    ...         print(get_filename(db))
    ...
    /.../test.sqlite
    """
    result = connection.execute('pragma database_list').fetchall()
    try:
        (_, _, filename), = result
    except ValueError:
        raise RuntimeError('Expected exactly one attached database')
    return filename
