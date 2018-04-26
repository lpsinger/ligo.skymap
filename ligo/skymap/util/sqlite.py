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

import sqlite3
_open = open


def _open_a(string):
    return sqlite3.connect(string)


def _open_r(string):
    return sqlite3.connect('file:{}?mode=ro'.format(string), uri=True)


def _open_w(string):
    with _open(string, 'wb'):
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
    connection : `sqlite3.Connection`

    Raises
    ------
    ValueError
        If the filename is invalid (e.g. ``/dev/stdin``), or if the requested
        mode is invalid
    OSError
        If the database could not be opened in the specified mode

    Examples
    --------

    >>> import tempfile
    >>> import os
    >>> with tempfile.TemporaryDirectory() as d:
    ...     open(os.path.join(d, 'test.sqlite'), 'w')
    ...
    <sqlite3.Connection object at 0x...>

    >>> with tempfile.TemporaryDirectory() as d:
    ...     open(os.path.join(d, 'test.sqlite'), 'r')
    ...
    Traceback (most recent call last):
      ...
    OSError: Failed to open database ...

    >>> open('/dev/stdin', 'r')
    Traceback (most recent call last):
      ...
    ValueError: Cannot open stdin/stdout as an SQLite database

    >>> open('test.sqlite', 'x')
    Traceback (most recent call last):
      ...
    ValueError: Invalid mode "x". Must be one of "arw".
    """
    if string in {'-', '/dev/stdin', '/dev/stdout'}:
        raise ValueError('Cannot open stdin/stdout as an SQLite database')
    try:
        opener = _openers[mode]
    except KeyError:
        raise ValueError('Invalid mode "{}". Must be one of "{}".'.format(
            mode, ''.join(sorted(_openers.keys()))))
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

    >>> import tempfile
    >>> import os
    >>> with tempfile.TemporaryDirectory() as d:
    ...     with sqlite3.connect(os.path.join(d, 'test.sqlite')) as db:
    ...         print(get_filename(db))
    ...
    /.../test.sqlite

    >>> with tempfile.TemporaryDirectory() as d:
    ...     with sqlite3.connect(os.path.join(d, 'test1.sqlite')) as db1, \\
    ...          sqlite3.connect(os.path.join(d, 'test2.sqlite')) as db2:
    ...         filename = get_filename(db1)
    ...         db2.execute('ATTACH DATABASE "{}" AS db2'.format(filename))
    ...         print(get_filename(db2))
    ...
    Traceback (most recent call last):
      ...
    RuntimeError: Expected exactly one attached database
    """
    result = connection.execute('pragma database_list').fetchall()
    try:
        (_, _, filename), = result
    except ValueError:
        raise RuntimeError('Expected exactly one attached database')
    return filename
