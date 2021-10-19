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
"""Read events from either HDF or LIGO-LW files."""
import builtins
import sqlite3

from ligo.lw.ligolw import Element
import h5py

from . import hdf, ligolw, sqlite

__all__ = ('MagicEventSource', 'open')


def _read_file_header(f, nbytes=16):
    """Read the first 16 bytes of a file

    This is presumed to include the characters that declare the
    file type.

    Parameters
    ----------
    f : file, str
        A file object or the path to a file

    Returns
    -------
    header : bytes
        A string (hopefully) describing the file type

    """
    try:
        pos = f.tell()
    except AttributeError:
        with builtins.open(f, "rb") as fobj:
            return fobj.read(nbytes)
    try:
        return f.read(nbytes)
    finally:
        f.seek(pos)


def MagicEventSource(f, *args, **kwargs):  # noqa: N802
    """Read events from LIGO-LW XML, LIGO-LW SQlite, or HDF5 files. The format
    is determined automatically using the :manpage:`file(1)` command, and then
    the file is opened using :obj:`.ligolw.open`, :obj:`.sqlite.open`, or
    :obj:`.hdf.open`, as appropriate.

    Returns
    -------
    `~ligo.skymap.io.events.EventSource`

    """
    if isinstance(f, h5py.File):
        opener = hdf.open
    elif isinstance(f, sqlite3.Connection):
        opener = sqlite.open
    elif isinstance(f, Element):
        opener = ligolw.open
    else:
        fileheader = _read_file_header(f)
        if fileheader.startswith(b'\x89HDF\r\n\x1a\n'):
            opener = hdf.open
        elif fileheader.startswith(b'SQLite format 3'):
            opener = sqlite.open
        elif fileheader.startswith((
                b'<?xml',  # XML
                b'\x1f\x8b\x08',  # GZIP
        )):
            opener = ligolw.open
        else:
            raise IOError('Unknown file format')
    return opener(f, *args, **kwargs)


open = MagicEventSource
