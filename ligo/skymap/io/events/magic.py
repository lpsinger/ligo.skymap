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
Read events from either HDF or LIGO-LW files.
"""
import os
import sqlite3
from subprocess import check_output

from glue.ligolw.ligolw import Element
import h5py

from . import hdf, ligolw, sqlite

__all__ = ('MagicEventSource', 'open')


def _get_file_type(f):
    """Determine the file type by calling the POSIX ``file`` utility.

    Parameters
    ----------
    f : file, str
        A file object or the path to a file

    Returns
    -------
    filetype : bytes
        A string describing the file type
    """
    try:
        f.read
    except AttributeError:
        filetype = check_output(
            ['file', f], env=dict(os.environ, POSIXLY_CORRECT='1'))
    else:
        filetype = check_output(
            ['file', '-'], env=dict(os.environ, POSIXLY_CORRECT='1'), stdin=f)
        f.seek(0)
    _, _, filetype = filetype.partition(b': ')
    return filetype.strip()


def MagicEventSource(f, *args, **kwargs):
    """
    Read events from either HDF or LIGO-LW files. The format is determined
    using the POSIX `file` command, which determines the file by looking for
    'magic' byte strings (hence the name of this module).

    Returns
    -------
    source : `EventSource`
    """
    if isinstance(f, h5py.File):
        opener = hdf.open
    elif isinstance(f, sqlite3.Connection):
        opener = sqlite.open
    elif isinstance(f, Element):
        opener = ligolw.open
    else:
        filetype = _get_file_type(f)
        if filetype == b'Hierarchical Data Format (version 5) data':
            opener = hdf.open
        elif filetype.startswith(b'SQLite 3.x database'):
            opener = sqlite.open
        elif filetype.startswith(b'XML') or filetype.startswith(b'gzip'):
            opener = ligolw.open
        else:
            raise IOError('Unknown file format')
    return opener(f, *args, **kwargs)


open = MagicEventSource
