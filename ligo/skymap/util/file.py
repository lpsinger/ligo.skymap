#
# Copyright (C) 2018-2019  Leo Singer
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
"""File tools."""
import errno
import os
import shutil
import tempfile


def rename(src, dst):
    """Like `os.rename`, but works across different devices because it
    catches and handles ``EXDEV`` (``Invalid cross-device link``) errors.
    """
    try:
        os.rename(src, dst)
    except OSError as e:
        if e.errno == errno.EXDEV:
            dir, suffix = os.path.split(dst)
            tmpfid, tmpdst = tempfile.mkstemp(dir=dir, suffix=suffix)
            try:
                os.close(tmpfid)
                shutil.copy2(src, tmpdst)
                os.rename(tmpdst, dst)
            except:  # noqa: E722
                os.remove(tmpdst)
                raise
        else:
            raise


def rm_f(filename):
    """Remove a file, or be silent if the file does not exist, like ``rm -f``.

    Examples
    --------
    >>> with tempfile.TemporaryDirectory() as d:
    ...     rm_f('test')
    ...     with open('test', 'w') as f:
    ...         print('Hello world', file=f)
    ...     rm_f('test')

    """
    try:
        os.remove(filename)
    except FileNotFoundError:
        pass
