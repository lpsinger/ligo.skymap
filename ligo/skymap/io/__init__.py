# Copyright (C) 2016-2018  Leo Singer
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
from __future__ import absolute_import
import os
import pkgutil
import six

__all__ = ()

# Import all symbols from all submodules of this module.
for _, module, _ in pkgutil.iter_modules([os.path.dirname(__file__)]):
    six.exec_('from . import {0};'
              '__all__ += getattr({0}, "__all__", ());'
              'from .{0} import *'.format(module))
    del module

# Clean up
del os, pkgutil, six
