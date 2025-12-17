#
# Copyright (C) 2025  Leo Singer
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


import logging

__all__ = ("basicConfig",)

_basic_config_kwargs = None


def basicConfig(**kwargs):
    """Configure root logger, but save kwargs to restore later in subprocesses."""
    global _basic_config_kwargs
    _basic_config_kwargs = kwargs
    logging.basicConfig(**kwargs)
