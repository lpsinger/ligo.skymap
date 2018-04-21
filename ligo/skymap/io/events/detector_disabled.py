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
Modify events by artificially disabling specified detectors.
"""
from .base import Event, EventSource

__all__ = ('DetectorDisabledEventSource', 'DetectorDisabledError')


class DetectorDisabledError(ValueError):
    pass


class DetectorDisabledEventSource(EventSource):

    def __init__(self, base_source, disabled_detectors, raises=True):
        self.base_source = base_source
        self.disabled_detectors = set(disabled_detectors)
        self.raises = raises

    def __iter__(self):
        return iter(self.base_source)

    def __getitem__(self, key):
        return DetectorDisabledEvent(self, self.base_source[key])

    def __len__(self):
        return len(self.base_source)


class DetectorDisabledEvent(Event):

    def __init__(self, source, base_event):
        self.source = source
        self.base_event = base_event

    @property
    def singles(self):
        disabled_detectors = self.source.disabled_detectors
        if self.source.raises:
            detectors = {s.detector for s in self.base_event.singles}
            if not detectors & disabled_detectors:
                raise DetectorDisabledError(
                    'Disabling detectors {{{}}} would have no effect on this '
                    'event with detectors {{{}}}'.format(
                        ' '.join(sorted(disabled_detectors)),
                        ' '.join(sorted(detectors))))
            if not detectors - disabled_detectors:
                raise DetectorDisabledError(
                    'Disabling detectors {{{}}} would exclude all data for '
                    'this event with detectors {{{}}}'.format(
                        ' '.join(sorted(disabled_detectors)),
                        ' '.join(sorted(detectors))))
        return tuple(s for s in self.base_event.singles
                     if s.detector not in disabled_detectors)

    @property
    def template_args(self):
        return self.base_event.template_args


open = DetectorDisabledEventSource
