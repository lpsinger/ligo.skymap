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
Read events from PyCBC-style HDF5 output.
"""
from .base import *
from operator import itemgetter
from itertools import groupby
import h5py
import numpy as np
import lal
from glue.segments import segment, segmentlist

__all__ = ('HDFEventSource',)


class _psd_segment(segment):

    def __new__(cls, psd, *args):
        return segment.__new__(cls, *args)

    def __init__(self, psd, *args):
        self.psd = psd


def _hdf_file(f):
    if isinstance(f, h5py.File):
        return f
    elif hasattr(f, 'read') and hasattr(f, 'name'):
        return h5py.File(f.name, 'r')
    else:
        return h5py.File(f, 'r')


def _classify_hdf_file(f, sample):
    if sample in f:
        return 'coincs'
    for key, value in f.items():
        if isinstance(value, h5py.Group):
            if 'psds' in value:
                return 'psds'
            if 'snr' in value and 'coa_phase' in value and 'end_time' in value:
                return 'triggers'
    if 'parameters' in f.attrs:
        return 'bank'
    raise ValueError('Unrecognized PyCBC file type')


class HDFEventSource(EventSource):

    def __init__(self, *files, **kwargs):
        sample = kwargs.get('sample', 'foreground')

        # Open the files and split them into coinc files, bank files, psds,
        # and triggers.
        key = itemgetter(0)
        files = [_hdf_file(f) for f in files]
        files = sorted(
            [(_classify_hdf_file(f, sample), f) for f in files], key=key)
        files = {key: list(v[1] for v in value)
                 for key, value in groupby(files, key)}

        try:
            coinc_file, = files['coincs']
        except (KeyError, ValueError):
            raise ValueError('You must provide exactly one coinc file.')
        try:
            bank_file, = files['bank']
        except (KeyError, ValueError):
            raise ValueError(
                'You must provide exactly one template bank file.')
        try:
            psd_files = files['psds']
        except KeyError:
            raise ValueError('You must provide PSD files.')
        try:
            trigger_files = files['triggers']
        except KeyError:
            raise ValueError('You must provide trigger files.')

        self._bank = bank_file

        key_prefix = 'detector_'
        detector_nums, self._ifos = zip(*sorted(
            (int(key[len(key_prefix):]), value)
            for key, value in coinc_file.attrs.items()
            if key.startswith(key_prefix)))

        coinc_group = coinc_file[sample]
        self._timeslide_interval = coinc_file.attrs.get(
            'timeslide_interval', 0)
        self._template_ids = coinc_group['template_id']
        self._timeslide_ids = coinc_group.get(
            'timeslide_id', np.zeros(len(self)))
        self._trigger_ids = [
            coinc_group['trigger_id{}'.format(detector_num)]
            for detector_num in detector_nums]

        triggers = {}
        for f in trigger_files:
            (ifo, group), = f.items()
            triggers[ifo] = [
                group['snr'], group['coa_phase'], group['end_time']]
        self._triggers = tuple(triggers[ifo] for ifo in self._ifos)

        psdseglistdict = {}
        for psd_file in psd_files:
            (ifo, group), = psd_file.items()
            psd = [group['psds'][str(i)] for i in range(len(group['psds']))]
            psdseglistdict[ifo] = segmentlist(
                _psd_segment(*segargs) for segargs in zip(
                    psd, group['start_time'], group['end_time']))
        self._psds = [psdseglistdict[ifo] for ifo in self._ifos]

    def __getitem__(self, id):
        return HDFEvent(self, id)

    def __iter__(self):
        return iter(range(len(self)))

    def __len__(self):
        return len(self._template_ids)


class HDFEvent(Event):

    def __init__(self, source, id):
        self._source = source
        self._id = id

    @property
    def singles(self):
        return tuple(
            HDFSingleEvent(
                ifo, self._id, i, trigger_ids[self._id],
                self._source._timeslide_interval, triggers,
                self._source._timeslide_ids, psds
            )
            for i, (ifo, trigger_ids, triggers, psds) in enumerate(zip(
                self._source._ifos, self._source._trigger_ids,
                self._source._triggers, self._source._psds
            ))
        )

    @property
    def template_args(self):
        bank = self._source._bank
        bank_id = self._source._template_ids[self._id]
        return {key: value[bank_id] for key, value in bank.items()}


class HDFSingleEvent(SingleEvent):

    def __init__(
            self, detector, _coinc_id, _detector_num, _trigger_id,
            _timeslide_interval, _triggers, _timeslide_ids, _psds):
        self._detector = detector
        self._coinc_id = _coinc_id
        self._detector_num = _detector_num
        self._trigger_id = _trigger_id
        self._timeslide_interval = _timeslide_interval
        self._triggers = _triggers
        self._timeslide_ids = _timeslide_ids
        self._psds = _psds

    @property
    def detector(self):
        return self._detector

    @property
    def snr(self):
        return self._triggers[0][self._trigger_id]

    @property
    def phase(self):
        return self._triggers[1][self._trigger_id]

    @property
    def time(self):
        value = self.zerolag_time

        # PyCBC does not specify which detector is time-shifted in time slides.
        # Since PyCBC's coincidence format currently supports only two
        # detectors, we arbitrarily apply half of the time slide to each
        # detector.
        shift = self._timeslide_ids[self._coinc_id] * self._timeslide_interval
        if self._detector_num == 0:
            value -= 0.5 * shift
        elif self._detector_num == 1:
            value += 0.5 * shift
        else:
            raise RuntimeError('This line should not be reached')
        return value

    @property
    def zerolag_time(self):
        return self._triggers[2][self._trigger_id]

    @property
    def psd(self):
        try:
            psd = self._psds[self._psds.find(self.zerolag_time)].psd
        except ValueError:
            raise ValueError(
                'No PSD found for detector {} at zero-lag GPS time {}'.format(
                    self.detector, self.zerolag_time))

        dyn_range_fac = psd.file.attrs['dynamic_range_factor']
        flow = psd.file.attrs['low_frequency_cutoff']
        df = psd.attrs['delta_f']
        kmin = int(flow / df)

        fseries = lal.CreateREAL8FrequencySeries(
            'psd', 0, kmin * df, df,
            lal.DimensionlessUnit, len(psd.value) - kmin)
        fseries.data.data = psd.value[kmin:] / np.square(dyn_range_fac)
        return fseries


open = HDFEventSource
