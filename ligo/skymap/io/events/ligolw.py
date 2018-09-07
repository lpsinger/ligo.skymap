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
Read events from pipedown/GstLal-style XML output.
"""
from collections import OrderedDict, defaultdict
import errno
from functools import lru_cache
import itertools
import logging
import operator
import os.path

from glue.ligolw import array, lsctables, param
from glue.ligolw.ligolw import Element, LIGOLWContentHandler, LIGO_LW
from glue.ligolw.lsctables import (
    CoincDefTable, CoincMapTable, CoincTable, ProcessTable, ProcessParamsTable,
    SnglInspiralID, SnglInspiralTable, TimeSlideID, TimeSlideTable)
from glue.ligolw.utils import load_filename, load_fileobj
import lal
import lal.series
from lalinspiral.thinca import InspiralCoincDef

from .base import Event, EventSource, SingleEvent

__all__ = ('LigoLWEventSource',)

log = logging.getLogger('BAYESTAR')


@array.use_in
@param.use_in
@lsctables.use_in
class ContentHandler(LIGOLWContentHandler):
    pass


def _read_xml(f, fallbackpath=None):
    if f is None:
        doc = filename = None
    elif isinstance(f, Element):
        doc = f
        filename = ''
    elif isinstance(f, str):
        try:
            doc = load_filename(f, contenthandler=ContentHandler)
        except IOError as e:
            if e.errno == errno.ENOENT and fallbackpath and \
                    not os.path.isabs(f):
                f = os.path.join(fallbackpath, f)
                doc = load_filename(f, contenthandler=ContentHandler)
            else:
                raise
        filename = f
    else:
        doc, _ = load_fileobj(f, contenthandler=ContentHandler)
        try:
            filename = f.name
        except AttributeError:
            filename = ''
    return doc, filename


class LigoLWEventSource(OrderedDict, EventSource):

    def __init__(self, f, psd_file=None, coinc_def=InspiralCoincDef,
                 fallbackpath=None, **kwargs):
        doc, filename = _read_xml(f)
        self._fallbackpath = (
            os.path.dirname(filename) if filename else fallbackpath)
        self._psds_for_file = lru_cache(maxsize=None)(self._psds_for_file)
        super(LigoLWEventSource, self).__init__(
            self._make_events(doc, psd_file, coinc_def))

    _template_keys = '''mass1 mass2
                        spin1x spin1y spin1z spin2x spin2y spin2z'''.split()

    _invert_phases = {
        'pycbc': False,
        'gstlal_inspiral': True,
        'gstlal_inspiral_postcohspiir_online': True,  # FIXME: wild guess
        'bayestar_realize_coincs': True,
        'bayestar-realize-coincs': True,
        'MBTAOnline': True
    }

    @classmethod
    def _phase_convention(cls, program):
        try:
            return cls._invert_phases[program]
        except KeyError:
            raise KeyError(
                ('The pipeline "{}" is unknown, so the phase '
                 'convention could not be deduced.').format(program))

    def _psds_for_file(self, f):
        doc, _ = _read_xml(f, self._fallbackpath)
        return lal.series.read_psd_xmldoc(doc, root_name=None)

    def _make_events(self, doc, psd_file, coinc_def):
        # Look up necessary tables.
        coinc_table = CoincTable.get_table(doc)
        coinc_map_table = CoincMapTable.get_table(doc)
        sngl_inspiral_table = SnglInspiralTable.get_table(doc)
        try:
            time_slide_table = TimeSlideTable.get_table(doc)
        except ValueError:
            offsets_by_time_slide_id = None
        else:
            offsets_by_time_slide_id = time_slide_table.as_dict()

        # Indices to speed up lookups by ID.
        key = operator.attrgetter('coinc_event_id')
        event_ids_by_coinc_event_id = {
            coinc_event_id:
                tuple(coinc_map.event_id for coinc_map in coinc_maps)
            for coinc_event_id, coinc_maps
            in itertools.groupby(sorted(coinc_map_table, key=key), key=key)}
        sngl_inspirals_by_event_id = {
            row.event_id: row for row in sngl_inspiral_table}

        # Filter rows by coinc_def if requested.
        if coinc_def is not None:
            coinc_def_table = CoincDefTable.get_table(doc)
            coinc_def_ids = {
                row.coinc_def_id for row in coinc_def_table
                if (row.search, row.search_coinc_type) ==
                (coinc_def.search, coinc_def.search_coinc_type)}
            coinc_table = [
                row for row in coinc_table
                if row.coinc_def_id in coinc_def_ids]

        snr_dict = dict(self._snr_series_by_sngl_inspiral(doc))

        process_table = ProcessTable.get_table(doc)
        program_for_process_id = {
            row.process_id: row.program for row in process_table}

        try:
            process_params_table = ProcessParamsTable.get_table(doc)
        except ValueError:
            psd_filenames_by_process_id = {}
        else:
            psd_filenames_by_process_id = {
                process_param.process_id: process_param.value
                for process_param in process_params_table
                if process_param.param == '--reference-psd'}

        ts0 = TimeSlideID(0)
        for time_slide_id in {coinc.time_slide_id for coinc in coinc_table}:
            if offsets_by_time_slide_id is None and time_slide_id == ts0:
                log.warning(
                    'Time slide record is missing for %s, '
                    'guessing that this is zero-lag', time_slide_id)

        for program in {program_for_process_id[coinc.process_id]
                        for coinc in coinc_table}:
            invert_phases = self._phase_convention(program)
            if invert_phases:
                log.warning(
                    'Using anti-FINDCHIRP phase convention; inverting phases. '
                    'This is currently the default and it is appropriate for '
                    'gstlal and MBTA but not pycbc as of observing run 1 '
                    '("O1"). The default setting is likely to change in the '
                    'future.')

        for coinc in coinc_table:
            coinc_event_id = coinc.coinc_event_id
            coinc_event_num = int(coinc_event_id)
            sngls = [sngl_inspirals_by_event_id[event_id] for event_id
                     in event_ids_by_coinc_event_id[coinc_event_id]]
            if offsets_by_time_slide_id is None and coinc.time_slide_id == ts0:
                offsets = defaultdict(float)
            else:
                offsets = offsets_by_time_slide_id[coinc.time_slide_id]

            template_args = [
                {key: getattr(sngl, key) for key in self._template_keys}
                for sngl in sngls]
            if any(d != template_args[0] for d in template_args[1:]):
                raise ValueError(
                    'Template arguments are not identical for all detectors!')
            template_args = template_args[0]

            invert_phases = self._phase_convention(
                program_for_process_id[coinc.process_id])

            singles = tuple(LigoLWSingleEvent(
                self, sngl.ifo, sngl.snr, sngl.coa_phase,
                float(sngl.end + offsets[sngl.ifo]), float(sngl.end),
                psd_file or psd_filenames_by_process_id.get(sngl.process_id),
                snr_dict.get(sngl.event_id), invert_phases)
                for sngl in sngls)

            event = LigoLWEvent(coinc_event_num, singles, template_args)

            yield coinc_event_num, event

    @classmethod
    def _snr_series_by_sngl_inspiral(cls, doc):
        for elem in doc.getElementsByTagName(LIGO_LW.tagName):
            try:
                if elem.Name != lal.COMPLEX8TimeSeries.__name__:
                    continue
                array.get_array(elem, 'snr')
                event_id = param.get_pyvalue(elem, 'event_id')
                if not isinstance(event_id, SnglInspiralID):
                    continue
            except (AttributeError, ValueError):
                continue
            else:
                yield event_id, lal.series.parse_COMPLEX8TimeSeries(elem)


class LigoLWEvent(Event):

    def __init__(self, id, singles, template_args):
        self._id = id
        self._singles = singles
        self._template_args = template_args

    @property
    def singles(self):
        return self._singles

    @property
    def template_args(self):
        return self._template_args


class LigoLWSingleEvent(SingleEvent):

    def __init__(self, source, detector, snr, phase, time, zerolag_time,
                 psd_file, snr_series, invert_phases):
        self._source = source
        self._detector = detector
        self._snr = snr
        self._phase = phase
        self._time = time
        self._zerolag_time = zerolag_time
        self._psd_file = psd_file
        self._snr_series = snr_series
        self._invert_phases = invert_phases

    @property
    def detector(self):
        return self._detector

    @property
    def snr(self):
        return self._snr

    @property
    def phase(self):
        value = self._phase
        if value is not None and self._invert_phases:
            value *= -1
        return value

    @property
    def time(self):
        return self._time

    @property
    def zerolag_time(self):
        return self._zerolag_time

    @property
    def psd(self):
        return self._source._psds_for_file(self._psd_file)[self._detector]

    @property
    def snr_series(self):
        value = self._snr_series
        if self._invert_phases and value is not None:
            value = lal.CutCOMPLEX8TimeSeries(value, 0, len(value.data.data))
            value.data.data = value.data.data.conj()
        return value


open = LigoLWEventSource
