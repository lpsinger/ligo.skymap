from contextlib import ExitStack
import errno
import gzip
import os
from pathlib import Path
import re
import subprocess

from ligo.lw.utils import load_filename
import h5py
import numpy as np
import pytest

from .. import events
from ..events.ligolw import ContentHandler
from ...util import sqlite

DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')


@pytest.fixture
def mock_gracedb(monkeypatch):
    stack = ExitStack()

    class MockGraceDb(ExitStack):

        def files(self, graceid, filename):
            path = os.path.join(DATA_PATH, '{}_{}'.format(graceid, filename))
            try:
                return stack.enter_context(open(path, 'rb'))
            except IOError as e:
                if e.errno != errno.ENOENT:
                    raise
                return stack.enter_context(gzip.GzipFile(path + '.gz', 'rb'))

    with stack:
        monkeypatch.setattr('ligo.gracedb.rest.GraceDb', MockGraceDb)
        yield


def raises(expected_exception, msg):
    return pytest.raises(expected_exception, match='^' + re.escape(msg) + '$')


def ligolw_assertions(source):
    """Test common assertions for test_ligolw and test_sqlite."""
    assert len(source) == 250
    event = source[821759]
    assert len(event.singles) == 2
    assert event.singles[0].snr_series is event.singles[1].snr_series is None
    assert event.singles[0].detector == 'H1'
    assert event.singles[0].snr == 8.5362396
    assert event.singles[0].phase == -0.81192881
    assert (event.singles[0].time == event.singles[0].zerolag_time ==
            967328914.866107842)
    psd = event.singles[0].psd
    assert psd.f0 == 0.0
    assert psd.deltaF == 0.125
    assert event.singles[1].detector == 'L1'
    assert event.singles[1].snr == 10.36818
    assert event.singles[1].phase == 1.9740163
    assert (event.singles[1].time == event.singles[1].zerolag_time ==
            967328914.866513726)
    psd = event.singles[0].psd
    assert psd.f0 == 0.0
    assert psd.deltaF == 0.125
    assert event.template_args == {
        'mass1': 1.56687,
        'mass2': 1.474779,
        'spin1x': 0.0,
        'spin1y': 0.0,
        'spin1z': 0.0,
        'spin2x': 0.0,
        'spin2y': 0.0,
        'spin2z': 0.0,
        'f_final': 2047.0}


def test_unknown():
    with pytest.raises(IOError, match='Unknown file format'):
        events.open(__file__)


def test_ligolw():
    """Test reading events from LIGO-LW XML files."""
    source = events.open(os.path.join(DATA_PATH, '2016_subset.xml.gz'))
    ligolw_assertions(source)

    assert str(source).startswith('<LigoLWEventSource({288172: <LigoLWEvent(')
    assert repr(source).startswith('LigoLWEventSource({288172: <LigoLWEvent(')
    event, *_ = source.values()
    assert str(event).startswith(
        "<LigoLWEvent(singles=(<LigoLWSingleEvent(detector='H1', snr=12.")
    assert repr(event).startswith(
        "<LigoLWEvent(singles=(<LigoLWSingleEvent(detector='H1', snr=12.")
    single_event, *_ = event.singles
    assert str(single_event) == "<LigoLWSingleEvent(detector='H1', snr=12.035994, phase=-1.0371021, time=970976257.4808338)>"  # noqa
    assert repr(single_event) == "<LigoLWSingleEvent(detector='H1', snr=12.035994, phase=-1.0371021, time=970976257.4808338)>"  # noqa


def test_ligolw_document():
    """Test reading events from LIGO-LW XML document."""
    xmldoc = load_filename(
        os.path.join(DATA_PATH, '2016_subset.xml.gz'),
        contenthandler=ContentHandler)
    source = events.open(xmldoc, fallbackpath=DATA_PATH)
    ligolw_assertions(source)


def test_ligolw_pathlib():
    """Test reading file provided with Pathlib instead of string"""
    source = events.open(Path(DATA_PATH, '2016_subset.xml.gz'))
    ligolw_assertions(source)


def test_sqlite(tmpdir):
    """Test reading events from SQLite files."""
    # Convert the test data to SQLite format in the temporary directory.
    xmlfilename = os.path.join(DATA_PATH, '2016_subset.xml.gz')
    dbfilename = str(tmpdir / '2016_subset.sqlite')
    subprocess.check_call(['ligolw_sqlite', '-p',
                           xmlfilename, '-d', dbfilename])

    # Symbolicly link the PSD directory into the temporary directory.
    os.symlink(os.path.join(DATA_PATH, 'gstlal_reference_psd'),
               str(tmpdir / 'gstlal_reference_psd'))

    source = events.open(dbfilename)
    ligolw_assertions(source)

    with open(dbfilename, 'rb') as f:
        source = events.open(f)
        ligolw_assertions(source)

    with sqlite.open(dbfilename, 'r') as db:
        source = events.open(db)
        ligolw_assertions(source)


def test_gracedb(mock_gracedb):
    """Test reading events from GraceDB records."""
    source = events.gracedb.open(['G211117', 'G197392'])
    assert len(source) == 2
    for i, (event_id, event) in enumerate(source.items()):
        if i == 0:
            assert event_id == 'G211117'
            assert (event.singles[0].snr_series is event.singles[1].snr_series
                    is None)
            assert event.singles[0].detector == 'H1'
            assert event.singles[0].snr == 9.0802174
            assert event.singles[0].phase == -0.13969257
            assert (event.singles[0].time == event.singles[0].zerolag_time ==
                    1135136350.647757924)
            psd = event.singles[0].psd
            assert psd.f0 == 0.0
            assert psd.deltaF == 0.125
            assert event.singles[1].detector == 'L1'
            assert event.singles[1].snr == 7.3947201
            assert event.singles[1].phase == -2.7356486
            assert (event.singles[1].time == event.singles[1].zerolag_time ==
                    1135136350.646883043)
            psd = event.singles[1].psd
            assert psd.f0 == 0.0
            assert psd.deltaF == 0.125
            assert event.template_args == {
                'mass1': 19.924686,
                'mass2': 6.4254546,
                'spin1x': 0.0,
                'spin1y': 0.0,
                'spin1z': 0.33962944,
                'spin2x': 0.0,
                'spin2y': 0.0,
                'spin2z': -0.1238557,
                'f_final': 1024.0}
        elif i == 1:
            assert event_id == 'G197392'
            assert (event.singles[0].snr_series is event.singles[1].snr_series
                    is None)
            assert event.singles[0].detector == 'H1'
            assert event.singles[0].snr == 6.9068823
            assert event.singles[0].phase == 1.8298783
            assert (event.singles[0].time == event.singles[0].zerolag_time ==
                    1128678900.444335938)
            psd = event.singles[0].psd
            assert psd.f0 == 30.0
            assert psd.deltaF == 0.125
            assert event.singles[1].detector == 'L1'
            assert event.singles[1].snr == 6.8389997
            assert event.singles[1].phase == -1.0297496
            assert (event.singles[1].time == event.singles[1].zerolag_time ==
                    1128678900.445068359)
            psd = event.singles[1].psd
            assert psd.f0 == 30.0
            assert psd.deltaF == 0.125
            assert event.template_args == {
                'mass1': 32.064007,
                'mass2': 14.607587,
                'spin1x': 0.0,
                'spin1y': 0.0,
                'spin1z': 0.34881824,
                'spin2x': 0.0,
                'spin2y': 0.0,
                'spin2z': -0.53029484,
                'f_final': 0.0}


def test_detector_disabled(mock_gracedb):
    """Test reading from event sources with certain detectors disabled."""
    graceids = ('G211117', 'G197392')
    base_source = events.gracedb.open(graceids)

    source = events.detector_disabled.open(base_source, ['H1'])
    assert len(source) == 2
    for graceid, (event_id, event) in zip(graceids, source.items()):
        assert event_id == graceid
        assert len(event.singles) == 1
        assert event.singles[0].detector == 'L1'

    for event, base_event in zip(source.values(), base_source.values()):
        assert event.template_args == base_event.template_args

    # Now test that exceptions are raised when they are called for.
    expected_message = ('Disabling detectors {H1, L1} would have no effect on '
                        'this event with detectors {H1 L1}')
    nonraising_source = events.detector_disabled.open(
        base_source, ['H1, L1'], raises=False)
    raising_source = events.detector_disabled.open(
        base_source, ['H1, L1'])
    for event in nonraising_source.values():
        event.singles
    for event in raising_source.values():
        with raises(events.DetectorDisabledError, expected_message):
            event.singles

    # Now test that exceptions are raised when they are called for.
    expected_message = ('Disabling detectors {H1 L1 V1} would exclude all '
                        'data for this event with detectors {H1 L1}')
    nonraising_source = events.detector_disabled.open(
        base_source, ['H1', 'L1', 'V1'], raises=False)
    raising_source = events.detector_disabled.open(
        base_source, ['H1', 'L1', 'V1'])
    for event in nonraising_source.values():
        event.singles
    for event in raising_source.values():
        with raises(events.DetectorDisabledError, expected_message):
            event.singles


def test_hdf(tmpdir):
    """Test reading events from HDF5 files."""
    # Create test input files
    ifos = ['L1', 'H1']
    filenames = []
    filename = str(tmpdir / 'coincs.hdf')
    filenames.append(filename)
    with h5py.File(filename, 'w') as coinc_file:
        coinc_file.attrs['timeslide_interval'] = np.e
        coinc_group = coinc_file.create_group('foreground')
        coinc_group['template_id'] = np.arange(5)
        coinc_group['timeslide_id'] = np.arange(0, 50, 10)

        filename = str(tmpdir / 'H1L1-BANK.hdf')
        filenames.append(filename)
        with h5py.File(filename, 'w') as bank_file:
            bank_file.attrs['parameters'] = []

        for i, ifo in enumerate(ifos):
            coinc_file.attrs['detector_{}'.format(i + 1)] = ifo
            coinc_group['trigger_id{}'.format(i + 1)] = np.arange(5)

            filename = str(tmpdir / (ifo + '_triggers.hdf'))
            filenames.append(filename)
            with h5py.File(filename, 'w') as trigger_file:
                trigger_group = trigger_file.create_group(ifo)
                trigger_group['snr'] = i + np.arange(5) * np.pi
                trigger_group['coa_phase'] = i + np.arange(5) * np.pi**2
                trigger_group['end_time'] = i + np.arange(5) * np.pi**3

            filename = str(tmpdir / (ifo + '_psds.hdf'))
            filenames.append(filename)
            with h5py.File(filename, 'w') as psd_file:
                psd_file.attrs['low_frequency_cutoff'] = 10.0
                psd_file.attrs['dynamic_range_factor'] = np.e
                psd_group = psd_file.create_group(ifo)
                psd_group['start_time'] = (np.arange(5) - 0.5) * np.pi**3
                psd_group['end_time'] = (np.arange(5) + 0.5) * np.pi**3
                psd_group = psd_group.create_group('psds')
                for j in range(5):
                    psd_group[str(j)] = np.concatenate(
                        (np.zeros(5), np.arange(5)**2)) * np.e**2
                    psd_group[str(j)].attrs['delta_f'] = 2.0

    with pytest.raises(
            ValueError,
            match='You must provide exactly one coinc file.'):
        source = events.open(*filenames[1:])

    with pytest.raises(
            ValueError,
            match='You must provide exactly one template bank file.'):
        source = events.open(*filenames[:1])

    with pytest.raises(
            ValueError,
            match='You must provide PSD files.'):
        source = events.open(*filenames[:2])

    with pytest.raises(
            ValueError,
            match='You must provide trigger files.'):
        source = events.open(*filenames[:2], filenames[3], filenames[5])

    # Test reading from filenames
    source = events.open(*filenames)
    assert len(source) == 5
    for coinc_id, coinc in source.items():
        for i, (ifo, single) in enumerate(zip(ifos, coinc.singles)):
            assert single.detector == ifo
            assert single.snr == i + coinc_id * np.pi
            assert single.phase == i + coinc_id * np.pi**2
            assert single.time == (
                single.zerolag_time +
                coinc_id * 10 * np.e * (-0.5 if i == 0 else +0.5))
            assert single.zerolag_time == i + coinc_id * np.pi**3
            assert single.psd.f0 == 10.0
            assert single.psd.deltaF == 2.0
            assert np.all(single.psd.data.data == np.arange(5)**2)
        assert coinc.template_args == {}

    # Test reading from h5py.File instances
    events.open(*(h5py.File(filename, 'r') for filename in filenames))

    # Test reading from file-like objects
    with ExitStack() as stack:
        events.open(*(
            stack.enter_context(open(filename, 'rb'))
            for filename in filenames))
