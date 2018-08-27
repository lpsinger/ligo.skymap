from collections import namedtuple

import lal
import pytest

from ...io.events.base import Event, SingleEvent
from .. import localize, rasterize


_MockSingleEvent = namedtuple(
    'MockSingleEvent', 'detector snr phase time psd snr_series')

_MockEvent = namedtuple(
    'MockEvent', 'singles template_args')


class MockSingleEvent(_MockSingleEvent, SingleEvent):
    pass


class MockEvent(_MockEvent, Event):
    pass


template_args = {'mass1': 1.414, 'mass2': 1.414}


def test_localize_0_detectors():
    """Running on an event with 0 detectors should raise an error. """
    test_event = MockEvent([], template_args)
    with pytest.raises(ValueError):
        localize(test_event)


def test_localize_1_detector():
    """Running on an event with 1 detector should produce a sky map that
    reflects the antenna pattern."""
    psd = lal.CreateREAL8FrequencySeries(
        None, 0, 0, 32, lal.DimensionlessUnit, 128)
    psd.data.data[:] = 1
    test_single_event = MockSingleEvent(
        'H1', 12.345, 0.6789, 0.1234, psd, None)
    test_event = MockEvent([test_single_event], template_args)
    skymap = localize(test_event)

    # Make sure that none of the extrinsic parameters are in the event history
    history = '\n'.join(skymap.meta['history'])
    for forbidden in ['snr=', '12.345', 'time=', '0.6789', 'phase=', '0.1234',
                      'mass1=', 'mass2=', '1.414']:
        assert forbidden not in history

    # FIXME: work out what this should be
    rasterized = rasterize(skymap)
    assert rasterized
