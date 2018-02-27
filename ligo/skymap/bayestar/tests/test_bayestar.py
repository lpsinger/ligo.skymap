from collections import namedtuple

import lal
import pytest

from ..sky_map import localize, rasterize


MockSingleEvent = namedtuple(
    'MockSingleEvent', 'detector snr phase time psd snr_series')

MockEvent = namedtuple(
    'MockEvent', 'singles template_args')

template_args = {'mass1': 1.4, 'mass2': 1.4}


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
    test_single_event = MockSingleEvent('H1', 10, 0, 0, psd, None)
    test_event = MockEvent([test_single_event], template_args)
    skymap = localize(test_event)
    rasterized = rasterize(skymap)
    # FIXME: work out what this should be
