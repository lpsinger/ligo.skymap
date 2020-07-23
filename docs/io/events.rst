Event Candidate I/O (`ligo.skymap.io.events`)
=============================================
.. module:: ligo.skymap.io.events

This module provides a high-level interface for reading event records from
LIGO/Virgo/KAGRA search pipelines.

* To read offline events from GstLAL, PyCBC, or
  :doc:`bayestar-realize-coincs </tool/bayestar_realize_coincs>`,
  use :obj:`ligo.skymap.io.events.open`. It can read GstLAL's LIGO-LW XML or
  SQLite files and PyCBC's HDF5 files produced by PyCBC.

* To read online events from GraceDB, use
  :obj:`ligo.skymap.io.events.gracedb.open`.

* To read events from any other source or format, write your own concrete
  event classes inheriting from the :class:`ligo.skymap.io.events.Event` and
  :class:`ligo.skymap.io.events.SingleEvent` base classes.

High-Level Openers
------------------

.. autofunction:: open
.. autofunction:: ligo.skymap.io.events.gracedb.open

Openers for Specific Formats
----------------------------

.. autofunction:: ligo.skymap.io.events.hdf.open
.. autofunction:: ligo.skymap.io.events.ligolw.open
.. autofunction:: ligo.skymap.io.events.sqlite.open

Abstract Base Classes
---------------------

To implement your own event format, write your own concrete classes using
:class:`~ligo.skymap.io.events.Event` and
:class:`~ligo.skymap.io.events.SingleEvent` as base classes. Here is an example
of a concrete class that produces events populated with random values::

    from lal import CreateREAL8FrequencySeries
    from lalsimulation import SimNoisePSD, SimNoisePSDaLIGOZeroDetHighPowerPtr
    from ligo.skymap.io.events import Event, SingleEvent
    import numpy as np

    class RandomEvent(Event):

        def __init__(self):
            self._singles = [RandomSingleEvent('H1'), RandomSingleEvent('L1')]
            self._template_args = {'mass1': np.random.uniform(1.0, 3.0),
                                   'mass1': np.random.uniform(1.0, 3.0)}

        @property
        def singles(self):
            return self._singles

        @property
        def template_args(self):
            return self._template_args


    class RandomSingleEvent(SingleEvent):

        def __init__(self, detector):
            self._detector = detector
            self._snr = np.random.uniform(1, 10)
            self._phase = np.random.uniform(-np.pi, np.pi)
            self._time = np.random.uniform(1e9, 1.1e9)
            self._psd = CreateREAL8FrequencySeries(name=None, epoch=None, f0=0.0,
                                                   deltaF=1.0, sampleUnits=None,
                                                   length=16384)
            SimNoisePSD(psd=self._psd, flow=10.0,
                        psdfunc=SimNoisePSDaLIGOZeroDetHighPowerPtr)
            self._psd.data.data[self._psd.data.data == 0] = np.inf

        @property
        def detector(self):
            return self._detector

        @property
        def snr(self):
            return self._snr

        @property
        def phase(self):
            return self._phase

        @property
        def time(self):
            return self._time

        @property
        def zerolag_time(self):
            return self._time

        @property
        def psd(self):
            return self._psd


    np.random.seed(0)
    event = RandomEvent()
    print(event)
    # output:
    # <RandomEvent(singles=[<RandomSingleEvent(detector='H1',
    # snr=5.939321535345923, phase=1.3520746650524709, time=1060276337.6071644)>,
    # <RandomSingleEvent(detector='L1', snr=5.903948646972072,
    # phase=-0.47969104306747123, time=1064589411.3066657)>])>


.. autoclass:: EventSource
    :show-inheritance:
.. autoclass:: Event
.. autoclass:: SingleEvent
