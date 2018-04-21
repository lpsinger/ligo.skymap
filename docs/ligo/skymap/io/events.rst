Event Candidate I/O (`ligo.skymap.io.events`)
=============================================
.. module:: ligo.skymap.io.events

This module provides a high-level interface for reading event records from
search pipelines.

Openers
-------

.. autofunction:: open
.. autofunction:: ligo.skymap.io.events.gracedb.open

Abstract Base Classes
---------------------

.. autoclass:: EventSource
    :show-inheritance:
.. autoclass:: Event
.. autoclass:: SingleEvent
