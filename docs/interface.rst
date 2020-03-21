###################################################
Interface definition for online detection pipelines
###################################################

The BAYESTAR rapid localization algorithm is designed as a post-processing
stage for compact binary merger search pipelines. The script
:doc:`bayestar-localize-lvalert </ligo/skymap/tool/bayestar_localize_lvalert>`
provide rapid sky localization as a service.

.. note::
    Other available modes of operation for BAYESTAR that are documented
    elsewhere include the script :doc:`bayestar-localize-coincs
    </ligo/skymap/tool/bayestar_localize_coincs>` for offline batch processing
    and the method :meth:`ligo.skymap.bayestar.localize` for directly calling
    BAYESTAR from Python code.

Sequence diagram
================

In online operation, search pipelines upload candidates to the
Gravitational-Wave Candidate Event Database (`GraceDB`_). The script
:doc:`bayestar-localize-lvalert </ligo/skymap/tool/bayestar_localize_lvalert>`
(or the equivalent :doc:`Celery <celery:index>` task in :doc:`GWCelery
<gwcelery:index>`) listens for and intercepts :doc:`LVAlert <gracedb:lvalert>`
pubsub messages. For each event created, the script downloads the pipeline's
output data products from GraceDB, performs rapid sky localization, and uploads
the resulting FITS file back to GraceDB.

The interactions between the search pipeline, GraceDB, and BAYESTAR are
illustrated in the `sequence diagram`_ below. Line styles have the following
meanings:

* Solid lines directed into GraceDB represent HTTP requests.
* Solid lines directed out of GraceDB represent HTTP responses.
* Dashed lines represent LVAlert pubsub messages.

.. mermaid::

    sequenceDiagram
    
        note over Search: New detection
        Search ->>+ GraceDB: Upload coinc.xml
        activate Search
        GraceDB ->>- Search: GraceDB ID: G123
        Search ->>+ GraceDB: Upload psd.xml.gz to G123
        deactivate Search
        GraceDB -->>- BAYESTAR: LVAlert: psd.xml.gz added to G123
        activate BAYESTAR
        BAYESTAR ->>+ GraceDB: Get coinc.xml from G123
        GraceDB ->>- BAYESTAR: coinc.xml
        BAYESTAR ->>+ GraceDB: Get psd.xml.gz from G123
        GraceDB ->>- BAYESTAR: psd.xml.gz
        note over BAYESTAR: Perform sky localization
        BAYESTAR ->> GraceDB: Upload bayestar.fits to G123
        deactivate BAYESTAR

Input files
===========

BAYESTAR requires the following two files attached to the GraceDB event:

* :file:`coinc.xml`: The initial event file.
* :file:`psd.xml.gz`: The power spectral density data file.
  This file *must* be uploaded to GraceDB with the :samp:`psd` tag applied.

The format of both files must be LIGO-LW (see :dcc:`T990023`). LIGO-LW is a
legacy XML-based format used by a variety of LIGO/Virgo/KAGRA software and
services for storing tabular datasets.

Unfortunately, LIGO-LW is a rather complicated format. We recommend using
either the :mod:`ligo.lw` module or GWPy's :ref:`tabular LIGO-LW I/O
<gwpy-table-io-ligolw>` feature to simplify reading and writing LIGO-LW files.

The :file:`coinc.xml` file
--------------------------

This file describes the search pipeline's matched filter output. It contains
the following information:

* Point estimates of the time, phase, and amplitude on arrival in each detector
* Intrinsic template parameters (masses and spins)
* Optionally, signal-to-noise time series for each detector

It must contain at least the
:class:`coinc <ligo.lw.lsctables.CoincTable>`,
:class:`coinc_event_map <ligo.lw.lsctables.CoincMapTable>`,
:class:`process <ligo.lw.lsctables.ProcessTable>`, and
:class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>`
LIGO-LW tables.

The :file:`psd.xml.gz` file
---------------------------

This file contains each analyzed detectors' estimated noise power spectral
density series.

.. _`GraceDB`: https://gracedb.ligo.org
.. _`sequence diagram`: https://en.wikipedia.org/wiki/Sequence_diagram
