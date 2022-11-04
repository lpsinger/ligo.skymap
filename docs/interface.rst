###################################################
Interface definition for online detection pipelines
###################################################

The BAYESTAR rapid localization algorithm is designed as a post-processing
stage for compact binary merger search pipelines. The script
:doc:`bayestar-localize-lvalert </tool/bayestar_localize_lvalert>`
provides rapid sky localization as a service.

.. note::
    Other available modes of operation for BAYESTAR that are documented
    elsewhere include the script :doc:`bayestar-localize-coincs
    </tool/bayestar_localize_coincs>` for offline batch processing
    and the method :meth:`ligo.skymap.bayestar.localize` for directly calling
    BAYESTAR from Python code.

Sequence diagram
================

In online operation, search pipelines upload candidates to the
Gravitational-Wave Candidate Event Database (`GraceDB`_). The script
:doc:`bayestar-localize-lvalert </tool/bayestar_localize_lvalert>`
(or the equivalent :doc:`Celery <celery:index>` task in :doc:`GWCelery
<gwcelery:index>`) listens for and intercepts
:doc:`IGWN Alert <gracedb:igwn_alert>` messages. For each event created, the
script downloads the pipeline's output data products from GraceDB, performs
rapid sky localization, and uploads the resulting FITS file back to GraceDB.

The interactions between the search pipeline, GraceDB, and BAYESTAR are
illustrated in the `sequence diagram`_ below. Line styles have the following
meanings:

* Solid lines directed into GraceDB represent HTTP requests.
* Solid lines directed out of GraceDB represent HTTP responses.
* Dashed lines represent IGWN Alert messages.

Sequence diagram for unified ``coinc.xml`` file
-----------------------------------------------

.. mermaid::

    sequenceDiagram

        note over Search: New detection
        Search ->>+ GraceDB: Upload coinc.xml
        note over GraceDB: Create event G123
        GraceDB -->>+ BAYESTAR: IGWN Alert: coinc.xml added to G123
        deactivate GraceDB
        BAYESTAR ->>+ GraceDB: Get coinc.xml from G123
        GraceDB ->>- BAYESTAR: coinc.xml
        note over BAYESTAR: Perform sky localization
        BAYESTAR ->> GraceDB: Upload bayestar.fits to G123
        deactivate BAYESTAR

Sequence diagram for separate ``coinc.xml`` and ``psd.xml.gz`` files
--------------------------------------------------------------------

.. mermaid::

    sequenceDiagram
    
        note over Search: New detection
        Search ->>+ GraceDB: Upload coinc.xml
        activate Search
        note over GraceDB: Create event G123
        GraceDB ->>- Search: GraceDB ID: G123
        Search ->>+ GraceDB: Upload psd.xml.gz to G123
        deactivate Search
        GraceDB -->>- BAYESTAR: IGWN Alert: psd.xml.gz added to G123
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

This section describes the interface between search pipelines and BAYESTAR. The
key words "MUST", "MUST NOT", "REQUIRED", "SHALL", "SHALL NOT", "SHOULD",
"SHOULD NOT", "RECOMMENDED", "MAY", and "OPTIONAL" in this document are to be
interpreted as described in :rfc:`2119`.

The following file MUST be uploaded to GraceDB:

* :file:`coinc.xml`: The event file, which SHOULD be the initial upload that
  creates the event.

The contents of the :file:`coinc.xml` file MUST conform to the
:ref:`event-data` section below. The :file:`coinc.xml` file SHOULD also contain
the data described in the :ref:`psd-data` section below. If the
:file:`coinc.xml` file does not include the PSD data, then the following
additional file MUST be uploaded to GraceDB:

* :file:`psd.xml.gz`: The power spectral density data file,
  which MUST be uploaded with the :samp:`psd` tag.

If the :file:`psd.xml.gz` is uploaded, then its contents MUST conform to the
:ref:`psd-data` section below.

The format of both files MUST be LIGO-LW (see :dcc:`T990023`). LIGO-LW is a
legacy XML-based format used by a variety of LIGO/Virgo/KAGRA software and
services for storing tabular datasets.

Unfortunately, LIGO-LW is a rather complicated format. We recommend using
either the :mod:`ligo.lw` module or GWPy's :ref:`tabular LIGO-LW I/O
<gwpy-table-io-ligolw>` feature to simplify reading and writing LIGO-LW files.

.. note::
    There are two variants of the LIGO-LW format, an old format implemented by
    :mod:`glue.ligolw` that uses string ("ilwdchar") row IDs, and a new format
    implemented by :mod:`ligo.lw` that uses integer row IDs. GraceDB and
    BAYESTAR can accept *either* format, but pipelines SHOULD upload files in
    the new format.

    The :program:`ligolw_no_ilwdchar` command-line tool provided by
    :mod:`ligo.lw` can convert from the new format to the old format.

.. _event-data:

Event data
----------

This event data describes the search pipeline's matched filter output. It MUST
include the point estimates of the time, phase, and amplitude on arrival in
each detector. It MUST provide the intrinsic template parameters (masses and
spins). It SHOULD include a signal-to-noise time series for each detector.

The event data MUST include at least the following LIGO-LW tables (in any
order):

:class:`process <ligo.lw.lsctables.ProcessTable>`
    * The :class:`process <ligo.lw.lsctables.ProcessTable>` table MUST contain
      at least one row with the :attr:`~ligo.lw.lsctables.Process.process_id`
      and :attr:`~ligo.lw.lsctables.Process.program` columns populated in order
      to identify the search pipeline.

    * The value of those rows' :attr:`~ligo.lw.lsctables.Process.program`
      column MUST be one of ``pycbc``, ``gstlal_inspiral``,
      ``gstlal_inspiral_postcohspiir_online``, ``MBTAOnline``,
      ``bayestar_realize_coincs``, or ``bayestar-realize-coincs``.

    * Additional valid columns of this table MAY be populated in order to
      identify the pipeline software version or include other metadata.
      Additional unrelated rows (e.g. to identify prior analysis steps such as
      template bank generation) MAY be included and will be ignored.

:class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>`
    * The :class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>` table
      MUST contain exactly one row per detector that the search analyzed.

    * The values of the :attr:`~ligo.lw.lsctables.SnglInspiral.event_id` column
      MUST be distinct across all rows.

    * The values of the following columns that specify the intrinsic template
      parameters MUST be identical across all
      rows: :attr:`~ligo.lw.lsctables.SnglInspiral.mass1`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.mass2`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.f_final`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.spin1x`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.spin1y`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.spin1z`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.spin2x`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.spin2y`, and
      :attr:`~ligo.lw.lsctables.SnglInspiral.spin2z`.

    * If the template has zero spin, then the spin columns MAY be left blank.
      If the template has aligned spins, then the _x_ and _y_ spin components
      MAY be left blank.

    * The :attr:`~ligo.lw.lsctables.SnglInspiral.end_time` and
      :attr:`~ligo.lw.lsctables.SnglInspiral.end_time_ns` columns MUST report
      the seconds and nanoseconds parts of the GPS time at which the same
      fiducial reference part of the signal (e.g., the time of merger, or the
      time at which the inspiral reaches reference frequency) is received in
      each detector. It SHOULD record the merger time. If the event is an
      "early warning" or pre-merger event, then it SHOULD record the predicted
      time of merger.

    * If the event is an early warning event, then the high-frequency cutoff
      frequency MUST be recorded in the
      :attr:`~ligo.lw.lsctables.SnglInspiral.f_final` column.

    * The :attr:`~ligo.lw.lsctables.SnglInspiral.snr` column MUST report the
      absolute value of the complex matched filter SNR of the best-matching
      template. It MUST NOT report a modified SNR-like quantity such as newSNR.

    * The :attr:`~ligo.lw.lsctables.SnglInspiral.coa_phase` column MUST report
      the argument of the complex matched filter SNR of the best-matching
      template.

    * If the search pipeline as identified by the
      :attr:`~ligo.lw.lsctables.Process.program` column in the :class:`process
      <ligo.lw.lsctables.ProcessTable>` table is ``pycbc``, then phase
      convention of the :attr:`~ligo.lw.lsctables.SnglInspiral.coa_phase`
      column MUST be that the matched filter output is linear in terms of the
      data. Otherwise, the phase convention MUST be that the matched filter
      output is antilinear in terms of the data.

    * The :attr:`~ligo.lw.lsctables.SnglInspiral.end_time`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.end_time_ns`,
      :attr:`~ligo.lw.lsctables.SnglInspiral.snr`, and
      :attr:`~ligo.lw.lsctables.SnglInspiral.coa_phase` columns MAY be blank
      for any row for which there is a corresponding SNR time series (see
      below).

    * Due to a `bug in GraceDB`_, *all* columns of the
      :class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>` table
      (including blank ones) must be present.

:class:`coinc <ligo.lw.lsctables.CoincTable>`
    * There MUST be exactly one row in the
      :class:`coinc <ligo.lw.lsctables.CoincTable>` table with at least the
      :attr:`~ligo.lw.lsctables.Coinc.coinc_event_id` column populated.

    * The value of the :attr:`~ligo.lw.lsctables.Coinc.process_id` column of
      the :class:`coinc <ligo.lw.lsctables.CoincTable>` tale MUST match the
      value of the :attr:`~ligo.lw.lsctables.Process.process_id` column in
      the:class:`process <ligo.lw.lsctables.ProcessTable>` table that
      identifies the search pipeline.

    * Note that due to `another bug in GraceDB`_, the
      :attr:`~ligo.lw.lsctables.Coinc.time_slide_id` column MUST be populated.
      It MAY have a legal dummy value such as ``time_slide:time_slide_id:0``.

:class:`coinc_event_map <ligo.lw.lsctables.CoincMapTable>`
    * There MUST be exactly one row in the
      :class:`coinc_event_map <ligo.lw.lsctables.CoincMapTable>` table for each
      row in the :class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>`
      table.

    * The value in each row's :attr:`~ligo.lw.lsctables.CoincMap.event_id`
      column must be set to the value of the
      :attr:`~ligo.lw.lsctables.SnglInspiral.event_id` column in the
      corresponding row of the
      :class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>` table.

    * The value in each row's :attr:`~ligo.lw.lsctables.CoincMap.table_name`
      column must be set ``sngl_inspiral``.

    * Each row MUST have the :attr:`~ligo.lw.lsctables.CoincMap.coinc_event_id`
      column set to the value of the
      :attr:`~ligo.lw.lsctables.Coinc.coinc_event_id` column in the one row of
      the :class:`coinc <ligo.lw.lsctables.CoincTable>` table.

:class:`coinc_inspiral <ligo.lw.lsctables.CoincInspiralTable>`
    * The :class:`coinc_inspiral <ligo.lw.lsctables.CoincInspiralTable>` table
      MUST be present because it is required by GraceDB (although it is ignored
      by BAYESTAR).

    * It MUST have exactly one row.

    * The value in the :attr:`~ligo.lw.lsctables.CoincInspiral.coinc_event_id`
      column MUST match the value in the corresponding column in the
      :class:`coinc <ligo.lw.lsctables.CoincTable>` table.

    * The following columns MUST be populated:
      :attr:`~ligo.lw.lsctables.CoincInspiral.coinc_event_id`,
      :attr:`~ligo.lw.lsctables.CoincInspiral.combined_far`,
      :attr:`~ligo.lw.lsctables.CoincInspiral.end_time`,
      :attr:`~ligo.lw.lsctables.CoincInspiral.end_time_ns`,
      :attr:`~ligo.lw.lsctables.CoincInspiral.ifos`, and
      :attr:`~ligo.lw.lsctables.CoincInspiral.snr`.

    * The :attr:`~ligo.lw.lsctables.CoincInspiral.mass` and
      :attr:`~ligo.lw.lsctables.CoincInspiral.mchirp` columns SHOULD be
      populated.

The :file:`coinc.xml` file SHOULD also provide SNR time series for each
detector.

* Each SNR time series MUST be stored inside a :class:`~ligo.lw.ligolw.LIGO_LW`
  element as a serialized :class:`~lal.COMPLEX8TimeSeries`. The function
  :func:`lal.sereries.build_COMPLEX8TimeSeries` can be used to serialize a
  :class:`~lal.COMPLEX8TimeSeries`.

* Each of the :class:`~ligo.lw.ligolw.LIGO_LW` elements for serialized SNR time
  series MUST contain a :class:`~ligo.lw.ligolw.Param` element to link it to a
  row in the :class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>`. The
  param name MUST be ``event_id:param`` and the param's type and value must
  match the :attr:`~ligo.lw.lsctables.SnglInspiral.event_id` column in the
  corresponding :class:`sngl_inspiral <ligo.lw.lsctables.SnglInspiralTable>`
  row.

* The SNR time series MUST have an odd number of samples, e.g., the length must
  be :math:`2 * n + 1` for some integer :math:`n`.

* The timestamp of the central sample (e.g. :math:`n` times the sample interval
  plus the epoch) MUST differ from the corresponding :class:`sngl_inspiral
  <ligo.lw.lsctables.SnglInspiralTable>` row's time (if present) by no more
  than one sample interval.

* The timestamps of the samples of the SNR time series MUST correspond to
  sample boundaries. The timestamps MUST NOT have any sub-sample time shift
  applied to them.

* For any detector that lacks an SNR time series, sub-sample interpolation
  SHOULD be applied by the search pipeline to obtain the values for the
  :attr:`~ligo.lw.lsctables.SnglInspiral.snr`,
  :attr:`~ligo.lw.lsctables.SnglInspiral.coa_phase`,
  :attr:`~ligo.lw.lsctables.SnglInspiral.end_time`, and
  :attr:`~ligo.lw.lsctables.SnglInspiral.end_time_ns` columns in the
  corresponding row of the :class:`sngl_inspiral
  <ligo.lw.lsctables.SnglInspiralTable>` table.

.. _psd-data:

PSD data
--------

The PSD data consists of each analyzed detectors' estimated noise power
spectral density (PSD) series.

* There MUST be exactly one PSD per detector analyzed.

* Each PSD MUST be stored inside a :class:`~ligo.lw.ligolw.LIGO_LW`
  element as a serialized :class:`~lal.REAL8FrequencySeries`. The
  :func:`lal.sereries.build_COMPLEX8TimeSeries` function or the
  :func:`lal.sereries.make_psd_xmldoc` function can be used to serialize
  :class:`~lal.REAL8FrequencySeries`.

* Each :class:`~ligo.lw.ligolw.LIGO_LW` element MUST contain a
  :class:`~ligo.lw.ligolw.Param` element to link it to a detector. The param's
  name MUST be ``instrument:param``, its type MUST be ``instrument:param``, and
  its value should be a detector prefix such (e.g. one of ``H1``, ``L1``,
  ``V1``, ``K1``, ``I1``, etc.)

* Any samples that are invalid because their frequencies are outside of the
  range analyzed by the search MUST be absent or have their values set to
  positive infinity. Invalid values MUST NOT be set to zero.

Example files
-------------

For a minimal example, see the mock :download:`coinc.xml <_static/coinc.xml>`
file.

.. _`GraceDB`: https://gracedb.ligo.org
.. _`sequence diagram`: https://en.wikipedia.org/wiki/Sequence_diagram
.. _`bug in GraceDB`: https://git.ligo.org/lscsoft/gracedb/-/merge_requests/44
.. _`another bug in GraceDB`: https://git.ligo.org/lscsoft/gracedb/-/issues/197
