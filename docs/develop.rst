.. highlight:: sh

Development Instructions
========================

Source dependencies
-------------------

If you are building `ligo.skymap` from source, then in addition to the
requirements in the :doc:`quick start section <quickstart/install>`, you will
also need:

*  `GSL <https://www.gnu.org/software/gsl>`_ ≥ 1.15
*  `chealpix <https://sourceforge.net/projects/healpix/files/Healpix_3.30/>`_
   (note: if missing, will be built automatically from `bundled sources
   <https://git.ligo.org/lscsoft/ligo.skymap/tree/master/cextern/chealpix>`_)

Building from source
--------------------

To build `ligo.skymap` from source, first clone the git repository::

    $ git clone https://git.ligo.org/lscsoft/ligo.skymap.git

Then install it with pip::

    $ pip install .

Environment variables that control the build
--------------------------------------------

There are several environment variables that control the build. To activate one
of these options, set it to any non-empty value when you run ``pip install``,
like this::

    $ env LIGO_SKYMAP_USE_SYSTEM_CHEALPIX=1 pip install .

Here is the full list of environment variables.

:envvar:`LIGO_SKYMAP_USE_SYSTEM_CHEALPIX`
    Use the system installation of chealpix rather than building chealpix from
    the bundled source code.

:envvar:`LIGO_SKYMAP_USE_ITTNOTIFY`
    Compile and link against the Intel® `Instrumentation and Tracing Technology
    (ITT)`_ API to add tracepoints for performance measurement using Intel®
    `VTune Profiler`_.

:envvar:`LIGO_SKYMAP_DISABLE_OPENMP`
    Disable OpenMP parallelization.

.. _`Instrumentation and Tracing Technology (ITT)`: https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/api-support/instrumentation-and-tracing-technology-apis.html
.. _`VTune Profiler`: https://software.intel.com/content/www/us/en/develop/tools/vtune-profiler.html`
