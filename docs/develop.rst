.. highlight:: sh

Development Instructions
========================

Source dependencies
-------------------

If you are building `ligo.skymap` from source, then in addition to the
requirements in the :doc:`quick start section <quickstart/install>`, you will
also need:

*  `GSL`_ ≥ 1.15
*  `chealpix`_
   (note: if missing, will be built automatically from `bundled sources`_)

You also need a C compiler with good support for the C11 standard. The
following compilers are known to work:

*  `GCC`_ ≥ 5 (≥ 8 recommended)
*  `Clang`_ ≥ 5.0
*  `Intel C/C++ Compiler`_

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

.. _`GSL`: https://www.gnu.org/software/gsl
.. _`chealpix`: https://sourceforge.net/projects/healpix/files/Healpix_3.30/
.. _`GCC`: https://gcc.gnu.org
.. _`Clang`: https://clang.llvm.org
.. _`Intel C/C++ Compiler`: https://software.intel.com/content/www/us/en/develop/tools/compilers/c-compilers.html
.. _`bundled sources`: https://git.ligo.org/lscsoft/ligo.skymap/tree/main/cextern/chealpix
.. _`Instrumentation and Tracing Technology (ITT)`: https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/api-support/instrumentation-and-tracing-technology-apis.html
.. _`VTune Profiler`: https://software.intel.com/content/www/us/en/develop/tools/vtune-profiler.html`
