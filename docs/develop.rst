.. highlight:: sh

Development Instructions
========================

Nightly builds
--------------

If you want to try out the latest, bleeding-edge, unreleased version of
ligo.skymap, then you can install the most recent nightly build by running the
following command::

    $ pip install --upgrade 'ligo.skymap>=0.0.0.dev0'

Source dependencies
-------------------

If you are building `ligo.skymap` from source, then in addition to the
requirements in the :doc:`quick start section <quickstart/install>`, you will
also need:

*  `GSL`_ ≥ 1.15
*  `chealpix`_
   (note: if missing, will be built automatically from `bundled sources`_)
*  `pkg-config`_

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

.. _python-version-policy:

Python version support policy
-----------------------------

Generally, the latest release of `ligo.skymap` supports the Python versions
specified by :doc:`neps:nep-0029-deprecation_policy`. We follow the lead of the
main Scientific Python packages that we depend upon (Numpy, Scipy, Astropy) to
determine when to drop support for older versions of Python and take advantage
of new Python language features.

.. _`GSL`: https://www.gnu.org/software/gsl
.. _`chealpix`: https://sourceforge.net/projects/healpix/files/Healpix_3.30/
.. _`pkg-config`: https://www.freedesktop.org/wiki/Software/pkg-config/
.. _`GCC`: https://gcc.gnu.org
.. _`Clang`: https://clang.llvm.org
.. _`Intel C/C++ Compiler`: https://software.intel.com/content/www/us/en/develop/tools/compilers/c-compilers.html
.. _`bundled sources`: https://git.ligo.org/lscsoft/ligo.skymap/tree/main/cextern/chealpix
.. _`Instrumentation and Tracing Technology (ITT)`: https://software.intel.com/content/www/us/en/develop/documentation/vtune-help/top/api-support/instrumentation-and-tracing-technology-apis.html
.. _`VTune Profiler`: https://software.intel.com/content/www/us/en/develop/tools/vtune-profiler.html`
