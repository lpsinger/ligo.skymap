.. highlight:: sh

Testing
=======

Unit tests
----------

This package has a unit test suite that is set up based on the `Astropy
Testing Guidelines`_. The unit tests are run automatically on every git push
using `GitLab Continuous Integration (CI)`_. See the repository's
`.gitlab-ci.yml`_ file for configuration details.

You can also run the unit tests manually by running these commands in the
source directory::

    $ pip install -e .[test]
    $ pytest

There are many options available to adjust what tests are run or how test
results are reported; see `Astropy's documentation on running tests`_ for
details.

Unit tests in ligo.skymap come in three forms:

-  **Test modules**: Generally, each Python subpackage of ligo.skymap has a
   `tests` directory containing :doc:`Pytest <pytest:index>`-style unit tests.
   For example, unit tests for :mod:`ligo.skymap.plot` are in
   `ligo/skymap/plot/tests`_.
-  **Doctests**: The documentation for many modules contains inline examples
   with expected console output called
   :doc:`doctests <python:library/doctest>`. When you run the test suite, it
   checks that the output matches. An example of a module with many doctests is
   :mod:`ligo.skymap.distance` (see `source`_).
-  **C unit tests**: The critical sections of C code have unit tests that are
   written in C using :doc:`GSL <gsl:index>`'s test framework (see
   `gsl_test.h`_). These tests are defined in the functions
   :c:func:`bayestar_test` in `src/bayestar_sky_map.c`_ and
   :c:func:`cubic_interp_test` in `src/cubic_interp_test.c`_.

Coverage analysis
-----------------

The CI pipeline does code coverage analysis when it runs the unit tests. See
the `coverage report`_ for the most recent build.

Acceptance tests
----------------

There is a suite of weekly `acceptance tests`_ for BAYESTAR that check that
the code reproduces localizations for past gravitational-wave events including
GW170814 and GW170817 as well as populations of simulated events.

.. _`Astropy Testing Guidelines`: https://docs.astropy.org/en/latest/development/testguide.html
.. _`GitLab Continuous Integration (CI)`: https://docs.gitlab.com/ee/ci/
.. _`.gitlab-ci.yml`: https://git.ligo.org/lscsoft/ligo.skymap/blob/main/.gitlab-ci.yml
.. _`Astropy's documentation on running tests`: https://docs.astropy.org/en/latest/development/testguide.html#running-tests
.. _`ligo/skymap/plot/tests`: https://git.ligo.org/lscsoft/ligo.skymap/-/blob/main/ligo/skymap/plot/tests
.. _`source`: https://git.ligo.org/lscsoft/ligo.skymap/-/blob/main/ligo/skymap/distance.py
.. _`gsl_test.h`: https://git.savannah.gnu.org/cgit/gsl.git/tree/test/gsl_test.h
.. _`src/bayestar_sky_map.c`: https://git.ligo.org/lscsoft/ligo.skymap/-/blob/main/src/bayestar_sky_map.c
.. _`src/cubic_interp_test.c`: https://git.ligo.org/lscsoft/ligo.skymap/-/blob/main/src/cubic_interp_test.c
.. _`coverage report`: https://lscsoft.docs.ligo.org/ligo.skymap/coverage.html
.. _`acceptance tests`: https://git.ligo.org/leo-singer/ligo-skymap-acceptance-tests-public
