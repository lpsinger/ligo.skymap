.. highlight:: sh

Testing
=======

Unit tests
----------

This package has a unit test suite that is set up following the `Astropy
Testing Guidelines`_. The unit tests are run automatically on every git push
using `GitLab Continuous Integration (CI)`_. See the repository's
`.gitlab-ci.yml`_ file for configuration details.

You can also run the unit tests manually by running these commands in the
source directory::

    $ pip install -e .[test]
    $ pytest

There are many options available to adjust what tests are run or how test
results are reported; see `Astropy's documentation on running-tests`_ for
details.

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
.. _`Astropy's documentation on running-tests`: https://docs.astropy.org/en/latest/development/testguide.html#running-tests
.. _`coverage report`: https://lscsoft.docs.ligo.org/ligo.skymap/coverage.html
.. _`acceptance tests`: https://git.ligo.org/lscsoft/ligo-skymap-acceptance-tests
