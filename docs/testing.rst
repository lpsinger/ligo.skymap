.. highlight:: sh

Testing
=======

Unit tests
----------

This package has a unit test suite that is set up following the Astropy
:doc:`astropy:development/testguide`. The unit tests are run automatically on
every git push using `GitLab Continuous Integration (CI)`_. See the
repository's `.gitlab-ci.yml`_ file for configuration details.

You can also run the unit tests manually by running this command in the source
directory::

    $ python setup.py test

There are many options available to adjust what tests are run or how test
results are reported; see Astropy's documentation on
:ref:`astropy:running-tests` for details.

Coverage analysis
-----------------

The CI pipeline does code coverage analysis when it runs the unit tests. See the `coverage report`_ for the most recent build.

Acceptance tests
----------------

There is a suite of nightly `acceptance tests`_ for BAYESTAR that check that
the code reproduces localizations for past gravitational-wave events including
GW170814 and GW170817 as well as populations of simulated events.

.. _`GitLab Continuous Integration (CI)`: https://docs.gitlab.com/ee/ci/
.. _`.gitlab-ci.yml`: https://git.ligo.org/lscsoft/ligo.skymap/blob/master/.gitlab-ci.yml
.. _`coverage report`: https://lscsoft.docs.ligo.org/ligo.skymap/coverage.html
.. _`acceptance tests`: https://git.ligo.org/lscsoft/ligo-skymap-acceptance-tests
