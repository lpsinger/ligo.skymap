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

Review
------

This package comprises an analysis pipeline of the Compact Binary Coalescence
(CBC) Working Group of the `LIGO Scientific Collaboration (LSC)`_. As such, it
is required to undergo an LSC `scientific review (private wiki link)`_. Minutes
of ligo.skymap's LSC review team are kept in the project's `private wiki`_.

Scope
~~~~~

The essential science capabilities that undergo LSC review are:

*   Rapid CBC localization (:doc:`BAYESTAR <bayestar/index>`)
*   Post-processing of posterior samples to 3D sky maps
    (:doc:`ligo-skymap-from-samples <tool/ligo_skymap_from_samples>`)
*   High-level 2D and 3D localization scripts
    (:doc:`ligo-skymap-plot <tool/ligo_skymap_plot>`,
    :doc:`ligo-skymap-plot-volume <tool/ligo_skymap_plot_volume>`)

Roles
~~~~~

All review tests (the unit tests and acceptance tests described above) should
be automated, take no more than 2 hours to run, and indicate success or failure
in a self-evident way. The task of the developers is to create and maintain the
automated tests. The task of the reviewers is to provide oversight to ensure
that the tests are necessary and sufficient to cover the scientific
functionality that is under review.

Changes
~~~~~~~

ligo.skymap is a stable, mature package. Most changes are conservative and
maintenance-oriented. The developers will not usually contact the reviewers
about these kinds of changes before merging into the main branch. Examples of
these kinds of change are adjustments to track API changes in Python, Numpy,
and Astropy.

However, the developers will flag potential changes that might need extra
scrutiny because they could science results by adding the `Review label`_ to
merge requests and requiring an approval from a reviewer before merging.
Examples of these changes are extracting new parameters from BAYESTAR (e.g.
inclination angles) or making significant changes to algorithm inner loops that
could affect floating point accuracy.

Releases
~~~~~~~~

All stable releases (versions that are triples of numbers of the form
``1.2.3``) of ligo.skymap have been approved by the review team. We create one
or more release candidates (versions that are of the form ``1.2.3rcN`` for some
number ``N``) until the latest release candidate satisfies all of the tests and
is verbally approved by the review team. Then we do a stable release. The
review team indicates its formal assent to the release by approving the
corresponding ticket in the LSC
`Software Change Control Board (SCCB) issue tracker (private link)`_.

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
.. _`LIGO Scientific Collaboration (LSC)`: https://www.ligo.org
.. _`scientific review (private wiki link)`: https://git.ligo.org/cbc-review/review/-/wikis/CBC-review-guidelines
.. _`private wiki`: https://git.ligo.org/lscsoft/ligo.skymap/-/wikis/home
.. _`Review label`: https://git.ligo.org/lscsoft/ligo.skymap/-/merge_requests?label_name%5B%5D=Review
.. _`Software Change Control Board (SCCB) issue tracker (private link)`: https://git.ligo.org/computing/sccb/-/issues
