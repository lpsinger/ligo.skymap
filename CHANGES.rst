#########
Changelog
#########

2.1.3 (unreleased)
==================

- No changes yet.

2.1.2 (2024-10-22)
==================

- Add sign convention for `gstlal_inspiral_coinc_extractor`, the program that
  is used to upload offline GstLAL events to GraceDB.

2.1.1 (2024-10-03)
==================

- Fix blank volume rendering plots for very well-localized events due to
  `changes in data type promotion rules in Numpy 2.x`__.

  __ https://numpy.org/devdocs/numpy_2_0_migration_guide.html#changes-to-numpy-data-type-promotion

- Restore support for Numpy 1.x while keeping support for Numpy 2.x.

2.1.0 (2024-08-23)
==================

- BAYESTAR now terminates gracefully if you interrupt it by typing control-C or
  sending the process in which it is running a SIGINT signal. Previously, there
  was a bug in the signal handling code that caused it to segfault on keyboard
  interrupt.

  Note that if you invoke BAYESTAR multiple times in different threads of the
  same process, then only one of the invocations will stop early due the
  interrupt, because signal handlers are process-wide.

- Enhancements to the command-line tool ``ligo-skymap-plot``:

  - Read in and reproject sky maps in multi-resolution format. This
    significantly decreases memory consumption when plotting high-resolution
    sky maps.

  - The algorithm that calculates credible areas printed in the plot legend
    now matches the algorithm in ``ligo.skymap.postprocess.crossmatch``: it
    does linear interpolation of area when the credible level falls between
    pixels.

- The ``ligo-skymap-plot-volume`` script now uses multi-resolution sky maps
  internally. This will significantly decrease its memory usage for
  high-resolution sky maps.

- Use a Wayback Machine URL to download the Mellinger sky panorama because the
  original URL is broken.

- Fixes related to the ``cut_dateline`` and ``cut_prime_meridian`` functions:

  - Adjust to API changes between Shapely 1.x and 2.x.

  - Add a dependency on Shapely.

  - Add test cases.

- Require numpy ≥ 2.0.0 due to changes in string representations of scalars
  that affected doctest test cases. See `astropy/astropy#15095`__.

  __ https://github.com/astropy/astropy/issues/15095

- There are now nightly builds of ligo.skymap. For instructions on installing
  the latest, unreleased, bleeding-edge version of ligo.skyamp, see
  https://lscsoft.docs.ligo.org/ligo.skymap/develop.html#nightly-builds.

2.0.1 (2024-05-30)
==================

- Drop support for Python 3.9 because Astropy 6.1.0 dropped Python 3.9 in
  accordance with the Numpy version support policy.

- Require astropy >= 6.0.0.

- Fix a memory alignment bug that caused ``ligo.skymap.moc.rasterize``
  and ``ligo-skymap-flatten`` to produce incorrect results for multiresolution
  sky maps that have an odd number of 32-bit columns. This impacted Swift BAT
  GUANO localizations which have a single 32-bit PROBDENSITY column.

- Ensure that ``ligo.skymap.distance.conditional_cdf`` and
  ``ligo.skymap.distance.marginal_cdf`` correctly handle corner cases for
  extreme values: both now return 0 for distances that are less than or equal
  to 0, and 1 for a distance of positive infinity (assuming that the
  ``distnorm`` argument does indeed normalize the distribution).

- Ensure that all ligo.skymap command-line tools close any files that they have
  opened. These command-line tools are used as functions in GWCelery, and were
  leaking unclosed file descriptors.

2.0.0 (2024-04-15)
==================

- Added options to center ``mollweide`` and ``aitoff`` projections. Thanks go
  to Sam Wyatt for this contribution.

- Added support for ``os.PathLike`` filenames when reading ligolw files. Thanks
  go to Thomas Sainrat for this contribution.

- Check for more invalid input corner cases in ``ligo.skymap.moc.rasterize``.

- Remove the ``--min-inclination`` and ``max-inclination`` options from
  ``bayestar-localize-coincs`` and ``bayestar-localize-lvalert``. These options
  are rarely used, and will be made obsolete by a future release that adds
  inclination posteriors to sky map output.

- Remove ``frameon=False`` in ``ligo-skymap-plot-volume`` so that it respects
  the (lack of the) ``--transparent`` option. This improves text and label 
  readability against dark backgrounds when transparent mode is not on. Thanks
  go to Geoffrey Mo for this contribution.

- Add documentation on the LIGO Scientific Collaboration (LSC) review process
  to the Testing section of the manual.

- Require Numpy >= 1.23.0. Rebuild for binary compatibility with Numpy 2.0.0.

- Add unit tests for Python 3.12.

1.1.2 (2023-10-03)
==================

- Update for compatibility with Matplotlib 3.8.0.

- Binary wheels for macOS x86_64 are now built against Big Sur (10.15), because
  Catalina (11) is past end of life.

- Fix deprecation warnings from importlib.resources.

1.1.1 (2023-07-08)
==================

- Fix a typo in the ``setup.cfg`` file that prevented correct interpretation of
  the minimum Python version. Contributed by
  `@ahnitz <https://github.com/ahnitz>`_.

1.1.0 (2023-07-07)
==================

- Add ``max_depth`` keyword argument to the call to
  ``MOC.from_valued_healpix_cells`` in ``ligo-skymap-contour-moc``.
  Contributed by `@parkma99 <https://github.com/parkma99>`_.

- Improve handling of the ``--output`` command line argument for
  ``ligo-skymap-contour-moc``:

  - Add ``-o`` as a short form.

  - Don't write to stdout by default; it does not make sense to write a binary
    FITS file to stdout.

  - Make the argument required.

- Drop dependency on distutils to prepare for its removal in Python 3.12.
  See `PEP 632 <https://peps.python.org/pep-0632/>`_.

- Drop support for Python 3.8.

- Vectorize ``find_ellipse`` over the ``cl`` argument.

- Tune compiler settings used to build wheels for PyPI:

  - Add the option ``-fvisibility=hidden`` to hide all symbols except for the
    Python entry point. This improves the efficiency of link-time optimization.
    On average, it speeds up BAYESTAR by about 5%.

  - Add the options ``-Ofast -fno-finite-math-only -flto`` on Linux aarch64
    and macOS, the targets on which we use gcc. These options approximate the
    configuration that we use for icc on Linux x86_64. On average, this change
    speeds up BAYESTAR on macOS by about 30%.

- Factor out the Python implementation of the BAYESTAR adaptive mesh refinement
  algorithm so that other libraries can use it. It is exposed as
  ``ligo.skymap.moc.bayestar_adaptive_grid``.

- Fix incorrectly rendered default values for some command line arguments in
  the documentation.

- Move coherence plots from GWCelery to ligo.skymap.

1.0.7 (2023-02-27)
==================

- Track an API change in Matplotlib 3.7.0. Update test baseline images.

- Update Linux wheels from manylinux2014 to manylinux_2_28.

- Require scipy ≠ 1.10.0 due to an unplanned API change in that version, which
  was fixed in 1.10.1.

- Add unit tests under Python 3.11 to the continuous integration pipeline.

1.0.6 (2023-02-03)
==================

- Fix an issue with OpenMP and Python multiprocessing that caused
  ``ligo-skymap-stats`` to parallelize inefficiently on Linux.

1.0.5 (2023-01-31)
==================

- Require scipy < 1.10.0 due to removal of ``multivariate_normal_gen.cov`` in
  that version. A future version of Scipy may add it back as a property; see
  `scipy/scipy#17896`__.

  __ https://github.com/scipy/scipy/issues/17896

1.0.4 (2022-12-06)
==================

- Change the default value of the ``origin`` card in FITS files generated by
  BAYESTAR and ``ligo-skymap-from-samples`` from ``LIGO/Virgo`` to
  ``LIGO/Virgo/KAGRA``.

- Build binary wheels for the aarch64 (Arm64) architecture on Linux.

1.0.3 (2022-10-11)
==================

- Update condor accounting group in ``bayestar-localize-coincs`` and
  ``bayestar-mcmc`` to ``ligo.dev.o4.cbc.pe.bayestar``.

- Track `pending deprecation of matplotlib.cm.register_cmap`__.
  Require matplotlib >= 3.5.0.

  __ https://matplotlib.org/stable/api/prev_api_changes/api_changes_3.6.0.html#pending-deprecation-top-level-cmap-registration-and-access-functions-in-mpl-cm

- The function ``ligo.skymap.postprocess.ellipse.find_ellipse`` will now return
  a tuple of the same length in all circumstances, even under error conditions.

1.0.2 (2022-08-18)
==================

- Add fast path for PowerPC and other architectures in ``uniq2order``.

1.0.1 (2022-08-17)
==================

- Replace deprecated
  ``astropy.cosmology.default_cosmology.get_cosmology_from_string``.

- Build wheels for arm64 on macOS.

- Add fast path for arm64 in ``uniq2order``.

1.0.0 (2022-06-01)
==================

- Run unit tests under Python 3.10.

- Update the `BAYESTAR interface definition document`_ to state that online CBC
  pipelines should now include their PSD files in the initial ``coinc.xml``
  upload, and should not upload a separate ``psd.xml.gz`` file.

  ``bayestar-localize-lvalert`` will now download ``psd.xml.gz`` (and log a
  warning) only if the PSD was not present in the ``coinc.xml`` file.

  .. _`BAYESTAR interface definition document`: https://lscsoft.docs.ligo.org/ligo.skymap/interface.html

- Several enhancements and bug fixes in ``bayestar-inject``:

  - Swap component masses if necessary so that mass1 >= mass2 always.

  - Rename the ``--min-snr`` option to ``--snr-threshold`` for consistency with
    the same option for ``bayestar-realize-coincs`. The old ``--min-snr``
    spelling is deprecated and will be
    removed in a future release.

  - Add the ``--min-triggers`` option to ``bayestar-inject`` to control the
    minimum number of triggers to form a coincidence, for consistency with
    ``bayestar-realize-coincs``.

  - Add the ``--distribution-samples`` option to load samples for the intrinsic
    mass and spin distribution from an external file.

- Linux wheels are now built against cfitsio 4.1.0. See
  https://github.com/lpsinger/ligo.skymap/issues/12.

- Add the ``request_disk`` flag when submitting ``bayestar-localize-coincs``
  jobs to HTCondor. This is now required on LIGO Data Grid clusters.

- Fix compatibility with Astropy 5.1.

0.6.1 (2022-01-18)
==================

- Skip Numpy 1.22.0 because of an issue with Astropy table aggregation.
  See `astropy#12706`_.

  .. _`astropy#12706`: https://github.com/astropy/astropy/issues/12706

- Skip lalsuite 7.2 due to an upstream regression. See `lalsuite!1757`_.

  .. _`lalsuite!1757`: https://git.ligo.org/lscsoft/lalsuite/-/merge_requests/1757

- Work around a regression in Numpy 1.22.0 that broke building third party
  packages using the limited Python C API. See `numpy#20818`_.

  .. _`numpy#20818`: https://github.com/numpy/numpy/pull/20818

- Update to python-ligo-lw >= 1.8.0.

0.6.0 (2021-12-01)
==================

- Rename ``master`` branch to ``main``.

- Add a ``max-distance`` option to ``bayestar-inject``.

- Increase verbosity of LAL error reporting so that the user gets more
  information for invalid waveform arguments.

- Wheels for macOS are now built against macOS 10.15 (Catalina) using GCC 11.

- Require Python >= 3.8 due Astropy and Numpy deprecation policy.
  See `APE 18`_ and `NEP 29`_.

  .. _`APE 18`: https://github.com/astropy/astropy-APEs/blob/main/APE18.rst
  .. _`NEP 29`: https://numpy.org/neps/nep-0029-deprecation_policy.html

- In ``bayestar_inject``, use the method ``vectorize_redshift_method`` instead
  of ``vectorize_if_needed`` from ``astropy.cosmology.utils``, because the
  latter was deprecated in Astropy 5.0 (see `astropy#12176`_).

  .. _`astropy#12176`: https://github.com/astropy/astropy/pull/12176

- Require astropy >= 5.0.

- Require python-ligo-lw <= 1.7.1 because of an API breaking change that will
  occur in the next version of python-ligo-lw. Support for new versions of
  python-ligo-lw will be added in an upcoming release. See `ligo.skymap#30`_.

  .. _`ligo.skymap#30`: https://git.ligo.org/lscsoft/ligo.skymap/-/issues/30

- Add support for all-sky projections in Galactic coordinates activated by
  creating Matplotlib axes with the keyword arguments like
  ``projection='galactic degrees mollweide'``.

- Add the ``mark_inset_circle`` and ``connect_inset_circle`` methods to
  ``AutoScaledWCSAxes`` in order to support circular insets (loupes).

- Determine input filetypes by reading the file header in Python rather than
  relying on a shell utility.

0.5.3 (2021-04-10)
==================

- Word-wrap the Python and command line arguments that are recorded in the
  ``HISTORY`` cards. This makes the arguments more legible, because Astropy's
  built-in FITS card wrapping behavior does not consider word breaks. It also
  works around a FITS validation regression in Astropy 4.2.1
  (see `astropy#11486`_).

  .. _`astropy#11486`: https://github.com/astropy/astropy/issues/11486

0.5.2 (2021-03-28)
==================

- Teach the ``astro zoom`` and ``astro globe`` projections to accept sky
  coordinates in any Astropy representation, including Cartesian coordinates.

- Enable SNR time series by default in ``bayestar-realize-coincs``.

- Update the required version of Matplotlib to >= 3.4.0, since it includes the
  bug fix for `matplotlib#18832`_.

- Update the required version of Astropy to >= 4.0.2 and != 4.2. Astropy 4.1
  now works with Matplotlib >= 3.4.0, but Astropy 4.2 introduced a bug
  affecting Numpy and sky coordinates that will be fixed in Astropy 4.2.1
  (see `astropy#11133`_).

  .. _`astropy#11133`: https://github.com/astropy/astropy/pull/11133

0.5.1 (2021-02-27)
==================

- This is the first release of ligo.skymap that is tested under and officially
  supports Python 3.9. (We were mostly waiting for LALSuite to be built for
  Python 3.9).

- Drop support for Python 3.6 because it is no longer supported by many other
  scientific Python packages like Matplotlib and Numpy.

- Update the required version of Astropy to >= 4.0.2 and < 4.1. Astropy 4.0.2
  includes a bug fix for cache handling on cluster filesystems (see
  `astropy#9970`_). Astropy 4.1 caused some issues with Matplotlib projection
  classes as a result of changes in hashing behavior of
  ``astropy.coordinates.SkyCoord`` (see `matplotlib#18832`_), which should be
  fixed in Matplotlib 3.4.0.

  .. _`astropy#9970`: https://github.com/astropy/astropy/issues/9970
  .. _`matplotlib#18832`: https://github.com/matplotlib/matplotlib/issues/18832

- Update the required version of LALSuite to >= 6.82 to work around an
  incompatibility between Numpy >= 1.20.0 and older versions of LALSuite
  (see `lalsuite#414`_).

  .. _`lalsuite#414`: https://git.ligo.org/lscsoft/lalsuite/-/issues/414

- Importing ligo.skymap no longer causes the
  ``astropy.coordinates.EarthLocation`` site registry to be populated with the
  locations of gravitational-wave observatories, because these sites are now
  included in Astropy's own data repository (see `astropy-data#89`_).

  .. _`astropy-data#89`: https://github.com/astropy/astropy-data/pull/89

- In the command line help for ``bayestar-localize-coincs`` and in the
  ``COMMENT`` card in the output FITS file, explain that the integer value in
  the ``OBJECT`` card in the FITS header is a row ID that refers to a
  coinc_event table row in the input LIGO-LW document.

- Add the ``--rescale-loglikelihood`` command line argument to expose
  BAYESTAR's log likelihood factor that accounts for excess technical sources
  of noise from the matched filter pipeline.

0.5.0 (2020-08-27)
==================

- Add ``--f-high`` option to ``bayestar-realize-coincs`` in order to simulate
  early warning triggers.

- In sky maps produced by ``bayestar-localize-coincs``, the FITS headers now
  contain ``OBJECT`` identifiers that are integer event IDs (such as ``1``)
  rather than strings (such as ``coinc_event:coinc_event_id:1``).

- The ``ligo-skymap-stats`` tool now recognizes FITS headers with either
  integer or string ``OBJECT`` identifiers.

- Use Astropy rather than LAL for GPS to UTC time conversion in FITS headers so
  that LALSuite is not a requirement for reading and writing FITS files.

- Refactor ``ligo-skymap-stats`` to unify its multiprocessing and progress bar
  implementation with other command line tools.

- Update the compiler version that is used to build Linux wheels to icc
  19.1.2.254 from Intel Parallel Studio XE 2020u2.

- Port the Python C extension to the limited stable Python API so that one
  binary wheel works for all supported Python versions for any given operating
  system. See `PEP 384 <https://www.python.org/dev/peps/pep-0384/>`_.

- Eliminate global static variables from the Python C extension to enable
  compatibility with Python subinterpreters. See
  `PEP 3121 <https://www.python.org/dev/peps/pep-3121/>`_.

- Improve the numerical stability of the method
  :meth:`ligo.skymap.distance.conditional_ppf` by reparametrizing the equation
  that is being solved. This method, which calculates the inverse of the
  distance CDF, works by solving the equation :math:`f(x) - p = 0` for
  :math:`x`, where :math:`f(x)` is the distance CDF, and :math:`p` is the
  desired probability.

  The reparametrized equation is :math:`log(1 - f(x)) - log(1 - p) = 0` if
  :math:`p > 1/2` and :math:`log(f(x)) - log(p) = 0` otherwise. This
  reparametrization is effective because it improves the dynamic range in the
  tails of the distribution. This same reparametrization had already proven
  effective in the related method :meth:`ligo.skymap.distance.marginal_ppf`.

  This change also fixes some rare corner cases where
  :meth:`~ligo.skymap.distance.marginal_ppf` returned silly values becauses it
  uses :meth:`~ligo.skymap.distance.conditional_ppf` internally to create its
  own initial guess. One example was the median distance for the binary neutron
  star candidate S191205ah. Before this patch, the result was negative and
  invalid::

      >>> from ligo.skymap.distance import marginal_ppf
      >>> from ligo.skymap.moc import uniq2pixarea
      >>> from ligo.skymap.io import read_sky_map
      >>> url = 'https://gracedb.ligo.org/apiweb/superevents/S191205ah/files/bayestar.multiorder.fits'
      >>> s = read_sky_map(url, moc=True)
      >>> marginal_ppf(0.5, s['PROBDENSITY'] * uniq2pixarea(s['UNIQ']),
      ...              s['DISTMU'], s['DISTSIGMA'], s['DISTNORM'])
      /Users/lpsinger/src/ligo.skymap/ligo/skymap/util/numpy.py:46: RuntimeWarning: invalid value encountered in marginal_ppf
        return func(*args, **kwargs)
      -223357.8508233767

  After this patch, the result is positive and sensible::

      >>> marginal_ppf(0.5, s['PROBDENSITY'] * uniq2pixarea(s['UNIQ']),
      ...              s['DISTMU'], s['DISTSIGMA'], s['DISTNORM'])
      362.7485740018039

- Increase the range of validity of the solver used in
  :meth:`ligo.skymap.distance.moments_to_parameters` for low-probability pixels
  that are very prior dominated. Sky maps that have many such pixels could have
  credible volumes repoted as infinity. The incidence of such cases should now
  be decreased.

- Correct the alignment of Numpy record arrays passed to
  :func:`ligo.skymap.moc.rasterize` in order to avoid possibly undefined
  behavior that was detected by UBSan.

0.4.0 (2020-07-26)
==================

- Normalize column names when an ASCII file is passed to
  ``ligo-skymap-from-samples``.

- Migrate LIGO-LW XML support from the ``glue.ligolw`` module to the newer and
  better maintained ``ligo.lw`` module.

- Teach BAYESTAR to accept either string row IDs (such as
  ``sngl_inspiral:event_id:1``) or integer row IDs (such as ``1``).

- The parallel ``map()`` implementation that is used by a number of the
  package's command line tools will now yield results in order as quickly as
  they arrive, rather than sorting all of the results at the end. This should
  provide a very modest speedup in some command line tools.

0.3.1 (2020-05-28)
==================

- Replace a call to the ``aligned_alloc`` function with the ``posix_memalign``
  function. The ``aligned_alloc`` function is part of the C11 standard library,
  but is missing on some platforms, particularly very old versions of macOS.

  This fixes an issue with building Conda packages.

0.3.0 (2020-05-26)
==================

- Fix an out of bounds access in the bicubic interpolation function that
  BAYESTAR uses to evaluate the integral over distance. Due to the relationship
  between the lookup table bounds and BAYESTAR's distance limits of
  integration, the corner case that caused out of bounds access was never
  triggered. This bug had no impact on localizations generated by BAYESTAR.

- More performance improvements in BAYESTAR providing a 2x speedup.
  For benchmark results, see the new `How fast is BAYESTAR?`_ section in the
  manual.

  - The function ``bicubic_interp_eval`` had not being effectively
    autovectorized by the compiler. Rewrite it in explicitly vector form using
    the `GCC vector extension`_ (which is also supported by clang and icc) and
    selected vector intrinsics. In x86_64 builds, gcc, clang, and icc will now
    emit SSE2, SSE4.1, and FMA instructions for this code.

  - Pre-evaluate the SNR=0 limit of the distance integral to move some
    conditionals and logarithms out of BAYESTAR's innermost loop.

  - Add loop count hints to improve the efficacy of loop unrolling.

  - Perform manual loop fission in ``bayestar_sky_map_toa_phoa_snr_pixel``.

- Update ligo.skymap to the latest version of the Astropy affiliated package
  template. Migrate package infrastructure from `APE 4`_ to `APE 17`_. The
  astropy-helpers submodule has been removed, and the package now includes a
  pyproject.toml file (see `PEP 517`_ and `PEP 518`_).

- As a consequence of migrating to `APE 17`_ and switching to
  `setuptools_scm`_, the version of ligo.skymap will be reported slightly
  differently. The ``ligo.skymap.__githash__`` variable has been removed, and
  instead the git hash will be part of the ``ligo.skymap.__version__`` version
  string for unreleased, local versions.

- Correspondingly, ``ligo.skymap`` tools that generate FITS files
  (``bayestar-localize-lvalert``, ``bayestar-localize-coincs``,
  ``ligo-skymap-from-samples``) will no longer populate the ``VCSREV`` and
  ``DATE-BLD`` keys in FITS headers.

  .. _`GCC vector extension`: https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html
  .. _`How fast is BAYESTAR?`: https://lscsoft.docs.ligo.org/ligo.skymap/performance.html
  .. _`APE 4`: https://github.com/astropy/astropy-APEs/blob/master/APE4.rst
  .. _`APE 17`: https://github.com/astropy/astropy-APEs/blob/master/APE17.rst
  .. _`PEP 517`: https://www.python.org/dev/peps/pep-0517/
  .. _`PEP 518`: https://www.python.org/dev/peps/pep-0518/
  .. _`setuptools_scm`: https://github.com/pypa/setuptools_scm

0.2.2 (2020-05-12)
==================

- Fix incorrect legends on histograms generated by ``ligo-skymap-plot-stats``.

- When the ``bayestar-localize-coincs`` or ``bayestar-localize-lvalert``
  scripts are called with ``--loglevel=info`` or higher, they will now output
  additional runtime measurements. Specifically, they will output the "real"
  time (wall clock time), "user" time (total time spent in userland across all
  threads), and "sys" time (total time spent in kernel land across all
  threads), similar to the UNIX :manpage:`time(1)` tool. Here is an example of
  the formatting::

      2020-05-12 18:57:12,024 INFO finished computationally-intensive section in real=0.918s, user=36.339s, sys=0.293s

0.2.1 (2020-05-04)
==================

- Speed up ``import ligo.skymap`` by up to a second by replacing uses of
  ``pkg_resources`` with the new Python standard library module
  ``importlib.resources`` (or, for Python < 3.7, the backport
  ``importlib_resources``). The old ``pkg_resources`` module is known to be
  slow because it does a lot of work on startup. (See, for example,
  https://github.com/pypa/setuptools/issues/926 and
  https://github.com/pypa/setuptools/issues/510.)

- Drop dependency on seaborn.

- Move some rarely used imports (``networkx`` and ``astropy.convolution``) from
  module scope to function scope to speed up imports by up to half a second on
  NFS filesystems.

0.2.0 (2020-04-21)
==================

- Update installation instructions to state that installation with pip requires
  pip 19.3 or newer. This has been the case since ligo.skymap 0.1.16.

- Teach BAYESTAR to respect the ``f_final`` column in the ``sngl_inspiral``
  table for pre-merger, early warning templates.

- Ensure that BAYESTAR's arrival time prior is long enough to contain at least
  half a cycle of the template autocorrelation sequence. Previously, the
  duration of the arrival time prior was calculated solely from the light
  travel times between the participating detectors. This fixes an issue where
  SNR time series for early-warning events could have been cropped to only 1-3
  samples.

- Change BAYESTAR's strategy for evaluating SNR time series from Catmull-Rom
  interpolation of the real and imaginary parts to Catmull-Rom interpolation of
  the amplitude and phase. The old interpolation method could produce
  oscillatory artifacts in the SNR amplitude if the data are nearly critically
  sampled, as is the case for early-warning BNS events. The new interpolation
  method is immune to this kind of artifact, and also has much faster
  convergence as a function of sample rate.

- Lift the code to apply time shifts to SNR series outside of BAYESTAR's inner
  loop because there are no data dependencies on the variables of integration.
  This is seen to speed up BAYESTAR by 30%.

- Add software version and command line arguments metadata to the output of
  ``ligo-skymap-plot-stats``.

- Fix a bug in the Lanczos sub-sample arrival time interpolant: the Lanczos
  kernel should be zero for ``abs(t) >= a``.

- Remove ``requirements.txt`` file and list dependencies in ``setup.cfg``
  instead.

- The ``bayestar-localize-coincs`` will no longer create HTCondor user log
  files because the large number of open log files could strain the filesystem
  if submitting from an NFS mount. This should reduce issues with held jobs on
  certain LIGO Data Grid clusters.

- Fix deprecation warning in ``ligo-skymap-stats``.

- Remove the deprecated ``ligo.skymap.postprocess.find_injection_moc`` method,
  which has been renamed to ``ligo.skymap.postprocess.crossmatch``.

0.1.16 (2020-02-26)
===================

- Update the compiler version that is used to build Linux wheels to icc
  19.1.0.166 from Intel Parallel Studio XE 2020u0. Due to C ABI requirements,
  the wheels are now built for the `manylinux2014
  <https://www.python.org/dev/peps/pep-0599/>`_ standard.

- Fix a unit test failure with astropy < 4.0.0.

- Add support for all combinations of map projection options, including
  ``geo degrees globe`` and ``geo degrees zoom``. Also, ``astro`` by itself is
  shorthand for ``astro hours``, and ``geo`` by itself is short for
  ``geo degrees``.

- ``ligo-skymap-plot`` now supports a variety of projections using the
  ``--projection`` option.

- Turn on continuous integration testing for Python 3.8.

- Change the license for the project as a whole to GPL 3.0 or later (GPLv3+).
  Previously, the source files had been a mix of GPLv2+ and GPLv3+.

- Add ``ligo-skymap-contour-moc`` command line to create a credible region 
  in a MOC (Multi Order Coverage) data structure. The input can be either a
  multiresolution or a flattened HEALPix probability map.

0.1.15 (2020-01-05)
===================

- Add support for the ``--detector-disabled`` command line option to the
  ``bayestar-localize-coincs`` tool, for consistency with
  ``bayestar-localize-lvalert`` tool.

- Remove installation dependency on astroquery, because it is only needed for
  the unit tests.

0.1.14 (2019-11-16)
===================

- Add a monkey patch to work around a regression in Astropy 3.2 that broke
  WCS transformations from ITRS to ICRS coordinates.
  See https://github.com/astropy/astropy/pull/9609.

- Fix a bug in the Python C extension code that could cause out-of-memory
  errors to be misreported as a SystemError with the message ``<built-in
  function rasterize> returned NULL without setting an error``, instead of as a
  MemoryError.

0.1.13 (2019-10-30)
===================

- The ``bayestar-inject`` script now assumes that the source distribution is
  specified per unit comoving volume per unit proper time, rather than per unit
  comoving volume per unit observer time. This is in agreement with the
  conventional definition for LIGO/Virgo astrophysical rates.

- The ``bayestar-inject`` and ``ligo-skymap-from-samples`` scripts now accept
  an optional integer value for the ``-j`` flag to set the number of
  subprocesses.

- ``ligo-skymap-from-samples`` will use all posterior samples if the value of
  the ``--maxpts`` argument is greater than or equal to the number of posterior
  samples.

- If the ``billiard`` package is present, then use it instead of the
  ``multiprocessing`` standard library module to parallelize
  ``ligo-skymap-from-samples`` so that the script's Python entry point can
  be called from daemon processes (for example, inside Celery tasks).

- Switch from WMAP9 to Planck15 cosmological parameters.

- ``ligo.skymap.kde.Clustered2DSkyKDE.as_healpix()`` has an optional
  ``top_nside`` to allow for better initial grid, before refinement.
  ``ligo-skymap-from-samples`` has an additional ``--top-nside`` argument,
  accordingly.

0.1.12 (2019-09-19)
===================

- Build macOS wheels with OpenMP.

- Record the command line with which ``ligo-skymap-stats`` was called by
  writing it to the ASCII table output as a comment line starting with ``#``.

0.1.11 (2019-08-28)
===================

- Fix a regression that caused ``ligo-skymap-flatten`` to fail for 2D sky maps.

0.1.10 (2019-08-28)
===================

- Add installation instructions for both pip and conda.

- Introduce the :mod:`ligo.skymap.postprocess.crossmatch` module for fast
  cross-matching of sky maps with galaxy redshift catalogs.

  This module used to be named :mod:`ligo.skymap.postprocess.find_injection`
  because it was originally designed for recovering injections (simulated
  signals) from sky localization simulations. We changed the name because
  galaxy cross matching is probably a more common use case than injection
  finding.

  The :func:`~ligo.skymap.postprocess.crossmatch.crossmatch` method also got
  some performance improvements for cross matching of large numbers of targets.
  Previously, to process :math:`n` targets, it took about :math:`(4 + 0.008 n)`
  seconds --- for a catalog of 300k targets, about 40 minutes. Now, it takes
  about 4 seconds total regardless of the number of targets.

  Note that the :mod:`ligo.skymap.postprocess.crossmatch` API is likely to
  change as documentation for it improves.

- Several performance improvements for BAYESTAR:

  - Add GCC branch prediction hints.

  - Exploit nested parallelism in radial integrator lookup table generation.

  - Calculate signal amplitudes using single-precision floating point.

  - Add tracepoints for Intel's Instrumentation and Tracing Technology (ITT)
    API, which can be enabled at build time by passing the ``--with-ittnotify``
    option to ``python setup.py build``.

0.1.9 (2019-08-02)
==================

- Switch from using the GNU Compiler Collection (gcc) to the Intel C Compiler
  (icc) for building optimized Linux binaries. On Intel Skylake machines, this
  can speed up BAYESTAR by 1.3x or more.

  Due to icc's C ABI requirements, Linux wheels now target the `manylinux2010
  <https://www.python.org/dev/peps/pep-0571/>`_ platform tag.

- In BAYESTAR, change the OpenMP scheduling kind from ``static`` (the default)
  to ``guided``. This improves CPU utilization by load-balancing work across
  threads more efficiently.

0.1.8 (2019-07-25)
==================

- Add ``ligo-skymap-constellations``, an easter egg program to list the most
  probable constellations for a localization, for fun and for public outreach
  purposes.

- Switch the implementation of the ``smooth`` option of ``imshow_hpx`` and
  ``contour_hpx`` from ``scipy.ndimage.gaussian_filter`` to
  ``astropy.convolution.convolve_fft`` in order to correctly handle points near
  the projection boundary where invalid values must be masked out.

- Register ``AutoScaledWCSAxes`` as a Matplotlib projection with the name
  ``astro wcs`` so that subclasses can be created using
  ``plt.axes(..., projection='astro wcs', header='...')``.

- Suppress Numpy warnings for HEALPix reprojection operations in WCS plots
  because it is normal for invalid values to occur when transforming pixels
  that lie outside of the projection.

- Add ``rotate`` option to ``astro globe``, ``geo globe``, and ``astro zoom``
  to rotate the plot in the plane of the screen about the center of the
  projection.

- Pass through keyword arguments from ``AutoScaledWCSAxes.scalebar()`` and
  ``AutoScaledWCSAxes.scalebar().label()`` to Matplotlib so that plot styles
  can be adjusted easily.

- Bump matplotlib version to >= 3.0.2 because of a bug that affected
  ``ligo-skymap-plot-stats``.

- The ``ligo-skymap-unflatten`` tool will now write multiresolution sky maps
  with pixels sorted by the ``UNIQ`` column, as required by the standard
  multi-order coverage map serialization in FITS.

- All functions in ``ligo.skymap.moc`` now assume that ``uniq`` is a signed
  integer. This makes it easier to call these functions with Numpy indexing
  routines, which work with signed integers. Also, saved multi-order sky maps
  will now be read correctly by tools such as ``fv`` from HEASOFT, which do not
  correctly handle unsigned integer columns.

- Add timestamps to the command line tools' default logging configuration in
  order to start characterizing the latency of BAYESTAR's data handling stages.

- Increase precision of BAYESTAR's run time measurement for the FITS headers.

0.1.7 (2019-04-24)
==================

- Add the ``ligo-skymap-plot-observability`` tool to plot observability windows
  for many sites at once. Conceptually, this tool is a variation of
  ``ligo-skymap-plot-airmass`` in which the sky position is integrated out.

- The ``ligo-skymap-plot-airmass`` tool will now use the color map's full
  dynamic range.

- Add ``order`` option to ``ligo.skymap.moc.rasterize`` and
  ``ligo.skymap.bayestar.rasterize`` and ``--nside`` option to
  ``ligo-skymap-flatten`` to support flattening multi-resolution HEALPix
  datasets to specified resolutions.

- ``ligo-skymap-stats`` now ignores skymaps with no corresponding entries in
  the inspinjfind database, instead of failing.

0.1.6 (2019-03-26)
==================

- Add options to ``ligo-skymap-plot-airmass`` to specify site coordinates
  explicitly rather than by a site nickname.

0.1.5 (2019-03-20)
==================

- Fix a bug caused by improper floating point comparison that caused some
  contours to be missing from the output of ``ligo-skymap-contour``.

- Speed up ``ligo-skymap-contour`` by skipping pixels that lie completely on
  the interior or exterior of the contour. For a typical LIGO/Virgo HEALPix map
  with a resolution of nside=512, the run time has decreased from about 42
  seconds to 3 seconds.

0.1.4 (2019-03-13)
==================

- The ``bayestar-localize-lvalert`` and ``ligo-skymap-from-samples`` tools will
  now generate multiresolution FITS files by default.

- Add ``--instrument`` option to ``ligo-skymap-from-samples`` to support
  storing metadata about which detectors contributed data.

0.1.3 (2019-03-04)
==================

- Fix a bug in ``ligo-skymap-plot-airmass`` that caused the airmass chart to be
  blank if the lower and upper credible levels were always in opposite
  hemispheres. The root cause was that ``plt.fill_between`` does not clip
  infinities to the plot's data range.

0.1.2 (2019-02-28)
==================

- Require lalsuite >6.53 and lscsoft-glue >=2.0.0 due to breaking changes in
  API and behavior for LIGO-LW XML reading.

0.1.1 (2019-02-20)
==================

- Pin lalsuite at <=6.52 and lscsoft-glue at <=1.60.0 due to breaking changes
  in API and behavior for LIGO-LW XML reading.

- Add the ``ligo-skymap-unflatten`` tool to convert flat, fixed resolution,
  implicitly indexed HEALPix files to multi-resolution HEALPix files. This
  tools is the inverse of ``ligo-skymap-flatten``.

0.1.0 (2019-02-01)
==================

- Migrate from glue.segments to ligo.segments.

- Add ``--min-inclination`` and ``max-inclination`` options to
  ``bayestar-localize-coincs`` and ``bayestar-localize-lvalert`` to control the
  limits of the isotropic prior over the inclination angle.

- Un-pin ligo-segments and require version >= 1.2.0 due to packaging
  bugfixes.

0.0.19 (2018-12-13)
===================

- Fix a bug that prevented the output of ligo-skymap-flatten from being
  gzip-compressed if the output filename ended in .gz.

- Require astropy >= 3.1 because some code that we previously had to
  monkeypatch went upstream. See
  https://github.com/astropy/astropy-healpix/pull/106.

- In the KDE clustering and ``ligo-skymap-from-samples``, disable OpenMP
  parallelism if Python mulitprocessing parallelism is enabled. This will
  prevent the program from spawning an excessive number of threads.

- ``ligo-skymap-plot`` no longer requires a DATE-OBS entry in the FITS header
  when plotting in astronomical coordinates.

0.0.18 (2018-11-19)
===================

- Fix a typo that caused ligo.skymap to always compile the bundled copy of
  chealpix instead of searching for a system version using pkgconfig.

- Un-pin Numpy version now that Numpy 1.15.4 is out.

- The ``bayestar-localize-lvalert`` and ``ligo-skymap-from-samples`` tools can
  now natively output multi-resolution HEALPix files, although they still
  natively output flat, fixed-resolution HEALPix files.

- Add the ``ligo-skymap-flatten`` tool to convert multi-resolution HEALPix
  files to flat, fixed-resolution, implicitly indexed HEALPix files.

- Bring back ``bayestar_samples_ppplot`` from LALInference as
  ``ligo-skymap-plot-pp-samples``, a tool for making P-P plots to compare a sky
  map with posterior samples.

- Add ``--cosmology`` feature to ``ligo-skymap-stats`` to calculate comoving
  volumes.

0.0.17 (2018-10-24)
===================

- In ``bayestar-mcmc``, correct a mistake in setting fixed parameters that
  undergo sampling transformations.

- By default, ``bayestar-realize-coincs`` will rewrite ``simulation_id`` values
  so that their integer values match the corresponding events'
  ``coinc_event_id`` values. The option ``--preserve-ids`` switches back to the
  old behavior of preserving the original ``simulation_id`` values.

- Track rename of ``ligo.gracedb.rest.GraceDb.service_url`` to
  ``ligo.gracedb.rest.GraceDb._service_url`` in ligo-gracedb >= 2.0.1.

- Update common files and submodules from the Astropy package template.

- Work around a change (possibly a regression?) in Numpy 1.15.3 that broke
  Astropy by requiring numpy <= 1.15.2. See
  <https://github.com/astropy/astropy/issues/7943>.

- Work around a bug introduced in ligo-segments 1.1.0 by requiring an earlier
  version of that package: its dependency on ligo-common, which does not
  correctly implement the namespace package ``ligo``, broke the continuous
  integration build.

- Depend on astropy-healpix >= 0.3 to pick up a bug fix related to HEALPix
  bilinear interpolation that affected ``ligo-skymap-plot``. See
  <https://github.com/astropy/astropy-healpix/pull/106>.

0.0.16 (2018-09-11)
===================

- Drop support for Python 3.5.

- The ``--condor-submit`` option of the ``bayestar-localize-coincs`` and
  ``bayestar-mcmc`` tools now passes the submit file directives to
  ``condor_submit`` via stdin rather than on the command line, so that the
  number of jobs is not limited by the operating system's maximum number of
  command line arguments.

- Print warnings from ``ligo.skymap.io.events.ligolw.open()`` only once per
  file to avoid excessive terminal output when reading large files.

- ``bayestar-realize-coincs`` now copies the process table from the injection
  file and fills in the SimInspiral table and associates coincidences with
  found injections. As a result, it is no longer necessary to run
  ``lalapps_inspinjfind`` on the output to find injections.

- ``bayestar-realize-coincs`` now prints a running count of the number of
  injections that have been found and saved.

0.0.15 (2018-09-04)
===================

- Parallelize ``bayestar-realize-coincs``.

- Add ``--min-distance`` and ``--max-distance`` options to
  ``bayestar-realize-coincs``.

- Add unit tests and binary wheels for Python 3.7.

0.0.14 (2018-08-28)
===================

- Increase lifetime of continuous integration artifacts. The unit tests take
  longer now because they are more complete.

0.0.13 (2018-08-27)
===================

- Add ``bayestar-mcmc`` tool for pure Markov Chain Monte Carlo parameter
  estimation, without sky map postprocessing but with options for holding
  parameters at fixed values.

- Fix a corner case in the initialization of the ``distance.marginal_ppf``
  solver that could cause NaN return values.

- Silence ``numpy.genfromtxt`` Unicode deprecation warning in
  ``ligo-skymap-plot-stats`` and update the minimum version of Numpy to 1.14.
  See the related `Numpy changelog entry
  <https://docs.scipy.org/doc/numpy/release.html#encoding-argument-for-text-io-functions>`_.

- Silence deprecation warning in ``ligo-skymap-plot-stats`` due to Matplotlib
  renaming the ``hist`` method's keyword argument from ``normed`` to
  ``density``.

- The ``bayestar-realize-coincs`` tool now copies over spins from the input
  ``sim_inspiral`` table to the output ``sngl_inspiral`` table.

- Switch the FFT implementation from LAL (which calls `FFTW
  <http://www.fftw.org>`_) to `scipy.fftpack
  <https://docs.scipy.org/doc/scipy/reference/tutorial/fftpack.html>`_, which
  is faster for small transform sizes (e.g. <= 1024).

- Add ``--seed`` option to ``bayestar-localize-coincs``,
  ``bayestar-localize-lvalert``, ``bayestar-mcmc``, and
  ``bayestar-realize-coincs``.

- Some reasonable sub-sample trigger interpolation schemes can return peak
  times that are almost a full sample away from the maximum sample if the SNR
  time series has a pronounced skew in one direction in the vicinity of the
  maximum. Such an example occurs for the ``catmull-rom`` interpolation method
  for the new unit tests in ``ligo.skymap.bayestar.tests.test_interpolation``.
  Because of this, relax the tolerance of BAYESTAR's sanity check on
  single-detector trigger times and SNR series timestamps to a full sample.

- Rewrite ``ligo-skymap-plot-stats`` to reduce code duplication.

- Add ``--measurement-error gaussian-noise`` option to
  ``bayestar-realize-coincs`` to simulate a matched filter in Gaussian noise.

- Remove deprecated module ``ligo.skymap.postprocess.detector_frame``.

0.0.12 (2018-07-18)
===================

- ``bayestar_localize_lvalert`` will now write the correct GraceDb URL
  to FITS file headers in the case that it is run with a non-default GraceDb
  server.

- BAYESTAR's SNR series time stamp assertions now include a bit more detail.

- Add phase convention for gstlal-spiir, which needs to be confirmed upstream.

- Fix datatype of simulated SNR time series produced by
  ``bayestar-realize-coincs``.

0.0.11 (2018-06-11)
===================

- Prebuilt binary wheels for macOS are now relocatable. See
  `delocate#38 <https://github.com/matthew-brett/delocate/pull/38>`_.

0.0.10 (2018-06-07)
===================

- Make lalsuite and lscsoft-glue required dependencies.

- The Python code is now required to pass linting by
  `Flake8 <http://flake8.pycqa.org/en/latest/>`_.

0.0.9 (2018-06-06)
==================

- On reading, rename columns from Fermi GBM HEALPix files to match the
  LIGO/Virgo convention. In particular, rename any column named `PROBABILITY`
  to `PROB`.

- Reduce the memory footprint of ``ligo-skymap-plot-airmass`` by transposing
  two nested loops.

- Make some cosmetic improvements to ``ligo-skymap-plot-airmass``:

  * Add altitude and local time axes.
  * Center plot on local solar midnight.
  * Adjust blending and z-order of twilight shading.

- ``ligo-skymap-plot-airmass`` will now write an airmass table to stdout.

- Rewrite the MCMC mode of BAYESTAR using ``ligo.skymap.ez_emcee``, a new
  reusable, fire-and-forget, parallel-tempering, MCMC sampler that features
  automated convergence testing and progress monitoring.

- Update common files from Astropy package template.

0.0.8 (2018-05-10)
==================

- Add ``ligo-skymap-combine``, a tool to combine sky localizations from
  different observations into a joint skymap.

0.0.7 (2018-04-27)
==================

- Move ``ligo.skymap.eigenframe.EigenFrame`` to
  ``ligo.skymap.coordinates.EigenFrame``.

- Add a new Astropy coordinate frame ``ligo.skymap.coordinates.DetectorFrame``
  to visualize triangulation rings with pairs of detectors.

- Deprecate all functions in ``ligo.skymap.postprocess.detector_frame``.

- Overhaul documentation so that all essential functionality is presented on
  the front page.

- Move ``ligo.skymap.command`` to top-level ``ligo.skymap.tool`` module.

- Require version 0.3.2 of the ``reproject`` package because of a regression
  that was caused by improper handling of nans in the ``astropy-healpix``
  package. See <https://github.com/astropy/astropy-healpix/pull/77>.

0.0.6 (2018-04-13)
==================

- Declare the top-level ``ligo`` module as a namespace package.

- Update common files from Astropy package template.

- Enable Python version check in ``setup.py`` and top-level namespace package.

0.0.5 (2018-04-12)
==================

- When running ``ligo-skymap-stats`` without injections, instead of writing
  ``nan`` values for irrelevant columns, don't write the columns in the first
  place.

- Start process of switching to tqdm for progress bars so that long-running
  operations show time estimates.

- In ``ligo-skymap-stats``, disable OpenMP parallelism if running with ``-j``
  to avoid creating a huge number of threads on machines with very many
  cores.

0.0.4 (2018-03-22)
==================

- Fix ``--condor-submit`` option for ``bayestar-localize-coincs``.

- Add ``--duty-cycle`` option to ``bayestar-realize-coincs``.

- Rename ``ligo-skymap-aggregate-found-injections`` to ``ligo-skymap-stats``
  and ``ligo-skymap-plot-found-injections`` to ``ligo-skymap-plot-stats``. The
  new ``ligo-skymap-stats`` program can generate summary statistics for
  skymaps, with or without injection-finding.

- This is the first version that has been tested and shown to reproduce the
  results in the "First Two Years" paper, which is the review benchmark.

0.0.3 (2018-03-21)
==================

- Bring back simulation tools from LALSuite.

- Add ``ligo-skymap-plot-airmass``, a tool for probabilistic airmass charts.

0.0.2 (2018-03-12)
==================

- Adjust CI configuration for uploading to PyPI.

0.0.1 (2018-03-12)
==================

- Initial release.
