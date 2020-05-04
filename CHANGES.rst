#########
Changelog
#########

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
