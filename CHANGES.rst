#########
Changelog
#########

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
