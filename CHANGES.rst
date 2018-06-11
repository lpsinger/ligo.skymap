#########
Changelog
#########

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

- Add ligo-skymap-plot-airmass, a tool for probabilistic airmass charts.

0.0.2 (2018-03-12)
==================

- Adjust CI configuration for uploading to PyPI.

0.0.1 (2018-03-12)
==================

- Initial release.
