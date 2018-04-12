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

0.0.4 (2017-03-22)
==================

- Fix ``--condor-submit`` option for ``bayestar-localize-coincs``.

- Add ``--duty-cycle`` option to ``bayestar-realize-coincs``.

- Rename ``ligo-skymap-aggregate-found-injections`` to ``ligo-skymap-stats``
  and ``ligo-skymap-plot-found-injections`` to ``ligo-skymap-plot-stats``. The
  new ``ligo-skymap-stats`` program can generate summary statistics for
  skymaps, with or without injection-finding.

- This is the first version that has been tested and shown to reproduce the
  results in the "First Two Years" paper, which is the review benchmark.

0.0.3 (2017-03-21)
==================

- Bring back simulation tools from LALSuite.

- Add ligo-skymap-plot-airmass, a tool for probabilistic airmass charts.

0.0.2 (2017-03-12)
==================

- Adjust CI configuration for uploading to PyPI.

0.0.1 (2017-03-12)
==================

- Initial release.
