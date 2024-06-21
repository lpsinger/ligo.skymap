.. highlight:: sh

Coherence vs. incoherence Bayes factor plots (`ligo-skymap-plot-coherence`)
===========================================================================

Example
-------

::

    $ curl -OL https://gracedb.ligo.org/api/superevents/S230707ai/files/bayestar.multiorder.fits,2
    $ ligo-skymap-plot-coherence bayestar.multiorder.fits,2 -o bayestar.coherence.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_plot_coherence import main
    from astropy.utils.data import download_file
    filename = download_file('https://gracedb.ligo.org/api/superevents/S230707ai/files/bayestar.multiorder.fits,2', cache=True)
    main([filename])

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_plot_coherence
    :func: parser
