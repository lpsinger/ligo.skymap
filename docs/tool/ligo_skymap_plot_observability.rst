.. highlight:: sh

Observability Plots (`ligo-skymap-plot-observability`)
======================================================

This is similar to :doc:`ligo-skymap-plot-airmass <ligo_skymap_plot_airmass>`,
but integrated over the sky in order to create a bar chart and show the
observability windows for several facilities at once.

Example
-------

To make an observability plot for GW170817::

    $ curl -OL https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz
    $ ligo-skymap-plot-observability --site 'Las Campanas Observatory' 'SAAO' 'Siding Spring Observatory' --time 2017-08-17 bayestar.fits.gz -o bayestar.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_plot_observability import main
    from astropy.utils.data import download_file
    filename = download_file('https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz', cache=True)
    main([filename, '--site', 'Las Campanas Observatory', 'SAAO', 'Siding Spring Observatory', '--time', '2017-08-17'])

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_plot_observability
    :func: parser
