.. highlight:: sh

Airmass Plots (`ligo-skymap-plot-airmass`)
==========================================

Example
-------

To make an airmass plot for GW170817::

    $ curl -OL https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz
    $ ligo-skymap-plot-airmass --site 'Las Campanas Observatory' --time 2017-08-17 bayestar.fits.gz -o bayestar.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_plot_airmass import main
    from astropy.utils.data import download_file
    filename = download_file('https://dcc.ligo.org/public/0146/G1701985/001/bayestar.fits.gz', cache=True)
    main([filename, '--site', 'Las Campanas Observatory', '--time', '2017-08-17'])

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_plot_airmass
    :func: parser
