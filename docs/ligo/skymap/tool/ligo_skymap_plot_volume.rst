.. highlight:: sh

3D Volume Rendering (`ligo-skymap-plot-volume`)
===============================================

Example
-------

To plot a 3D volume rendering for GW170104::

    $ curl -O https://losc.ligo.org/s/events/GW170104/P1500227/bayestar.fits.gz
    $ ligo-skymap-plot-volume bayestar.fits.gz -o bayestar.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_plot_volume import main
    from astropy.utils.data import download_file
    filename = download_file('https://losc.ligo.org/s/events/GW170104/P1500227/bayestar.fits.gz', cache=True)
    main([filename])

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_plot_volume
    :func: parser
