.. highlight:: sh

3D Volume Rendering (`ligo-skymap-plot-volume`)
===============================================

Example
-------

To plot a 3D volume rendering for GW170104 and optionally mark a position
(here, RA=135° Dec=+37° luminosity distance=500 Mpc)::

    $ curl -OL https://gwosc.org/s/events/GW170104/P1500227/bayestar.fits.gz
    $ ligo-skymap-plot-volume bayestar.fits.gz --radecdist 135 37 500 -o bayestar.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_plot_volume import main
    from astropy.utils.data import download_file
    filename = download_file('https://gwosc.org/s/events/GW170104/P1500227/bayestar.fits.gz', cache=True)
    main([filename, '--radecdist', '135', '37', '500'])

To plot a single projection plane::

    $ curl -OL https://gwosc.org/s/events/GW170104/P1500227/bayestar.fits.gz
    $ ligo-skymap-plot-volume bayestar.fits.gz --projection 1 -o bayestar-plane1.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_plot_volume import main
    from astropy.utils.data import download_file
    filename = download_file('https://gwosc.org/s/events/GW170104/P1500227/bayestar.fits.gz', cache=True)
    main([filename, '--projection', '1'])

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_plot_volume
    :func: parser
