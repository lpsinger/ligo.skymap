.. highlight:: sh

Joint Sky Localization (`ligo-skymap-combine`)
==============================================

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_combine
    :func: parser

Example
-------

Combine the LIGO-only localization of GW170817 and the Fermi/GBM localization
of GRB 170817A into a more precise localization::

    $ curl -OL https://dcc.ligo.org/public/0146/G1701985/001/bayestar_no_virgo.fits.gz
    $ curl -OL https://gammaray.nsstc.nasa.gov/gbm/science/grbs/grb170817a/gbuts_healpix_systematic.fit
    $ ligo-skymap-combine bayestar_no_virgo.fits.gz gbuts_healpix_systematic.fit hanford_livingston_gbm.fits.gz
    $ ligo-skymap-plot hanford_livingston_gbm.fits.gz -o hanford_livingston_gbm.png

.. plot::
   :context: reset
   :align: center

    from ligo.skymap.tool.ligo_skymap_combine import main as combine
    from ligo.skymap.tool.ligo_skymap_plot import main as plot
    from astropy.utils.data import download_file
    from tempfile import NamedTemporaryFile
    gw_map = download_file('https://dcc.ligo.org/public/0146/G1701985/001/bayestar_no_virgo.fits.gz', cache=True)
    gbm_map = download_file('https://gammaray.nsstc.nasa.gov/gbm/science/grbs/grb170817a/gbuts_healpix_systematic.fit', cache=True)
    with NamedTemporaryFile(suffix='.fits.gz') as combined_map:
        combine([gw_map, gbm_map, combined_map.name])
        plot([combined_map.name])
