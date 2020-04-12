###########
ligo.skymap
###########

The `ligo.skymap` package provides tools for reading, writing, generating, and
visualizing gravitational-wave probability maps from LIGO and Virgo. Some of
the key features of this package are:

*  Command line tool :doc:`bayestar-localize-coincs
   <ligo/skymap/tool/bayestar_localize_coincs>`: BAYESTAR, providing rapid,
   coherent, Bayesian, 3D position reconstruction for compact binary
   coalescence events

*  Command line tool :doc:`ligo-skymap-from-samples
   <ligo/skymap/tool/ligo_skymap_from_samples>`: Create 3D sky maps from
   posterior sample chains using kernel density estimation

*  Command line tool :doc:`ligo-skymap-plot
   <ligo/skymap/tool/ligo_skymap_plot>`: An everyday tool for plotting HEALPix
   maps

*  Module :mod:`ligo.skymap.plot.allsky`: Astronomical mapmaking tools for
   perfectionists and figure connoisseurs

.. figure:: _static/localization.svg
   :alt: GW170817 localization

   This illustration of the position of GW170817, prepared using `ligo.skymap`,
   illustrates the main features of this package: rapid localization of compact
   binaries with BAYESTAR [#BAYESTAR]_, three-dimensional density estimation
   [#GoingTheDistance]_ [#GoingTheDistanceSupplement]_, cross-matching with
   galaxy catalogs, and visualization of gravitational-wave sky maps.
   Reproduced from [#IlluminatingGravitationalWaves]_.

**********************
Quick Start, Tutorials
**********************

.. toctree::
   :maxdepth: 1

   quickstart/install
   help
   contributing
   quickstart/bayestar-injections

*******
Modules
*******

Coordinate Frames (`ligo.skymap.coordinates`)
---------------------------------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/coordinates/detector
   ligo/skymap/coordinates/eigenframe

I/O and Data Format Support (`ligo.skymap.io`)
----------------------------------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/io/events
   ligo/skymap/io/fits
   ligo/skymap/io/hdf5

Plotting and Visualization (`ligo.skymap.plot`)
-----------------------------------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/plot/allsky
   ligo/skymap/plot/backdrop
   ligo/skymap/plot/marker
   ligo/skymap/plot/pp

Sky Map Postprocessing (`ligo.skymap.postprocess`)
--------------------------------------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/postprocess/contour
   ligo/skymap/postprocess/cosmology
   ligo/skymap/postprocess/crossmatch
   ligo/skymap/postprocess/ellipse
   ligo/skymap/postprocess/util

Localization
------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/bayestar
   ligo/skymap/bayestar/ez_emcee
   ligo/skymap/distance
   ligo/skymap/healpix_tree
   ligo/skymap/kde
   ligo/skymap/moc

Utilities
---------

.. toctree::
   :maxdepth: 1

   ligo/skymap/util/file
   ligo/skymap/util/numpy
   ligo/skymap/util/sqlite
   ligo/skymap/tool/index

******************
Command Line Tools
******************

BAYESTAR Rapid Sky Localization
-------------------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/tool/bayestar_localize_coincs
   ligo/skymap/tool/bayestar_localize_lvalert
   ligo/skymap/tool/bayestar_mcmc
   ligo/skymap/tool/bayestar_realize_coincs
   ligo/skymap/tool/bayestar_sample_model_psd
   ligo/skymap/tool/bayestar_inject

Sky Map Visualization
---------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/tool/ligo_skymap_contour
   ligo/skymap/tool/ligo_skymap_contour_moc
   ligo/skymap/tool/ligo_skymap_plot
   ligo/skymap/tool/ligo_skymap_plot_airmass
   ligo/skymap/tool/ligo_skymap_plot_observability
   ligo/skymap/tool/ligo_skymap_plot_volume

Postprocessing
--------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/tool/ligo_skymap_combine
   ligo/skymap/tool/ligo_skymap_constellations
   ligo/skymap/tool/ligo_skymap_from_samples
   ligo/skymap/tool/ligo_skymap_plot_stats
   ligo/skymap/tool/ligo_skymap_stats
   ligo/skymap/tool/ligo_skymap_flatten
   ligo/skymap/tool/ligo_skymap_unflatten

**************
Developer Info
**************

.. toctree::
   :maxdepth: 1

   changes
   develop
   testing
   interface
   GitLab Project Page <https://git.ligo.org/lscsoft/ligo.skymap>
   GitHub Mirror <https://github.com/lpsinger/ligo.skymap>
   Python Package Index <https://pypi.org/project/ligo.skymap>

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

**********
References
**********

.. [#BAYESTAR]
   Singer, L. P., & Price, L. R. 2016, *Phys. Rev. D*, 93, 024013.
   https://doi.org/10.1103/PhysRevD.93.024013

.. [#GoingTheDistance]
   Singer, L. P., Chen, H.-Y., Holz, D. E., et al. 2016, *Astropys. J. Lett.*,
   829, L15. https://doi.org/10.3847/2041-8205/829/1/L15

.. [#GoingTheDistanceSupplement]
   Singer, L. P., Chen, H.-Y., Holz, D. E., et al. 2016, *Astropys. J. Supp.*,
   226, 10. https://doi.org/10.3847/0067-0049/226/1/10

.. [#IlluminatingGravitationalWaves]
   Kasliwal, M. M., Nakar, E., Singer, L. P. et al. 2019, *Science*, 358, 1559.
   https://10.1126/science.aap9455
