###########
ligo.skymap
###########

.. image:: https://vignette.wikia.nocookie.net/memoryalpha/images/c/cf/Picard_and_Data_in_stellar_cartography.jpg/revision/latest/scale-to-width-down/640?cb=20100527083827&path-prefix=en
   :alt: Picard and Data in stellar cartography
   :width: 640px
   :height: 269px

*LIGO's stellar cartography department.*

The `ligo.skymap` package provides tools for reading, writing, generating, and
visualizing gravitational-wave probability maps from LIGO and Virgo. Some of
the key features of this package are:

*  :doc:`ligo/skymap/tool/bayestar_localize_coincs`: BAYESTAR, providing rapid,
   coherent, Bayesian, 3D position reconstruction for compact binary
   coalescence events

*  :doc:`ligo/skymap/tool/ligo_skymap_from_samples`: Create 3D sky maps from
   posterior sample chains using kernel density estimation

*  :doc:`ligo/skymap/tool/ligo_skymap_plot`: An everyday tool for plotting
   HEALPix maps

*  Module :mod:`ligo.skymap.plot`: Astronomical mapmaking tools for
   perfectionists and figure connoisseurs

**********************
Quick Start, Tutorials
**********************

.. toctree::
   :maxdepth: 1

   quickstart/install
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
   ligo/skymap/postprocess/ellipse
   ligo/skymap/postprocess/find_injection
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

Sky Map Visualization
---------------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/tool/ligo_skymap_contour
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
   GitLab Project Page <http://git.ligo.org/lscsoft/ligo.skymap>
   Python Package Index <https://pypi.org/project/ligo.skymap>

*****
Index
*****

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
