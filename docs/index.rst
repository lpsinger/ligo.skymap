###########
ligo.skymap
###########

.. image:: https://vignette.wikia.nocookie.net/memoryalpha/images/c/cf/Picard_and_Data_in_stellar_cartography.jpg/revision/latest/scale-to-width-down/640?cb=20100527083827&path-prefix=en
   :alt: Picard and Data in stellar cartography
   :width: 640px
   :height: 269px

*LIGO's stellar cartography department.*

The `ligo.skymap` package provides tools for reading, writing, generating, and
visualizing gravitational-wave probability maps from LIGO and Virgo. It
provides several tools that used to live in `LALSuite
<http://git.ligo.org/lscsoft/lalsuite>`_ and elsewhere, but in the form of a
tiny Python package that is easier to install. Some of the key features of this
package are:

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

*****************
API Documentation
*****************

.. toctree::
   :maxdepth: 1

   ligo/skymap/bayestar
   ligo/skymap/bayestar/ez_emcee
   ligo/skymap/coordinates/index
   ligo/skymap/distance
   ligo/skymap/healpix_tree
   ligo/skymap/io/index
   ligo/skymap/kde
   ligo/skymap/moc
   ligo/skymap/plot/index
   ligo/skymap/postprocess/index
   ligo/skymap/tool/index
   ligo/skymap/util/index

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
   ligo/skymap/tool/ligo_skymap_plot_volume

Postprocessing
--------------

.. toctree::
   :maxdepth: 1

   ligo/skymap/tool/ligo_skymap_combine
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
