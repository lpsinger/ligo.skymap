###########
ligo.skymap
###########

.. figure:: https://vignette.wikia.nocookie.net/memoryalpha/images/c/cf/Picard_and_Data_in_stellar_cartography.jpg/revision/latest/scale-to-width-down/640?cb=20100527083827&path-prefix=en
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

*  `Script bayestar-localize-coincs`_: BAYESTAR, providing rapid,
   coherent, Bayesian, 3D position reconstruction for compact binary
   coalescence events

*  `Script ligo-skymap-from-samples`_: Create 3D sky maps from
   posterior sample chains using kernel density estimation

*  `Script ligo-skymap-plot`_: An everyday tool for plotting
   HEALPix maps

*  `Module ligo.skymap.plot`_: Astronomical mapmaking tools for
   perfectionists and figure connoisseurs

See the `installation instructions`_ or the `full documentation`_.

.. _`Script bayestar-localize-coincs`: https://lscsoft.docs.ligo.org/ligo.skymap/ligo/skymap/tool/bayestar_localize_coincs.html
.. _`Script ligo-skymap-from-samples`: https://lscsoft.docs.ligo.org/ligo.skymap/ligo/skymap/tool/ligo_skymap_from_samples.html
.. _`Script ligo-skymap-plot`: https://lscsoft.docs.ligo.org/ligo.skymap/ligo/skymap/tool/ligo_skymap_plot.html
.. _`Module ligo.skymap.plot`: https://lscsoft.docs.ligo.org/ligo.skymap/ligo/skymap/plot
.. _`installation instructions`: https://lscsoft.docs.ligo.org/ligo.skymap/quickstart/install.html
.. _`full documentation`: https://lscsoft.docs.ligo.org/ligo.skymap
