###########
ligo.skymap
###########

The `ligo.skymap` package provides tools for reading, writing, generating, and
visualizing gravitational-wave probability maps from LIGO, Virgo, and KAGRA.
Some of the key features of this package are:

*  `Command line tool bayestar-localize-coincs`_ and
   `bayestar-localize-lvalert`_: BAYESTAR, providing rapid, coherent, Bayesian,
   3D position reconstruction for compact binary coalescence events

*  `Command line tool ligo-skymap-from-samples`_: Create 3D sky maps from
   posterior sample chains using kernel density estimation

*  `Command line tool ligo-skymap-plot`_: An everyday tool for plotting
   HEALPix maps

*  `Module ligo.skymap.plot`_: Astronomical mapmaking tools for
   perfectionists and figure connoisseurs

To get started, see the `installation instructions`_ or the `full
documentation`_.

.. figure:: https://lscsoft.docs.ligo.org/ligo.skymap/_images/localization.svg
   :alt: GW170817 localization

   This illustration of the position of GW170817, prepared using `ligo.skymap`,
   illustrates the main features of this package: rapid localization of compact
   binaries with BAYESTAR [#BAYESTAR]_, three-dimensional density estimation
   [#GoingTheDistance]_ [#GoingTheDistanceSupplement]_, cross-matching with
   galaxy catalogs, and visualization of gravitational-wave sky maps.
   Reproduced from [#IlluminatingGravitationalWaves]_.

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
   https://doi.org/10.1126/science.aap9455

.. _`Command line tool bayestar-localize-coincs`: https://lscsoft.docs.ligo.org/ligo.skymap/tool/bayestar_localize_coincs.html
.. _`bayestar-localize-lvalert`: https://lscsoft.docs.ligo.org/ligo.skymap/tool/bayestar_localize_lvalert.html
.. _`Command line tool ligo-skymap-from-samples`: https://lscsoft.docs.ligo.org/ligo.skymap/tool/ligo_skymap_from_samples.html
.. _`Command line tool ligo-skymap-plot`: https://lscsoft.docs.ligo.org/ligo.skymap/tool/ligo_skymap_plot.html
.. _`Module ligo.skymap.plot`: https://lscsoft.docs.ligo.org/ligo.skymap/#plotting-and-visualization-ligo-skymap-plot
.. _`installation instructions`: https://lscsoft.docs.ligo.org/ligo.skymap/quickstart/install.html
.. _`full documentation`: https://lscsoft.docs.ligo.org/ligo.skymap
