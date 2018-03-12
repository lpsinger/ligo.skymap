###########
ligo.skymap
###########

.. image:: https://vignette.wikia.nocookie.net/memoryalpha/images/c/cf/Picard_and_Data_in_stellar_cartography.jpg/revision/latest/scale-to-width-down/640?cb=20100527083827&path-prefix=en
   :alt: Picard and Data in stellar cartography
   :width: 640px
   :height: 269px

*LIGO's stellar cartography department.*

=========  ===============  ==========  ========
**Build**  **Python code**  **C code**  **Docs**
|_build_|  |_python_code_|  |_c_code_|  |docs|__
=========  ===============  ==========  ========

.. |_build_| image:: https://git.ligo.org/leo-singer/ligo.skymap/badges/master/pipeline.svg
   :alt: pipeline status
   :target: https://git.ligo.org/leo-singer/ligo.skymap/pipelines

.. |_python_code_| image:: https://git.ligo.org/leo-singer/ligo.skymap/badges/master/coverage.svg?job=coverage:py
   :alt: coverage report
   :target: https://leo-singer.docs.ligo.org/ligo.skymap/cov/py

.. |_c_code_| image:: https://git.ligo.org/leo-singer/ligo.skymap/badges/master/coverage.svg?job=coverage:c
   :alt: coverage report
   :target: https://leo-singer.docs.ligo.org/ligo.skymap/cov/c

.. |docs| replace:: Latest
__ https://leo-singer.docs.ligo.org/ligo.skymap/

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

.. note:: Since the migration of these tools from LALSuite is very recent, the
   modules, scripts, and documentation in this package are subject to change.
