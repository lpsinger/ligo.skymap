.. highlight:: sh

Development Instructions
========================

Source dependencies
-------------------

If you are building `ligo.skymap` from source, then in addition to the
requirements in the :doc:`quick start section <install>`, you will also need:

*  `GSL <https://www.gnu.org/software/gsl>`_ ≥ 1.15
*  `chealpix <https://sourceforge.net/projects/healpix/files/Healpix_3.30/>`_
   (note: if missing, will be built automatically from `bundled sources
   <https://git.ligo.org/lscsoft/ligo.skymap/tree/master/cextern/chealpix>`_)

Building from source
--------------------

To build `ligo.skymap` from source, first clone the git repository::

    $ git clone https://git.ligo.org/lscsoft/ligo.skymap.git

Then install it with pip::

    $ pip install -e .
