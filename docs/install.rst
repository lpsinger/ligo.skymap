Installation
============

Required dependencies
---------------------

`ligo.skymap` has the following Python requirements:

*  `Python <https://www.python.org>`_ ≥ 3.5
*  `Astropy <http://www.astropy.org>`_ ≥ 1.3.1
*  `Numpy <http://www.numpy.org>`_ ≥ 1.10
*  `Healpy <http://healpy.readthedocs.io>`_ ≥ 1.9.1
   (note, considering a transition to
   `astropy-healpix <http://astropy-healpix.readthedocs.io>`_,
   which is easier to install)
*  `h5py <https://www.h5py.org>`_
*  `Matplotlib <https://matplotlib.org>`_ ≥ 2.1.0
*  `Pillow <http://pillow.readthedocs.io>`_ ≥ 2.5.0
*  `Reproject <https://reproject.readthedocs.io>`_ ≥ 0.3.2
*  `Scipy <https://www.scipy.org>`_ ≥ 0.14

Source dependencies
-------------------

If you are building `ligo.skymap` from source, you will also need:

*  `GSL <https://www.gnu.org/software/gsl>`_ ≥ 1.15
*  `chealpix <https://sourceforge.net/projects/healpix/files/Healpix_3.30/>`_
   (note, source code bundled with `ligo.skymap`)

Optional dependencies
---------------------

The following packages are optional for specific features:

*  `LALSuite <https://pypi.python.org/pypi/lalsuite>`_ for using BAYESTAR
*  `pytest <https://docs.pytest.org>`_ for running the tests suite

Quick start
-----------

The recommended way to install `ligo.skymap` is with
`pip <https://pip.pypa.io>`_::

    $ pip install ligo.skymap

Note, if you want to run BAYESTAR, you will also need to install LALSuite::

    $ pip install lalsuite ligo.skymap
