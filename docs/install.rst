.. highlight:: sh

Installation
============

The recommended way to install the latest stable release of `ligo.skymap` is
with `pip <https://pip.pypa.io>`_. This method will install a pre-built binary
from the `Python package index <https://pypi.org/project/ligo.skymap/>`_ and
will install all required Python packages automatically.

Basic requirements
------------------

*  Linux or macOS
*  `Python <https://www.python.org>`_ ≥ 3.5
*  `pip <https://pip.pypa.io>`_

Python dependencies
-------------------

When you use pip to install `ligo.skymap` with pip, it will automatically
install the following required Python packages:

*  Astropy_ ≥ 3.0
*  `Numpy <http://www.numpy.org>`_ ≥ 1.13
*  `Healpy <http://healpy.readthedocs.io>`_ ≥ 1.9.1
   (note, considering a transition to
   `astropy-healpix <http://astropy-healpix.readthedocs.io>`_,
   which is easier to install)
*  `h5py <https://www.h5py.org>`_
*  `Matplotlib <https://matplotlib.org>`_ ≥ 2.1.0
*  `Pillow <http://pillow.readthedocs.io>`_ ≥ 2.5.0
*  `Reproject <https://reproject.readthedocs.io>`_ = 0.3.2
*  `Scipy <https://www.scipy.org>`_ ≥ 0.14

Optional dependencies
---------------------

The following packages are optional for specific features:

*  `LALSuite <https://pypi.python.org/pypi/lalsuite>`_ for using BAYESTAR
*  `pytest <https://docs.pytest.org>`_ for running the test suite
*  `Astroplan <http://astroplan.readthedocs.io/>`_ for airmass charts

Quick start
-----------

Just run this command::

    $ pip install ligo.skymap

If you want to run BAYESTAR, you will also need to install LALSuite::

    $ pip install lalsuite ligo.skymap
