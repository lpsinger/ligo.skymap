.. highlight:: sh

Installation
============

.. important:: The `ligo.skymap` package requires `Python`_ 3.8 or later.

On Linux or macOS x86_64 systems, we recommend installing `ligo.skymap` using
`pip`_ or `conda`_, either of which will automatically install all of the
additional :ref:`Python dependencies <python-dependencies>`.

(On other operating systems and architectures, you can :doc:`install from
source <../develop>`.)

Option 1: pip
-------------

To install `ligo.skymap` using `pip`_, you will need pip 19.3 or later. You can
check what version of pip you have by running this command::

    $ pip --version
    pip 20.0.2 from /usr/local/lib/python3.8/site-packages/pip (python 3.8)

If your version of pip is too old, then you can update pip to the most recent
version by running this command::

    $ pip install --upgrade pip

Then, just run this command::

    $ pip install ligo.skymap

You are now ready to get started using `ligo.skymap`.

Option 2: conda
---------------

If you are using the Anaconda Python distribution or the lightweight Miniconda
version, you can install `ligo.skymap` using `conda`_. First, enable the
`conda-forge`_ repository by running these commands::

    $ conda config --add channels conda-forge
    $ conda config --set channel_priority strict

Then, run this command::

    $ conda install ligo.skymap

You are now ready to get started using `ligo.skymap`.

.. _Python: https://www.python.org
.. _`pip`: https://pip.pypa.io
.. _`Python package index`: https://pypi.org/project/ligo.skymap/
.. _`conda`: https://conda.io
.. _`conda-forge`: https://conda-forge.org

.. _python-dependencies:
.. note:: When you use pip to install `ligo.skymap` with pip or conda, it will
          automatically install the following required Python packages:

          *  `Astroplan <http://astroplan.readthedocs.io>`_ ≥ 0.7
          *  `Astropy`_ ≥ 5.0
          *  `astropy-healpix <https://astropy-healpix.readthedocs.io>`_ ≥ 0.3
          *  `Healpy <http://healpy.readthedocs.io>`_
          *  `h5py <https://www.h5py.org>`_
          *  `LALSuite <https://pypi.python.org/pypi/lalsuite>`_ ≥ 6.53, ≠ 7.2
          *  `ligo-gracedb <https://pypi.org/project/ligo-gracedb/>`_ ≥ 2.0.1
          *  `ligo-segments <https://pypi.org/project/ligo-segments/>`_ ≥ 1.2.0
          *  `Matplotlib <https://matplotlib.org>`_ ≥ 3.4.0
          *  `NetworkX <https://networkx.github.io>`_
          *  `Numpy <http://www.numpy.org>`_ ≥ 1.19.3, ≠ 1.22.0
          *  `Pillow <http://pillow.readthedocs.io>`_ ≥ 2.5.0
          *  `ptemcee <https://github.com/willvousden/ptemcee>`_
          *  `python-ligo-lw <https://pypi.org/project/python-ligo-lw/>`_ ≥ 1.8.0
          *  `Reproject <https://reproject.readthedocs.io>`_ ≥ 0.3.2
          *  `Scipy <https://www.scipy.org>`_ ≥ 0.14
          *  `tqdm <https://tqdm.github.io>`_ ≥ 4.27.0
          *  `pytz <http://pytz.sourceforge.net>`_

          The following packages are optional for specific features:

          *  `pytest <https://docs.pytest.org>`_ for running the test suite

Optional LALSimulation Data
---------------------------

The following instructions are only relevant if you are installing ligo.skymap
for the purpose of generating localizations with BAYESTAR (e.g., for analysis
of LIGO/Virgo data or for simulations) and you are **not** using a LIGO Data
Grid cluster.

Some gravitational waveform approximants in LALSuite (notably, reduced order
models) rely on external data files. These data files are part of
`lalsuite-extra`_, which must be installed separately. To install these data
files, run the following commands::

    $ curl -O https://software.igwn.org/lscsoft/source/lalsuite-extra-1.3.0.tar.gz
    $ tar xf lalsuite-extra-1.3.0.tar.gz
    $ cd lalsuite-extra-1.3.0
    $ ./configure --prefix=$HOME/.local
    $ make install

Then, add the following line to your shell profile script (``~/.profile``,
``~/.bashrc``, or similar)::

    export LAL_DATA_PATH=$HOME/.local/share/lalsimulation

Then log out and log back in.

.. _`lalsuite-extra`: https://git.ligo.org/lscsoft/lalsuite-extra
