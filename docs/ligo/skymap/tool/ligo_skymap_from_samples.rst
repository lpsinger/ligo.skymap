Posterior Samples to FITS Files (`ligo-skymap-from-samples`)
============================================================

Example
-------

This example shows how to extract the conditional distance posterior for a
given sky location from the kernel density estimator.

First, run ``ligo-skymap-from-samples`` to create the KDE and the 3D sky map
from the posterior sample chain:

.. code-block:: sh

    $ ligo-skymap-from-samples --maxpts 1000 posterior_samples.dat

Then you can load the pickle in a Python interpreter and evaluate the KDE
at any point:

.. code-block:: pycon

    >>> from astropy.coordinates import SkyCoord
    >>> from matplotlib import pyplot as plt
    >>> import numpy as np
    >>> import pickle
    >>> with open('skypost.obj', 'rb') as f:
    ...     skypost = pickle.load(f)
    ...
    >>> coord = SkyCoord.from_name('NGC 4993')
    >>> distance = np.arange(1, 100)
    >>> coords = np.column_stack((np.tile(coord.ra.rad, len(distance)),
    ...                           np.tile(coord.dec.rad, len(distance)),
    ...                           distance))
    >>> post = skypost.posterior_spherical(coords)
    >>> plt.plot(distance, post)
    >>> plt.show()

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_from_samples
    :func: parser
