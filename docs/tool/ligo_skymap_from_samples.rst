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
    >>> with open('skypost.obj', 'rb') as f:  # doctest: +SKIP
    ...     skypost = pickle.load(f)  # doctest: +SKIP
    ...
    >>> coord = SkyCoord.from_name('NGC 4993')
    >>> distance = np.arange(1, 100)
    >>> coords = np.column_stack((np.tile(coord.ra.rad, len(distance)),
    ...                           np.tile(coord.dec.rad, len(distance)),
    ...                           distance))  # doctest: +SKIP
    >>> post = skypost.posterior_spherical(coords)  # doctest: +SKIP
    >>> plt.plot(distance, post)  # doctest: +SKIP
    >>> plt.show()  # doctest: +SKIP

.. argparse::
    :module: ligo.skymap.tool.ligo_skymap_from_samples
    :func: parser
