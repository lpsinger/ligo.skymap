Fire-and-Forget MCMC Sampling (`ligo.skymap.bayestar.ez_emcee`)
===============================================================

.. autofunction:: ligo.skymap.bayestar.ez_emcee.ez_emcee

    **Examples**

    .. code-block:: pycon

        >>> from ligo.skymap.bayestar.ez_emcee import ez_emcee
        >>> from matplotlib import pyplot as plt
        >>> import numpy as np
        >>>
        >>> def log_prob(params):
        ...     """Eggbox function"""
        ...     return 5 * np.log((2 + np.cos(0.5 * params).prod(-1)))
        ...
        >>> lo = [-3*np.pi, -3*np.pi]
        >>> hi = [+3*np.pi, +3*np.pi]
        >>> chain = ez_emcee(log_prob, lo, hi, vectorize=True)
        Sampling:  51%|██  | 8628/16820 [00:04<00:04, 1966.74it/s, accept=0.535, acl=62]
        >>> plt.plot(chain[:, 0], chain[:, 1], '.')

    .. image:: eggbox.png
