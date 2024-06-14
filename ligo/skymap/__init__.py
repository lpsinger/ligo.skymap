# This file is adapted from the Astropy package template, which is licensed
# under a 3-clause BSD style license - see licenses/TEMPLATE_LICENSE.rst

# Packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *   # noqa
# ----------------------------------------------------------------------------

__all__ = ('omp',)


class Omp:
    """OpenMP runtime settings.

    Attributes
    ----------
    num_threads : int
        Adjust the number of OpenMP threads. Getting and setting this attribute
        call :man:`omp_get_num_threads` and :man:`omp_set_num_threads`
        respectively.

    """

    @property
    def num_threads(self):
        from .core import get_num_threads
        return get_num_threads()

    @num_threads.setter
    def num_threads(self, value):
        from .core import set_num_threads
        set_num_threads(value)


class Itt():
    """Context manager to run code with Intel ITT tracing temporarily enabled.

    If you use this context manager, then you should start profiling with
    collection paused.

    See https://www.intel.com/content/www/us/en/docs/vtune-profiler/user-guide/2023-1/instrumentation-and-tracing-technology-apis.html
    """  # noqa: E501

    def __enter__(self):
        from .core import itt_resume
        itt_resume()
        return self

    def __exit__(self, *args, **kwargs):
        from .core import itt_pause
        itt_pause()


omp = Omp()
del Omp

itt = Itt()
del Itt
