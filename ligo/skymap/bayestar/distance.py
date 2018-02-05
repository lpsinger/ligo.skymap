import warnings
from .. import distance
from ..distance import *
__all__ = distance.__all__

warnings.warn('lalinference.bayestar.distance is deprecated, use lalinference.distance instead', DeprecationWarning)
