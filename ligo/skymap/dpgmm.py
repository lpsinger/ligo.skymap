from figaro.mixture import DPGMM
from figaro.utils import get_priors
from .coordinates import EigenFrame
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from tqdm import tqdm
from . import moc
import copyreg

class _SkyDPGMMMeta(type):  # noqa: N802
    """Metaclass to make dynamically created subclasses of Clustered2DSkyKDE
    picklable.
    """


def _SkyDPGMM_pickle(cls):  # noqa: N802
    """Pickle dynamically created subclasses of Clustered2DSkyKDE."""
    return type, (cls.__name__, cls.__bases__, {'frame': cls.frame})


# Register function to pickle subclasses of Clustered2DSkyKDE.
copyreg.pickle(_SkyDPGMMMeta, _SkyDPGMM_pickle)


def _SkyDPGMM_factory(name, frame):  # noqa: N802
    """Unpickle instances of dynamically created subclasses of
    Clustered2DSkyKDE.

    FIXME: In Python 3, we could make this a class method of Clustered2DSkyKDE.
    Unfortunately, Python 2 is picky about pickling bound class methods.
    """
    new_cls = type(name, (SkyDPGMM,), {'frame': frame})
    return super(SkyDPGMM, SkyDPGMM).__new__(new_cls)

class SkyDPGMM:

    frame = None

    @classmethod
    def transform(cls, pts):
        pts = SkyCoord(*pts.T, unit=u.rad).transform_to(cls.frame).spherical
        return np.column_stack((pts.lon.rad, np.sin(pts.lat.rad)))
    
    def __new__(cls, pts, *args, **kwargs):
        frame = EigenFrame.for_coords(SkyCoord(*pts.T, unit=u.rad))
        name = '{:s}_{:x}'.format(cls.__name__, id(frame))
        new_cls = type(name, (cls,), {'frame': frame})
        return super().__new__(new_cls)

    def __init__(self, pts, **kwargs):
        # transform to eigenframe
        pts = self.transform(pts)

        # build DPGMM model
        prior_pars = None #get_priors(self.bounds, pts)
        model = DPGMM(self.bounds, prior_pars=prior_pars)
        for s in tqdm(pts):
            model.add_new_point(s)

        self.model = model

    @property
    def bounds(self):
        """Bounds on ra and dec"""
        return [[0, 2*np.pi], [-np.pi/2, np.pi/2]]

    def __call__(self, pts):
        return self.model.pdf(self.transform(pts))
    
    def __reduce__(self):
        """Pickle instances of dynamically created subclasses of
        Clustered2DSkyKDE.
        """
        factory_args = self.__class__.__name__, self.frame
        return _SkyDPGMM_factory, factory_args, self.__dict__
    
    def as_healpix(self, top_nside=16, rounds=8):
        return moc.bayestar_adaptive_grid(self, top_nside=top_nside,
                                          rounds=rounds)

