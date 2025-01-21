from figaro.mixture import DPGMM
from figaro.utils import get_priors
from .coordinates import EigenFrame
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from tqdm import tqdm
from . import moc


class SkyDPGMM:


    def __init__(self, pts, **kwargs):
        self.frame = EigenFrame.for_coords(SkyCoord(*pts.T, unit=u.rad))
        
        # transform to eigenframe
        pts = self.transform(pts)

        # build DPGMM model
        bounds = [[0, 2*np.pi], [-np.pi/2, np.pi/2]]
        prior_pars = get_priors(bounds, pts)
       
        model = DPGMM(bounds, prior_pars=prior_pars)
        for s in tqdm(pts):
            model.add_new_point(s)

        self.model = model

    def transform(self, pts):
        pts = SkyCoord(*pts.T, unit=u.rad).transform_to(self.frame).spherical
        return np.column_stack((pts.lon.rad, np.sin(pts.lat.rad)))

    def __call__(self, pts):
        return self.model.pdf(self.transform(pts))

    def as_healpix(self, top_nside=16, rounds=8):
        return moc.bayestar_adaptive_grid(self, top_nside=top_nside,
                                          rounds=rounds)
    
