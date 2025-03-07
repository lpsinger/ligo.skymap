#
# Copyright (C) 2012-2025  Ethan Marx <emarx@mit.edu>
#                          Leo P. Singer <leo.singer@ligo.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#

from astropy.coordinates import SkyCoord
from astropy import units as u
try:
    from figaro.mixture import DPGMM
    from figaro.utils import get_priors
except ModuleNotFoundError as e:
    raise RuntimeError('In order to use the DPGMM feature'
                       'you must install `figaro`') from e

import numpy as np
from tqdm.auto import tqdm

from .coordinates import EigenFrame
from . import moc


class SkyDPGMM:

    def __init__(self, pts, **kwargs):
        self.frame = EigenFrame.for_coords(SkyCoord(*pts.T, unit=u.rad))

        # transform to eigenframe
        pts = self.transform(pts)

        # build DPGMM model
        bounds = [[0, 2*np.pi], [-1, 1]]
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
