from astropy import units as u
import numpy as np
import pytest

from ..cosmology import cosmo, dVC_dVL_for_DL, dVC_dVL_for_z, dDL_dz_for_z
from ...util.math import derivative

redshift_range = pytest.mark.parametrize('z', np.logspace(-6, 2))


@redshift_range
def test_dDL_dz(z):
    expected = derivative(
        cosmo.luminosity_distance, z, dx=0.001 * z
    ).to_value(u.Mpc)
    result = dDL_dz_for_z(z).to_value(u.Mpc)
    assert expected == pytest.approx(result)


@redshift_range
def test_dVC_dVL(z):
    dVC_dz = cosmo.differential_comoving_volume(z)

    DH = cosmo.hubble_distance
    DL = cosmo.luminosity_distance(z)
    DM = cosmo.comoving_transverse_distance(z)

    dDC_dz = DH * cosmo.inv_efunc(z)  # DC = integral of inv_efunc

    assert cosmo.Ok0 == 0  # Flat universe
    dDM_dz = dDC_dz  # otherwise this expression is more complicated

    dDL_dz = DM + (1 + z) * dDM_dz  # DL = (1 + z) DM

    dVL_dz = DL**2 * dDL_dz * u.sr**-1  # by definition

    expected = (dVC_dz / dVL_dz).to_value(u.dimensionless_unscaled)

    result = dVC_dVL_for_z(z)
    assert expected == pytest.approx(result)

    result = dVC_dVL_for_DL(DL.to_value(u.Mpc))
    assert expected == pytest.approx(result)
