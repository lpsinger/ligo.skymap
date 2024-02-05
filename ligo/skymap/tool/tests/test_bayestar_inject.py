from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
import lal
import lalsimulation
import numpy as np
import pytest

from ...util.math import derivative
from ..bayestar_inject import get_decisive_snr, cell_max, GWCosmo


def test_get_decisive_snr():
    snrs = [1.0, 3.0, 2.0, 4.0]
    assert get_decisive_snr(snrs, 1) == 4.0
    assert get_decisive_snr(snrs, 2) == 3.0
    assert get_decisive_snr(snrs, 3) == 2.0
    assert get_decisive_snr(snrs, 4) == 1.0
    with pytest.raises(IndexError):
        get_decisive_snr(snrs, 5)


def get_snr_at_z_lalsimulation(cosmo, z, mass1, mass2, f_low, f_high, psd):
    """Calculate optimal SNR the LALSimulation way."""
    params = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertRedshift(params, z)

    # "Signal" waveform with requested inclination angle.
    # Take (arbitrarily) only the plus polarization.
    h, _ = lalsimulation.SimInspiralFD(
        mass1 * lal.MSUN_SI, mass2 * lal.MSUN_SI,
        0, 0, 0, 0, 0, 0, cosmo.comoving_distance(z).to_value(u.m),
        0, 0, 0, 0, 0, 1e-3 * (f_high - f_low), f_low, f_high, f_low,
        params, lalsimulation.IMRPhenomPv2)
    return lalsimulation.MeasureSNRFD(h, psd, f_low, f_high)


@pytest.mark.parametrize('mtotal', [2.8, 10.0, 50.0, 100.0])
@pytest.mark.parametrize('z', [0.001, 0.01, 0.1, 1.0, 2.0])
def test_z_at_snr(mtotal, z):
    gwcosmo = GWCosmo(cosmo)
    f_low = 10
    f_high = 4096
    df = 0.1
    mass1 = mass2 = 0.5 * mtotal

    psd = lal.CreateREAL8FrequencySeries(
        '', 0, f_low, df, lal.DimensionlessUnit, int((f_high - f_low) // df))
    lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087(psd, f_low)

    snr = get_snr_at_z_lalsimulation(
        gwcosmo.cosmo, z, mass1, mass2, f_low, f_high, psd)
    z_solution = gwcosmo.z_at_snr(
        [psd], 'IMRPhenomPv2', f_low, snr, 1, mass1, mass2, 0, 0)

    assert z_solution == pytest.approx(z, rel=1e-2)


def test_get_max_z():
    gwcosmo = GWCosmo(cosmo)
    f_low = 10
    f_high = 4096
    df = 0.1
    waveform = 'IMRPhenomPv2'
    snr = 8
    m1 = np.asarray([50.0])
    m2 = np.asarray([30.0, 50.0])
    x1 = np.asarray([-1.0, 0.0, 1.0])
    x2 = np.asarray([-1.0, -0.5, 0.5, 1.0])

    psd = lal.CreateREAL8FrequencySeries(
        '', 0, f_low, df, lal.DimensionlessUnit, int((f_high - f_low) // df))
    lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087(psd, f_low)

    result = gwcosmo.get_max_z(
        [psd], waveform, f_low, snr, 1,
        *np.meshgrid(m1, m2, x1, x2, indexing='ij'))
    # Check that shape matches
    assert result.shape == (1, 2, 3, 4)
    # Spot check some individual cells
    for im1, m1_ in enumerate(m1):
        for im2, m2_ in enumerate(m2):
            for ix1, x1_ in enumerate(x1):
                for ix2, x2_ in enumerate(x2):
                    expected = gwcosmo.z_at_snr(
                        [psd], waveform, f_low, snr, 1, m1_, m2_, x1_, x2_)
                    assert result[im1, im2, ix1, ix2] == expected


def test_sensitive_volume_0():
    gwcosmo = GWCosmo(cosmo)
    assert gwcosmo.sensitive_volume(0) == 0


@pytest.mark.parametrize('z', [0.001, 0.01, 0.1, 1.0, 2.0])
def test_sensitive_volume(z):
    """Test that the sensitive volume has the correct derivative with z."""
    gwcosmo = GWCosmo(cosmo)
    dVC_dz = gwcosmo.cosmo.differential_comoving_volume(z)
    expected = (dVC_dz / (1 + z) * 4 * np.pi * u.sr).to_value(u.Mpc**3)
    actual = derivative(
        gwcosmo.sensitive_volume, x0=z, dx=1e-6 * z).to_value(u.Mpc**3)
    assert expected == pytest.approx(actual)


@pytest.mark.parametrize('z', [0.0, 0.001, 0.01, 0.1, 1.0, 2.0])
def test_sensitive_distance(z):
    gwcosmo = GWCosmo(cosmo)
    expected = gwcosmo.sensitive_volume(z).to_value(u.Mpc**3)
    actual = (
        4/3 * np.pi * gwcosmo.sensitive_distance(z)**3).to_value(u.Mpc**3)
    assert expected == pytest.approx(actual)


def test_cell_max():
    values = np.random.uniform(size=(5, 6, 7))
    maxima = cell_max(values)
    np.testing.assert_array_equal(maxima.shape, np.asarray(values.shape) - 1)
    assert np.all(maxima >= values[+1:, +1:, +1:])
    assert np.all(maxima >= values[+1:, +1:, :-1])
    assert np.all(maxima >= values[+1:, :-1, +1:])
    assert np.all(maxima >= values[+1:, :-1, :-1])
    assert np.all(maxima >= values[:-1, +1:, +1:])
    assert np.all(maxima >= values[:-1, +1:, :-1])
    assert np.all(maxima >= values[:-1, :-1, +1:])
    assert np.all(maxima >= values[:-1, :-1, :-1])
