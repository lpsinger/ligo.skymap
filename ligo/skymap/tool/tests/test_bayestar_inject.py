from astropy.cosmology import default_cosmology
from astropy import units as u
import lal
import lalsimulation
import numpy as np
import pytest

from ...bayestar.filter import (abs2, abscissa, get_f_lso, InterpolatedPSD,
                                signal_psd_series)
from ..bayestar_inject import (get_decisive_snr, z_at_snr,
                               z_at_comoving_distance, cell_max)


def test_get_decisive_snr():
    assert get_decisive_snr([]) == 0.0
    assert get_decisive_snr([1.0, 3.0, 2.0, 4.0]) == 3.0
    assert get_decisive_snr([4.0]) == 4.0


def get_snr_at_z_lalsimulation(cosmo, z, mass1, mass2, f_low, f_high, psd):
    params = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertRedshift(
        params, z)

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
    cosmo = default_cosmology.get_cosmology_from_string('WMAP9')
    f_low = 10
    f_high = 4096
    df = 0.1
    mass1 = mass2 = 0.5 * mtotal

    psd = lal.CreateREAL8FrequencySeries(
        '', 0, f_low, df, lal.DimensionlessUnit, int((f_high - f_low) // df))
    lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087(psd, f_low)

    snr = get_snr_at_z_lalsimulation(cosmo, z, mass1, mass2, f_low, f_high, psd)
    z_solution = z_at_snr(
        cosmo, [psd], 'IMRPhenomPv2', f_low, snr, (mass1, mass2, 0, 0))

    assert z_solution == pytest.approx(z, rel=1e-2)


def test_z_at_comoving_distance():
    cosmo = default_cosmology.get_cosmology_from_string('WMAP9')
    expected_redshift = np.concatenate(([0], np.logspace(-8, 3), [np.inf]))
    comoving_distance = cosmo.comoving_distance(expected_redshift)
    redshift = z_at_comoving_distance(cosmo, comoving_distance)
    np.testing.assert_allclose(redshift, expected_redshift)


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
