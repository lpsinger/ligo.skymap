from astropy.cosmology import default_cosmology
from astropy import units as u
import lal
import lalsimulation
import numpy as np
import pytest

from ...bayestar.filter import (abs2, abscissa, get_f_lso, InterpolatedPSD,
                                signal_psd_series, sngl_inspiral_psd)
from ..bayestar_inject import (get_decisive_snr, get_max_comoving_distance,
                               get_snr_at_z, z_at_comoving_distance)


def test_get_decisive_snr():
    assert get_decisive_snr([]) == 0.0
    assert get_decisive_snr([1.0, 3.0, 2.0, 4.0]) == 3.0
    assert get_decisive_snr([4.0]) == 4.0


def get_snr_at_z_lalsimulation(cosmo, z, mass1, mass2, f_low, f_high, df, S):
    params = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(
        params, lalsimulation.PNORDER_NEWTONIAN)
    lalsimulation.SimInspiralWaveformParamsInsertPNAmplitudeOrder(
        params, lalsimulation.PNORDER_NEWTONIAN)
    lalsimulation.SimInspiralWaveformParamsInsertRedshift(
        params, z)

    # "Signal" waveform with requested inclination angle
    Hplus, Hcross = lalsimulation.SimInspiralFD(
        mass1 * lal.MSUN_SI, mass2 * lal.MSUN_SI,
        0, 0, 0, 0, 0, 0, cosmo.comoving_distance(z).to_value(u.m),
        0, 0, 0, 0, 0, df, f_low, f_high, f_low,
        params, lalsimulation.TaylorF2)

    # Force `plus' and `cross' waveform to be in quadrature.
    H = 0.5 * (Hplus.data.data + 1j * Hcross.data.data)

    # Throw away signal above merger.
    H[abscissa(Hplus) >= get_f_lso(mass1, mass2) / (1 + z)] = 0

    # Create output frequency series.
    signal_psd = lal.CreateREAL8FrequencySeries(
        'signal PSD', 0, Hplus.f0, Hplus.deltaF, Hplus.sampleUnits**2, len(H))
    signal_psd.data.data = abs2(H)
    HS = signal_psd_series(signal_psd, S)
    assert np.all(np.isfinite(HS.data.data))
    return np.sqrt(4 * np.trapz(HS.data.data, dx=HS.deltaF))


@pytest.mark.parametrize('z', [0.01, 0.1, 1.0, 2.0])
def test_get_snr_at_z(z):
    cosmo = default_cosmology.get_cosmology_from_string('WMAP9')
    f_low = 10
    f_high = 4096
    df = 1
    mass1 = mass2 = 1.4

    psd = lal.CreateREAL8FrequencySeries(
        '', 0, f_low, df, lal.DimensionlessUnit, (f_high - f_low) // df)
    lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087(psd, f_low)
    # Strip last sample because it's zero
    S = InterpolatedPSD(abscissa(psd)[:-1], psd.data.data[:-1])
    H = sngl_inspiral_psd('TaylorF2zeroPN', mass1, mass2)

    snr = get_snr_at_z(cosmo, [S], H.f0, H.deltaF, len(H.data.data),
                       InterpolatedPSD(abscissa(H), H.data.data, fill_value=0),
                       z)
    expected_snr = get_snr_at_z_lalsimulation(
        cosmo, z, mass1, mass2, f_low, f_high, df, S)

    assert snr == pytest.approx(expected_snr, rel=1e-3)


def test_get_max_comoving_distance():
    cosmo = default_cosmology.get_cosmology_from_string('WMAP9')
    f_low = 10
    f_high = 4096
    df = 1
    mass1 = mass2 = 1.4
    z = 1.0

    params = (mass1, mass2, 0.0, 0.0)

    psd = lal.CreateREAL8FrequencySeries(
        '', 0, f_low, df, lal.DimensionlessUnit, (f_high - f_low) // df)
    lalsimulation.SimNoisePSDaLIGODesignSensitivityP1200087(psd, f_low)
    # Strip last sample because it's zero
    S = InterpolatedPSD(abscissa(psd)[:-1], psd.data.data[:-1])

    H = sngl_inspiral_psd('TaylorF2zeroPN', mass1, mass2)
    snr = get_snr_at_z(cosmo, [S], H.f0, H.deltaF, len(H.data.data),
                       InterpolatedPSD(abscissa(H), H.data.data, fill_value=0),
                       z)

    comoving_distance = get_max_comoving_distance(
        cosmo, [S], 'TaylorF2zeroPN', f_low, snr, params)
    expected_comoving_distance = cosmo.comoving_distance(z).value

    assert comoving_distance == pytest.approx(expected_comoving_distance)


def test_z_at_comoving_distance():
    cosmo = default_cosmology.get_cosmology_from_string('WMAP9')
    expected_redshift = np.concatenate(([0], np.logspace(-8, 3), [np.inf]))
    comoving_distance = cosmo.comoving_distance(expected_redshift)
    redshift = z_at_comoving_distance(cosmo, comoving_distance)
    np.testing.assert_allclose(redshift, expected_redshift)
