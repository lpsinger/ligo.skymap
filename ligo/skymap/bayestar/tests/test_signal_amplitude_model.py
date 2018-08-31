import lal
import lalinspiral
import lalmetaio
import lalsimulation
import numpy as np
from pytest import approx, mark

from ...core import signal_amplitude_model


def get_eff_dist(detector, ra, dec, inclination, polarization, epoch, gmst):
    sim_inspiral = lalmetaio.SimInspiralTable()
    sim_inspiral.distance = 1
    sim_inspiral.longitude = ra
    sim_inspiral.latitude = dec
    sim_inspiral.inclination = inclination
    sim_inspiral.polarization = polarization
    sim_inspiral.geocent_end_time = epoch
    sim_inspiral.end_time_gmst = gmst
    _ = lal.LIGOTimeGPS()
    _, eff_dist = lalinspiral.InspiralSiteTimeAndDist(
        sim_inspiral, detector, _)
    return eff_dist


def get_complex_antenna(response, ra, dec, gmst):
    Fplus, Fcross = lal.ComputeDetAMResponse(response, ra, dec, 0, gmst)
    return Fplus + 1j * Fcross


@mark.parametrize('ra', np.arange(-np.pi, 1.1 * np.pi, 0.8 * np.pi))
@mark.parametrize('dec', np.arange(-np.pi, 1.1 * np.pi, 0.8 * np.pi))
@mark.parametrize('inclination', np.arange(-np.pi, 1.1 * np.pi, 0.8 * np.pi))
@mark.parametrize('polarization', np.arange(-np.pi, 1.1 * np.pi, 0.8 * np.pi))
@mark.parametrize('epoch', np.arange(1000000000.0, 1000086400.0, 14400.0))
@mark.parametrize('instrument', ['H1'])
def test_bayestar_signal_amplitude_model(ra, dec, inclination, polarization,
                                         epoch, instrument):
    """Test BAYESTAR signal amplitude model against LAL injection code."""
    detector = lalsimulation.DetectorPrefixToLALDetector(instrument)
    epoch = lal.LIGOTimeGPS(epoch)
    gmst = lal.GreenwichMeanSiderealTime(epoch)

    exp_i_twopsi = np.exp(2j * polarization)
    u = np.cos(inclination)
    u2 = np.square(u)
    F = get_complex_antenna(detector.response, ra, dec, gmst)
    result = signal_amplitude_model(F, exp_i_twopsi, u, u2)

    abs_expected = 1 / get_eff_dist(
        detector, ra, dec, inclination, polarization, epoch, gmst)

    # This is the *really* slow way of working out the signal amplitude:
    # generate a frequency-domain waveform and inject it.

    params = lal.CreateDict()
    lalsimulation.SimInspiralWaveformParamsInsertPNPhaseOrder(
        params, lalsimulation.PNORDER_NEWTONIAN)
    lalsimulation.SimInspiralWaveformParamsInsertPNAmplitudeOrder(
        params, lalsimulation.PNORDER_NEWTONIAN)

    # Calculate antenna factors
    Fplus, Fcross = lal.ComputeDetAMResponse(
        detector.response, ra, dec, polarization, gmst)

    # Check that the way I compute the antenna factors matches
    F = get_complex_antenna(detector.response, ra, dec, gmst)
    F *= np.exp(-2j * polarization)

    assert F.real == approx(Fplus, abs=4 * np.finfo(np.float64).eps)
    assert F.imag == approx(Fcross, abs=4 * np.finfo(np.float64).eps)

    # "Template" waveform with inclination angle of zero
    Htemplate, Hcross = lalsimulation.SimInspiralFD(
        1.4 * lal.MSUN_SI, 1.4 * lal.MSUN_SI,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 100, 101, 100, params,
        lalsimulation.TaylorF2)

    # Discard any non-quadrature phase component of "template"
    Htemplate.data.data += 1j * Hcross.data.data

    # Normalize "template"
    h = np.sum(np.square(Htemplate.data.data.real) +
               np.square(Htemplate.data.data.imag))
    h = 2 / h
    Htemplate.data.data *= h

    # "Signal" waveform with requested inclination angle
    Hsignal, Hcross = lalsimulation.SimInspiralFD(
        1.4 * lal.MSUN_SI, 1.4 * lal.MSUN_SI,
        0, 0, 0, 0, 0, 0, 1, inclination, 0, 0, 0, 0, 1, 100, 101, 100,
        params, lalsimulation.TaylorF2)

    # Project "signal" using antenna factors
    Hsignal.data.data = Fplus * Hsignal.data.data + Fcross * Hcross.data.data

    # Work out complex amplitude by comparing with "template" waveform
    expected = np.sum(Htemplate.data.data.conj() * Hsignal.data.data)

    # Test to nearly float (32-bit) precision because
    # lalinspiral.InspiralSiteTimeAndDist returns result as float.
    assert abs(expected) == approx(
        abs_expected, abs=1.5 * np.finfo(np.float32).eps)
    assert abs(result) == approx(
        abs_expected, abs=1.5 * np.finfo(np.float32).eps)
    assert result.real == approx(
        expected.real, abs=4 * np.finfo(np.float64).eps)
    assert result.imag == approx(
        expected.imag, abs=4 * np.finfo(np.float64).eps)
