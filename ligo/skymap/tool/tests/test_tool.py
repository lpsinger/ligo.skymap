import subprocess

import pytest

from . import run_entry_point


@pytest.fixture
def psd(tmpdir):
    filename = str(tmpdir / 'psd.xml')
    run_entry_point('bayestar-sample-model-psd',
                    '-o', filename,
                    '--H1=aLIGOZeroDetHighPower',
                    '--L1=aLIGOZeroDetHighPower',
                    '--V1=AdvVirgo')
    return filename


@pytest.fixture
def inj(tmpdir):
    filename = str(tmpdir / 'inj.xml')
    subprocess.check_call(['lalapps_inspinj',
                           '-o', filename,
                           '--m-distr=fixMasses',
                           '--fixed-mass1=1.4',
                           '--fixed-mass2=1.4',
                           '--gps-start-time=1000000000',
                           '--gps-end-time=1000000001',
                           '--time-step=1',
                           '--l-distr=random',
                           '--i-distr=uniform',
                           '--d-distr=uniform',
                           '--min-distance=1',
                           '--max-distance=200e3',
                           '--waveform=TaylorF2threePointFivePN',
                           '--disable-spin',
                           '--f-lower=10'])
    return filename


@pytest.fixture
def coinc(inj, psd, tmpdir):
    filename = str(tmpdir / 'coinc.xml')
    run_entry_point('bayestar-realize-coincs',
                    '--reference-psd', psd, '-o', filename, inj,
                    '--enable-snr-series', '--detector', 'H1', 'L1', 'V1')
    return filename


def test_bayestar_localize_coincs(coinc, psd, tmpdir):
    run_entry_point('bayestar-localize-coincs', coinc, '-o', str(tmpdir))
