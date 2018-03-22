import os
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


@pytest.fixture
def inj_coinc(inj, coinc, tmpdir):
    filename = str(tmpdir / 'inj_coinc.xml')
    subprocess.check_call(['ligolw_add', inj, coinc, '-o', filename])
    subprocess.check_call(['lalapps_inspinjfind', filename])
    return filename


@pytest.fixture
def inj_coinc_sqlite(inj_coinc, tmpdir):
    filename = str(tmpdir / 'inj_coinc.sqlite')
    subprocess.check_call(['ligolw_sqlite', inj_coinc, '-p', '-d', filename])
    return filename


@pytest.fixture
def localize_coincs(coinc, psd, tmpdir):
    run_entry_point('bayestar-localize-coincs', coinc, '-o', str(tmpdir))
    return tmpdir


# FIXME: Skip until
# https://git.ligo.org/lscsoft/lalsuite/merge_requests/192
# is fixed in a release
@pytest.mark.skip
def test_stats(inj_coinc_sqlite, localize_coincs, tmpdir):
    filename = str(tmpdir / 'bayestar.out')
    run_entry_point('ligo-skymap-stats', '-o', filename,
                    inj_coinc_sqlite, os.path.join(localize_coincs, '*.fits'))
