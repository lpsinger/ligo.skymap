import os
import subprocess

import pytest

import numpy as np
import healpy as hp

from . import run_entry_point
from ... import io
from ... import distance


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


def test_combine(tmpdir):
    fn1 = str(tmpdir / 'skymap1.fits.gz')
    fn2 = str(tmpdir / 'skymap2.fits.gz')
    fn3 = str(tmpdir / 'joint_skymap.fits.gz')

    # generate a hemisphere of constant probability
    nside1 = 32
    npix1 = hp.nside2npix(nside1)
    m1 = np.zeros(npix1)
    disc_idx = hp.query_disc(nside1, (1, 0, 0), np.pi/2)
    m1[disc_idx] = 1
    m1 /= m1.sum()
    hp.write_map(fn1, m1, column_names=['PROBABILITY'],
                 extra_header=[('INSTRUME', 'X1')])

    # generate another hemisphere of constant probability
    # but with higher resolution and rotated 90 degrees
    nside2 = 64
    npix2 = hp.nside2npix(nside2)
    m2 = np.zeros(npix2)
    disc_idx = hp.query_disc(nside2, (0, 1, 0), np.pi/2)
    m2[disc_idx] = 1
    m2 /= m2.sum()
    hp.write_map(fn2, m2, column_names=['PROBABILITY'],
                 extra_header=[('INSTRUME', 'Y1')])

    run_entry_point('ligo-skymap-combine', fn1, fn2, fn3)

    m3 = hp.read_map(fn3, nest=True)
    npix3 = len(m3)
    nside3 = hp.npix2nside(npix3)
    pix_area3 = hp.nside2pixarea(nside3)

    # resolution must match the highest original resolution
    assert npix3 == npix2
    # probability must be normalized to 1
    assert m3.sum() == pytest.approx(1)
    # support must be ¼ of the sphere
    tolerance = 10 * hp.nside2pixarea(nside1)
    assert sum(m3 > 0) * pix_area3 == pytest.approx(np.pi, abs=tolerance)

    # generate a BAYESTAR-like map with mock distance information
    d_mu = np.zeros_like(m1)
    d_sigma = np.ones_like(m1)
    d_norm = np.ones_like(m1)
    io.write_sky_map(fn1, np.vstack((m1, d_mu, d_sigma, d_norm)).T)

    run_entry_point('ligo-skymap-combine', fn1, fn2, fn3)

    m3, meta3 = io.read_sky_map(fn3, nest=True, distances=True)

    # check that marginal distance moments match what was simulated
    mean, std, _ = distance.parameters_to_moments(d_mu[0], d_sigma[0])
    assert meta3['distmean'] == pytest.approx(mean)
    assert meta3['diststd'] == pytest.approx(std)
