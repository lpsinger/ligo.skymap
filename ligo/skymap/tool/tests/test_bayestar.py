from glue.ligolw import lsctables as lsctables
from glue.ligolw import utils as ligolw_utils
from lalinspiral.thinca import InspiralCoincDef
import numpy as np
import pytest

from ... import io
from ...io.events.ligolw import ContentHandler
from . import run_entry_point, run_glue, run_lalsuite


@pytest.fixture
def psd(tmpdir):
    filename = str(tmpdir / 'psd.xml')
    run_entry_point('bayestar-sample-model-psd',
                    '-o', filename,
                    '--H1=aLIGOZeroDetHighPower',
                    '--L1=aLIGOZeroDetHighPower',
                    '--V1=AdVDesignSensitivityP1200087')
    return filename


@pytest.fixture
def inj(tmpdir):
    filename = str(tmpdir / 'inj.xml')
    run_lalsuite('lalapps_inspinj',
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
                 '--f-lower=10')
    return filename


@pytest.fixture
def coinc(inj, psd, tmpdir):
    filename = str(tmpdir / 'coinc.xml')
    run_entry_point('bayestar-realize-coincs',
                    '--measurement-error', 'gaussian-noise', '--f-low', '25',
                    '--reference-psd', psd, '-o', filename, inj,
                    '--enable-snr-series', '--detector', 'H1', 'L1', 'V1')
    return filename


@pytest.fixture
def coinc_without_inj(coinc, tmpdir):
    """Produce a coinc.xml file with the found coincs stripped out."""
    filename = str(tmpdir / 'coinc_without_inj.xml')
    xmldoc = ligolw_utils.load_filename(coinc, contenthandler=ContentHandler)

    # Prune coinc_def table
    coinc_def_table = lsctables.CoincDefTable.get_table(xmldoc)
    included = [row for row in coinc_def_table
                if row.search_coinc_type == InspiralCoincDef.search_coinc_type
                and row.search == InspiralCoincDef.search]
    included_coinc_def_ids = {row.coinc_def_id for row in included}
    coinc_def_table[:] = included

    # Prune coinc table
    coinc_table = lsctables.CoincTable.get_table(xmldoc)
    included = [row for row in coinc_table
                if row.coinc_def_id in included_coinc_def_ids]
    included_coinc_ids = {row.coinc_event_id for row in included}
    coinc_table[:] = included

    # Prune coinc_map table
    coinc_map_table = lsctables.CoincMapTable.get_table(xmldoc)
    coinc_map_table[:] = [row for row in coinc_map_table
                          if row.coinc_event_id in included_coinc_ids]

    ligolw_utils.write_filename(xmldoc, filename)
    return filename


@pytest.fixture
def coinc_sqlite(coinc, tmpdir):
    filename = str(tmpdir / 'coinc.sqlite')
    run_glue('ligolw_sqlite', coinc, '--ilwdchar-compat',
             '-p', '-r', '-d', filename)
    return filename


@pytest.fixture
def localize_coincs(coinc, psd, tmpdir):
    run_entry_point('bayestar-localize-coincs', coinc, '-o', str(tmpdir))
    return str(tmpdir / '0.fits')


# Note: any test that uses this fixture should be marked with
# @pytest.mark.internet_off to make sure that it does not actually
# contact GraceDb.
@pytest.fixture
def localize_lvalert(coinc_without_inj, psd, tmpdir, monkeypatch):

    class MockGraceDb:

        def __init__(self, service_url, *args, **kwargs):
            self._service_url = service_url

        def files(self, graceid, filename):
            assert graceid == 'G1234'
            mock_filename = {'coinc.xml': coinc_without_inj,
                             'psd.xml.gz': psd}[filename]
            return open(mock_filename, 'rb')

        def writeLog(self, *args, **kwargs):
            pass

    monkeypatch.setattr('ligo.gracedb.rest.GraceDb', MockGraceDb)

    filename = str(tmpdir / 'bayestar.fits')
    run_entry_point('bayestar-localize-lvalert', 'G1234', '-N', '-o', filename)
    return filename


@pytest.mark.internet_off
def test_bayestar(localize_coincs, localize_lvalert, coinc_sqlite, tmpdir):
    """Test bayestar-realize-coincs, bayestar-localize-coincs,
    bayestar-localize-lvalert, and ligo-skymap-stats."""
    # Check that bayestar-localize-coincs and bayestar-localize-lvalert
    # produce the same output.
    skymap1, meta1 = io.read_sky_map(localize_coincs, distances=True)
    skymap2, meta2 = io.read_sky_map(localize_lvalert, distances=True)
    for col1, col2 in zip(skymap1, skymap2):
        np.testing.assert_allclose(col1, col2)
    for key in 'gps_time origin vcs_version vcs_revision build_date'.split():
        assert meta1[key] == meta2[key]
    for key in 'distmean diststd log_bci log_bsn'.split():
        np.testing.assert_allclose(meta1[key], meta2[key])

    # Test ligo-skymap-stats.
    out1 = str(tmpdir / 'stats1.out')
    out2 = str(tmpdir / 'stats2.out')
    args = ('ligo-skymap-stats', '--modes', '-p', '90', '-a', '100', '-o')
    run_entry_point(*args, out1, localize_coincs, '-d', coinc_sqlite)
    run_entry_point(*args, out2, localize_lvalert, '-j', '2')
