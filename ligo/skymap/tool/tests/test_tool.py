import pytest

from . import run_entry_point


def test_bayestar_sample_model_psd(tmpdir):
    filename = str(tmpdir / 'psd.xml')
    run_entry_point('bayestar-sample-model-psd',
                    '--H1', 'aLIGOZeroDetHighPower',
                    '--L1', 'aLIGOZeroDetHighPower',
                    '--V1', 'AdvVirgo', '-o', filename)
