from astropy.table import Table
import numpy as np
import pytest

from . import run_entry_point


@pytest.mark.parametrize('ndatasets', [1, 2])
@pytest.mark.parametrize('groupby', [None, 'far', 'snr'])
def test_plot_stats(tmpdir, groupby, ndatasets):
    """Test ligo-skymap-plot-stats."""
    filenames = [str(tmpdir / '{}.out'.format(i)) for i in range(ndatasets)]
    n = 250

    for filename in filenames:
        Table({
            'searched_area': np.random.uniform(0, 1e4, size=n),
            'p_value': np.random.uniform(0, 1, size=n),
            'snr': np.random.uniform(10, 12, size=n),
            'far': 10**np.random.uniform(-12, -10, size=n)
        }).write(
            filename, format='ascii.tab'
        )

    args = ['ligo-skymap-plot-stats',
            '--format', 'png', '-o', str(tmpdir)] + filenames
    if groupby is not None:
        args.extend(['--group-by', groupby])
    run_entry_point(*args)


def test_plot_stats_missing(tmpdir):
    """Test ligo-skymap-plot-stats, with group by value missing."""
    filename = str(tmpdir / '0.out')
    n = 250

    Table({
        'p_value': np.random.uniform(0, 1, size=n),
    }).write(
        filename, format='ascii.tab'
    )

    args = ['ligo-skymap-plot-stats',
            '--group-by', 'snr', '-o', str(tmpdir), filename]
    with pytest.raises(RuntimeError) as excinfo:
        run_entry_point(*args)
    assert 'The following files had no "snr" column' in str(excinfo.value)


def test_plot_stats_invalid(tmpdir):
    """Test ligo-skymap-plot-stats, with group by value invalid."""
    filename = str(tmpdir / '0.out')
    n = 250

    Table({
        'p_value': np.random.uniform(0, 1, size=n),
        'snr': np.tile(np.nan, n)
    }).write(
        filename, format='ascii.tab'
    )

    args = ['ligo-skymap-plot-stats',
            '--group-by', 'snr', '-o', str(tmpdir), filename]
    with pytest.raises(RuntimeError) as excinfo:
        run_entry_point(*args)
    assert ('The following files had invalid values in the "snr" column'
            in str(excinfo.value))
