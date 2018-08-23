import numpy as np
import pytest

from ... import io
from . import run_entry_point


@pytest.mark.parametrize('geo', [False, True])
def test_plot(tmpdir, geo):
    """Test ligo-skymap-plot."""
    fitsfilename = str(tmpdir / 'skymap.fits')
    pngfilename = str(tmpdir / 'skymap.png')
    io.write_sky_map(fitsfilename, np.arange(12) / 12,
                     gps_time=1e9, objid='foobar')
    args = ['ligo-skymap-plot', fitsfilename, '-o', pngfilename, '--annotate',
            '--colorbar', '--radec', '320', '45']
    if geo:
        args.append('--geo')
    args.extend(['--contour', '90'])
    run_entry_point(*args)
