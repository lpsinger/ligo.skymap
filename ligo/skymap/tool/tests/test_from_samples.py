from astropy.coordinates import (CartesianRepresentation, SkyCoord,
                                 SphericalRepresentation)
from astropy.table import Table
from astropy.utils.misc import NumpyRNGContext
from astropy import units as u
import numpy as np
from scipy import stats
import pytest

from ...io.hdf5 import write_samples
from . import run_entry_point


@pytest.fixture
def seed():
    with NumpyRNGContext(0):
        yield


@pytest.fixture
def samples(seed, tmpdir):
    mean = SkyCoord(ra=stats.uniform(0, 360).rvs() * u.hourangle,
                    dec=np.arcsin(stats.uniform(-1, 1).rvs()) * u.radian,
                    distance=stats.uniform(100, 200).rvs()).cartesian.xyz.value
    eigvals = stats.uniform(0, 1).rvs(3)
    eigvals *= len(eigvals) / eigvals.sum()
    cov = stats.random_correlation.rvs(eigvals) * 100
    pts = stats.multivariate_normal(mean, cov).rvs(200)
    pts = SkyCoord(pts, representation_type=CartesianRepresentation)
    pts.representation_type = SphericalRepresentation
    time = stats.uniform(-0.01, 0.01).rvs(200) + 1e9
    table = Table({
        'ra': pts.ra.rad,
        'dec': pts.dec.rad,
        'distance': pts.distance.value,
        'time': time
    })
    filename = str(tmpdir / 'samples.hdf5')
    write_samples(table, filename, path='/posterior_samples')
    return filename


def test_from_samples(samples, tmpdir):
    """Test ligo-skyamp-from-samples."""
    run_entry_point('ligo-skymap-from-samples', '--seed', '150914',
                    samples, '-o', str(tmpdir))
