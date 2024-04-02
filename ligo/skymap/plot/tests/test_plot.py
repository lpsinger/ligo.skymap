from itertools import chain, combinations, product
import platform

from astropy.coordinates import CartesianRepresentation, SkyCoord
import astropy_healpix as ah
import matplotlib
matplotlib.use('agg')
from astropy import units as u  # noqa: E402
import numpy as np  # noqa: E402
import healpy as hp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from ..bayes_factor import plot_bayes_factor  # noqa: E402
from ..marker import reticle  # noqa: E402

skip_if_macos_arm64 = pytest.mark.skipif(
    platform.system() == 'Darwin' and platform.machine() == 'arm64',
    reason='Tick labels vary on macOS arm64')


def pp_plot():
    # Re-initialize the random seed to make the unit test repeatable
    np.random.seed(0)
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111, projection='pp_plot')
    p_values = np.arange(1, 20) / 20
    return fig, ax, p_values


@pytest.fixture
def rcparams():
    with plt.rc_context({'figure.dpi': 72, 'savefig.dpi': 72}):
        yield


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=10)
def test_pp_plot_steps(rcparams):
    """Test P--P plot with drawstyle='steps'."""
    fig, ax, p_values = pp_plot()
    ax.add_confidence_band(len(p_values))
    ax.add_diagonal()
    ax.add_lightning(len(p_values), 20, drawstyle='steps')
    ax.add_series(p_values, drawstyle='steps')
    ax.add_worst(p_values)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=10)
def test_pp_plot_lines(rcparams):
    """Test P--P plot with drawstyle='steps'."""
    fig, ax, p_values = pp_plot()
    ax.add_confidence_band(len(p_values))
    ax.add_diagonal()
    ax.add_lightning(len(p_values), 20, drawstyle='lines')
    ax.add_series(p_values, drawstyle='lines')
    ax.add_diagonal()
    ax.add_series(p_values)
    ax.add_worst(p_values)
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=10)
def test_pp_plot_default(rcparams):
    """Test P--P plot with drawstyle='steps'."""
    fig, ax, p_values = pp_plot()
    ax.add_confidence_band(len(p_values))
    ax.add_diagonal()
    ax.add_lightning(len(p_values), 20)
    ax.add_series(p_values)
    ax.add_worst(p_values)
    return fig


@pytest.mark.parametrize('proj', ['aitoff', 'mollweide'])
@pytest.mark.parametrize('units', ['degrees', 'hours'])
@pytest.mark.parametrize('coordsys', [
    'astro',
    pytest.param('geo', marks=skip_if_macos_arm64),
    pytest.param('galactic', marks=skip_if_macos_arm64)])
@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_allsky_axes(rcparams, coordsys, units, proj):
    """Test projection of a HEALPix image onto allsky axes, either
    in celestial or earth-fixed coordinates.
    """
    # Set up axes. (The obstime has an effect only for geographic axes.)
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111, projection=coordsys + ' ' + units + ' ' + proj,
                         obstime='2017-08-17T12:41:04.444458')

    # Build a low-resolution example HEALPix sky map:
    # the value is equal to the right ascension.
    nside = 8
    npix = ah.nside_to_npix(nside)
    ra, dec = hp.pix2ang(nside, np.arange(npix), lonlat=True)
    img = np.sin(np.deg2rad(ra))

    # Plot, show grid, and return figure.
    ax.imshow_hpx((img, 'ICRS'))
    ax.grid()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_allsky_obstime(rcparams):
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111, projection='geo degrees mollweide',
                         obstime='2017-08-17T12:41:04.444458')
    ax.grid()
    return fig


@pytest.mark.skipif(
    not (platform.system() == 'Linux' and platform.machine == 'x86_64'),
    reason='Tick label positions vary on different operating systems')
@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_globe_axes(rcparams):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], projection='astro globe',
                      center='197.45d -23.38d')
    ax.grid()
    return fig


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_zoom_axes(rcparams):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.6], projection='astro zoom',
                      center='197.45d -23.38d', radius='90 arcmin')
    ax.scalebar((0.1, 0.1), 30 * u.arcmin)
    ax.grid()
    for key in ['ra', 'dec']:
        ax.coords[key].set_auto_axislabel(False)
    return fig


@pytest.mark.parametrize(
    'projection', [
        'astro zoom',
        'astro globe',
        'astro mollweide',
        'astro aitoff'
    ]
)
def test_center_cartesian(projection):
    """Test that zoom axes accept coordinates in other representations."""
    fig = plt.figure()
    center = SkyCoord(0, 0, 1, representation_type=CartesianRepresentation)
    fig.add_axes([0.2, 0.2, 0.6, 0.6], projection=projection, center=center)


@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_reticle():
    which_list = [''.join(d) for d in
                  chain.from_iterable(combinations('lrtb', n)
                  for n in range(2, 5))]
    inners = [0.0, 0.2, 0.4]
    outers = [0.8, 0.9, 1.0]
    angles = [0.0, 22.5, 45.0]
    args_list = list(product(inners, outers, angles))

    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    for args, x in zip(args_list, np.linspace(0.1, 0.9, len(args_list))):
        for which, y in zip(which_list,
                            np.linspace(0.1, 0.9, len(which_list))):
            ax.plot(x, y, marker=reticle(*args, which=which))
    return fig


@pytest.mark.parametrize('proj', ['aitoff', 'mollweide'])
@pytest.mark.parametrize('cent', ['197.45d -23.38d', None])
@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_center_projections(rcparams, proj, cent):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_axes(111, projection=f'astro {proj}',
                      center=cent)
    ax.grid()
    return fig


@pytest.mark.mpl_image_compare(tolerance=1.5)
def test_plot_bayes_factor():
    fig, ax = plot_bayes_factor(6.3, title='BAYESTAR is awesome')
    return fig
