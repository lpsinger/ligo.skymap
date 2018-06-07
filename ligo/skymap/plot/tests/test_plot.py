from itertools import chain, combinations, product

import matplotlib
matplotlib.use('agg')
from astropy import units as u  # noqa: E402
import numpy as np  # noqa: E402
import healpy as hp  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from ..marker import reticle  # noqa: E402


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
@pytest.mark.parametrize('coordsys', ['astro', 'geo'])
@pytest.mark.mpl_image_compare(remove_text=True, tolerance=1.5)
def test_allsky_axes(rcparams, coordsys, units, proj):
    """Test projection of a HEALPix image onto allsky axes, either
    in celestial or earth-fixed coordinates."""
    # Set up axes. (The obstime has an effect only for geographic axes.)
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111, projection=coordsys + ' ' + units + ' ' + proj,
                         obstime='2017-08-17T12:41:04.444458')

    # Build a low-resolution example HEALPix sky map:
    # the value is equal to the right ascension.
    nside = 8
    npix = hp.nside2npix(nside)
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
    return fig


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
