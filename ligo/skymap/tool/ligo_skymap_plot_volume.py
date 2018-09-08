# -*- coding: utf-8 -*-
#
# Copyright (C) 2015-2018  Leo Singer
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Plot a volumetric posterior in three-projection view.
"""


from . import ArgumentParser, FileType, figure_parser


def parser():
    parser = ArgumentParser(parents=[figure_parser])
    parser.add_argument(
        '--annotate', default=False, action='store_true',
        help='annotate plot with information about the event')
    parser.add_argument(
        '--max-distance', metavar='Mpc', type=float,
        help='maximum distance of plot in Mpc')
    parser.add_argument(
        '--contour', metavar='PERCENT', type=float, nargs='+',
        help='plot contour enclosing this percentage of'
        ' probability mass')
    parser.add_argument(
        '--radecdist', nargs=3, type=float, action='append', default=[],
        help='right ascension (deg), declination (deg), and distance to mark')
    parser.add_argument(
        '--chain', metavar='CHAIN.hdf5', type=FileType('rb'),
        help='optionally plot a posterior sample chain')
    parser.add_argument(
        '--projection', type=int, choices=list(range(4)), default=0,
        help='Plot one specific projection, or 0 for all projections')
    parser.add_argument(
        'input', metavar='INPUT.fits[.gz]', type=FileType('rb'),
        default='-', nargs='?', help='Input FITS file')
    parser.add_argument(
        '--align-to', metavar='SKYMAP.fits[.gz]', type=FileType('rb'),
        help='Align to the principal axes of this sky map')
    parser.set_defaults(figure_width='3.5', figure_height='3.5')
    return parser


def main(args=None):
    opts = parser().parse_args(args)

    # Create progress bar.
    from glue.text_progress_bar import ProgressBar
    progress = ProgressBar()
    progress.update(-1, 'Starting up')

    # Late imports
    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    from matplotlib import transforms
    from .. import io
    from ..plot import marker
    from ..distance import (parameters_to_marginal_moments, principal_axes,
                            volume_render, marginal_pdf)
    import healpy as hp
    import numpy as np
    import scipy.stats
    import seaborn

    # Read input, determine input resolution.
    progress.update(-1, 'Loading FITS file')
    (prob, mu, sigma, norm), metadata = io.read_sky_map(
        opts.input.name, distances=True)
    npix = len(prob)
    nside = hp.npix2nside(npix)

    progress.update(-1, 'Preparing projection')

    if opts.align_to is None or opts.input.name == opts.align_to.name:
        prob2, mu2, sigma2, norm2 = prob, mu, sigma, norm
    else:
        (prob2, mu2, sigma2, norm2), _ = io.read_sky_map(
            opts.align_to.name, distances=True)
    if opts.max_distance is None:
        mean, std = parameters_to_marginal_moments(prob2, mu2, sigma2)
        max_distance = mean + 2.5 * std
    else:
        max_distance = opts.max_distance
    R = np.ascontiguousarray(principal_axes(prob2, mu2, sigma2))

    if opts.chain:
        chain = io.read_samples(opts.chain.name)
        chain = np.dot(R.T, (hp.ang2vec(
            0.5 * np.pi - chain['dec'], chain['ra']) *
            np.atleast_2d(chain['dist']).T).T)

    fig = plt.figure(frameon=False)
    n = 1 if opts.projection else 2
    gs = gridspec.GridSpec(
        n, n, left=0.01, right=0.99, bottom=0.01, top=0.99,
        wspace=0.05, hspace=0.05)

    imgwidth = int(opts.dpi * opts.figure_width / n)
    s = np.linspace(-max_distance, max_distance, imgwidth)
    xx, yy = np.meshgrid(s, s)

    # Color palette for markers
    colors = seaborn.color_palette(n_colors=len(opts.radecdist) + 1)

    truth_marker = marker.reticle(
        inner=0.5 * np.sqrt(2), outer=1.5 * np.sqrt(2), angle=45)

    for iface, (axis0, axis1, (sp0, sp1)) in enumerate((
            (1, 0, [0, 0]),
            (0, 2, [1, 1]),
            (1, 2, [1, 0]),)):

        if opts.projection and opts.projection != iface + 1:
            continue

        progress.update(text='Plotting projection {0}'.format(iface + 1))

        # Marginalize onto the given face
        density = volume_render(
            xx.ravel(), yy.ravel(), max_distance, axis0, axis1, R, False,
            prob, mu, sigma, norm).reshape(xx.shape)

        # Plot heat map
        ax = fig.add_subplot(
            gs[0, 0] if opts.projection else gs[sp0, sp1], aspect=1)
        ax.imshow(
            density, origin='lower',
            extent=[-max_distance, max_distance, -max_distance, max_distance],
            cmap=opts.colormap)

        # Add contours if requested
        if opts.contour:
            flattened_density = density.ravel()
            indices = np.argsort(flattened_density)[::-1]
            cumsum = np.empty_like(flattened_density)
            cs = np.cumsum(flattened_density[indices])
            cumsum[indices] = cs / cs[-1] * 100
            cumsum = np.reshape(cumsum, density.shape)
            u, v = np.meshgrid(s, s)
            contourset = ax.contour(
                u, v, cumsum, levels=opts.contour, linewidths=0.5)

        # Mark locations
        for (ra, dec, dist), color in zip(opts.radecdist, colors[1:]):
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            xyz = np.dot(R.T, hp.ang2vec(theta, phi) * dist)
            ax.plot(
                xyz[axis0], xyz[axis1], marker=truth_marker,
                markeredgecolor=color, markerfacecolor='none',
                markeredgewidth=1)

        # Plot chain
        if opts.chain:
            ax.plot(chain[axis0], chain[axis1], '.k', markersize=0.5)

        # Hide axes ticks
        ax.set_xticks([])
        ax.set_yticks([])

        # Set axis limits
        ax.set_xlim([-max_distance, max_distance])
        ax.set_ylim([-max_distance, max_distance])

        # Mark origin (Earth)
        ax.plot(
            [0], [0], marker=marker.earth, markersize=5,
            markerfacecolor='none', markeredgecolor='black',
            markeredgewidth=0.75)

        if iface == 2:
            ax.invert_xaxis()

    # Add contour labels if contours requested
    if opts.contour:
        ax.clabel(contourset, fmt='%d%%', fontsize=7)

    if not opts.projection:
        # Add scale bar, 1/4 width of the plot
        ax.plot(
            [0.0625, 0.3125], [0.0625, 0.0625],
            color='black', linewidth=1, transform=ax.transAxes)
        ax.text(
            0.0625, 0.0625,
            '{0:d} Mpc'.format(int(np.round(0.5 * max_distance))),
            fontsize=8, transform=ax.transAxes, verticalalignment='bottom')

        # Create marginal distance plot.
        progress.update(-1, 'Plotting distance')
        gs1 = gridspec.GridSpecFromSubplotSpec(5, 5, gs[0, 1])
        ax = fig.add_subplot(gs1[1:-1, 1:-1])

        # Plot marginal distance distribution, integrated over the whole sky.
        d = np.linspace(0, max_distance)
        ax.fill_between(d, marginal_pdf(d, prob, mu, sigma, norm),
                        alpha=0.5, color=colors[0])

        # Plot conditional distance distribution at true position
        # and mark true distance.
        for (ra, dec, dist), color in zip(opts.radecdist, colors[1:]):
            theta = 0.5 * np.pi - np.deg2rad(dec)
            phi = np.deg2rad(ra)
            ipix = hp.ang2pix(nside, theta, phi)
            ax.fill_between(d, scipy.stats.norm(
                mu[ipix], sigma[ipix]).pdf(d) * norm[ipix] * np.square(d),
                alpha=0.5, color=color)
            ax.axvline(dist, color='black', linewidth=0.5)
            ax.plot(
                [dist], [-0.15], marker=truth_marker, markeredgecolor=color,
                markerfacecolor='none', markeredgewidth=1, clip_on=False,
                transform=transforms.blended_transform_factory(
                    ax.transData, ax.transAxes))
            ax.axvline(dist, color='black', linewidth=0.5)

        # Scale axes
        ax.set_xticks([0, max_distance])
        ax.set_xticklabels(
            ['0', "{0:d}\nMpc".format(int(np.round(max_distance)))],
            fontsize=9)
        ax.set_yticks([])
        ax.set_xlim(0, max_distance)
        ax.set_ylim(0, ax.get_ylim()[1])

        if opts.annotate:
            text = []
            try:
                objid = metadata['objid']
            except KeyError:
                pass
            else:
                text.append('event ID: {}'.format(objid))
            try:
                distmean = metadata['distmean']
                diststd = metadata['diststd']
            except KeyError:
                pass
            else:
                text.append(u'distance: {}±{} Mpc'.format(
                            int(np.round(distmean)), int(np.round(diststd))))
            ax.text(0, 1, '\n'.join(text), transform=ax.transAxes, fontsize=7,
                    ha='left', va='bottom', clip_on=False)

    progress.update(-1, 'Saving')
    opts.output()
