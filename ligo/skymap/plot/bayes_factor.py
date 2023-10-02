#
# Copyright (C) 2019-2023  Leo Singer
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
"""Bullet charts for Bayes factors."""

from matplotlib import pyplot as plt
import numpy as np

__all__ = ('plot_bayes_factor',)


def plot_bayes_factor(logb,
                      values=(1, 3, 5),
                      labels=('', 'strong', 'very strong'),
                      xlim=7, title=None, palette='RdYlBu',
                      var_label="B"):
    """Visualize a Bayes factor as a `bullet graph`_.

    Make a bar chart of a log Bayes factor as compared to a set of subjective
    threshold values. By default, use the thresholds from
    Kass & Raftery (1995).

    .. _`bullet graph`: https://en.wikipedia.org/wiki/Bullet_graph

    Parameters
    ----------
    logb : float
        The natural logarithm of the Bayes factor.
    values : list
        A list of floating point values for human-friendly confidence levels.
    labels : list
        A list of string labels for human-friendly confidence levels.
    xlim : float
        Limits of plot (`-xlim` to `+xlim`).
    title : str
        Title for plot.
    palette : str
        Color palette.
    var_label : str
        The variable symbol used in plotting

    Returns
    -------
    fig : Matplotlib figure
    ax : Matplotlib axes

    Examples
    --------

    .. plot::
       :include-source:

        from ligo.skymap.plot.bayes_factor import plot_bayes_factor
        plot_bayes_factor(6.3, title='BAYESTAR is awesome')

    """
    with plt.style.context('seaborn-v0_8-notebook'):
        fig, ax = plt.subplots(figsize=(6, 1.7), tight_layout=True)
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(-0.5, 0.5)
        ax.set_yticks([])
        ax.set_title(title)
        ax.set_ylabel(r'$\ln\,{}$'.format(var_label), rotation=0,
                      rotation_mode='anchor',
                      ha='right', va='center')

        # Add human-friendly labels
        ticks = (*(-x for x in reversed(values)), 0, *values)
        ticklabels = (
            *(f'{s}\nevidence\nagainst'.strip() for s in reversed(labels)), '',
            *(f'{s}\nevidence\nfor'.strip() for s in labels))
        ax.set_xticks(ticks)
        ax.set_xticklabels(ticklabels)
        plt.setp(ax.get_xticklines(), visible=False)
        plt.setp(ax.get_xticklabels()[:len(ticks) // 2], ha='right')
        plt.setp(ax.get_xticklabels()[len(ticks) // 2:], ha='left')

        # Plot colored bands for confidence thresholds
        fmt = plt.FuncFormatter(lambda x, _: f'{x:+g}'.replace('+0', '0'))
        ax2 = ax.twiny()
        ax2.set_xlim(*ax.get_xlim())
        ax2.set_xticks(ticks)
        ax2.xaxis.set_major_formatter(fmt)
        levels = (-xlim, *ticks, xlim)
        colors = plt.get_cmap(palette)(np.arange(1, len(levels)) / len(levels))
        ax.barh(0, np.diff(levels), 1, levels[:-1],
                linewidth=plt.rcParams['xtick.major.width'],
                color=colors, edgecolor='white')

        # Plot bar for log Bayes factor value
        ax.barh(0, logb, 0.5, color='black',
                linewidth=plt.rcParams['xtick.major.width'],
                edgecolor='white')

        for ax_ in fig.axes:
            ax_.grid(False)
            for spine in ax_.spines.values():
                spine.set_visible(False)

    return fig, ax
