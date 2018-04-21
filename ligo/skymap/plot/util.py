#
# Copyright (C) 2012-2018  Leo Singer
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
Miscellaneous plotting utilities
"""
import matplotlib
from matplotlib import text
from matplotlib import ticker
from matplotlib import patheffects

__all__ = ('colorbar', 'outline_text')


def colorbar(*args):
    from matplotlib import pyplot as plt

    usetex = matplotlib.rcParams['text.usetex']
    locator = ticker.AutoLocator()
    formatter = ticker.ScalarFormatter(useMathText=not usetex)
    formatter.set_scientific(True)
    formatter.set_powerlimits((1e-1, 100))

    # Plot colorbar
    cb = plt.colorbar(*args,
                      orientation='horizontal', shrink=0.4,
                      ticks=locator, format=formatter)

    if cb.orientation == 'vertical':
        axis = cb.ax.yaxis
    else:
        axis = cb.ax.xaxis

    # Move order of magnitude text into last label.
    ticklabels = [label.get_text() for label in axis.get_ticklabels()]
    # Avoid putting two '$' next to each other if we are in tex mode.
    if usetex:
        fmt = '{{{0}}}{{{1}}}'
    else:
        fmt = u'{0}{1}'
    ticklabels[-1] = fmt.format(ticklabels[-1], formatter.get_offset())
    axis.set_ticklabels(ticklabels)
    last_ticklabel = axis.get_ticklabels()[-1]
    last_ticklabel.set_horizontalalignment('left')

    # Draw edges in colorbar bands to correct thin white bands that
    # appear in buggy PDF viewers. See:
    # https://github.com/matplotlib/matplotlib/pull/1301
    cb.solids.set_edgecolor("face")

    # Done.
    return cb


def outline_text(ax):
    """Add a white outline to all text to make it stand out from the
    background."""
    effects = [patheffects.withStroke(linewidth=2, foreground='w')]
    for artist in ax.findobj(text.Text):
        artist.set_path_effects(effects)
