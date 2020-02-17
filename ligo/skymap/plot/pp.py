#
# Copyright (C) 2012-2020  Leo Singer
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
"""Axes subclass for making probability--probability (P--P) plots.

Example
-------
You can create new P--P plot axes by passing the keyword argument
``projection='pp_plot'`` when creating new Matplotlib axes.

.. plot::
   :context: reset
   :include-source:
   :align: center

    import ligo.skymap.plot
    from matplotlib import pyplot as plt
    import numpy as np

    n = 100
    p_values_1 = np.random.uniform(size=n) # One experiment
    p_values_2 = np.random.uniform(size=n) # Another experiment
    p_values_3 = np.random.uniform(size=n) # Yet another experiment

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111, projection='pp_plot')
    ax.add_confidence_band(n, alpha=0.95) # Add 95% confidence band
    ax.add_diagonal() # Add diagonal line
    ax.add_lightning(n, 20) # Add some random realizations of n samples
    ax.add_series(p_values_1, p_values_2, p_values_3) # Add our data

Or, you can call the constructor of `PPPlot` directly.

.. plot::
   :context: reset
   :include-source:
   :align: center

    from ligo.skymap.plot import PPPlot
    from matplotlib import pyplot as plt
    import numpy as np

    n = 100

    rect = [0.1, 0.1, 0.8, 0.8] # Where to place axes in figure
    fig = plt.figure(figsize=(5, 5))
    ax = PPPlot(fig, rect)
    fig.add_axes(ax)
    ax.add_confidence_band(n, alpha=0.95)
    ax.add_lightning(n, 20)
    ax.add_diagonal()

"""
import matplotlib
from matplotlib import axes
from matplotlib.projections import projection_registry
import scipy.stats
import numpy as np

__all__ = ('PPPlot',)


class PPPlot(axes.Axes):
    """Construct a probability--probability (P--P) plot."""

    name = 'pp_plot'

    def __init__(self, *args, **kwargs):
        # Call parent constructor
        super().__init__(*args, **kwargs)

        # Square axes, limits from 0 to 1
        self.set_aspect(1.0)
        self.set_xlim(0.0, 1.0)
        self.set_ylim(0.0, 1.0)

    @staticmethod
    def _make_series(p_values):
        for ps in p_values:
            if np.ndim(ps) == 1:
                ps = np.sort(np.atleast_1d(ps))
                n = len(ps)
                xs = np.concatenate(([0.], ps, [1.]))
                ys = np.concatenate(([0.], np.arange(1, n + 1) / n, [1.]))
            elif np.ndim(ps) == 2:
                xs = np.concatenate(([0.], ps[0], [1.]))
                ys = np.concatenate(([0.], ps[1], [1.]))
            else:
                raise ValueError('All series must be 1- or 2-dimensional')
            yield xs
            yield ys

    def add_series(self, *p_values, **kwargs):
        """Add a series of P-values to the plot.

        Parameters
        ----------
        p_values : `numpy.ndarray`
            One or more lists of P-values.

            If an entry in the list is one-dimensional, then it is interpreted
            as an unordered list of P-values. The ranked values will be plotted
            on the horizontal axis, and the cumulative fraction will be plotted
            on the vertical axis.

            If an entry in the list is two-dimensional, then the first subarray
            is plotted on the horizontal axis and the second subarray is
            plotted on the vertical axis.

        drawstyle : {'steps', 'lines', 'default'}
            Plotting style. If ``steps``, then plot steps to represent a
            piecewise constant function. If ``lines``, then connect points with
            straight lines. If ``default`` then use steps if there are more
            than 2 pixels per data point, or else lines.

        Other parameters
        ----------------
        kwargs :
            optional extra arguments to `matplotlib.axes.Axes.plot`

        """
        # Construct sequence of x, y pairs to pass to plot()
        args = list(self._make_series(p_values))
        min_n = min(len(ps) for ps in p_values)

        # Make copy of kwargs to pass to plot()
        kwargs = dict(kwargs)
        ds = kwargs.pop('drawstyle', 'default')
        if (ds == 'default' and 2 * min_n > self.bbox.width) or ds == 'lines':
            kwargs['drawstyle'] = 'default'
        else:
            kwargs['drawstyle'] = 'steps-post'

        return self.plot(*args, **kwargs)

    def add_worst(self, *p_values):
        """Mark the point at which the deviation is largest.

        Parameters
        ----------
        p_values : `numpy.ndarray`
            Same as in `add_series`.

        """
        series = list(self._make_series(p_values))
        for xs, ys in zip(series[0::2], series[1::2]):
            i = np.argmax(np.abs(ys - xs))
            x = xs[i]
            y = ys[i]
            if y == x:
                continue
            self.plot([x, x, 0], [0, y, y], '--', color='black', linewidth=0.5)
            if y < x:
                self.plot([x, y], [y, y], '-', color='black', linewidth=1)
                self.text(
                    x, y, ' {0:.02f} '.format(np.around(x - y, 2)),
                    ha='left', va='top')
            else:
                self.plot([x, x], [x, y], '-', color='black', linewidth=1)
                self.text(
                    x, y, ' {0:.02f} '.format(np.around(y - x, 2)),
                    ha='right', va='bottom')

    def add_diagonal(self, *args, **kwargs):
        """Add a diagonal line to the plot, running from (0, 0) to (1, 1).

        Other parameters
        ----------------
        kwargs :
            optional extra arguments to `matplotlib.axes.Axes.plot`

        """
        # Make copy of kwargs to pass to plot()
        kwargs = dict(kwargs)
        kwargs.setdefault('color', 'black')
        kwargs.setdefault('linestyle', 'dashed')
        kwargs.setdefault('linewidth', 0.5)

        # Plot diagonal line
        return self.plot([0, 1], [0, 1], *args, **kwargs)

    def add_lightning(self, nsamples, ntrials, **kwargs):
        """Add P-values drawn from a random uniform distribution, as a visual
        representation of the acceptable scatter about the diagonal.

        Parameters
        ----------
        nsamples : int
            Number of P-values in each trial
        ntrials : int
            Number of line series to draw.

        Other parameters
        ----------------
        kwargs :
            optional extra arguments to `matplotlib.axes.Axes.plot`

        """
        # Draw random samples
        args = np.random.uniform(size=(ntrials, nsamples))

        # Make copy of kwargs to pass to plot()
        kwargs = dict(kwargs)
        kwargs.setdefault('color', 'black')
        kwargs.setdefault('alpha', 0.5)
        kwargs.setdefault('linewidth', 0.25)

        # Plot series
        return self.add_series(*args, **kwargs)

    def add_confidence_band(
            self, nsamples, alpha=0.95, annotate=True, **kwargs):
        """Add a target confidence band.

        Parameters
        ----------
        nsamples : int
            Number of P-values
        alpha : float, default: 0.95
            Confidence level
        annotate : bool, optional, default: True
            If True, then label the confidence band.

        Other parameters
        ----------------
        **kwargs :
            optional extra arguments to `matplotlib.axes.Axes.fill_betweenx`

        """
        n = nsamples
        k = np.arange(0, n + 1)
        p = k / n
        ci_lo, ci_hi = scipy.stats.beta.interval(alpha, k + 1, n - k + 1)

        # Make copy of kwargs to pass to fill_betweenx()
        kwargs = dict(kwargs)
        kwargs.setdefault('color', 'lightgray')
        kwargs.setdefault('edgecolor', 'gray')
        kwargs.setdefault('linewidth', 0.5)
        fontsize = kwargs.pop('fontsize', 'x-small')

        if annotate:
            percent_sign = r'\%' if matplotlib.rcParams['text.usetex'] else '%'
            label = 'target {0:g}{1:s}\nconfidence band'.format(
                100 * alpha, percent_sign)

            self.annotate(
                label,
                xy=(1, 1),
                xytext=(0, 0),
                xycoords='axes fraction',
                textcoords='offset points',
                annotation_clip=False,
                horizontalalignment='right',
                verticalalignment='bottom',
                fontsize=fontsize,
                arrowprops=dict(
                    arrowstyle="->",
                    shrinkA=0, shrinkB=2, linewidth=0.5,
                    connectionstyle="angle,angleA=0,angleB=45,rad=0"))

        return self.fill_betweenx(p, ci_lo, ci_hi, **kwargs)

    @classmethod
    def _as_mpl_axes(cls):
        """Support placement in figure using the `projection` keyword argument.

        See http://matplotlib.org/devel/add_new_projection.html.
        """
        return cls, {}


projection_registry.register(PPPlot)
