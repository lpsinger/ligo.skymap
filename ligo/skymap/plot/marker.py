#
# Copyright (C) 2016-2019  Leo Singer
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
"""Specialized markers."""

from matplotlib.path import Path
import numpy as np

__all__ = ('earth', 'reticle')


earth = Path.unit_circle()
verts = np.concatenate([earth.vertices, [[-1, 0], [1, 0], [0, -1], [0, 1]]])
codes = np.concatenate([earth.codes, [Path.MOVETO, Path.LINETO] * 2])
earth = Path(verts, codes)
del verts, codes
earth.__doc__ = """
The Earth symbol (circle and cross).

Examples
--------
.. plot::
   :context: reset
   :include-source:
   :align: center

    from matplotlib import pyplot as plt
    from ligo.skymap.plot.marker import earth

    plt.plot(0, 0, markersize=20, markeredgewidth=2,
             markerfacecolor='none', marker=earth)

"""


def reticle(inner=0.5, outer=1.0, angle=0.0, which='lrtb'):
    """Create a reticle or crosshairs marker.

    Parameters
    ----------
    inner : float
        Distance from the origin to the inside of the crosshairs.
    outer : float
        Distance from the origin to the outside of the crosshairs.
    angle : float
        Rotation in degrees; 0 for a '+' orientation and 45 for 'x'.

    Returns
    -------
    path : `matplotlib.path.Path`
        The new marker path, suitable for passing to Matplotlib functions
        (e.g., `plt.plot(..., marker=reticle())`)

    Examples
    --------
    .. plot::
       :context: reset
       :include-source:
       :align: center

        from matplotlib import pyplot as plt
        from ligo.skymap.plot.marker import reticle

        markers = [reticle(inner=0),
                   reticle(which='lt'),
                   reticle(which='lt', angle=45)]

        fig, ax = plt.subplots(figsize=(6, 2))
        ax.set_xlim(-0.5, 2.5)
        ax.set_ylim(-0.5, 0.5)
        for x, marker in enumerate(markers):
            ax.plot(x, 0, markersize=20, markeredgewidth=2, marker=marker)

    """
    angle = np.deg2rad(angle)
    x = np.cos(angle)
    y = np.sin(angle)
    rotation = [[x, y], [-y, x]]
    vertdict = {'l': [-1, 0], 'r': [1, 0], 'b': [0, -1], 't': [0, 1]}
    verts = [vertdict[direction] for direction in which]
    codes = [Path.MOVETO, Path.LINETO] * len(verts)
    verts = np.dot(verts, rotation)
    verts = np.swapaxes([inner * verts, outer * verts], 0, 1).reshape(-1, 2)
    return Path(verts, codes)
