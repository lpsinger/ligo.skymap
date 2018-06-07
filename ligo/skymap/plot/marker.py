#
# Copyright (C) 2016  Leo Singer
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
Specialized markers
"""

from matplotlib.path import Path
import numpy as np

__all__ = ('earth', 'reticle')


earth = Path.unit_circle()
verts = np.concatenate([earth.vertices, [[-1, 0], [1, 0], [0, -1], [0, 1]]])
codes = np.concatenate([earth.codes, [Path.MOVETO, Path.LINETO] * 2])
earth = Path(verts, codes)
del verts, codes
earth.__doc__ = """The Earth symbol (circle and cross)."""


def reticle(inner=0.5, outer=1.0, angle=0.0, which='lrtb'):
    """
    Create a reticle (crosshairs) marker.

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
    """
    angle = np.deg2rad(angle)
    x = np.cos(angle)
    y = np.sin(angle)
    R = [[x, y], [-y, x]]
    vertdict = {'l': [-1, 0], 'r': [1, 0], 'b': [0, -1], 't': [0, 1]}
    verts = [vertdict[direction] for direction in which]
    codes = [Path.MOVETO, Path.LINETO] * len(verts)
    verts = np.dot(verts, R)
    verts = np.swapaxes([inner * verts, outer * verts], 0, 1).reshape(-1, 2)
    return Path(verts, codes)
