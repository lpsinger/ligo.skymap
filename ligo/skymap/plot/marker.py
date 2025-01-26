#
# Copyright (C) 2016-2025  Leo Singer
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

import numpy as np
from matplotlib.markers import MarkerStyle
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

__all__ = ('earth', 'sun', 'moon', 'reticle')


earth = Path.unit_circle()
earth = MarkerStyle(
    Path(
        np.concatenate((earth.vertices, [[-1, 0], [1, 0], [0, -1], [0, 1]])),
        np.concatenate((earth.codes, [Path.MOVETO, Path.LINETO] * 2))
    ),
    fillstyle='none'
)
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

    plt.plot(0, 0, marker=earth, markersize=20, markeredgewidth=2)

"""


sun = Path.unit_circle()
sun = MarkerStyle(
    Path(
        np.concatenate((sun.vertices, [[0, 0], [1e-3, 0]])),
        np.concatenate((sun.codes, [Path.MOVETO, Path.LINETO]))
    ),
    capstyle='round',
    joinstyle='round',
    fillstyle='none'
)
sun.__doc__ = """
The Sun symbol (circle and dot).

Examples
--------
.. plot::
   :context: reset
   :include-source:
   :align: center

    from matplotlib import pyplot as plt
    from ligo.skymap.plot.marker import sun

    plt.plot(0, 0, marker=sun, markersize=20, markeredgewidth=2)

"""


def moon(phase, shadow=False):
    """Create a marker in the shape of the Moon.

    Parameters
    ----------
    phase : float
        Lunar phase in degrees between -180 and 180.
    shadow : bool
        If set, then the shadowed portion of the Moon is included in the
        marker, and its fill color can be set independently using the
        ``markerfacecoloralt`` keyword argument for
        :meth:`~matplotlib.axes.Axes.plot` (see
        :ref:`matplotlib:marker_fill_styles`).

    Returns
    -------
    markerstyle : matplotlib.markers.MarkerStyle

    Examples
    --------

    .. plot::
       :context: reset
       :include-source:
       :align: center

        from matplotlib import pyplot as plt
        from matplotlib.ticker import MultipleLocator
        import numpy as np
        from ligo.skymap.plot.marker import moon

        d_phase = 30
        phases = np.arange(-180, 180 + d_phase, d_phase)

        fig, ax = plt.subplots(figsize=(8, 3), tight_layout=True)
        ax.xaxis.set_major_locator(MultipleLocator(d_phase))
        for phase in phases:
            ax.plot(phase, 4, ms=20, marker=moon(phase, shadow=False), mfc="none", mec="black")
            ax.plot(phase, 3, ms=20, marker=moon(phase, shadow=False), mfc="goldenrod", mec="none")
            ax.plot(phase, 2, ms=20, marker=moon(phase, shadow=False), mfc="goldenrod", mec="k")
            ax.plot(phase, 1, ms=20, marker=moon(phase, shadow=True), mfc="goldenrod", mfcalt="gray", mec="none")
            ax.plot(phase, 0, ms=20, marker=moon(phase, shadow=True), mfc="goldenrod", mfcalt="gray", mec="black")
        ax.set_yticks(
            [0, 1, 2, 3, 4],
            ["shadow, fill, stroke", "shadow, fill", "fill, stroke", "fill", "stroke"],
        )
        ax.set_ylim(-0.5, 4.5)

    """  # noqa: E501
    angle = np.deg2rad(90 - phase)
    sign = np.sign(np.cos(angle))
    arc = Path.arc(90, 270, 9)

    path1 = arc.transformed(Affine2D().scale(sign * np.sin(angle), 1))
    path2 = arc.transformed(Affine2D().scale(-sign, 1))
    path3 = arc.transformed(Affine2D().scale(sign, 1))

    light_path = Path(
        np.concatenate((path1.vertices, path2.vertices[::-1])),
        np.concatenate((path1.codes, path2.codes[:0:-1], [Path.CLOSEPOLY])),
    )
    dark_path = Path(
        np.concatenate((path1.vertices, path3.vertices[::-1])),
        np.concatenate((path1.codes, path3.codes[:0:-1], [Path.CLOSEPOLY])),
    )

    markerstyle = MarkerStyle(light_path, joinstyle='miter')
    if shadow:
        markerstyle._alt_path = dark_path
        markerstyle._alt_transform = markerstyle._transform
    return markerstyle


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
