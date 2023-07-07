#
# Copyright (C) 2013-2023  Leo Singer
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
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
"""Functions that support the command line interface for plotting tools."""

import argparse
import os
import sys

import matplotlib

from ..plot import cmap  # noqa
from . import FileType, HelpChoicesAction, type_with_sideeffect, version_string

# Set no-op Matplotlib backend to defer importing anything that requires a GUI
# until we have determined that it is necessary based on the command line
# arguments.
if 'matplotlib.pyplot' in sys.modules:
    from matplotlib import pyplot as plt
    plt.switch_backend('Template')
else:
    matplotlib.use('Template', warn=False, force=True)
    from matplotlib import pyplot as plt

__all__ = ('get_figure_parser',)


class MatplotlibFigureType(FileType):

    def __init__(self):
        super().__init__('wb')

    @staticmethod
    def __show():
        from matplotlib import pyplot as plt
        return plt.show()

    @staticmethod
    def get_savefig_metadata(format):
        program, _ = os.path.splitext(os.path.basename(sys.argv[0]))
        cmdline = ' '.join([program] + sys.argv[1:])
        metadata = {'Title': cmdline}
        if format == 'png':
            metadata['Software'] = version_string
        elif format in {'pdf', 'ps', 'eps'}:
            metadata['Creator'] = version_string
        return metadata

    def __save(self):
        from matplotlib import pyplot as plt
        _, ext = os.path.splitext(self.string)
        format = ext.lower().lstrip('.')
        metadata = self.get_savefig_metadata(format)
        return plt.savefig(self.string, metadata=metadata)

    def __call__(self, string):
        from matplotlib import pyplot as plt
        if string == '-':
            plt.switch_backend(matplotlib.rcParamsOrig['backend'])
            return self.__show
        else:
            with super().__call__(string):
                pass
            plt.switch_backend('agg')
            self.string = string
            return self.__save


@type_with_sideeffect(str)
def colormap(value):
    from matplotlib import rcParams
    rcParams['image.cmap'] = value


@type_with_sideeffect(float)
def figwidth(value):
    from matplotlib import rcParams
    rcParams['figure.figsize'][0] = float(value)


@type_with_sideeffect(float)
def figheight(value):
    from matplotlib import rcParams
    rcParams['figure.figsize'][1] = float(value)


@type_with_sideeffect(int)
def dpi(value):
    from matplotlib import rcParams
    rcParams['figure.dpi'] = rcParams['savefig.dpi'] = float(value)


@type_with_sideeffect(int)
def transparent(value):
    from matplotlib import rcParams
    rcParams['savefig.transparent'] = bool(value)


def get_figure_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(
        'figure options', 'Options that affect figure output format')
    group.add_argument(
        '-o', '--output', metavar='FILE.{pdf,png}',
        default='-', type=MatplotlibFigureType(),
        help='output file, or - to plot to screen')
    group.add_argument(
        '--colormap', default='cylon', choices=plt.colormaps(), type=colormap,
        metavar='CMAP', help='matplotlib colormap')
    group.add_argument(
        '--help-colormap', action=HelpChoicesAction, choices=plt.colormaps())
    group.add_argument(
        '--figure-width', metavar='INCHES', type=figwidth, default='8',
        help='width of figure in inches')
    group.add_argument(
        '--figure-height', metavar='INCHES', type=figheight, default='6',
        help='height of figure in inches')
    group.add_argument(
        '--dpi', metavar='PIXELS', type=dpi, default=300,
        help='resolution of figure in dots per inch')
    group.add_argument(
        '--transparent', const='1', default='0', nargs='?', type=transparent,
        help='Save image with transparent background')
    return parser
