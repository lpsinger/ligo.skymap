#
# Copyright (C) 2016-2018  Leo Singer
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
"""Version control friendly Matplotlib PDF backend"""

import inspect

from matplotlib.backend_bases import register_backend
from matplotlib.backends import backend_pdf
from matplotlib.backends.backend_pdf import (
    FigureCanvasBase, MixedModeRenderer, RendererPdf)

__all__ = ()


class PdfFile(backend_pdf.PdfFile):

    def __init__(self, filename):
        super(PdfFile, self).__init__(filename)
        del self.infoDict['CreationDate']


class PdfPages(backend_pdf.PdfPages):

    # Copied verbatim from
    # `matplotlib.backends.backend_pdf.PdfPages.__init__`
    # in order to pick up local definition of `PdfFile`.
    def __init__(self, filename, keep_empty=True):
        self._file = PdfFile(filename)
        self.keep_empty = keep_empty
    __init__.__doc__ = backend_pdf.PdfPages.__init__.__doc__


# Execute verbatim definition from
# `matplotlib.backends.backend_pdf.FigureCanvasPdf`
# in order to pick up local definition of `PdfFile`.
exec(inspect.getsource(backend_pdf.FigureCanvasPdf))


class FigureManagerPdf(backend_pdf.FigureManagerPdf):
    pass


FigureCanvas = FigureCanvasPdf
FigureManager = FigureManagerPdf


def register():
    """Enable saving version-control friendly PDF files with Matplotlib.

    After calling, all subsequently saved PDF files will have the CreationDate
    metadata field stripped out in order to generate reproducible and
    deterministic output."""
    description = 'Portable Document Format (version control friendly)'
    register_backend('pdf', __name__, description=description)
