#
# Copyright (C) 2020-2022  Leo Singer
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
"""Tools for adapting LIGO-LW row ID formats."""
import re

from ligo.lw.ligolw import Param
from ligo.lw.lsctables import TableByName
from ligo.lw.table import Column, TableStream
from ligo.lw.types import FormatFunc, FromPyType

__all__ = ('use_in',)

IDTypes = {'ilwd:char', 'ilwd:char_u'}

ROWID_PYTYPE = int
ROWID_TYPE = FromPyType[ROWID_PYTYPE]
ROWID_FORMATFUNC = FormatFunc[ROWID_TYPE]


_ilwd_regex = re.compile(r'\s*\w+:\w+:(\d+)\s*')


def ilwd_to_int(ilwd):
    match = _ilwd_regex.fullmatch(ilwd)
    if not match:
        raise ValueError(f'"{ilwd}" is not formatt like an ilwd')
    return int(match[1])


def use_in(ContentHandler):
    """Convert from old-style to new-style row IDs on the fly.

    This is loosely adapted from :func:`ligo.lw.utils.ilwd.strip_ilwdchar`.

    Notes
    -----
    When building a ContentHandler, this must be the _outermost_ decorator,
    outside of :func:`ligo.lw.lsctables.use_in`, :func:`ligo.lw.param.use_in`,
    or :func:`ligo.lw.table.use_in`.

    Examples
    --------
    >>> from importlib.resources import as_file, files
    >>> from ligo.lw import array, ligolw, lsctables, param, table, utils
    >>> from ligo.skymap.util import ilwd
    >>> @ilwd.use_in
    ... @lsctables.use_in
    ... @param.use_in
    ... @table.use_in
    ... class ContentHandler(ligolw.LIGOLWContentHandler):
    ...     pass
    >>> with as_file(files('ligo.skymap.io.tests.data').joinpath(
    ...         'G197392_coinc.xml.gz')) as f:
    ...     xmldoc = utils.load_filename(f, contenthandler=ContentHandler)
    >>> table = lsctables.SnglInspiralTable.get_table(xmldoc)
    >>> table[0].process_id
    0

    """

    def endElementNS(self, uri_localname, qname,
                     __orig_endElementNS=ContentHandler.endElementNS):
        """Convert values of <Param> elements from ilwdchar to int."""
        if isinstance(self.current, Param) and self.current.Type in IDTypes:
            new_value = ilwd_to_int(self.current.pcdata)
            self.current.Type = ROWID_TYPE
            self.current.pcdata = ROWID_FORMATFUNC(new_value)
        __orig_endElementNS(self, uri_localname, qname)

    remapped = {}

    def startColumn(self, parent, attrs,
                    __orig_startColumn=ContentHandler.startColumn):
        """Convert types in <Column> elements from ilwdchar to int.

        Notes
        -----
        This method is adapted from
        :func:`ligo.lw.utils.ilwd.strip_ilwdchar`.

        """
        result = __orig_startColumn(self, parent, attrs)

        # If this is an ilwdchar column, then create a function to convert its
        # rows' values for use in the startStream method below.
        if result.Type in IDTypes:
            remapped[(id(parent), result.Name)] = ilwd_to_int
            result.Type = ROWID_TYPE

        # If this is an ilwdchar column, then normalize the column name.
        if parent.Name in TableByName:
            validcolumns = TableByName[parent.Name].validcolumns
            if result.Name not in validcolumns:
                stripped_column_to_valid_column = {
                    Column.ColumnName(name): name for name in validcolumns}
                if result.Name in stripped_column_to_valid_column:
                    result.setAttribute(
                        'Name', stripped_column_to_valid_column[result.Name])

        return result

    def startStream(self, parent, attrs,
                    __orig_startStream=ContentHandler.startStream):
        """Convert values in table <Stream> elements from ilwdchar to int.

        Notes
        -----
        This method is adapted from
        :meth:`ligo.lw.table.TableStream.config`.

        """
        result = __orig_startStream(self, parent, attrs)
        if isinstance(result, TableStream):
            loadcolumns = set(parent.columnnames)
            if parent.loadcolumns is not None:
                # FIXME:  convert loadcolumns attributes to sets to
                # avoid the conversion.
                loadcolumns &= set(parent.loadcolumns)
            result._tokenizer.set_types([
                (remapped.pop((id(parent), colname), pytype)
                 if colname in loadcolumns else None)
                for pytype, colname
                in zip(parent.columnpytypes, parent.columnnames)])
        return result

    ContentHandler.endElementNS = endElementNS
    ContentHandler.startColumn = startColumn
    ContentHandler.startStream = startStream

    return ContentHandler
