# Copyright (C) 2016-2022  Leo Singer, John Veitch
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
"""Read HDF5 posterior sample chain HDF5 files."""

import numpy as np
import h5py
from astropy.table import Column, Table

# Constants from lalinference module
POSTERIOR_SAMPLES = 'posterior_samples'
LINEAR = 0
CIRCULAR = 1
FIXED = 2
OUTPUT = 3

__all__ = ('read_samples', 'write_samples')


def _identity(x):
    return x


_colname_map = (('rightascension', 'ra', _identity),
                ('right_ascension', 'ra', _identity),
                ('declination', 'dec', _identity),
                ('logdistance', 'dist', np.exp),
                ('distance', 'dist', _identity),
                ('luminosity_distance', 'dist', _identity),
                ('polarisation', 'psi', _identity),
                ('chirpmass', 'mc', _identity),
                ('chirp_mass', 'mc', _identity),
                ('a_spin1', 'a1', _identity),
                ('a_1', 'a1', _identity),
                ('a_spin2', 'a2', _identity),
                ('a_2', 'a2', _identity),
                ('tilt_spin1', 'tilt1', _identity),
                ('tilt_1', 'tilt1', _identity),
                ('tilt_spin2', 'tilt2', _identity),
                ('tilt_2', 'tilt2', _identity),
                ('geocent_time', 'time', _identity))


def _remap_colnames(table):
    for old_name, new_name, func in _colname_map:
        if old_name in table.colnames:
            table[new_name] = func(table.columns.pop(old_name))


def _find_table(group, tablename):
    """Recursively search an HDF5 group or file for a dataset by name.

    Parameters
    ----------
    group : `h5py.File` or `h5py.Group`
        The file or group to search
    tablename : str
        The name of the table to search for

    Returns
    -------
    dataset : `h5py.Dataset`
        The dataset whose name is `tablename`

    Raises
    ------
    KeyError
        If the table is not found or if multiple matching tables are found

    Examples
    --------
    Check that we can find a file by name:

    >>> import os.path
    >>> import tempfile
    >>> table = Table(np.eye(3), names=['a', 'b', 'c'])
    >>> with tempfile.TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     table.write(filename, path='foo/bar', append=True)
    ...     table.write(filename, path='foo/bat', append=True)
    ...     table.write(filename, path='foo/xyzzy/bat', append=True)
    ...     with h5py.File(filename, 'r') as f:
    ...         _find_table(f, 'bar')
    <HDF5 dataset "bar": shape (3,), type "|V24">

    Check that an exception is raised if the table is not found:

    >>> with tempfile.TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     table.write(filename, path='foo/bar', append=True)
    ...     table.write(filename, path='foo/bat', append=True)
    ...     table.write(filename, path='foo/xyzzy/bat', append=True)
    ...     with h5py.File(filename, 'r') as f:
    ...         _find_table(f, 'plugh')
    Traceback (most recent call last):
        ...
    KeyError: 'Table not found: plugh'

    Check that an exception is raised if multiple tables are found:

    >>> with tempfile.TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     table.write(filename, path='foo/bar', append=True)
    ...     table.write(filename, path='foo/bat', append=True)
    ...     table.write(filename, path='foo/xyzzy/bat', append=True)
    ...     with h5py.File(filename, 'r') as f:
    ...         _find_table(f, 'bat')
    Traceback (most recent call last):
        ...
    KeyError: 'Multiple tables called bat exist: foo/bat, foo/xyzzy/bat'

    """
    results = {}

    def visitor(key, value):
        _, _, name = key.rpartition('/')
        if name == tablename:
            results[key] = value

    group.visititems(visitor)

    if len(results) == 0:
        raise KeyError('Table not found: {0}'.format(tablename))

    if len(results) > 1:
        raise KeyError('Multiple tables called {0} exist: {1}'.format(
            tablename, ', '.join(sorted(results.keys()))))

    table, = results.values()
    return table


def read_samples(filename, path=None, tablename=POSTERIOR_SAMPLES):
    """Read an HDF5 sample chain file.

    Parameters
    ----------
    filename : str
        The path of the HDF5 file on the filesystem.
    path : str, optional
        The path of the dataset within the HDF5 file.
    tablename : str, optional
        The name of table to search for recursively within the HDF5 file.
        By default, search for 'posterior_samples'.

    Returns
    -------
    chain : `astropy.table.Table`
        The sample chain as an Astropy table.

    Examples
    --------
    Test reading a file written using the Python API:

    >>> import os.path
    >>> import tempfile
    >>> table = Table([
    ...     Column(np.ones(10), name='foo', meta={'vary': FIXED}),
    ...     Column(np.arange(10), name='bar', meta={'vary': LINEAR}),
    ...     Column(np.arange(10) * np.pi, name='bat', meta={'vary': CIRCULAR}),
    ...     Column(np.arange(10), name='baz', meta={'vary': OUTPUT})
    ... ])
    >>> with tempfile.TemporaryDirectory() as dir:
    ...     filename = os.path.join(dir, 'test.hdf5')
    ...     write_samples(table, filename, path='foo/bar/posterior_samples')
    ...     len(read_samples(filename))
    10

    Test reading a file that was written using the LAL HDF5 C API:

    >>> from importlib.resources import files
    >>> with files('ligo.skymap.io.tests.data').joinpath(
    ...         'test.hdf5').open('rb') as f:
    ...     table = read_samples(f)
    >>> table.colnames
    ['uvw', 'opq', 'lmn', 'ijk', 'def', 'abc', 'ghi', 'rst']

    """
    with h5py.File(filename, 'r') as f:
        if path is not None:  # Look for a given path
            table = f[path]
        else:  # Look for a given table name
            table = _find_table(f, tablename)
        table = Table.read(table)

    # Restore vary types.
    for i, column in enumerate(table.columns.values()):
        column.meta['vary'] = table.meta.get(
            'FIELD_{0}_VARY'.format(i), OUTPUT)

    # Restore fixed columns from table attributes.
    for key, value in table.meta.items():
        # Skip attributes from H5TB interface
        # (https://www.hdfgroup.org/HDF5/doc/HL/H5TB_Spec.html).
        if key == 'CLASS' or key == 'VERSION' or key == 'TITLE' or \
                key.startswith('FIELD_'):
            continue
        table.add_column(Column([value] * len(table), name=key,
                         meta={'vary': FIXED}))

    # Delete remaining table attributes.
    table.meta.clear()

    # Normalize column names.
    _remap_colnames(table)

    # Done!
    return table


def write_samples(table, filename, metadata=None, **kwargs):
    """Write an HDF5 sample chain file.

    Parameters
    ----------
    table : `astropy.table.Table`
        The sample chain as an Astropy table.
    filename : str
        The path of the HDF5 file on the filesystem.
    metadata: dict (optional)
        Dictionary of (path, value) pairs of metadata attributes
        to add to the output file
    kwargs: dict
        Any keyword arguments for `astropy.table.Table.write`.

    Examples
    --------
    Check that we catch columns that are supposed to be FIXED but are not:

    >>> table = Table([
    ...     Column(np.arange(10), name='foo', meta={'vary': FIXED})
    ... ])
    >>> write_samples(table, 'bar.hdf5', 'bat/baz')
    Traceback (most recent call last):
        ...
    AssertionError: 
    Arrays are not equal
    Column foo is a fixed column, but its values are not identical
    ...

    And now try writing an arbitrary example to a temporary file.

    >>> import os.path
    >>> import tempfile
    >>> table = Table([
    ...     Column(np.ones(10), name='foo', meta={'vary': FIXED}),
    ...     Column(np.arange(10), name='bar', meta={'vary': LINEAR}),
    ...     Column(np.arange(10) * np.pi, name='bat', meta={'vary': CIRCULAR}),
    ...     Column(np.arange(10), name='baz', meta={'vary': OUTPUT}),
    ...     Column(np.ones(10), name='plugh'),
    ...     Column(np.arange(10), name='xyzzy')
    ... ])
    >>> with tempfile.TemporaryDirectory() as dir:
    ...     write_samples(
    ...         table, os.path.join(dir, 'test.hdf5'), path='bat/baz',
    ...         metadata={'bat/baz': {'widget': 'shoephone'}})

    """  # noqa: W291
    # Copy the table so that we do not modify the original.
    table = table.copy()

    # Make sure that all tables have a 'vary' type.
    for column in table.columns.values():
        if 'vary' not in column.meta:
            if np.all(column[0] == column[1:]):
                column.meta['vary'] = FIXED
            else:
                column.meta['vary'] = OUTPUT
    # Reconstruct table attributes.
    for colname, column in tuple(table.columns.items()):
        if column.meta['vary'] == FIXED:
            np.testing.assert_array_equal(column[1:], column[0],
                                          'Column {0} is a fixed column, but '
                                          'its values are not identical'
                                          .format(column.name))
            table.meta[colname] = column[0]
            del table[colname]
    for i, column in enumerate(table.columns.values()):
        table.meta['FIELD_{0}_VARY'.format(i)] = column.meta.pop('vary')
    table.write(filename, format='hdf5', **kwargs)
    if metadata:
        with h5py.File(filename, 'r+') as hdf:
            for internal_path, attributes in metadata.items():
                for key, value in attributes.items():
                    try:
                        hdf[internal_path].attrs[key] = value
                    except KeyError:
                        raise KeyError(
                            'Unable to set metadata {0}[{1}] = {2}'.format(
                                internal_path, key, value))
