#
# Copyright (C) 2013-2024  Leo Singer
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
"""Functions that support the command line interface."""

import argparse
import glob
import inspect
import itertools
import logging
import os
import sys

import numpy as np

from ..util import sqlite
from .. import version

version_string = version.__package__ + ' ' + version.version


def _try_close(obj):
    try:
        obj.close()
    except AttributeError:
        pass


class Namespace(argparse.Namespace):
    """A Namespace that can be used in a context to close all open files."""

    def __enter__(self):
        return self

    def __exit__(self, *args, **kwargs):
        for value in vars(self).values():
            if isinstance(value, list):
                for item in value:
                    _try_close(item)
            _try_close(value)


class FileType(argparse.FileType):
    """Inherit from :class:`argparse.FileType` to enable opening stdin or
    stdout in binary mode.

    This is a workaround for https://bugs.python.org/issue14156, which is fixed
    in Python 3.9.12, 3.10.3, and >=3.11.0.

    FIXME: Remove this when we drop support for Python 3.10.
    """

    def __call__(self, string):
        if string == '-' and 'b' in self._mode:
            if 'r' in self._mode:
                return sys.stdin.buffer
            elif 'w' in self._mode:
                return sys.stdout.buffer
        return super().__call__(string)


class EnableAction(argparse.Action):

    def __init__(self,
                 option_strings,
                 dest,
                 default=True,
                 required=False,
                 help=None):
        opt, = option_strings
        if not opt.startswith('--enable-'):
            raise ValueError('Option string must start with --enable-')
        option_strings = [opt, opt.replace('--enable-', '--disable-')]
        super().__init__(
            option_strings,
            dest=dest,
            nargs=0,
            default=default,
            required=required,
            help=help)

    def __call__(self, parser, namespace, values, option_string):
        if option_string.startswith('--enable-'):
            setattr(namespace, self.dest, True)
        elif option_string.startswith('--disable-'):
            setattr(namespace, self.dest, False)
        else:
            raise RuntimeError('This code cannot be reached')


class GlobAction(argparse._StoreAction):
    """Generate a list of filenames from a list of filenames and globs."""

    def __call__(self, parser, namespace, values, *args, **kwargs):
        values = list(
            itertools.chain.from_iterable(glob.iglob(s) for s in values))
        if values:
            super().__call__(parser, namespace, values, *args, **kwargs)
        nvalues = getattr(namespace, self.dest)
        nvalues = 0 if nvalues is None else len(nvalues)
        if self.nargs == argparse.OPTIONAL:
            if nvalues > 1:
                msg = 'expected at most one file'
            else:
                msg = None
        elif self.nargs == argparse.ONE_OR_MORE:
            if nvalues < 1:
                msg = 'expected at least one file'
            else:
                msg = None
        elif self.nargs == argparse.ZERO_OR_MORE:
            msg = None
        elif int(self.nargs) != nvalues:
            msg = 'expected exactly %s file' % self.nargs
            if self.nargs != 1:
                msg += 's'
        else:
            msg = None
        if msg is not None:
            msg += ', but found '
            msg += '{} file'.format(nvalues)
            if nvalues != 1:
                msg += 's'
            raise argparse.ArgumentError(self, msg)


def get_waveform_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(
        'waveform options', 'Options that affect template waveform generation')
    # FIXME: The O1 uberbank high-mass template, SEOBNRv2_ROM_DoubleSpin, does
    # not support frequencies less than 30 Hz.
    group.add_argument(
        '--f-low', type=float, metavar='Hz', default=30,
        help='Low frequency cutoff')
    group.add_argument(
        '--f-high-truncate', type=float, default=0.95,
        help='Truncate waveform at this fraction of the maximum frequency of '
        'the PSD')
    group.add_argument(
        '--waveform', default='o2-uberbank',
        help='Template waveform approximant: e.g., TaylorF2threePointFivePN')
    return parser


def get_posterior_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(
        'posterior options', 'Options that affect the BAYESTAR posterior')
    group.add_argument(
        '--min-distance', type=float, metavar='Mpc',
        help='Minimum distance of prior in megaparsecs')
    group.add_argument(
        '--max-distance', type=float, metavar='Mpc',
        help='Maximum distance of prior in megaparsecs')
    group.add_argument(
        '--prior-distance-power', type=int, metavar='-1|2', default=2,
        help='Distance prior: -1 for uniform in log, 2 for uniform in volume')
    group.add_argument(
        '--cosmology', action='store_true',
        help='Use cosmological comoving volume prior')
    group.add_argument(
        '--enable-snr-series', action=EnableAction,
        help='Enable input of SNR time series')
    group.add_argument(
        '--rescale-loglikelihood', type=float, default=0.83,
        help='Rescale log likelihood by the square of this factor to account '
        'for excess technical noise from search pipeline')
    return parser


def get_mcmc_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(
        'BAYESTAR MCMC options', 'BAYESTAR options for MCMC sampling')
    group.add_argument(
        '--mcmc', action='store_true',
        help='Use MCMC sampling instead of Gaussian quadrature')
    group.add_argument(
        '--chain-dump', action='store_true',
        help='For MCMC methods, dump the sample chain to disk')
    return parser


class HelpChoicesAction(argparse.Action):

    def __init__(self,
                 option_strings,
                 choices=(),
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS):
        name = option_strings[0].replace('--help-', '')
        super().__init__(
            option_strings=option_strings,
            dest=dest,
            default=default,
            nargs=0,
            help='show supported values for --' + name + ' and exit')
        self._name = name
        self._choices = choices

    def __call__(self, parser, namespace, values, option_string=None):
        print('Supported values for --' + self._name + ':')
        for choice in self._choices:
            print(choice)
        parser.exit()


def type_with_sideeffect(type):
    def decorator(sideeffect):
        def func(value):
            ret = type(value)
            sideeffect(ret)
            return ret
        return func
    return decorator


@type_with_sideeffect(str)
def loglevel_type(value):
    try:
        value = int(value)
    except ValueError:
        value = value.upper()
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        level=value)


class LogLevelAction(argparse._StoreAction):

    def __init__(
            self, option_strings, dest, nargs=None, const=None, default=None,
            type=None, choices=None, required=False, help=None, metavar=None):
        metavar = '|'.join(logging._levelToName.values())
        type = loglevel_type
        super().__init__(
            option_strings, dest, nargs=nargs, const=const, default=default,
            type=type, choices=choices, required=required, help=help,
            metavar=metavar)


@type_with_sideeffect(int)
def seed(value):
    np.random.seed(value)


def get_random_parser():
    parser = argparse.ArgumentParser(add_help=False)
    group = parser.add_argument_group(
        'random number generator options',
        'Options that affect the Numpy pseudo-random number genrator')
    group.add_argument(
        '--seed', type=seed, help='Pseudo-random number generator seed '
        '[default: initialized from /dev/urandom or clock]')
    return parser


class HelpFormatter(argparse.RawDescriptionHelpFormatter,
                    argparse.ArgumentDefaultsHelpFormatter):
    pass


class ArgumentParser(argparse.ArgumentParser):
    """An ArgumentParser subclass with some sensible defaults.

    - Any ``.py`` suffix is stripped from the program name, because the
      program is probably being invoked from the stub shell script.

    - The description is taken from the docstring of the file in which the
      ArgumentParser is created.

    - If the description is taken from the docstring, then whitespace in
      the description is preserved.

    - A ``--version`` option is added that prints the version of ligo.skymap.

    - The Namespace object returned by ``parser.parse_args()`` can be used as a
      context manager, in a ``with::`` statement, to close all opened files.
    """

    def __init__(self,
                 prog=None,
                 usage=None,
                 description=None,
                 epilog=None,
                 parents=[],
                 prefix_chars='-',
                 fromfile_prefix_chars=None,
                 argument_default=None,
                 conflict_handler='error',
                 add_help=True):
        parent_frame = inspect.currentframe().f_back
        if prog is None:
            prog = parent_frame.f_code.co_filename
            prog = os.path.basename(prog)
            prog = prog.replace('_', '-').replace('.py', '')
        if description is None:
            description = parent_frame.f_globals.get('__doc__', None)
        super().__init__(
            prog=prog,
            usage=usage,
            description=description,
            epilog=epilog,
            parents=parents,
            formatter_class=HelpFormatter,
            prefix_chars=prefix_chars,
            fromfile_prefix_chars=fromfile_prefix_chars,
            argument_default=argument_default,
            conflict_handler=conflict_handler,
            add_help=add_help)
        self.register('action', 'glob', GlobAction)
        self.register('action', 'loglevel', LogLevelAction)
        self.add_argument(
            '--version', action='version', version=version_string)
        self.add_argument(
            '-l', '--loglevel', action='loglevel', default='INFO')

    def parse_known_args(self, args=None, namespace=None):
        if namespace is None:
            namespace = Namespace()
        return super().parse_known_args(args, namespace)


class DirType:
    """Factory for directory arguments."""

    def __init__(self, create=False):
        self._create = create

    def __call__(self, string):
        if self._create:
            try:
                os.makedirs(string, exist_ok=True)
            except OSError as e:
                raise argparse.ArgumentTypeError(e.message)
        else:
            try:
                os.listdir(string)
            except OSError as e:
                raise argparse.ArgumentTypeError(e)
        return string


class SQLiteType(FileType):
    """Open an SQLite database, or fail if it does not exist.

    Here is an example of trying to open a file that does not exist for
    reading (mode='r'). It should raise an exception:

    >>> import tempfile
    >>> filetype = SQLiteType('r')
    >>> filename = tempfile.mktemp()
    >>> # Note, simply check or a FileNotFound error in Python 3.
    >>> filetype(filename)
    Traceback (most recent call last):
      ...
    argparse.ArgumentTypeError: ...

    If the file already exists, then it's fine:

    >>> import sqlite3
    >>> filetype = SQLiteType('r')
    >>> with tempfile.NamedTemporaryFile() as f:
    ...     with sqlite3.connect(f.name) as db:
    ...         _ = db.execute('create table foo (bar char)')
    ...     filetype(f.name)
    <sqlite3.Connection object at ...>

    Here is an example of opening a file for writing (mode='w'), which should
    overwrite the file if it exists. Even if the file was not an SQLite
    database beforehand, this should work:

    >>> filetype = SQLiteType('w')
    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     print('This is definitely not an SQLite file.', file=f)
    ...     f.flush()
    ...     with filetype(f.name) as db:
    ...         db.execute('create table foo (bar char)')
    <sqlite3.Cursor object at ...>

    Here is an example of opening a file for appending (mode='a'), which should
    NOT overwrite the file if it exists. If the file was not an SQLite database
    beforehand, this should raise an exception.

    >>> filetype = SQLiteType('a')
    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     print('This is definitely not an SQLite file.', file=f)
    ...     f.flush()
    ...     with filetype(f.name) as db:
    ...         db.execute('create table foo (bar char)')
    Traceback (most recent call last):
      ...
    sqlite3.DatabaseError: ...

    And if the database did exist beforehand, then opening for appending
    (mode='a') should not clobber existing tables.

    >>> filetype = SQLiteType('a')
    >>> with tempfile.NamedTemporaryFile() as f:
    ...     with sqlite3.connect(f.name) as db:
    ...         _ = db.execute('create table foo (bar char)')
    ...     with filetype(f.name) as db:
    ...         db.execute('select count(*) from foo').fetchone()
    (0,)
    """

    def __init__(self, mode):
        if mode not in 'arw':
            raise ValueError('Unknown file mode: {}'.format(mode))
        self.mode = mode

    def __call__(self, string):
        try:
            return sqlite.open(string, self.mode)
        except OSError as e:
            raise argparse.ArgumentTypeError(e)


def _sanitize_arg_value_for_xmldoc(value):
    if hasattr(value, 'read'):
        return value.name
    elif isinstance(value, tuple):
        return tuple(_sanitize_arg_value_for_xmldoc(v) for v in value)
    elif isinstance(value, list):
        return [_sanitize_arg_value_for_xmldoc(v) for v in value]
    else:
        return value


def register_to_xmldoc(xmldoc, parser, opts, **kwargs):
    from ligo.lw.utils import process
    params = {key: _sanitize_arg_value_for_xmldoc(value)
              for key, value in opts.__dict__.items()}
    return process.register_to_xmldoc(
        xmldoc, parser.prog, params, **kwargs, version=version_string)


start_msg = '\
Waiting for input on stdin. Type control-D followed by a newline to terminate.'
stop_msg = 'Reached end of file. Exiting.'


def iterlines(file, start_message=start_msg, stop_message=stop_msg):
    """Iterate over non-emtpy lines in a file."""
    is_tty = os.isatty(file.fileno())

    if is_tty:
        print(start_message, file=sys.stderr)

    while True:
        # Read a line.
        line = file.readline()

        if not line:
            # If we reached EOF, then exit.
            break

        # Strip off the trailing newline and any whitespace.
        line = line.strip()

        # Emit the line if it is not empty.
        if line:
            yield line

    if is_tty:
        print(stop_message, file=sys.stderr)


_compress_arg_map = {
    '.bz2': 'bz2',
    '.gz': 'gz',
    '.xz': 'xz',
    '.zst': 'zst'
}


def write_fileobj(xmldoc, f):
    import ligo.lw.utils

    _, ext = os.path.splitext(f.name)
    compress = _compress_arg_map.get(ext.lower())

    with ligo.lw.utils.SignalsTrap():
        ligo.lw.utils.write_fileobj(xmldoc, f, compress=compress)
