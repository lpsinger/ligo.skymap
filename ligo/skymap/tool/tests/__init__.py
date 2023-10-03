from importlib import metadata
from importlib.resources import as_file, files
import multiprocessing
import os
import shutil
import stat
import subprocess
import sys
import tempfile

dist = 'ligo.skymap'
group = 'console_scripts'

__all__ = ('entry_points', 'run_entry_point', 'run_ligolw', 'run_lalsuite')

entry_points = {entry_point.name: entry_point
                for entry_point in metadata.distribution(dist).entry_points
                if entry_point.group == group}


def run_entry_point(name, *args):
    main = entry_points[name].load()
    try:
        main(args)
    except SystemExit as e:
        if e.code != 0:
            raise subprocess.CalledProcessError(e.code, [name, *args])


def exec_ligolw(name, *args):
    main, = (entry_point.load() for entry_point
             in metadata.distribution('python-ligo-lw').entry_points
             if entry_point.name == name)
    sys.argv = [name, *args]
    main()


def run_ligolw(name, *args):
    """Run an external tool that is provided by python-ligo-lw.

    This is trivial if python-ligo-lw has actually been installed and LALSuite
    is in the PATH. If python-ligo-lw is not installed and is only present as
    an egg, then things get more complicated.
    """
    path = shutil.which(name)
    if path:
        # The tool has been installed, so we can just call it.
        subprocess.check_call([path, *args])
    else:
        # The tool has not been installed, so we have to try to run it using
        # importlib.metadata.
        process = multiprocessing.Process(
            target=exec_ligolw, args=[name, *args])
        process.start()
        process.join()
        if process.exitcode != 0:
            raise subprocess.CalledProcessError(
                process.exitcode, [name, *args])


def run_lalsuite(name, *args):
    """Run an external tool that is provided by LALSuite.

    This is trivial if lalsuite has actually been installed it is in the PATH.
    If LALSuite is not installed and is only present as an egg, then things get
    more complicated because of how the LALSuite wheel packs binaries as
    package data.
    """
    path = shutil.which(name)
    if path:
        # The tool has been installed, so we can just call it.
        subprocess.check_call([path, *args])
    else:
        # The tool has not been installed, so we have to find the underlying
        # binary inside the lalapps module.
        with as_file(files('lalapps.bin').joinpath(name)) as path:
            # Copy to a temporary file so that we can make it executable.
            # For some reason, when eggs are extracted, permissions are not
            # preserved.
            with tempfile.NamedTemporaryFile(dir=os.path.dirname(path)) as tmp:
                with open(path, 'rb') as orig:
                    shutil.copyfileobj(orig, tmp)
                tmp.flush()
                fd = tmp.fileno()
                stat_result = os.fstat(fd)
                mode = stat_result.st_mode | stat.S_IXUSR
                os.fchmod(fd, mode)
                subprocess.check_call([tmp.name, *args])
