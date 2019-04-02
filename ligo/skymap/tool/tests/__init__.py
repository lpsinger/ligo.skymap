import distutils.spawn
import multiprocessing
import os
import shutil
import stat
import subprocess
import sys
import tempfile

import pkg_resources

dist = 'ligo.skymap'
group = 'console_scripts'

__all__ = ('entry_points', 'run_entry_point', 'run_glue', 'run_lalsuite')

entry_points = sorted(pkg_resources.get_entry_map(dist, group).keys())


def run_entry_point(name, *args):
    main = pkg_resources.load_entry_point(dist, group, name)
    try:
        main(args)
    except SystemExit as e:
        if e.code != 0:
            raise subprocess.CalledProcessError(e.code, [name, *args])


def exec_glue(name, *args):
    provider = pkg_resources.get_provider('glue')
    sys.argv = [name, *args]
    provider.run_script(name, {'__name__': '__main__'})


def run_glue(name, *args):
    """Run an external tool that is provided by Glue.

    This is trivial if glue has actually been installed and LALSuite is
    in the PATH. If glue is not installed and is only present as an egg,
    then things get more complicated.
    """
    path = distutils.spawn.find_executable(name)
    if path:
        # The tool has been installed, so we can just call it.
        subprocess.check_call([path, *args])
    else:
        # The tool has not been installed, so we have to try to run it using
        # pkg_resources.
        process = multiprocessing.Process(target=exec_glue, args=[name, *args])
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
    path = distutils.spawn.find_executable(name)
    if path:
        # The tool has been installed, so we can just call it.
        subprocess.check_call([path, *args])
    else:
        # The tool has not been installed, so we have to find the underlying
        # binary inside the lalapps module.
        path = pkg_resources.resource_filename(
            'lalapps', os.path.join('bin', name))

        # Copy to a temporary file so that we can make it executable.
        # For some reason, when eggs are extracted, permissions are not
        # preserved.
        with tempfile.NamedTemporaryFile(dir=os.path.dirname(path)) as temp:
            with open(path, 'rb') as orig:
                shutil.copyfileobj(orig, temp)
            temp.flush()
            fd = temp.fileno()
            stat_result = os.fstat(fd)
            mode = stat_result.st_mode | stat.S_IXUSR
            os.fchmod(fd, mode)
            subprocess.check_call([temp.name, *args])
