#!/usr/bin/env python
# This file is adapted from the Astropy affiliated package template, which is
# licensed under a 3-clause BSD style license - see astropy_helpers/LICENSE.rst

import builtins

# Ensure that astropy-helpers is available
import ah_bootstrap  # noqa

from setuptools import setup
from setuptools.config import read_configuration

from astropy_helpers.setup_helpers import register_commands, get_package_info
from astropy_helpers.version_helpers import generate_version_py

import pkg_resources

# Monkey patch find_packages to only locate our namespace package
from setuptools import PEP420PackageFinder
from astropy_helpers import setup_helpers

def find_packages(where='.', exclude=(), include=('*',)):
    return PEP420PackageFinder.find(
        where=where, exclude=exclude, include=('ligo.*',))

setup_helpers._find_packages = find_packages

# Store the package name in a built-in variable so it's easy
# to get from other parts of the setup infrastructure
builtins._ASTROPY_PACKAGE_NAME_ = read_configuration('setup.cfg')['metadata']['name']

# Create a dictionary with setup command overrides. Note that this gets
# information about the package (name and version) from the setup.cfg file.
cmdclass = register_commands()

# Freeze build information in version.py. Note that this gets information
# about the package (name and version) from the setup.cfg file.
version = generate_version_py()

# Get configuration information from all of the various subpackages.
# See the docstring for setup_helpers.update_package_files for more
# details.
package_info = get_package_info()

# Note that requires and provides should not be included in the call to
# ``setup``, since these are now deprecated. See this link for more details:
# https://groups.google.com/forum/#!topic/astropy-dev/urYO8ckB2uM

with open('requirements.txt', 'r') as f:
    install_requires = [str(r) for r in pkg_resources.parse_requirements(f)]

setup(version=version, cmdclass=cmdclass, install_requires=install_requires,
      **package_info)
