import os

# If this package has tests data in the tests/data directory, add them to
# the paths here, see commented example
paths = [os.path.join('baseline', '*.png')]

def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.plot.tests': paths}
