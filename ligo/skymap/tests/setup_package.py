import os

# If this package has tests data in the tests/data directory, add them to
# the paths here, see commented example
paths = ['coveragerc',
        os.path.join('data', '*.hdf5')
         ]

def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.tests': paths}
