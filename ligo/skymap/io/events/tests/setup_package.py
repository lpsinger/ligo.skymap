import os

# If this package has tests data in the tests/data directory, add them to
# the paths here, see commented example
paths = ['coveragerc',
        os.path.join('data', '*.xml.gz'),
        os.path.join('data', 'gstlal_reference_psd', '*.xml.gz')
         ]

def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.io.events.tests': paths}
