paths = ['data/*.hdf5',
         'data/*.xml.gz',
         'data/gstlal_reference_psd/*.xml.gz']

def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.io.tests': paths}
