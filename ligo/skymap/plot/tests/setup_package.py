paths = ['baseline/*.png']

def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.plot.tests': paths}
