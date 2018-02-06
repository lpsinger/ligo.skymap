paths = ['cylon.csv', 'ne_simplified_coastline.json']

def get_package_data():
    return {
        _ASTROPY_PACKAGE_NAME_ + '.plot': paths}
