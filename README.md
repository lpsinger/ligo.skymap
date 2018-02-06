# ligo.skymap

Includes sky map I/O, plotting, and postprocessing codes; BAYESTAR; and KDE
posterior samples to FITS workflow (skyarea).

Structured as an [Astropy affiliated package](http://www.astropy.org/affiliated/) and based on the [Astropy package template](https://github.com/astropy/package-template).

## Python dependencies

* Python >= 3.5
* astropy
* healpy (note, transition to astropy_healpix, which is easier to install)
* h5py
* matplotlib
* pillow
* reproject
* scipy
* shapely

## LALSuite dependencies

These should become optional because they are only needed by BAYESTAR.

* glue
* Currently depends on SWIG bindings for lal, lalsimulation and lalinference.
* Can eliminate dependencies on lalinference.
* lal and lalsimulation should be available through `minlal` wrapper.

## Build dependencies

* chealpix
* gsl

These should both be bundled in the `cextern/` directory like astropy does with
`cfitsio`, `wcslib`, and `erfa` so that standalone builds without dependencies
are possible.
