/*
 * Copyright (C) 2013-2017  Leo Singer
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 */

#include "config.h"
#include <Python.h>
/* Ignore warnings in Numpy API itself */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>
#pragma GCC diagnostic pop
#include <chealpix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_nan.h>
#include <lal/bayestar_sky_map.h>
#include <stddef.h>
#include "six.h"


static void capsule_free(PyObject *self)
{
    free(PyCapsule_GetPointer(self, NULL));
}


#define INPUT_LIST_OF_ARRAYS(NAME, NPYTYPE, DEPTH, CHECK) \
{ \
    const Py_ssize_t n = PySequence_Length(NAME##_obj); \
    if (n < 0) \
        goto fail; \
    else if (n != nifos) { \
        PyErr_SetString(PyExc_ValueError, #NAME \
        " appears to be the wrong length for the number of detectors"); \
        goto fail; \
    } \
    for (unsigned int iifo = 0; iifo < nifos; iifo ++) \
    { \
        PyObject *obj = PySequence_ITEM(NAME##_obj, iifo); \
        if (!obj) goto fail; \
        PyArrayObject *npy = (PyArrayObject *) \
            PyArray_ContiguousFromAny(obj, NPYTYPE, DEPTH, DEPTH); \
        Py_XDECREF(obj); \
        if (!npy) goto fail; \
        NAME##_npy[iifo] = npy; \
        { CHECK } \
        NAME[iifo] = PyArray_DATA(npy); \
    } \
}

#define FREE_INPUT_LIST_OF_ARRAYS(NAME) \
{ \
    for (unsigned int iifo = 0; iifo < nifos; iifo ++) \
        Py_XDECREF(NAME##_npy[iifo]); \
}

#define INPUT_VECTOR_NIFOS(CTYPE, NAME, NPYTYPE) \
    NAME##_npy = (PyArrayObject *) \
        PyArray_ContiguousFromAny(NAME##_obj, NPYTYPE, 1, 1); \
    if (!NAME##_npy) goto fail; \
    if (PyArray_DIM(NAME##_npy, 0) != nifos) \
    { \
        PyErr_SetString(PyExc_ValueError, #NAME \
            " appears to be the wrong length for the number of detectors"); \
        goto fail; \
    } \
    const CTYPE *NAME = PyArray_DATA(NAME##_npy);

#define INPUT_VECTOR_DOUBLE_NIFOS(NAME) \
    INPUT_VECTOR_NIFOS(double, NAME, NPY_DOUBLE)


static PyArray_Descr *sky_map_descr;


static PyArray_Descr *sky_map_create_descr(void)
{
    PyArray_Descr *dtype = NULL;
    PyObject *dtype_dict = Py_BuildValue("{s(ssss)s(cccc)s(IIII)}",
        "names", "UNIQ", "PROBDENSITY", "DISTMEAN", "DISTSTD",
        "formats", NPY_ULONGLONGLTR, NPY_DOUBLELTR, NPY_DOUBLELTR, NPY_DOUBLELTR,
        "offsets",
        (unsigned int) offsetof(bayestar_pixel, uniq),
        (unsigned int) offsetof(bayestar_pixel, value[0]),
        (unsigned int) offsetof(bayestar_pixel, value[1]),
        (unsigned int) offsetof(bayestar_pixel, value[2]));

    if (dtype_dict)
    {
        PyArray_DescrConverter(dtype_dict, &dtype);
        Py_DECREF(dtype_dict);
    }

    return dtype;
}


static PyObject *sky_map_toa_phoa_snr(
    PyObject *NPY_UNUSED(module), PyObject *args, PyObject *kwargs)
{
    /* Input arguments */
    double min_distance;
    double max_distance;
    int prior_distance_power;
    int cosmology;
    double gmst;
    unsigned int nifos;
    unsigned long nsamples = 0;
    double sample_rate;
    PyObject *epochs_obj;
    PyObject *snrs_obj;
    PyObject *responses_obj;
    PyObject *locations_obj;
    PyObject *horizons_obj;

    /* Names of arguments */
    static const char *keywords[] = {"min_distance", "max_distance",
        "prior_distance_power", "cosmology", "gmst", "sample_rate", "epochs",
        "snrs", "responses", "locations", "horizons", NULL};

    /* Parse arguments */
    /* FIXME: PyArg_ParseTupleAndKeywords should expect keywords to be const */
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddiiddOOOOO",
        keywords, &min_distance, &max_distance, &prior_distance_power,
        &cosmology, &gmst, &sample_rate, &epochs_obj, &snrs_obj,
        &responses_obj, &locations_obj, &horizons_obj)) return NULL;
    #pragma GCC diagnostic pop

    /* Determine number of detectors */
    {
        Py_ssize_t n = PySequence_Length(epochs_obj);
        if (n < 0) return NULL;
        nifos = n;
    }

    /* Return value */
    PyObject *out = NULL;
    double log_bci, log_bsn;

    /* Numpy array objects */
    PyArrayObject *epochs_npy = NULL, *snrs_npy[nifos], *responses_npy[nifos],
        *locations_npy[nifos], *horizons_npy = NULL;
    memset(snrs_npy, 0, sizeof(snrs_npy));
    memset(responses_npy, 0, sizeof(responses_npy));
    memset(locations_npy, 0, sizeof(locations_npy));

    /* Arrays of pointers for inputs with multiple dimensions */
    const float complex *snrs[nifos];
    const float (*responses[nifos])[3];
    const double *locations[nifos];

    /* Gather C-aligned arrays from Numpy types */
    INPUT_VECTOR_DOUBLE_NIFOS(epochs)
    INPUT_LIST_OF_ARRAYS(snrs, NPY_CFLOAT, 1,
        npy_intp dim = PyArray_DIM(npy, 0);
        if (iifo == 0)
            nsamples = dim;
        else if ((unsigned long)dim != nsamples)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of snrs to be vectors of the same length");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(responses, NPY_FLOAT, 2,
        if (PyArray_DIM(npy, 0) != 3 || PyArray_DIM(npy, 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of responses to be 3x3 arrays");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(locations, NPY_DOUBLE, 1,
        if (PyArray_DIM(npy, 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of locations to be vectors of length 3");
            goto fail;
        }
    )
    INPUT_VECTOR_DOUBLE_NIFOS(horizons)

    /* Call function */
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    size_t len;
    bayestar_pixel *pixels;
    Py_BEGIN_ALLOW_THREADS
    pixels = bayestar_sky_map_toa_phoa_snr(&len, &log_bci, &log_bsn,
        min_distance, max_distance, prior_distance_power, cosmology, gmst,
        nifos, nsamples, sample_rate, epochs, snrs, responses, locations,
        horizons);
    Py_END_ALLOW_THREADS
    gsl_set_error_handler(old_handler);

    PyErr_CheckSignals();

    if (!pixels)
        goto fail;

    /* Prepare output object */
    PyObject *capsule = PyCapsule_New(pixels, NULL, capsule_free);
    if (!capsule)
        goto fail;

    npy_intp dims[] = {len};
    Py_INCREF(sky_map_descr);
    out = PyArray_NewFromDescr(&PyArray_Type,
        sky_map_descr, 1, dims, NULL, pixels, NPY_ARRAY_DEFAULT, NULL);
    if (!out)
    {
        Py_DECREF(capsule);
        goto fail;
    }

    if (PyArray_SetBaseObject((PyArrayObject *) out, capsule))
    {
        Py_DECREF(out);
        out = NULL;
        goto fail;
    }

fail: /* Cleanup */
    Py_XDECREF(epochs_npy);
    FREE_INPUT_LIST_OF_ARRAYS(snrs)
    FREE_INPUT_LIST_OF_ARRAYS(responses)
    FREE_INPUT_LIST_OF_ARRAYS(locations)
    Py_XDECREF(horizons_npy);
    if (out) {
        out = Py_BuildValue("Ndd", out, log_bci, log_bsn);
    }
    return out;
};


static PyObject *log_likelihood_toa_phoa_snr(
    PyObject *NPY_UNUSED(module), PyObject *args, PyObject *kwargs)
{
    /* Input arguments */
    double ra;
    double sin_dec;
    double distance;
    double u;
    double twopsi;
    double t;
    double gmst;
    unsigned int nifos;
    unsigned long nsamples = 0;
    double sample_rate;
    PyObject *epochs_obj;
    PyObject *snrs_obj;
    PyObject *responses_obj;
    PyObject *locations_obj;
    PyObject *horizons_obj;

    /* Names of arguments */
    static const char *keywords[] = {"params", "gmst", "sample_rate", "epochs",
        "snrs", "responses", "locations", "horizons", NULL};

    /* Parse arguments */
    /* FIXME: PyArg_ParseTupleAndKeywords should expect keywords to be const */
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(dddddd)ddOOOOO",
        keywords, &ra, &sin_dec, &distance, &u, &twopsi, &t, &gmst,
        &sample_rate, &epochs_obj, &snrs_obj, &responses_obj, &locations_obj,
        &horizons_obj)) return NULL;
    #pragma GCC diagnostic pop

    /* Determine number of detectors */
    {
        Py_ssize_t n = PySequence_Length(epochs_obj);
        if (n < 0) return NULL;
        nifos = n;
    }

    /* Return value */
    PyObject *out = NULL;

    /* Numpy array objects */
    PyArrayObject *epochs_npy = NULL, *snrs_npy[nifos], *responses_npy[nifos],
        *locations_npy[nifos], *horizons_npy = NULL;
    memset(snrs_npy, 0, sizeof(snrs_npy));
    memset(responses_npy, 0, sizeof(responses_npy));
    memset(locations_npy, 0, sizeof(locations_npy));

    /* Arrays of pointers for inputs with multiple dimensions */
    const float complex *snrs[nifos];
    const float (*responses[nifos])[3];
    const double *locations[nifos];

    /* Gather C-aligned arrays from Numpy types */
    INPUT_VECTOR_DOUBLE_NIFOS(epochs)
    INPUT_LIST_OF_ARRAYS(snrs, NPY_CFLOAT, 1,
        npy_intp dim = PyArray_DIM(npy, 0);
        if (iifo == 0)
            nsamples = dim;
        else if ((unsigned long)dim != nsamples)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of snrs to be vectors of the same length");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(responses, NPY_FLOAT, 2,
        if (PyArray_DIM(npy, 0) != 3 || PyArray_DIM(npy, 1) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of responses to be 3x3 arrays");
            goto fail;
        }
    )
    INPUT_LIST_OF_ARRAYS(locations, NPY_DOUBLE, 1,
        if (PyArray_DIM(npy, 0) != 3)
        {
            PyErr_SetString(PyExc_ValueError,
                "expected elements of locations to be vectors of length 3");
            goto fail;
        }
    )
    INPUT_VECTOR_DOUBLE_NIFOS(horizons)

    /* Call function */
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const double ret = bayestar_log_likelihood_toa_phoa_snr(ra, sin_dec,
        distance, u, twopsi, t, gmst, nifos, nsamples, sample_rate, epochs,
        snrs, responses, locations, horizons);
    gsl_set_error_handler(old_handler);

    /* Prepare output object */
    out = PyFloat_FromDouble(ret);

fail: /* Cleanup */
    Py_XDECREF(epochs_npy);
    FREE_INPUT_LIST_OF_ARRAYS(snrs)
    FREE_INPUT_LIST_OF_ARRAYS(responses)
    FREE_INPUT_LIST_OF_ARRAYS(locations)
    Py_XDECREF(horizons_npy);
    return out;
};


static void signal_amplitude_model_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double complex *) &args[4][i * steps[4]] =
            bayestar_signal_amplitude_model(
            *(double complex *) &args[0][i * steps[0]],
            *(double complex *) &args[1][i * steps[1]],
            *(double *)         &args[2][i * steps[2]],
            *(double *)         &args[3][i * steps[3]]);
        #pragma GCC diagnostic pop
    }
}


static PyObject *test(
    PyObject *NPY_UNUSED(module), PyObject *NPY_UNUSED(arg))
{
    int ret;
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    Py_BEGIN_ALLOW_THREADS
    ret = bayestar_test();
    Py_END_ALLOW_THREADS
    gsl_set_error_handler(old_handler);
    return PyLong_FromLong(ret);
}


static PyMethodDef methods[] = {
    {"toa_phoa_snr", (PyCFunction)sky_map_toa_phoa_snr,
        METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"log_likelihood_toa_phoa_snr", (PyCFunction)log_likelihood_toa_phoa_snr,
        METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"test", (PyCFunction)test,
        METH_NOARGS, "fill me in"},
    {NULL, NULL, 0, NULL}
};


static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_sky_map", NULL, -1, methods,
    NULL, NULL, NULL, NULL
};


static const PyUFuncGenericFunction
    signal_amplitude_model_loops[] = {signal_amplitude_model_loop};

static const char signal_amplitude_model_ufunc_types[] = {
    NPY_CDOUBLE, NPY_CDOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_CDOUBLE};

static void *const no_ufunc_data[] = {NULL};


PyMODINIT_FUNC PyInit__sky_map(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC PyInit__sky_map(void)
{
    PyObject *module = NULL;

    gsl_set_error_handler_off();
    import_array();
    import_umath();

    sky_map_descr = sky_map_create_descr();
    if (!sky_map_descr)
        goto done;

    module = PyModule_Create(&moduledef);
    if (!module)
        goto done;

    /* Ignore warnings in Numpy API */
    #pragma GCC diagnostic push
    #ifdef __clang__
    #pragma GCC diagnostic ignored "-Wincompatible-pointer-types-discards-qualifiers"
    #else
    #pragma GCC diagnostic ignored "-Wdiscarded-qualifiers"
    #endif

    PyModule_AddObject(
        module, "signal_amplitude_model", PyUFunc_FromFuncAndData(
            signal_amplitude_model_loops, no_ufunc_data,
            signal_amplitude_model_ufunc_types, 1, 4, 1, PyUFunc_None,
            "signal_amplitude_model", NULL, 0));

    #pragma GCC diagnostic pop

done:
    return module;
}


SIX_COMPAT_MODULE(_sky_map)
