/*
 * Copyright (C) 2013-2018  Leo Singer
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

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#ifdef _OPENMP
#include <omp.h>
#endif
#include <limits.h>
#include <stddef.h>
#include <chealpix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_nan.h>
#include <Python.h>
/* Ignore warnings in Numpy API itself */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>
#pragma GCC diagnostic pop
#include "bayestar_distance.h"
#include "bayestar_moc.h"
#include "bayestar_sky_map.h"
#include "cubic_interp.h"
#include "omp_interruptible.h"


typedef struct {
    PyObject_HEAD
} Omp;


static PyObject *
Omp_get_num_threads(Omp *NPY_UNUSED(self), void *NPY_UNUSED(closure))
{
    int ret;
#ifdef _OPENMP
    /*
     * omp_get_num_threads() always returns 1 in single-threaded code, so
     * wrap it inside an itty bitty parallel section.
     * See https://stackoverflow.com/questions/11071116
     */
    #pragma omp parallel
    {
        #pragma omp single
        ret = omp_get_num_threads();
    }
#else
    ret = 1;
#endif
    return PyLong_FromLong(ret);
}


static int Omp_set_num_threads(Omp *NPY_UNUSED(self), PyObject *value,
                               void *NPY_UNUSED(closure))
{
    unsigned long value_ulong = PyLong_AsUnsignedLong(value);
    if (value_ulong > INT_MAX)
    {
        PyErr_SetString(
            PyExc_OverflowError,
            "omp.num_threads must be less than or equal to INT_MAX");
        return -1;
    }
#ifdef _OPENMP
    omp_set_num_threads((int) value_ulong);
#endif
    return 0;
}


static PyGetSetDef Omp_getsetdefs[] = {
    {"num_threads", (getter)Omp_get_num_threads, (setter)Omp_set_num_threads,
     "Number of OpenMP threads", NULL},
    {NULL}
};


static PyTypeObject OmpType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "ligo.skymap.omp",         /* tp_name */
    sizeof(Omp),               /* tp_basicsize */
    0,                         /* tp_itemsize */
    0,                         /* tp_dealloc */
    0,                         /* tp_print */
    0,                         /* tp_getattr */
    0,                         /* tp_setattr */
    0,                         /* tp_reserved */
    0,                         /* tp_repr */
    0,                         /* tp_as_number */
    0,                         /* tp_as_sequence */
    0,                         /* tp_as_mapping */
    0,                         /* tp_hash  */
    0,                         /* tp_call */
    0,                         /* tp_str */
    0,                         /* tp_getattro */
    0,                         /* tp_setattro */
    0,                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,        /* tp_flags */
    "Global OpenMP settings",  /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    0,                         /* tp_methods */
    0,                         /* tp_members */
    Omp_getsetdefs,            /* tp_getset */
};



/*****************************************************************************/


static void conditional_pdf_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[4][i * steps[4]] = bayestar_distance_conditional_pdf(
        *(double *) &args[0][i * steps[0]],
        *(double *) &args[1][i * steps[1]],
        *(double *) &args[2][i * steps[2]],
        *(double *) &args[3][i * steps[3]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void conditional_cdf_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[4][i * steps[4]] = bayestar_distance_conditional_cdf(
        *(double *) &args[0][i * steps[0]],
        *(double *) &args[1][i * steps[1]],
        *(double *) &args[2][i * steps[2]],
        *(double *) &args[3][i * steps[3]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void conditional_ppf_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[4][i * steps[4]] = bayestar_distance_conditional_ppf(
        *(double *) &args[0][i * steps[0]],
        *(double *) &args[1][i * steps[1]],
        *(double *) &args[2][i * steps[2]],
        *(double *) &args[3][i * steps[3]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void moments_to_parameters_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        bayestar_distance_moments_to_parameters(
            *(double *) &args[0][i * steps[0]],
            *(double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void parameters_to_moments_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        bayestar_distance_parameters_to_moments(
            *(double *) &args[0][i * steps[0]],
            *(double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void volume_render_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    if (dimensions[1] != 3)
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid dimension");
        return;
    }

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const long long nside = npix2nside64(dimensions[2]);

    /* FIXME: Check that array arguments are stored contiguously */

    OMP_BEGIN_INTERRUPTIBLE
    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        if (OMP_WAS_INTERRUPTED)
            OMP_EXIT_LOOP_EARLY
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[11][i * steps[11]] = bayestar_volume_render(
            *(double *)   &args[0][i * steps[0]],
            *(double *)   &args[1][i * steps[1]],
            *(double *)   &args[2][i * steps[2]],
            *(int *)      &args[3][i * steps[3]],
            *(int *)      &args[4][i * steps[4]],
             (double *)   &args[5][i * steps[5]],
            nside,
            *(npy_bool *) &args[6][i * steps[6]],
             (double *)   &args[7][i * steps[7]],
             (double *)   &args[8][i * steps[8]],
             (double *)   &args[9][i * steps[9]],
             (double *)   &args[10][i * steps[10]]);
        #pragma GCC diagnostic pop
    }
    OMP_END_INTERRUPTIBLE

    gsl_set_error_handler(old_handler);
}


static void marginal_pdf_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    /* Assert that array arguments are stored contiguously */
    if (steps[6] != sizeof(double))
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid dimension");
        return;
    }

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const npy_intp npix = dimensions[1];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[5][i * steps[5]] =
            bayestar_distance_marginal_pdf(
            *(double *) &args[0][i * steps[0]], npix,
             (double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void marginal_cdf_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    /* Assert that array arguments are stored contiguously */
    if (steps[6] != sizeof(double))
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid dimension");
        return;
    }

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const npy_intp npix = dimensions[1];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[5][i * steps[5]] =
            bayestar_distance_marginal_cdf(
            *(double *) &args[0][i * steps[0]], npix,
             (double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


static void marginal_ppf_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    /* Assert that array arguments are stored contiguously */
    if (steps[6] != sizeof(double))
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid dimension");
        return;
    }

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const npy_intp npix = dimensions[1];

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *) &args[5][i * steps[5]] =
            bayestar_distance_marginal_ppf(
            *(double *) &args[0][i * steps[0]], npix,
             (double *) &args[1][i * steps[1]],
             (double *) &args[2][i * steps[2]],
             (double *) &args[3][i * steps[3]],
             (double *) &args[4][i * steps[4]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


/*****************************************************************************/


static void capsule_free(PyObject *self)
{
    free(PyCapsule_GetPointer(self, NULL));
}


static PyObject *rasterize(PyObject *NPY_UNUSED(module), PyObject *arg)
{
    PyArrayObject *arr = (PyArrayObject *) PyArray_FromAny(
        arg, NULL, 1, 1, NPY_ARRAY_CARRAY_RO, NULL);
    PyObject *uniq_key = PyUnicode_FromString("UNIQ");
    PyObject *new_fields = PyDict_New();
    PyObject *capsule = NULL;
    PyArrayObject *ret = NULL;

    if (!arr || !uniq_key || !new_fields)
        goto done;

    PyObject *fields = PyArray_DESCR(arr)->fields;
    if (!fields)
    {
        PyErr_SetString(PyExc_ValueError, "expected record array");
        goto done;
    }

    PyObject *uniq_field = PyDict_GetItem(fields, uniq_key);
    if (!uniq_field)
    {
        PyErr_SetString(PyExc_ValueError, "array does not have UNIQ column");
        goto done;
    }

    PyObject *uniq_dtype;
    int uniq_offset;
    if (!PyArg_ParseTuple(uniq_field, "Oi", &uniq_dtype, &uniq_offset))
        return NULL;

    if (!PyArray_DescrCheck(uniq_dtype))
    {
        PyErr_SetString(PyExc_ValueError, "not a dtype");
        goto done;
    }

    int uniq_typenum = ((PyArray_Descr *) uniq_dtype)->type_num;
    if (!PyArray_EquivTypenums(uniq_typenum, NPY_UINT64))
    {
        PyErr_SetString(PyExc_ValueError, "'uniq' field must be uint64");
        goto done;
    }

    if (uniq_offset != 0)
    {
        PyErr_SetString(PyExc_ValueError, "'uniq' field must be at offset 0");
        goto done;
    }

    const void *pixels = PyArray_DATA(arr);
    const size_t offset = ((PyArray_Descr *) uniq_dtype)->elsize;
    const size_t itemsize = PyArray_ITEMSIZE(arr) - offset;
    const size_t len = PyArray_SIZE(arr);
    size_t npix;

    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(fields, &pos, &key, &value))
    {
        PyObject *new_descr;
        int new_offset;

        if (PyObject_RichCompareBool(uniq_key, key, Py_EQ))
            continue;

        if (!PyArg_ParseTuple(value, "Oi", &new_descr, &new_offset))
            goto done;

        new_offset -= offset;
        PyObject *new_value = Py_BuildValue("Oi", new_descr, new_offset);
        if (!new_value)
            goto done;
        if (PyDict_SetItem(new_fields, key, new_value))
            goto done;
    }

    void *out;
    Py_BEGIN_ALLOW_THREADS
    out = moc_rasterize64(pixels, offset, itemsize, len, &npix);
    Py_END_ALLOW_THREADS
    if (!out)
        goto done;

    /* Prepare output object */
    capsule = PyCapsule_New(out, NULL, capsule_free);
    if (!capsule)
        goto done;

    PyArray_Descr *descr;
    if (PyArray_DescrConverter(new_fields, &descr) != NPY_SUCCEED)
        goto done;

    npy_intp dims[] = {npix};
    ret = (PyArrayObject *) PyArray_NewFromDescr(&PyArray_Type, descr, 1, dims, NULL,
        out, NPY_ARRAY_DEFAULT, NULL);
    if (!ret)
        goto done;

    if (PyArray_SetBaseObject(ret, capsule))
    {
        Py_DECREF(ret);
        ret = NULL;
    } else {
        capsule = NULL;
    }

done:
    Py_XDECREF(arr);
    Py_XDECREF(uniq_key);
    Py_XDECREF(new_fields);
    Py_XDECREF(capsule);
    return (PyObject *) ret;
}


static void nest2uniq_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(uint64_t *) &args[2][i * steps[2]] = nest2uniq64(
        *(int8_t *)   &args[0][i * steps[0]],
        *(uint64_t *) &args[1][i * steps[1]]);
        #pragma GCC diagnostic pop
    }
}


static void uniq2nest_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(int8_t *)   &args[1][i * steps[1]] = uniq2nest64(
        *(uint64_t *) &args[0][i * steps[0]],
         (uint64_t *) &args[2][i * steps[2]]);
        #pragma GCC diagnostic pop
    }
}


static void uniq2order_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(int8_t *)   &args[1][i * steps[1]] = uniq2order64(
        *(uint64_t *) &args[0][i * steps[0]]);
        #pragma GCC diagnostic pop
    }
}


static void uniq2pixarea_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *)  &args[1][i * steps[1]] = uniq2pixarea64(
        *(int64_t *) &args[0][i * steps[0]]);
        #pragma GCC diagnostic pop
    }
}


static void uniq2ang_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0];

    for (npy_intp i = 0; i < n; i ++)
    {
        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        uniq2ang64(
            *(uint64_t *) &args[0][i * steps[0]],
             (double *)   &args[1][i * steps[1]],
             (double *)   &args[2][i * steps[2]]);
        #pragma GCC diagnostic pop
    }
}


/*****************************************************************************/


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
    double min_inclination;
    double max_inclination;
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
    static const char *keywords[] = {"min_inclination", "max_inclination",
        "min_distance", "max_distance", "prior_distance_power", "cosmology",
        "gmst", "sample_rate", "epochs", "snrs", "responses", "locations",
        "horizons", NULL};

    /* Parse arguments */
    /* FIXME: PyArg_ParseTupleAndKeywords should expect keywords to be const */
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wincompatible-pointer-types"
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ddddiiddOOOOO",
        keywords, &min_inclination, &max_inclination, &min_distance,
        &max_distance, &prior_distance_power, &cosmology, &gmst, &sample_rate,
        &epochs_obj, &snrs_obj, &responses_obj, &locations_obj, &horizons_obj))
        return NULL;
    #pragma GCC diagnostic pop

    if (cosmology && prior_distance_power != 2)
    {
        PyErr_SetString(PyExc_ValueError,
            "BAYESTAR supports cosmological priors only for for prior_distance_power=2");
        return NULL;
    }

    if (min_inclination < 0)
    {
        PyErr_SetString(PyExc_ValueError,
            "min_inclination must be >= 0");
        return NULL;
    }

    if (max_inclination > M_PI_2)
    {
        PyErr_SetString(PyExc_ValueError,
            "max_inclination must be <= 90 (deg)");
        return NULL;
    }

    if (min_inclination > max_inclination)
    {
        PyErr_SetString(PyExc_ValueError,
            "min_inclination must be <= max_inclination");
        return NULL;
    }

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
        min_inclination, max_inclination, min_distance, max_distance,
        prior_distance_power, cosmology, gmst, nifos, nsamples, sample_rate,
        epochs, snrs, responses, locations, horizons);
    Py_END_ALLOW_THREADS
    gsl_set_error_handler(old_handler);

    /* Give the Python interpreter a chance to check for pending signals
     * and raise a KeyboardInterrupt if necessary.
     *
     * Why do we need to call free(pixels) below? Here's why. Chances are that
     * any interrupt occurred while bayestar_sky_map_toa_phoa_snr was in one of
     * its OpenMP parallel sections, in which case the function intercepted
     * the signal, terminated the parallel section early, and returned NULL.
     * However, it is possible (though unlikely) that the interrupt could have
     * occurred during the brief serial sections and that the function could
     * have run successfully to completion. In that case, we still need to
     * return NULL to Python and report the exception! */
    if (PyErr_CheckSignals())
    {
        free(pixels);
        pixels = NULL;
        goto fail;
    }

    if (!pixels)
    {
        /* The only remaining way for bayestar_sky_map_toa_phoa_snr to
         * return NULL is if a call to malloc failed. */
        PyErr_SetString(PyExc_MemoryError, "Out of memory");
        goto fail;
    }

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


static void log_posterior_toa_phoa_snr_loop(
    char **args, npy_intp *dimensions, npy_intp *steps, void *NPY_UNUSED(data))
{
    const npy_intp n = dimensions[0],
               nifos = dimensions[1],
            nsamples = dimensions[2],
                ndim = dimensions[3];

    if (ndim != 3)
    {
        PyErr_SetString(PyExc_RuntimeError, "Invalid dimension");
        return;
    }

    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

    #pragma omp parallel for
    for (npy_intp i = 0; i < n; i ++)
    {
        const float complex *snrs[nifos];
        const float (*responses[nifos])[3];
        const double *locations[nifos];

        for (npy_intp j = 0; j < nifos; j ++)
        {
            snrs[j] = (const float complex *)
                &args[13][i * steps[13] + j * steps[19]];
            responses[j] = (const float (*)[3])
                &args[14][i * steps[14] + j * steps[21]];
            locations[j] = (const double *)
                &args[15][i * steps[15] + j * steps[24]];
        }

        /* FIXME: args must be void ** to avoid alignment warnings */
        #pragma GCC diagnostic push
        #pragma GCC diagnostic ignored "-Wcast-align"
        *(double *)   &args[17][i * steps[17]] = bayestar_log_posterior_toa_phoa_snr(
        *(double *)   &args[0][i * steps[0]],
        *(double *)   &args[1][i * steps[1]],
        *(double *)   &args[2][i * steps[2]],
        *(double *)   &args[3][i * steps[3]],
        *(double *)   &args[4][i * steps[4]],
        *(double *)   &args[5][i * steps[5]],
        *(double *)   &args[6][i * steps[6]],
        *(double *)   &args[7][i * steps[7]],
        *(int *)      &args[8][i * steps[8]],
        *(npy_bool *) &args[9][i * steps[9]],
        *(double *)   &args[10][i * steps[10]],
        nifos, nsamples,
        *(double *) &args[11][i * steps[11]],
         (const double *) &args[12][i * steps[12]],
         snrs, responses, locations,
         (const double *) &args[16][i * steps[16]]);
        #pragma GCC diagnostic pop
    }

    gsl_set_error_handler(old_handler);
}


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
    ret = bayestar_test() + cubic_interp_test();
    Py_END_ALLOW_THREADS
    gsl_set_error_handler(old_handler);
    return PyLong_FromLong(ret);
}


/*****************************************************************************/


static const PyUFuncGenericFunction
    conditional_pdf_loops[] = {conditional_pdf_loop},
    conditional_cdf_loops[] = {conditional_cdf_loop},
    conditional_ppf_loops[] = {conditional_ppf_loop},
    moments_to_parameters_loops[] = {moments_to_parameters_loop},
    parameters_to_moments_loops[] = {parameters_to_moments_loop},
    volume_render_loops[] = {volume_render_loop},
    marginal_pdf_loops[] = {marginal_pdf_loop},
    marginal_cdf_loops[] = {marginal_cdf_loop},
    marginal_ppf_loops[] = {marginal_ppf_loop},
    nest2uniq_loops[] = {nest2uniq_loop},
    uniq2nest_loops[] = {uniq2nest_loop},
    uniq2order_loops[] = {uniq2order_loop},
    uniq2pixarea_loops[] = {uniq2pixarea_loop},
    uniq2ang_loops[] = {uniq2ang_loop},
    log_posterior_toa_phoa_snr_loops[] = {log_posterior_toa_phoa_snr_loop},
    signal_amplitude_model_loops[] = {signal_amplitude_model_loop};

static const char log_posterior_toa_phoa_snr_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
    NPY_DOUBLE, NPY_DOUBLE, NPY_INT, NPY_BOOL, NPY_DOUBLE, NPY_DOUBLE,
    NPY_DOUBLE, NPY_CFLOAT, NPY_FLOAT, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};

static const char volume_render_ufunc_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INTP, NPY_INT, NPY_DOUBLE, NPY_BOOL,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};

static const char double_ufunc_types[] = {
                      NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE,
                      NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE},
                  nest2uniq_types[] = {NPY_INT8, NPY_UINT64, NPY_UINT64},
                  uniq2nest_types[] = {NPY_UINT64, NPY_INT8, NPY_UINT64},
                  uniq2order_types[] = {NPY_UINT64, NPY_INT8},
                  uniq2pixarea_types[] = {NPY_UINT64, NPY_DOUBLE},
                  uniq2ang_types[] = {NPY_UINT64, NPY_DOUBLE, NPY_DOUBLE},
                  signal_amplitude_model_ufunc_types[] = {
                      NPY_CDOUBLE, NPY_CDOUBLE, NPY_DOUBLE,
                      NPY_DOUBLE, NPY_CDOUBLE};

static void *const no_ufunc_data[] = {NULL};

static const char modulename[] = "core";

static PyMethodDef methods[] = {
    {"rasterize", rasterize, METH_O, "fill me in"},
    {"toa_phoa_snr", (PyCFunction)sky_map_toa_phoa_snr,
        METH_VARARGS | METH_KEYWORDS, "fill me in"},
    {"test", (PyCFunction)test,
        METH_NOARGS, "fill me in"},
    {NULL, NULL, 0, NULL}
};

static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    modulename, NULL, -1, methods,
    NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC PyInit_core(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC PyInit_core(void)
{
    PyObject *module = NULL;

    gsl_set_error_handler_off();
    import_array();
    import_umath();

    OmpType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&OmpType) < 0)
        return NULL;

    sky_map_descr = sky_map_create_descr();
    if (!sky_map_descr)
        goto done;

    module = PyModule_Create(&moduledef);
    if (!module)
        goto done;

    PyModule_AddObject(module, "omp",
        PyObject_CallFunctionObjArgs((PyObject *) &OmpType, NULL));

    /* Ignore warnings in Numpy API */
    #pragma GCC diagnostic push
    #ifdef __clang__
    #pragma GCC diagnostic ignored "-Wincompatible-pointer-types-discards-qualifiers"
    #else
    #pragma GCC diagnostic ignored "-Wdiscarded-qualifiers"
    #endif

    PyModule_AddObject(
        module, "log_posterior_toa_phoa_snr", PyUFunc_FromFuncAndDataAndSignature(
            log_posterior_toa_phoa_snr_loops, no_ufunc_data,
            log_posterior_toa_phoa_snr_types, 1, 17, 1, PyUFunc_None,
            "log_posterior_toa_phoa_snr", NULL, 0,
            "(),(),(),(),(),(),(),(),(),(),(),(),(nifos),(nifos,nsamples),(nifos,ndim,ndim),(nifos,ndim),(nifos)->()"));

    PyModule_AddObject(
        module, "conditional_pdf", PyUFunc_FromFuncAndData(
            conditional_pdf_loops, no_ufunc_data,
            double_ufunc_types, 1, 4, 1, PyUFunc_None,
            "conditional_pdf", NULL, 0));

    PyModule_AddObject(
        module, "conditional_cdf", PyUFunc_FromFuncAndData(
            conditional_cdf_loops, no_ufunc_data,
            double_ufunc_types, 1, 4, 1, PyUFunc_None,
            "conditional_cdf", NULL, 0));

    PyModule_AddObject(
        module, "conditional_ppf", PyUFunc_FromFuncAndData(
            conditional_ppf_loops, no_ufunc_data,
            double_ufunc_types, 1, 4, 1, PyUFunc_None,
            "conditional_ppf", NULL, 0));

    PyModule_AddObject(
        module, "moments_to_parameters", PyUFunc_FromFuncAndData(
            moments_to_parameters_loops, no_ufunc_data,
            double_ufunc_types, 1, 2, 3, PyUFunc_None,
            "moments_to_parameters", NULL, 0));

    PyModule_AddObject(
        module, "parameters_to_moments", PyUFunc_FromFuncAndData(
            parameters_to_moments_loops, no_ufunc_data,
            double_ufunc_types, 1, 2, 3, PyUFunc_None,
            "parameters_to_moments", NULL, 0));

    PyModule_AddObject(
        module, "volume_render", PyUFunc_FromFuncAndDataAndSignature(
            volume_render_loops, no_ufunc_data,
            volume_render_ufunc_types, 1, 11, 1, PyUFunc_None,
            "volume_render", NULL, 0,
            "(),(),(),(),(),(i,i),(),(n),(n),(n),(n)->()"));

    PyModule_AddObject(
        module, "marginal_pdf", PyUFunc_FromFuncAndDataAndSignature(
            marginal_pdf_loops, no_ufunc_data,
            double_ufunc_types, 1, 5, 1, PyUFunc_None,
            "marginal_pdf", NULL, 0,
            "(),(n),(n),(n),(n)->()"));

    PyModule_AddObject(
        module, "marginal_cdf", PyUFunc_FromFuncAndDataAndSignature(
            marginal_cdf_loops, no_ufunc_data,
            double_ufunc_types, 1, 5, 1, PyUFunc_None,
            "marginal_cdf", NULL, 0,
            "(),(n),(n),(n),(n)->()"));

    PyModule_AddObject(
        module, "marginal_ppf", PyUFunc_FromFuncAndDataAndSignature(
            marginal_ppf_loops, no_ufunc_data,
            double_ufunc_types, 1, 5, 1, PyUFunc_None,
            "marginal_ppf", NULL, 0,
            "(),(n),(n),(n),(n)->()"));

    PyModule_AddObject(
        module, "nest2uniq", PyUFunc_FromFuncAndData(
            nest2uniq_loops, no_ufunc_data,
            nest2uniq_types, 1, 2, 1, PyUFunc_None,
            "nest2uniq", NULL, 0));

    PyModule_AddObject(
        module, "uniq2nest", PyUFunc_FromFuncAndData(
            uniq2nest_loops, no_ufunc_data,
            uniq2nest_types, 1, 1, 2, PyUFunc_None,
            "uniq2nest", NULL, 0));

    PyModule_AddObject(
        module, "uniq2order", PyUFunc_FromFuncAndData(
            uniq2order_loops, no_ufunc_data,
            uniq2order_types, 1, 1, 1, PyUFunc_None,
            "uniq2order", NULL, 0));

    PyModule_AddObject(
        module, "uniq2pixarea", PyUFunc_FromFuncAndData(
            uniq2pixarea_loops, no_ufunc_data,
            uniq2pixarea_types, 1, 1, 1, PyUFunc_None,
            "uniq2pixarea", NULL, 0));

    PyModule_AddObject(
        module, "uniq2ang", PyUFunc_FromFuncAndData(
            uniq2ang_loops, no_ufunc_data,
            uniq2ang_types, 1, 1, 2, PyUFunc_None,
            "uniq2ang", NULL, 0));

    PyModule_AddObject(
        module, "signal_amplitude_model", PyUFunc_FromFuncAndData(
            signal_amplitude_model_loops, no_ufunc_data,
            signal_amplitude_model_ufunc_types, 1, 4, 1, PyUFunc_None,
            "signal_amplitude_model", NULL, 0));
    #pragma GCC diagnostic pop

done:
    return module;
}
