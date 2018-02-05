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
#include <chealpix.h>
#include <gsl/gsl_errno.h>
#include <Python.h>
/* Ignore warnings in Numpy API itself */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#include <numpy/arrayobject.h>
#include <numpy/ufuncobject.h>
#pragma GCC diagnostic pop
#include <lal/bayestar_distance.h>
#include "six.h"
#include "omp_interruptible.h"


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
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const long long nside = npix2nside64(dimensions[2]);

    /* Assert that rotation matrix is 3x3 */
    assert(dimensions[1] == 3);

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
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const npy_intp npix = dimensions[1];

    /* Assert that array arguments are stored contiguously */
    assert(steps[6] == sizeof(double));

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
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const npy_intp npix = dimensions[1];

    /* Assert that array arguments are stored contiguously */
    assert(steps[6] == sizeof(double));

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
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();
    const npy_intp n = dimensions[0];
    const npy_intp npix = dimensions[1];

    /* Assert that array arguments are stored contiguously */
    assert(steps[6] == sizeof(double));

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


static const PyUFuncGenericFunction
    conditional_pdf_loops[] = {conditional_pdf_loop},
    conditional_cdf_loops[] = {conditional_cdf_loop},
    conditional_ppf_loops[] = {conditional_ppf_loop},
    moments_to_parameters_loops[] = {moments_to_parameters_loop},
    parameters_to_moments_loops[] = {parameters_to_moments_loop},
    volume_render_loops[] = {volume_render_loop},
    marginal_pdf_loops[] = {marginal_pdf_loop},
    marginal_cdf_loops[] = {marginal_cdf_loop},
    marginal_ppf_loops[] = {marginal_ppf_loop};

static const char volume_render_ufunc_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_INTP, NPY_INT, NPY_DOUBLE, NPY_BOOL,
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};

static const char double_ufunc_types[] = {
    NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE, NPY_DOUBLE};

static void *const no_ufunc_data[] = {NULL};

static const char modulename[] = "_distance";


static PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    modulename, NULL, -1, NULL,
    NULL, NULL, NULL, NULL
};


PyMODINIT_FUNC PyInit__distance(void); /* Silence -Wmissing-prototypes */
PyMODINIT_FUNC PyInit__distance(void)
{
    PyObject *module = NULL;

    gsl_set_error_handler_off();
    import_array();
    import_umath();

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

    #pragma GCC diagnostic pop

done:
    return module;
}


SIX_COMPAT_MODULE(_distance)
