/*
 * Copyright (C) 2017-2024  Leo Singer
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#include "bayestar_moc.h"
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <string.h>
#include <chealpix.h>

#include "branch_prediction.h"


int64_t nest2uniq64(uint8_t order, int64_t nest)
{
    if (nest < 0)
        return -1;
    else
        return nest + ((int64_t) 1 << 2 * (order + 1));
}


int8_t uniq2order64(int64_t uniq)
{
    if (uniq < 4)
        return -1;

    int8_t order;
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    int64_t o;
    asm("bsrq %1, %0\n\t"
        : "=r" (o)
        : "rm" (uniq));
    order = o;
#else
    order = 63 - __builtin_clzll(uniq);
#endif
    return (order >> 1) - 1;
}


double uniq2pixarea64(int64_t uniq)
{
    int8_t order = uniq2order64(uniq);
    if (order < 0)
        return GSL_NAN;
    else
        return ldexp(M_PI / 3, -2 * order);
}


int8_t uniq2nest64(int64_t uniq, int64_t *nest)
{
    int8_t order = uniq2order64(uniq);
    if (order < 0)
        *nest = -1;
    else
        *nest = uniq - ((int64_t) 1 << 2 * (order + 1));
    return order;
}


void uniq2ang64(int64_t uniq, double *theta, double *phi)
{
    int64_t nest;
    int8_t order = uniq2nest64(uniq, &nest);
    if (order < 0) {
        *theta = *phi = GSL_NAN;
    } else {
        int64_t nside = (int64_t) 1 << order;
        pix2ang_nest64(nside, nest, theta, phi);
    }
}


void *moc_rasterize64(
    const void *pixels, size_t offset, size_t in_stride, size_t out_stride,
    size_t len, size_t *npix, int8_t order)
{
    /* If the parameter order >= 0, then rasterize at that order.
     * Otherwise, find maximum order. Note: normally MOC datasets are stored in
     * order of ascending MOC index, so the last pixel should have the highest
     * order. However, our rasterization algorithm doesn't depend on this
     * sorting, so let's just do a linear search for the maximum order. */
    int8_t max_order;
    {
        int64_t max_uniq = 0;
        for (size_t i = 0; i < len; i ++)
        {
            const void *pixel = (const char *) pixels + i * in_stride;
            const int64_t uniq = *(const int64_t *) pixel;
            if (uniq > max_uniq)
                max_uniq = uniq;
        }
        max_order = uniq2order64(max_uniq);
    }
    if (UNLIKELY(max_order < 0)) {
        GSL_ERROR_NULL("invalid UNIQ value", GSL_EINVAL);
    }

    /* Don't handle downsampling here, because we don't know how to do
     * reduction across pixels without more knowledge of the pixel datatype and
     * contents. */
    if (order >= max_order)
        max_order = order;
    else if (order >= 0)
        GSL_ERROR_NULL("downsampling not implemented", GSL_EUNIMPL);

    /* Allocate output. */
    *npix = 12 * ((size_t) 1 << 2 * max_order);
    void *ret = calloc(*npix, out_stride);
    if (!ret)
        GSL_ERROR_NULL("not enough memory to allocate image", GSL_ENOMEM);

    /* Paint pixels into output. */
    for (size_t i = 0; i < len; i ++)
    {
        const void *pixel = (const char *) pixels + i * in_stride;
        int64_t nest;
        order = uniq2nest64(*(const int64_t *) pixel, &nest);
        if (UNLIKELY(order < 0)) {
            free(ret);
            GSL_ERROR_NULL("invalid UNIQ value", GSL_EINVAL);
        }
        const size_t reps = (size_t) 1 << 2 * (max_order - order);
        for (size_t j = 0; j < reps; j ++)
            memcpy((char *) ret + (nest * reps + j) * out_stride,
                (const char *) pixel + offset, out_stride);
    }

    return ret;
}
