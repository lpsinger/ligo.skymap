/*
 * Copyright (C) 2017  Leo Singer
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


#include "bayestar_moc.h"
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <stdlib.h>
#include <string.h>
#include <chealpix.h>


uint64_t nest2uniq64(uint8_t order, uint64_t nest)
{
    return nest + ((int64_t) 1 << 2 * (order + 1));
}


int8_t uniq2order64(uint64_t uniq)
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
    for (order = -1; uniq; uniq >>= 1, order ++)
        /* noop */;
#endif
    return (order >> 1) - 1;
}


double uniq2pixarea64(uint64_t uniq)
{
    return ldexp(M_PI / 3, -2 * uniq2order64(uniq));
}


int8_t uniq2nest64(uint64_t uniq, uint64_t *nest)
{
    int8_t order = uniq2order64(uniq);
    *nest = uniq - ((uint64_t) 1 << 2 * (order + 1));
    return order;
}


void uniq2ang64(uint64_t uniq, double *theta, double *phi)
{
    uint64_t nest;
    int8_t order = uniq2nest64(uniq, &nest);
    int64_t nside = (int64_t) 1 << order;
    pix2ang_nest64(nside, nest, theta, phi);
}


void *moc_rasterize64(const void *pixels, size_t offset, size_t itemsize, size_t len, size_t *npix)
{
    /* Calculate pixel size. */
    const size_t pixelsize = offset + itemsize;

    /* Find maximum order. Note: normally MOC datasets are stored in order
     * of ascending MOC index, so the last pixel should have the highest order.
     * However, our rasterization algorithm doesn't depend on this sorting,
     * so let's just do a linear search for the maximum order. */
    int8_t max_order;
    {
        uint64_t max_uniq = 0;
        for (size_t i = 0; i < len; i ++)
        {
            const void *pixel = (const char *) pixels + i * pixelsize;
            const uint64_t uniq = *(const uint64_t *) pixel;
            if (uniq > max_uniq)
                max_uniq = uniq;
        }
        max_order = uniq2order64(max_uniq);
    }

    /* Allocate output. */
    *npix = 12 * ((size_t) 1 << 2 * max_order);
    void *ret = calloc(*npix, itemsize);
    if (!ret)
        GSL_ERROR_NULL("not enough memory to allocate image", GSL_ENOMEM);

    /* Paint pixels into output. */
    for (size_t i = 0; i < len; i ++)
    {
        const void *pixel = (const char *) pixels + i * pixelsize;
        uint64_t nest;
        const int8_t order = uniq2nest64(*(const uint64_t *) pixel, &nest);
        const size_t reps = (size_t) 1 << 2 * (max_order - order);
        for (size_t j = 0; j < reps; j ++)
            memcpy((char *) ret + (nest * reps + j) * itemsize, (const char *) pixel + offset, itemsize);
    }

    return ret;
}
