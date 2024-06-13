/*
 * Copyright (C) 2015-2024  Leo Singer
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


#include "cubic_interp.h"
#include "branch_prediction.h"
#include "vmath.h"
#include <math.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>

/* Allow contraction of a * b + c to a faster fused multiply-add operation.
 * This pragma is supposedly standard C, but only clang seems to support it.
 * On other compilers, floating point contraction is ON by default at -O3. */
#if defined(__clang__) || defined(__llvm__)
#pragma STDC FP_CONTRACT ON
#endif

#define VCLIP(x, a, b) VMIN(VMAX((x), (a)), (b))
#define VCUBIC(a, t) (t * (t * (t * a[0] + a[1]) + a[2]) + a[3])


struct cubic_interp {
    double f, t0, length;
    double a[][4];
};


struct bicubic_interp {
    v2df fx, x0, xlength;
    v4df a[][4];
};


/*
 * Calculate coefficients of the interpolating polynomial in the form
 *      a[0] * t^3 + a[1] * t^2 + a[2] * t + a[3]
 */
static void cubic_interp_init_coefficients(
    double *a, const double *z, const double *z1)
{
    if (UNLIKELY(!isfinite(z1[1] + z1[2])))
    {
        /* If either of the inner grid points are NaN or infinite,
         * then fall back to nearest-neighbor interpolation. */
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;
        a[3] = z[1];
    } else if (UNLIKELY(!isfinite(z1[0] + z1[3]))) {
        /* If either of the outer grid points are NaN or infinite,
         * then fall back to linear interpolation. */
        a[0] = 0;
        a[1] = 0;
        a[2] = z[2] - z[1];
        a[3] = z[1];
    } else {
        /* Otherwise, all of the grid points are finite.
         * Use cubic interpolation. */
        a[0] = 1.5 * (z[1] - z[2]) + 0.5 * (z[3] - z[0]);
        a[1] = z[0] - 2.5 * z[1] + 2 * z[2] - 0.5 * z[3];
        a[2] = 0.5 * (z[2] - z[0]);
        a[3] = z[1];
    }
}


cubic_interp *cubic_interp_init(
    const double *data, int n, double tmin, double dt)
{
    const int length = n + 6;
    cubic_interp *interp = malloc(sizeof(*interp) + length * sizeof(*interp->a));
    if (LIKELY(interp))
    {
        interp->f = 1 / dt;
        interp->t0 = 3 - interp->f * tmin;
        interp->length = length;
        for (int i = 0; i < length; i ++)
        {
            double z[4];
            for (int j = 0; j < 4; j ++)
            {
                z[j] = data[VCLIP(i + j - 4, 0, n - 1)];
            }
            cubic_interp_init_coefficients(interp->a[i], z, z);
        }
    }
    return interp;
}


void cubic_interp_free(cubic_interp *interp)
{
    free(interp);
}


double cubic_interp_eval(const cubic_interp *interp, double t)
{
    if (UNLIKELY(isnan(t)))
        return t;

    double x = t, xmin = 0.0, xmax = interp->length - 1.0;
    x *= interp->f;
    x += interp->t0;
    x = VCLIP(x, xmin, xmax);

    double ix = VFLOOR(x);
    x -= ix;

    const double *a = interp->a[(int) ix];
    return VCUBIC(a, x);
}


bicubic_interp *bicubic_interp_init(
    const double *data, int ns, int nt,
    double smin, double tmin, double ds, double dt)
{
    const int slength = ns + 6;
    const int tlength = nt + 6;
    bicubic_interp *interp = aligned_alloc(
        alignof(bicubic_interp),
        sizeof(*interp) + slength * tlength * sizeof(*interp->a));
    if (LIKELY(interp))
    {
        interp->fx[0] = 1 / ds;
        interp->fx[1] = 1 / dt;
        interp->x0[0] = 3 - interp->fx[0] * smin;
        interp->x0[1] = 3 - interp->fx[1] * tmin;
        interp->xlength[0] = slength;
        interp->xlength[1] = tlength;

        for (int is = 0; is < slength; is ++)
        {
            for (int it = 0; it < tlength; it ++)
            {
                double a[4][4], a1[4][4];
                for (int js = 0; js < 4; js ++)
                {
                    double z[4];
                    int ks = VCLIP(is + js - 4, 0, ns - 1);
                    for (int jt = 0; jt < 4; jt ++)
                    {
                        int kt = VCLIP(it + jt - 4, 0, nt - 1);
                        z[jt] = data[ks * ns + kt];
                    }
                    cubic_interp_init_coefficients(a[js], z, z);
                }
                for (int js = 0; js < 4; js ++)
                {
                    for (int jt = 0; jt < 4; jt ++)
                    {
                        a1[js][jt] = a[jt][js];
                    }
                }
                for (int js = 0; js < 4; js ++)
                {
                    cubic_interp_init_coefficients(a[js], a1[js], a1[3]);
                }
                memcpy(interp->a[is * slength + it], a, sizeof(a));
            }
        }
    }
    return interp;
}


void bicubic_interp_free(bicubic_interp *interp)
{
    free(interp);
}


double bicubic_interp_eval(const bicubic_interp *interp, double s, double t)
{
    if (UNLIKELY(isnan(s) || isnan(t)))
        return s + t;

    v2df x = {s, t}, xmin = {0.0, 0.0}, xmax = interp->xlength - 1.0;
    x *= interp->fx;
    x += interp->x0;
    x = VCLIP(x, xmin, xmax);

    v2df ix = VFLOOR(x);
    x -= ix;

    const v4df *a = interp->a[(int) (ix[0] * interp->xlength[0] + ix[1])];
    v4df b = VCUBIC(a, x[1]);
    return VCUBIC(b, x[0]);
}
