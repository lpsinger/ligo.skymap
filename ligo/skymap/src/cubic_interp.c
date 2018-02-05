/*
 * Copyright (C) 2015-2017  Leo Singer
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


#include "cubic_interp.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>


static int clip_int(int t, int min, int max)
{
    if (t < min)
        return min;
    else if (t > max)
        return max;
    else
        return t;
}


static double clip_double(double t, double min, double max)
{
    if (t < min)
        return min;
    else if (t > max)
        return max;
    else
        return t;
}


/*
 * Calculate coefficients of the interpolating polynomial in the form
 *      a[0] * t^3 + a[1] * t^2 + a[2] * t + a[3]
 */
static void cubic_interp_init_coefficients(
    double *a, const double *z, const double *z1)
{
    if (!isfinite(z1[1] + z1[2]))
    {
        /* If either of the inner grid points are NaN or infinite,
         * then fall back to nearest-neighbor interpolation. */
        a[0] = 0;
        a[1] = 0;
        a[2] = 0;
        a[3] = z[1];
    } else if (!isfinite(z1[0] + z1[3])) {
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


static double cubic_eval(const double *a, double t)
{
    return ((a[0] * t + a[1]) * t + a[2]) * t + a[3];
}


static void cubic_interp_index(
    double f, double t0, double length, double *t, double *i)
{
    *t = modf(clip_double(*t * f + t0, 0, length), i);
}


cubic_interp *cubic_interp_init(
    const double *data, int n, double tmin, double dt)
{
    cubic_interp *interp;
    const int length = n + 6;
    interp = malloc(sizeof(*interp) + length * sizeof(*interp->a));
    if (interp)
    {
        interp->f = 1 / dt;
        interp->t0 = 3 - interp->f * tmin;
        interp->length = length;
        for (int i = 0; i < length; i ++)
        {
            double z[4];
            for (int j = 0; j < 4; j ++)
            {
                z[j] = data[clip_int(i + j - 4, 0, n - 1)];
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
    double i;
    if (isnan(t))
        return t;
    cubic_interp_index(interp->f, interp->t0, interp->length, &t, &i);
    return cubic_eval(interp->a[(int) i], t);
}


bicubic_interp *bicubic_interp_init(
    const double *data, int ns, int nt,
    double smin, double tmin, double ds, double dt)
{
    bicubic_interp *interp;
    const int slength = ns + 6;
    const int tlength = nt + 6;
    interp = malloc(sizeof(*interp) + slength * tlength * sizeof(*interp->a));
    if (interp)
    {
        interp->fs = 1 / ds;
        interp->ft = 1 / dt;
        interp->s0 = 3 - interp->fs * smin;
        interp->t0 = 3 - interp->ft * tmin;
        interp->slength = slength;
        interp->tlength = tlength;
        for (int is = 0; is < slength; is ++)
        {
            for (int it = 0; it < tlength; it ++)
            {
                double a[4][4], a1[4][4];
                for (int js = 0; js < 4; js ++)
                {
                    double z[4];
                    int ks = clip_int(is + js - 4, 0, ns - 1);
                    for (int jt = 0; jt < 4; jt ++)
                    {
                        int kt = clip_int(it + jt - 4, 0, nt - 1);
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
    const double (*a)[4];
    double b[4];
    double is, it;

    if (isnan(s) || isnan(t))
        return s + t;
    cubic_interp_index(interp->fs, interp->s0, interp->slength, &s, &is);
    cubic_interp_index(interp->ft, interp->t0, interp->tlength, &t, &it);
    a = interp->a[(int) (is * interp->slength + it)];
    for (int i = 0; i < 4; i ++)
        b[i] = cubic_eval(a[i], s);
    return cubic_eval(b, t);
}
