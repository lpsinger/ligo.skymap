/*
 * Copyright (C) 2020  Leo Singer
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

#ifndef VMATH_H
#define VMATH_H

#ifndef __cplusplus

#include <math.h>


/* Vector types (gcc/clang/icc vector extension to the C language) */

typedef double v2df __attribute__ ((vector_size (2 * sizeof(double))));
typedef double v4df __attribute__ ((vector_size (4 * sizeof(double))));


/* Vectorized math functions using x86-64 intrinsics if available */

#ifdef __x86_64__
#include <immintrin.h>
#endif

#define V2DF_BINARY_OP(func, scalarfunc) \
static v2df v2df_ ## func(v2df a, v2df b) \
{ \
    v2df result; \
    for (int i = 0; i < 2; i ++) \
        result[i] = scalarfunc(a[i], b[i]); \
    return result; \
}

#define V2DF_UNARY_OP(func, scalarfunc) \
static v2df v2df_ ## func(v2df a) \
{ \
    v2df result; \
    for (int i = 0; i < 2; i ++) \
        result[i] = scalarfunc(a[i]); \
    return result; \
}

#ifdef __SSE2__
static v2df v2df_min(v2df a, v2df b) { return _mm_min_pd(a, b); }
static v2df v2df_max(v2df a, v2df b) { return _mm_max_pd(a, b); }
#else
V2DF_BINARY_OP(min, fmin)
V2DF_BINARY_OP(max, fmax)
#endif

#ifdef __SSE4_1__
static v2df v2df_floor(v2df a) { return _mm_floor_pd(a); }
#else
V2DF_UNARY_OP(floor, floor)
#endif


/* C11 generics for selected math functions */

static int int_min(int a, int b)
{
    return a < b ? a : b;
}

static int int_max(int a, int b)
{
    return a > b ? a : b;
}

#define VMIN(a, b) _Generic((a), \
    v2df: v2df_min, \
    int: int_min, \
    double: fmin \
)((a), (b))

#define VMAX(a, b) _Generic((a), \
    v2df: v2df_max, \
    int: int_max, \
    double: fmax \
)((a), (b))

#define VFLOOR(a) _Generic((a), \
    v2df: v2df_floor, \
    double: floor \
)(a)


#endif /* __cplusplus */

#endif /* VMATH_H */
