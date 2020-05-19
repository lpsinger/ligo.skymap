/*
 * Copyright (C) 2015-2017  Leo Singer
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <assert.h>


int cubic_interp_test(void)
{
    {
        static const double data[] = {0, 0, 0, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = -10; t <= 10; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = 0;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for zero input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {1, 1, 1, 1};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = -10; t <= 10; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = 1;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for unit input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {1, 0, 1, 4};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = gsl_pow_2(t);
            gsl_test_abs(result, expected, 10 * GSL_DBL_EPSILON,
                "testing cubic interpolant for quadratic input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            GSL_POSINF, GSL_POSINF, GSL_POSINF, GSL_POSINF};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_POSINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for +inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, GSL_POSINF, GSL_POSINF, GSL_POSINF};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_POSINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for +inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            GSL_POSINF, GSL_POSINF, GSL_POSINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_POSINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for +inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, GSL_POSINF, GSL_POSINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_POSINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for +inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, 0, GSL_POSINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0.01; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = 0;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for +inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, GSL_NEGINF, GSL_POSINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        const double result = cubic_interp_eval(interp, 1);
        cubic_interp_free(interp);
        const double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing cubic interpolant for +inf input");
    }

    {
        static const double data[] = {
            0, GSL_POSINF, GSL_NEGINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        const double result = cubic_interp_eval(interp, 0);
        cubic_interp_free(interp);
        const double expected = GSL_POSINF;
        gsl_test_abs(result, expected, 0,
            "testing cubic interpolant for +inf input");
    }

    {
        static const double data[] = {
            0, GSL_NEGINF, GSL_NEGINF, GSL_NEGINF};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_NEGINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for -inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            GSL_NEGINF, GSL_NEGINF, GSL_NEGINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_NEGINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for -inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, GSL_NEGINF, GSL_NEGINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_NEGINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for -inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, 0, GSL_NEGINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0.01; t <= 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = 0;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for -inf input");
        }
        cubic_interp_free(interp);
    }

    {
        static const double data[] = {
            0, GSL_NEGINF, GSL_POSINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        const double result = cubic_interp_eval(interp, 0);
        cubic_interp_free(interp);
        const double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing cubic interpolant for -inf input");
    }

    {
        static const double data[] = {
            0, GSL_POSINF, GSL_NEGINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        const double result = cubic_interp_eval(interp, 1);
        cubic_interp_free(interp);
        const double expected = GSL_NEGINF;
        gsl_test_abs(result, expected, 0,
            "testing cubic interpolant for -inf input");
    }

    {
        static const double data[] = {
            0, GSL_NEGINF, GSL_POSINF, 0};
        cubic_interp *interp = cubic_interp_init(data, 4, -1, 1);
        assert(interp);
        for (double t = 0.01; t < 1; t += 0.01)
        {
            const double result = cubic_interp_eval(interp, t);
            const double expected = GSL_NEGINF;
            gsl_test_abs(result, expected, 0,
                "testing cubic interpolant for indeterminate input");
        }
        cubic_interp_free(interp);
    }


    {
        static const double constants[] = {
            0, 1, GSL_POSINF, GSL_NEGINF, GSL_NAN};
        for (unsigned k = 0; k < sizeof(constants) / sizeof(*constants); k ++)
        {
            double data[4][4];
            for (int i = 0; i < 4; i ++)
                for (int j = 0; j < 4; j ++)
                    data[i][j] = constants[k];
            bicubic_interp *interp = bicubic_interp_init(
                *data, 4, 4, -1, -1, 1, 1);
            for (double s = -5; s <= 2; s += 0.1)
            {
                for (double t = -5; t <= 1; t += 0.1)
                {
                    const double result = bicubic_interp_eval(interp, s, t);
                    const double expected = constants[k];
                    gsl_test_abs(result, expected, 0,
                        "testing bicubic interpolant for constant %g input",
                        constants[k]);
                }
            }
            assert(interp);
            bicubic_interp_free(interp);
        }
    }

    for (int k = 1; k < 3; k ++)
    {
        {
            double data[4][4];
            for (int i = 0; i < 4; i ++)
                for (int j = 0; j < 4; j ++)
                    data[i][j] = gsl_pow_int(i - 1, k);
            bicubic_interp *interp = bicubic_interp_init(
                *data, 4, 4, -1, -1, 1, 1);
            for (double s = 0; s <= 1; s += 0.1)
            {
                for (double t = 0; t <= 1; t += 0.1)
                {
                    const double result = bicubic_interp_eval(interp, s, t);
                    const double expected = gsl_pow_int(s, k);
                    gsl_test_abs(result, expected, 10 * GSL_DBL_EPSILON,
                        "testing bicubic interpolant for s^%d input", k);
                }
            }
            assert(interp);
            bicubic_interp_free(interp);
        }

        {
            double data[4][4];
            for (int i = 0; i < 4; i ++)
                for (int j = 0; j < 4; j ++)
                    data[i][j] = gsl_pow_int(j - 1, k);
            bicubic_interp *interp = bicubic_interp_init(
                *data, 4, 4, -1, -1, 1, 1);
            for (double s = 0; s <= 1; s += 0.1)
            {
                for (double t = 0; t <= 1; t += 0.1)
                {
                    const double result = bicubic_interp_eval(interp, s, t);
                    const double expected = gsl_pow_int(t, k);
                    gsl_test_abs(result, expected, 10 * GSL_DBL_EPSILON,
                        "testing bicubic interpolant for t^%d input", k);
                }
            }
            assert(interp);
            bicubic_interp_free(interp);
        }

        {
            double data[4][4];
            for (int i = 0; i < 4; i ++)
                for (int j = 0; j < 4; j ++)
                    data[i][j] = gsl_pow_int(i - 1, k) + gsl_pow_int(j - 1, k);
            bicubic_interp *interp = bicubic_interp_init(
                *data, 4, 4, -1, -1, 1, 1);
            for (double s = 0; s <= 1; s += 0.1)
            {
                for (double t = 0; t <= 1; t += 0.1)
                {
                    const double result = bicubic_interp_eval(interp, s, t);
                    const double expected = gsl_pow_int(s, k)
                                          + gsl_pow_int(t, k);
                    gsl_test_abs(result, expected, 10 * GSL_DBL_EPSILON,
                        "testing bicubic interpolant for s^%d + t^%d input",
                        k, k);
                }
            }
            gsl_test(!gsl_isnan(bicubic_interp_eval(interp, 0, GSL_NAN)),
                "testing that bicubic interpolant for nan input returns nan");
            gsl_test(!gsl_isnan(bicubic_interp_eval(interp, GSL_NAN, 0)),
                "testing that bicubic interpolant for nan input returns nan");
            assert(interp);
            bicubic_interp_free(interp);
        }
    }

    return gsl_test_summary();
}
