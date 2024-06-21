/*
 * Copyright (C) 2024  Leo Singer
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

#include "find_floor.h"

#include <stddef.h>
#include <gsl/gsl_test.h>


long int find_floor(const long int *a, long int x, long int n) {
    long int base = 0;
    while (n > 0) {
        long int mid = base + n / 2;
        if (x < a[mid]) {
            n /= 2;
        } else {
            base = mid + 1;
            n -= n / 2 + 1;
        }
    }
    return base - 1;
}


int find_floor_test(void)
{
    const long int a[] = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18};
    long int n = sizeof(a) / sizeof(a[0]);

    gsl_test_int(
        find_floor(NULL, -1, 0), -1,
        "find_floor returns -1 for empty array");

    for (long int x = -2; x < 0; x ++) {
        gsl_test_int(
            find_floor(a, x, n), -1,
            "find_floor(range(0, 20, 2), %d)", x);
    }

    for (long int x = 0; x < 20; x ++)
        gsl_test_int(
            find_floor(a, x, n), x / 2,
            "find_floor(range(0, 20, 2), %d)", x);

    for (long int x = 0; x < 20; x ++)
        gsl_test_int(
            find_floor(a, x, n), x / 2,
            "find_floor(range(0, 18, 2), %d)", x);

    return gsl_test_summary();
}
