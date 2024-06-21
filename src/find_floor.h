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


#ifndef FIND_FLOOR_H
#define FIND_FLOOR_H

#ifndef __cplusplus

/**
 * Find the index of the greatest item in a that is less than or equal to x.
 * a must be a sorted array of length n.
 *
 * If the array has length 0, or if all of its elements are greater than x,
 * then return -1.
 */
long int find_floor(const long int *a, long int x, long int n);

int find_floor_test(void);

#endif /* __cplusplus */

#endif /* FIND_FLOOR_H */
