/*
 * Copyright (C) 2017-2019  Leo Singer
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


#ifndef BAYESTAR_MOC_H
#define BAYESTAR_MOC_H

#ifndef __cplusplus

#include <stdint.h>
#include <stddef.h>

/* All of these functions should eventually be contributed to HEALPix. */

/* Convert a NESTED pixel index to NUNIQ ordering. */
__attribute__ ((const))
int64_t nest2uniq64(uint8_t order, int64_t nest);

__attribute__ ((const))
int8_t uniq2order64(int64_t uniq);

__attribute__ ((const))
double uniq2pixarea64(int64_t uniq);

/* Convert a NUNIQ pixel index to NESTED ordering. */
int8_t uniq2nest64(int64_t uniq, int64_t *nest);

void uniq2ang64(int64_t uniq, double *theta, double *phi);

__attribute__ ((malloc))
void *moc_rasterize64(const void *pixels, size_t offset, size_t in_stride, size_t out_stride, size_t len, size_t *npix, int8_t order);

#endif /* __cplusplus */

#endif /* BAYESTAR_MOC_H */
