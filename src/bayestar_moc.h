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


#ifndef BAYESTAR_MOC_H
#define BAYESTAR_MOC_H

/* exclude from SWIG interface and C++ code */
/* FIXME: maybe use GSL vector and matrix types that SWIG can wrap easily */
#if !defined(SWIG) && !defined(__cplusplus)

#include <stdint.h>
#include <stddef.h>

/* All of these functions should eventually be contributed to HEALPix. */

/* Convert a NESTED pixel index to NUNIQ ordering. */
uint64_t nest2uniq64(uint8_t order, uint64_t nest);

int8_t uniq2order64(uint64_t uniq);

double uniq2pixarea64(uint64_t uniq);

/* Convert a NUNIQ pixel index to NESTED ordering. */
int8_t uniq2nest64(uint64_t uniq, uint64_t *nest);

void pix2ang_uniq64(uint64_t uniq, double *theta, double *phi);

void *moc_rasterize64(const void *pixels, size_t offset, size_t itemsize, size_t len, size_t *npix);

#endif /* !defined(SWIG) && !defined(__cplusplus) */

#endif /* BAYESTAR_MOC_H */
