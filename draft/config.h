/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Brice Boyer <brice.boyer@lip6.fr>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */



#ifndef __GB_config_H
#define __GB_config_H

#undef ONLY_FFLAS

#define AVX
#undef SSE
#if defined(AVX) || defined(SSE)
#define SIMD
#endif /* SIMD */

#include <malloc.h>
#include <immintrin.h>

#define BLOCK_CSR

#ifndef _OPENMP
#ifdef BLOCK_CSR
#define MAT_SUB_BLK 8               /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*8) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#else
#define MAT_SUB_BLK 16              /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*8) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#endif /* BLOCK_CSR */
#else /* OPENMP present */
#ifdef BLOCK_CSR
#define MAT_SUB_BLK 16              /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*32) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#else
#define MAT_SUB_BLK 16              /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*32) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#endif /* BLOCK_CSR */
#endif /* OPENMP */

#ifdef SSE
#define UNRL 2
#endif

#ifdef AVX
#define UNRL 4
#endif

#ifndef SIMD
#define UNRL 4
#endif

#define CONV_A
#define CONV_C
#define USE_B_SPARSE
#define DEROULE
#undef USE_SAXPY
#define USE_SAXPY2
#undef USE_SAXPYn
#define STATS

#ifndef SIMD
#warning "you should enable simd operations"
#endif



#endif /* __GB_config_H */

/* vim: set ft=c: */
