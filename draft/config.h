#ifndef __GB_config_H
#define __GB_config_H

#define BLOCK_CSR
#define AVX
#define CONV_A
#define CONV_C

#include "malloc.h"

/* #ifdef AVX */
#include <immintrin.h>
/* #endif */

#ifndef _OPENMP
#ifdef BLOCK_CSR
#define MAT_SUB_BLK 8               /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*8) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#else
#define MAT_SUB_BLK 16              /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*8) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#endif
#else
#ifdef BLOCK_CSR
#define MAT_SUB_BLK 16              /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*32) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#else
#define MAT_SUB_BLK 16              /* process submatrix MAT_SUB_BLK by MAT_SUB_BLK */
#define MAT_ROW_BLK (MAT_SUB_BLK*32) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#endif
#endif
#define USE_B_SPARSE
#define DEROULE
#define UNRL 4
#ifdef AVX
#if (UNRL != 4)
#error "avx requires UNRL = 4"
#endif
#endif
#undef USE_SAXPY
#define USE_SAXPY2
#undef USE_SAXPYn
#undef STATS
#define ALIGN(x) DIVIDE_INTO((x),UNRL)*UNRL


#endif /* __GB_config_H */

/* vim: set ft=c: */
