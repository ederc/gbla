/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
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

#include "elimination.h"

#ifndef min
#define min(a,b) \
          ({ __typeof__ (a) _a = (a); \
                    __typeof__ (b) _b = (b); \
                    _a < _b ? _a : _b; })
#endif

#define DDEBUG 0
#define DDEBUG_C 0
#define DDEBUG_D 0
#define DDEBUG_DD 0
#define DDEBUG_DONE 0
#define DDEBUG_D_ONE 0
#define DENSE_MULTILINE_C_AFTER_ELIM  0

#define TASKS 0
#define GBLA_WITH_FFLAS 0

#ifdef GBLA_USE_AVX
#include <immintrin.h>

#define NOTALIGNED(a) \
  ((long)(a) % ALIGNT)

#if (!defined(NDEBUG) || defined(DEBUG) || defined(DDEBUG))
#define CHECK_ALIGN(a) \
  if ( NOTALIGNED(a) ) printf("not aligned %d\n", __LINE__)
#else
#define CHECK_ALIGN(a) (a)
#endif

#ifdef HAVE_AVX
#define ELEM_4  __m256d
#define ELEM_2  __m128d


#define STORE_L_2(a,b) \
  _mm_storel_pd(a,b)

#define STORE_H_2(a,b) \
  _mm_storeh_pd(a,b)


#define COPY_2(y,x) \
  STORE_2(y, LOAD_2(x))

#define SET1_4(a) \
  _mm256_set1_pd(a)
#define SET_4(a,b,c,d) \
  _mm256_set_pd((d), (c), (b), (a))

/* can do better if FMA */

#define STORE_4(a,b) \
  _mm256_store_pd((a), (b))
#define LOAD_4(a) \
  _mm256_load_pd((a))
#define COPY_4(y,x) \
  STORE_4(y, LOAD_4(x))
#define GREATER_4(a,b) \
  _mm256_cmp_pd((a), (b), _CMP_GT_OS)
#define LESSER_4(a,b) \
  _mm256_cmp_pd((a), (b),    _CMP_LT_OS)
#define VAND_4(a,b) \
  _mm256_and_pd((a), (b))
#define VOR_4(a,b) \
  _mm256_or_pd((a), (b))
#define ADD_4(a,b) \
  _mm256_add_pd((a), (b))
#define SUB_4(a,b) \
  _mm256_sub_pd((a), (b))
#define MUL_4(a,b) \
  _mm256_mul_pd((a), (b))
#define FLOOR_4(a) \
  _mm256_floor_pd((a))
#define ZERO_4() \
  _mm256_setzero_pd()

/* can do better if FMA */
#define FNMADD_4(c,a,b) \
  SUB_4((c), MUL_4((a), (b)))
#define FMADD_4(c,a,b) \
  ADD_4((c), MUL_4((a), (b)))

#define NORML_MOD_4(C, P, NEGP, MIN, MAX, Q, T) \
{ \
  Q = GREATER_4(C,MAX); \
  T = LESSER_4(C,MIN); \
  Q = VAND_4(Q,  NEGP); \
  T = VAND_4(T, P); \
  Q = VOR_4(Q, T); \
  C = ADD_4(C, Q); \
}

#define FLOAT_MOD_4(C, P, INVP,Q)\
{\
  Q = MUL_4(C,INVP);\
  Q = FLOOR_4(Q);\
  C = FNMADD_4(C, Q,P);\
}

#else
#error "you want avx..."
#endif /* HAVE_AVX */


#endif /* SIMD */

#define DEBUG_NEW_ELIM 0

int elim_fl_A_hybrid_block(hbm_fl_t **A_in, hbm_fl_t *B, mod_t modulus, int nthrds)
{
  hbm_fl_t *A = *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ci_t clB  = get_number_hybrid_col_blocks(B);
  const ci_t clA  = get_number_hybrid_col_blocks(A);
  const ri_t rlA  = get_number_hybrid_row_blocks(A);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<clB; ++i) {
#pragma omp task
      {
        rc  = elim_fl_A_hybrid_blocks_task(A, B, i, rlA, modulus);
      }
    }
#pragma omp taskwait
  }
  // free A
  free_hybrid_submatrix(&A, nthrds);
  return 0;
}

int elim_fl_A_hybrid_blocks_task(hbm_fl_t *A, hbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus)
{
  bi_t i, ctr;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_A; ++j) {
    ctr = 0;
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (B->blocks[j][block_col_idx_B] != NULL) {
      ctr = 1;
      copy_hybrid_to_wide_block(&B->blocks[j][block_col_idx_B], wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<j; ++k) {
      /*
      printf("k %d | lci %d\n",k,block_col_idx_B);
      if (k==22 && j==23) {
        int bi = k;
        int bj= block_col_idx_B;
        printf("B[%d][%d] ----------------------------------------\n",bi,bj);
        if (B->blocks[bi][bj] != NULL) {
          for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
            if (B->blocks[bi][bj][ii] == NULL) {
              for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
                printf("0 | ");
              }
            } else {
              for (int jj=0; jj<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++jj) {
                if (B->blocks[bi][bj][ii][jj].val == NULL) {
                  for (int kk=0; kk<__GBLA_SIMD_INNER_SIZE; ++kk) {
                    printf("0 | ");
                  }
                } else {
                  for (int kk=0; kk<__GBLA_SIMD_INNER_SIZE; ++kk) {
                    printf("%ld | ",B->blocks[bi][bj][ii][jj].val[kk]);
                  }
                }
              }
            }
            printf("\n");
          }
        }
      }
      */
      if ((A->blocks[j][k] != NULL) && (B->blocks[k][block_col_idx_B] != NULL)) {
        ctr = 1;
        red_hybrid_rectangular(A->blocks[j][k], B->blocks[k][block_col_idx_B],
            wide_block);
        // do the diagonal block from A
      }
      /*
      if (j==23) {
        printf("RECT %d < %d----------------------------------------\n",k,j);
        for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
          for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
            printf("%ld | ", wide_block[ii][jj]);
          }
          printf("\n");
        }
      }
      */
    }
    if (ctr == 1)
      red_hybrid_triangular(A->blocks[j][j], wide_block, modulus);
    /*
    printf("ctr %d -- TRIANGULAR %d----------------------------------------\n",ctr,j);
    for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
      for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
        printf("%ld | ", wide_block[ii][jj]);
      }
      printf("\n");
    }
    */
    copy_wide_to_hybrid_block(wide_block, &B->blocks[j][block_col_idx_B]);

#if DDDEBUG
    printf("after copying\n");
    if (B->blocks[j][block_col_idx_B] != NULL) {
      for (int kk=0; kk<B->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (B->blocks[j][block_col_idx_B][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<B->blocks[j][block_col_idx_B][kk].sz; ++ll) {
            printf("%d %d ",B->blocks[j][block_col_idx_B][kk].val[2*ll], B->blocks[j][block_col_idx_B][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }
  free_wide_block(&wide_block);

  return 0;
}

int elim_fl_A_sparse_dense_block(sb_fl_t **A_in, dbm_fl_t *B, mod_t modulus, int nthrds)
{
  sb_fl_t *A  = *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ci_t clB  = get_number_dense_col_blocks(B);
  const ci_t clA  = get_number_sparse_col_blocks(A);
  const ri_t rlA  = get_number_sparse_row_blocks(A);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=clB; i>0; --i) {
#pragma omp task
      {
        rc  = elim_fl_A_sparse_dense_blocks_task(A, B, i-1, rlA, modulus);
      }
    }
#pragma omp taskwait
  }
  // free A
  free_sparse_submatrix(&A, nthrds);
  *A_in = A;
  return 0;
}

int elim_fl_A_sparse_dense_blocks_task(const sb_fl_t *A, dbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus)
{
  bi_t i, ctr;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_A; ++j) {
    ctr = 0;
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (B->blocks[j][block_col_idx_B].val != NULL) {
      ctr = 1;
      copy_dense_to_wide_block(B->blocks[j][block_col_idx_B].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<j; ++k) {
      /*
         printf("k %d | lci %d\n",k,block_col_idx_B);
         if (k==22 && j==23) {
         int bi = k;
         int bj= block_col_idx_B;
         printf("B[%d][%d] ----------------------------------------\n",bi,bj);
         if (B->blocks[bi][bj].val != NULL) {
         for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
         for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
         printf("%ld | ", B->blocks[bi][bj].val[ii*__GBLA_SIMD_BLOCK_SIZE+jj]);
         }
         printf("\n");
         }
         }
         }
         */
      if ((A->blocks[j][k].val != NULL) && (B->blocks[k][block_col_idx_B].val != NULL)) {
        ctr = 1;
        red_sparse_dense_rectangular(&A->blocks[j][k],
            B->blocks[k][block_col_idx_B].val, wide_block);
        // do the diagonal block from A
        /*
        printf("RECT %d < %d----------------------------------------\n",k,j);
        for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
          for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
            printf("%ld | ", wide_block[ii][jj]);
          }
          printf("\n");
        }
        */
      }
    }
    if (ctr == 1)
      red_sparse_triangular(&A->blocks[j][j], wide_block, modulus);
    /*
    printf("ctr %d -- TRIANGULAR %d----------------------------------------\n",ctr,j);
    for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
      for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
        printf("%ld | ", wide_block[ii][jj]);
      }
      printf("\n");
    }
    */
    copy_wide_to_dense_block(wide_block, &B->blocks[j][block_col_idx_B].val);

#if DDDEBUG
    printf("after copying\n");
    if (B->blocks[j][block_col_idx_B] != NULL) {
      for (int kk=0; kk<B->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (B->blocks[j][block_col_idx_B][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<B->blocks[j][block_col_idx_B][kk].sz; ++ll) {
            printf("%d %d ",B->blocks[j][block_col_idx_B][kk].val[2*ll], B->blocks[j][block_col_idx_B][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }
  free_wide_block(&wide_block);

  return 0;
}

int elim_fl_A_hybrid_dense_block(hbm_fl_t **A_in, dbm_fl_t *B, mod_t modulus, int nthrds)
{
  hbm_fl_t *A = *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ci_t clB  = get_number_dense_col_blocks(B);
  const ci_t clA  = get_number_hybrid_col_blocks(A);
  const ri_t rlA  = get_number_hybrid_row_blocks(A);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<clB; ++i) {
#pragma omp task
      {
        rc  = elim_fl_A_hybrid_dense_blocks_task(A, B, i, rlA, modulus);
      }
    }
#pragma omp taskwait
  }
  // free A
  free_hybrid_submatrix(&A, nthrds);
  return 0;
}

int elim_fl_A_hybrid_dense_blocks_task(hbm_fl_t *A, dbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus)
{
  bi_t i, ctr;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_A; ++j) {
    ctr = 0;
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (B->blocks[j][block_col_idx_B].val != NULL) {
      ctr = 1;
      copy_dense_to_wide_block(B->blocks[j][block_col_idx_B].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<j; ++k) {
      /*
         printf("k %d | lci %d\n",k,block_col_idx_B);
         if (k==22 && j==23) {
         int bi = k;
         int bj= block_col_idx_B;
         printf("B[%d][%d] ----------------------------------------\n",bi,bj);
         if (B->blocks[bi][bj].val != NULL) {
         for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
         for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
         printf("%ld | ", B->blocks[bi][bj].val[ii*__GBLA_SIMD_BLOCK_SIZE+jj]);
         }
         printf("\n");
         }
         }
         }
         */
      if ((A->blocks[j][k] != NULL) && (B->blocks[k][block_col_idx_B].val != NULL)) {
        ctr = 1;
        red_hybrid_dense_rectangular(A->blocks[j][k],
            B->blocks[k][block_col_idx_B].val, wide_block);
        // do the diagonal block from A
      }
      /*
        printf("RECT %d < %d----------------------------------------\n",k,j);
        for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
          for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
            printf("%ld | ", wide_block[ii][jj]);
          }
          printf("\n");
        }
        */
    }
    if (ctr == 1)
      red_hybrid_triangular(A->blocks[j][j], wide_block, modulus);
    /*
    printf("ctr %d -- TRIANGULAR %d----------------------------------------\n",ctr,j);
    for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
      for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
        printf("%ld | ", wide_block[ii][jj]);
      }
      printf("\n");
    }
    */
    copy_wide_to_dense_block(wide_block, &B->blocks[j][block_col_idx_B].val);

#if DDDEBUG
    printf("after copying\n");
    if (B->blocks[j][block_col_idx_B] != NULL) {
      for (int kk=0; kk<B->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (B->blocks[j][block_col_idx_B][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<B->blocks[j][block_col_idx_B][kk].sz; ++ll) {
            printf("%d %d ",B->blocks[j][block_col_idx_B][kk].val[2*ll], B->blocks[j][block_col_idx_B][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }
  free_wide_block(&wide_block);

  return 0;
}

int elim_fl_A_dense_block(dbm_fl_t **A_in, dbm_fl_t *B, mod_t modulus, int nthrds)
{
  dbm_fl_t *A = *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ci_t clB  = get_number_dense_col_blocks(B);
  const ci_t clA  = get_number_dense_col_blocks(A);
  const ri_t rlA  = get_number_dense_row_blocks(A);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<clB; ++i) {
#pragma omp task
      {
        rc  = elim_fl_A_dense_blocks_task(A, B, i, rlA, modulus);
      }
    }
#pragma omp taskwait
  }
  // free A
  free_dense_submatrix(&A, nthrds);
  return 0;
}

int elim_fl_A_dense_blocks_task(const dbm_fl_t *A, dbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus)
{
  bi_t i, ctr;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_A; ++j) {
    ctr = 0;
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (B->blocks[j][block_col_idx_B].val != NULL) {
      ctr = 1;
      copy_dense_to_wide_block(B->blocks[j][block_col_idx_B].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<j; ++k) {
      /*
         printf("k %d | lci %d\n",k,block_col_idx_B);
         if (k==22 && j==23) {
         int bi = k;
         int bj= block_col_idx_B;
         printf("B[%d][%d] ----------------------------------------\n",bi,bj);
         if (B->blocks[bi][bj].val != NULL) {
         for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
         for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
         printf("%ld | ", B->blocks[bi][bj].val[ii*__GBLA_SIMD_BLOCK_SIZE+jj]);
         }
         printf("\n");
         }
         }
         }
         */
      if ((A->blocks[j][k].val != NULL) && (B->blocks[k][block_col_idx_B].val != NULL)) {
        ctr = 1;
        red_dense_rectangular(A->blocks[j][k].val,
            B->blocks[k][block_col_idx_B].val, wide_block);
        // do the diagonal block from A
        /*
        printf("RECT %d < %d----------------------------------------\n",k,j);
        for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
          for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
            printf("%ld | ", wide_block[ii][jj]);
          }
          printf("\n");
        }
        */
      }
    }
    if (ctr == 1)
      red_dense_triangular(A->blocks[j][j].val, wide_block, modulus);
    /*
    printf("ctr %d -- TRIANGULAR %d----------------------------------------\n",ctr,j);
    for (int ii=0; ii<__GBLA_SIMD_BLOCK_SIZE; ++ii) {
      for (int jj=0; jj<__GBLA_SIMD_BLOCK_SIZE; ++jj) {
        printf("%ld | ", wide_block[ii][jj]);
      }
      printf("\n");
    }
    */
    copy_wide_to_dense_block(wide_block, &B->blocks[j][block_col_idx_B].val);

#if DDDEBUG
    printf("after copying\n");
    if (B->blocks[j][block_col_idx_B] != NULL) {
      for (int kk=0; kk<B->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (B->blocks[j][block_col_idx_B][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<B->blocks[j][block_col_idx_B][kk].sz; ++ll) {
            printf("%d %d ",B->blocks[j][block_col_idx_B][kk].val[2*ll], B->blocks[j][block_col_idx_B][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }
  free_wide_block(&wide_block);

  return 0;
}

//#if __GBLA_COLUMN_B
int elim_fl_C_intermediate_block(sb_fl_t *B, ibm_fl_t **C_in, dbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  ibm_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = get_number_dense_col_blocks(D);
  const ci_t clC  = get_number_intermediate_col_blocks(C);
  const ri_t rlC  = get_number_intermediate_row_blocks(C);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=clD; i>0; --i) {
#pragma omp task
      {
        rc  = elim_fl_C_intermediate_blocks_task(B, C, D, i-1, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  // free C
  free_intermediate_submatrix(&C, nthrds);
  return 0;
}
int elim_fl_C_dense_sparse_block(sb_fl_t *B, dbm_fl_t **C_in, dbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  dbm_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = get_number_dense_col_blocks(D);
  const ci_t clC  = get_number_dense_col_blocks(C);
  const ri_t rlC  = get_number_dense_row_blocks(C);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=clD; i>0; --i) {
#pragma omp task
      {
        rc  = elim_fl_C_dense_sparse_blocks_task(B, C, D, i-1, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  // free C
  free_dense_submatrix(&C, nthrds);
  return 0;
}
//#endif

int elim_fl_C_sparse_sparse_block(sb_fl_t *B, sb_fl_t **C_in, dbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  sb_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = get_number_dense_col_blocks(D);
  const ci_t clC  = get_number_sparse_col_blocks(C);
  const ri_t rlC  = get_number_sparse_row_blocks(C);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=clD; i>0; --i) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_sparse_blocks_task(B, C, D, i-1, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  // free C
  free_sparse_submatrix(&C, nthrds);
  return 0;
}

int elim_fl_C_sparse_dense_block(dbm_fl_t *B, sb_fl_t **C_in, dbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  sb_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = get_number_dense_col_blocks(D);
  const ci_t clC  = get_number_sparse_col_blocks(C);
  const ri_t rlC  = get_number_sparse_row_blocks(C);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=clD; i>0; --i) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_blocks_task(B, C, D, i-1, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  // free C
  free_sparse_submatrix(&C, nthrds);
  return 0;
}

int elim_fl_C_dense_block(dbm_fl_t *B, dbm_fl_t **C_in, dbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  dbm_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = get_number_dense_col_blocks(D);
  const ci_t clC  = get_number_dense_col_blocks(C);
  const ri_t rlC  = get_number_dense_row_blocks(C);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<clD; ++i) {
#pragma omp task
      {
        rc  = elim_fl_C_dense_blocks_task(B, C, D, i, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  // free C
  free_dense_submatrix(&C, nthrds);
  return 0;
}

//#if __GBLA_COLUMN_B
int elim_fl_C_intermediate_blocks_task(sb_fl_t *B, ibm_fl_t *C, dbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_C; ++j) {
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (D->blocks[j][block_col_idx_D].val != NULL) {
      copy_dense_to_wide_block(D->blocks[j][block_col_idx_D].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<nbcols_C; ++k) {
      if (B->blocks[k][block_col_idx_D].val != NULL) {
        if (C->blocks[j][k].val != NULL) {
        //printf("%u | %u | %u\n",j,k,block_col_idx_D);
          red_dense_sparse_rectangular(C->blocks[j][k].val,
              &B->blocks[k][block_col_idx_D], wide_block);
        }
        if (C->blocks[j][k].sz != NULL) {
        //printf("%u | %u | %u\n",j,k,block_col_idx_D);
          red_intermediate_sparse_sparse_rectangular(&C->blocks[j][k],
              &B->blocks[k][block_col_idx_D], wide_block);
        }
      }
    }
    // cut down elements mod field characteristic
    modulo_wide_block(&wide_block, modulus);

    copy_wide_to_dense_block(wide_block, &D->blocks[j][block_col_idx_D].val);
  }
  free_wide_block(&wide_block);

  return 0;
}
int elim_fl_C_dense_sparse_blocks_task(sb_fl_t *B, dbm_fl_t *C, dbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_C; ++j) {
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (D->blocks[j][block_col_idx_D].val != NULL) {
      copy_dense_to_wide_block(D->blocks[j][block_col_idx_D].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<nbcols_C; ++k) {
      if ((C->blocks[j][k].val != NULL) && (B->blocks[k][block_col_idx_D].val != NULL)) {
        //printf("%u | %u | %u\n",j,k,block_col_idx_D);
        red_dense_sparse_rectangular(C->blocks[j][k].val,
            &B->blocks[k][block_col_idx_D], wide_block);
      }
    }
    // cut down elements mod field characteristic
    modulo_wide_block(&wide_block, modulus);

    copy_wide_to_dense_block(wide_block, &D->blocks[j][block_col_idx_D].val);
  }
  free_wide_block(&wide_block);

  return 0;
}
//#endif

int elim_fl_C_sparse_sparse_blocks_task(sb_fl_t *B, sb_fl_t *C, dbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_C; ++j) {
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (D->blocks[j][block_col_idx_D].val != NULL) {
      copy_dense_to_wide_block(D->blocks[j][block_col_idx_D].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<nbcols_C; ++k) {
      if ((C->blocks[j][k].val != NULL) && (B->blocks[k][block_col_idx_D].val != NULL)) {
        //printf("%u | %u | %u\n",j,k,block_col_idx_D);
        red_sparse_sparse_rectangular(&C->blocks[j][k],
            &B->blocks[k][block_col_idx_D], wide_block);
      }
    }
    // cut down elements mod field characteristic
    modulo_wide_block(&wide_block, modulus);

    copy_wide_to_dense_block(wide_block, &D->blocks[j][block_col_idx_D].val);
  }
  free_wide_block(&wide_block);

  return 0;
}

int elim_fl_C_sparse_dense_blocks_task(dbm_fl_t *B, sb_fl_t *C, dbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_C; ++j) {
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (D->blocks[j][block_col_idx_D].val != NULL) {
      copy_dense_to_wide_block(D->blocks[j][block_col_idx_D].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<nbcols_C; ++k) {
      if ((C->blocks[j][k].val != NULL) && (B->blocks[k][block_col_idx_D].val != NULL)) {
        red_sparse_dense_rectangular(&C->blocks[j][k],
            B->blocks[k][block_col_idx_D].val, wide_block);
      }
    }
    // cut down elements mod field characteristic
    modulo_wide_block(&wide_block, modulus);

    copy_wide_to_dense_block(wide_block, &D->blocks[j][block_col_idx_D].val);
  }
  free_wide_block(&wide_block);

  return 0;
}

int elim_fl_C_dense_blocks_task(dbm_fl_t *B, dbm_fl_t *C, dbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t **wide_block;
  
  init_wide_blocks(&wide_block);
  for (j=0; j<nbrows_C; ++j) {
    const ri_t first_block_idx  = 0;

    set_wide_block_to_zero(wide_block, __GBLA_SIMD_BLOCK_SIZE);

    // copy sparse block data to dense representation
    if (D->blocks[j][block_col_idx_D].val != NULL) {
      copy_dense_to_wide_block(D->blocks[j][block_col_idx_D].val, wide_block);
    }
    // do all rectangular blocks
    for (k=0; k<nbcols_C; ++k) {
      if ((C->blocks[j][k].val != NULL) && (B->blocks[k][block_col_idx_D].val != NULL)) {
        red_dense_rectangular(C->blocks[j][k].val,
            B->blocks[k][block_col_idx_D].val, wide_block);
      }
    }
    // cut down elements mod field characteristic
    modulo_wide_block(&wide_block, modulus);

    copy_wide_to_dense_block(wide_block, &D->blocks[j][block_col_idx_D].val);
  }
  free_wide_block(&wide_block);

  return 0;
}
#if TASKS
int elim_fl_A_block(sbm_fl_t **A_in, sbm_fl_t *B, mod_t modulus, int nthrds) {
  sbm_fl_t *A = *A_in;
  ci_t i/*, rc*/;
  ri_t j, k;
  const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth);
  const ri_t rlA  = (ri_t) ceil((float) A->nrows / A->bheight);
  const ci_t clA  = (ci_t) ceil((float) A->ncols / A->bwidth);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    /*  each task takes one block column of B */
    for (i=0; i<clB; ++i) {
#pragma omp task
      {
        /*rc  =*/ elim_fl_A_blocks_task(A, B, i, rlA, modulus);
      }
    }
#pragma omp taskwait
  }
  /*  free A */
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i] != NULL) {
          for (k=0; k<A->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(A->blocks[j][i][k].idx);
            A->blocks[j][i][k].idx  = NULL;
            free(A->blocks[j][i][k].val);
            A->blocks[j][i][k].val  = NULL;
          }
          free(A->blocks[j][i]);
          A->blocks[j][i] = NULL;
        }
      }
      free(A->blocks[j]);
      A->blocks[j]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
  return 0;
}

int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A, mod_t modulus) {
  int ret;
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[B->bheight] __attribute__((aligned(0x1000)));
  /* re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *)); */
  uint64_t size = B->bwidth * sizeof(re_l_t);
  for (i=0; i<B->bheight; ++i) {
    ret = posix_memalign((void **)&dense_block[i], ALIGNT, size);
  }
  for (j=0; j<nbrows_A; ++j) {
    /* const ri_t first_block_idx  = 0; */

    /*  set dense block entries to zero */
    for (k=0; k<B->bheight; ++k)
      memset(dense_block[k], 0, size);

    /*  copy sparse block data to dense representation */
    if (B->blocks[j][block_col_idx_B] != NULL)
      copy_sparse_to_dense_block(B->blocks[j][block_col_idx_B], dense_block,
          B->bheight, B->bwidth);

    for (k=0; k<j; ++k) {
      red_with_rectangular_block(A->blocks[j][k], B->blocks[k][block_col_idx_B],
          dense_block, B->bheight, 1, modulus);
      /*
         printf("RECTANGULAR DONE\n");
         for (int kk=0; kk<B->bheight; ++kk) {
         for (int ll=0; ll<B->bheight; ++ll) {
         printf("(%d,%d) %ld ",kk,ll,dense_block[kk][ll]);
         }
         printf("\n");
         }
         */
    }

    red_with_triangular_block(A->blocks[j][j], dense_block,
        B->bheight, 1, modulus);
    /*printf("TRIANGULAR DONE\n");
      for (int kk=0; kk<B->bheight; ++kk) {
      for (int ll=0; ll<B->bheight; ++ll) {
      printf("%ld ",dense_block[kk][ll]);
      }
      printf("\n");
      }
      */


    /* printf("OUT BEFORE %p\n",B->blocks[j][block_col_idx_B]); */
    copy_dense_block_to_sparse(
        dense_block, &B->blocks[j][block_col_idx_B],
        B->bheight, B->bwidth, modulus);
    /* printf("OUT AFTERWARDS %p\n",B->blocks[j][block_col_idx_B]); */
#if DDDEBUG
    printf("after copying\n");
    if (B->blocks[j][block_col_idx_B] != NULL) {
      for (int kk=0; kk<B->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (B->blocks[j][block_col_idx_B][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<B->blocks[j][block_col_idx_B][kk].sz; ++ll) {
            printf("%d %d ",B->blocks[j][block_col_idx_B][kk].val[2*ll], B->blocks[j][block_col_idx_B][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }
  for (i=0; i<B->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#else
int elim_fl_A_block(sbm_fl_t **A_in, sbm_fl_t *B, mod_t modulus, int nthrds) {
  sbm_fl_t *A = *A_in;
  ci_t i/*, rc*/;
  ri_t j, k;
  int tid;
  /* const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth); */
  const ri_t rlA  = (ri_t) ceil((float) A->nrows / A->bheight);
  const ci_t clA  = (ci_t) ceil((float) A->ncols / A->bwidth);

  reduce_A_next_col_to_reduce = 0;
  omp_init_lock(&reduce_A_lock);

#pragma omp parallel shared(reduce_A_next_col_to_reduce) num_threads(nthrds)
  {
#pragma omp for nowait
    for (tid=0; tid<nthrds; ++tid) {
      /*rc  =*/ elim_fl_A_blocks_task(A, B, tid, rlA, modulus);
    }
  }
  omp_destroy_lock(&reduce_A_lock);
  /*  free A */
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i] != NULL) {
          for (k=0; k<A->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(A->blocks[j][i][k].idx);
            A->blocks[j][i][k].idx  = NULL;
            free(A->blocks[j][i][k].val);
            A->blocks[j][i][k].val  = NULL;
          }
          free(A->blocks[j][i]);
          A->blocks[j][i] = NULL;
        }
      }
      free(A->blocks[j]);
      A->blocks[j]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
  return 0;
}

int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A, mod_t modulus) {
  int ret;
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[B->bheight] __attribute__((aligned(0x1000)));
  /* re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *)); */
  uint64_t size = B->bwidth * sizeof(re_l_t);
  for (i=0; i<B->bheight; ++i) {
    ret = posix_memalign((void **)&dense_block[i], ALIGNT, size);
  }

  const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth);
  uint32_t lci; /*  local column index */
  while (1) {
    omp_set_lock(&reduce_A_lock);
    if (reduce_A_next_col_to_reduce < clB) {
      lci = reduce_A_next_col_to_reduce;
      ++reduce_A_next_col_to_reduce;
    } else {
      omp_unset_lock(&reduce_A_lock);
      break;
    }
    omp_unset_lock(&reduce_A_lock);
    for (j=0; j<nbrows_A; ++j) {
      /* const ri_t first_block_idx  = 0; */

      /*  set dense block entries to zero */
      for (k=0; k<B->bheight; ++k) {
        memset(dense_block[k], 0, size);
      }
      /*  copy sparse block data to dense representation */
      if (B->blocks[j][lci] != NULL)
        copy_sparse_to_dense_block(B->blocks[j][lci], dense_block,
            B->bheight, B->bwidth);
      
      for (k=0; k<j; ++k) {
        /*
        printf("k %d | lci %d\n",k,lci);
        if (lci==8 && k==34) {  
          int bi=k;
          int bj=lci;
          printf("B[%d][%d] ----------------------------------------\n",bi,bj);
          for (int ii=0; ii<B->bheight/2; ++ii) {
            if (B->blocks[bi][bj] != NULL) {
            if (B->blocks[bi][bj][ii].dense == 0) {
              for (int jj=0; jj<B->blocks[bi][bj][ii].sz; ++jj) {
                printf("%d -- ", B->blocks[bi][bj][ii].idx[jj]);
                printf("%d %d ", B->blocks[bi][bj][ii].val[2*jj], B->blocks[bi][bj][ii].val[2*jj+1]);
              }
              printf("\n");
            } else {
              for (int jj=0; jj<B->blocks[bi][bj][ii].sz; ++jj) {
                printf("%d -- ", jj);
                printf("%d %d ", B->blocks[bi][bj][ii].val[2*jj], B->blocks[bi][bj][ii].val[2*jj+1]);
              }
              printf("\n");
            }
          }
          }
        }
        */
        red_with_rectangular_block(A->blocks[j][k], B->blocks[k][lci],
            dense_block, B->bheight, 1, modulus);
        /*
        if (j == 23) {
          printf("RECT %d < %d----------------------------------------\n",k,j);
          for (int ii=0; ii<B->bheight; ++ii) {
            for (int jj=0; jj<B->bheight; ++jj) {
              printf("%ld | ", dense_block[ii][jj]);
            }
            printf("\n");
          }
        }
        */
      }

      red_with_triangular_block(A->blocks[j][j], dense_block,
          B->bheight, 1, modulus);
      /*
      printf("TRIANGULAR %d----------------------------------------\n",j);
      for (int ii=0; ii<B->bheight; ++ii) {
        for (int jj=0; jj<B->bheight; ++jj) {
          printf("%ld | ", dense_block[ii][jj]);
        }
        printf("\n");
      }
      */
      copy_dense_block_to_sparse(
          dense_block, &B->blocks[j][lci],
          B->bheight, B->bwidth, modulus);
    }
  }
  for (i=0; i<B->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#endif
#if TASKS
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C_in, sbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  sbm_fl_t *C = *C_in;

  ci_t i/*, rc*/;
  ri_t j, k;
  const ci_t clD  = (ci_t) ceil((float) D->ncols / D->bwidth);
  const ri_t rlC  = (ri_t) ceil((float) C->nrows / C->bheight);
  const ri_t clC  = (ci_t) ceil((float) C->ncols / C->bwidth);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    /*  each task takes one block column of B */
    for (i=0; i<clD; ++i) {
#pragma omp task
      {
        /*rc  =*/ elim_fl_C_blocks_task(B, C, D, i, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  /*  free C */
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlC; ++j) {
      for (i=0; i<clC; ++i) {
        if (C->blocks[j][i] != NULL) {
          for (k=0; k<C->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(C->blocks[j][i][k].idx);
            C->blocks[j][i][k].idx  = NULL;
            free(C->blocks[j][i][k].val);
            C->blocks[j][i][k].val  = NULL;
          }
          free(C->blocks[j][i]);
          C->blocks[j][i] = NULL;
        }
      }
      free(C->blocks[j]);
      C->blocks[j]  = NULL;
    }
  }
  free(C->blocks);
  C->blocks = NULL;
  free(C);
  C = NULL;
  *C_in = C;
  return 0;
}

int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const int inv_scalars, const mod_t modulus) {
  int ret;
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[D->bheight] __attribute__((aligned(0x1000)));
  /* re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *)); */
  uint64_t size = D->bwidth * sizeof(re_l_t);
  for (i=0; i<D->bheight; ++i) {
    ret = posix_memalign((void **)&dense_block[i], ALIGNT, size);
  }

  /*  take maximum number and check if blocks are NULL in loop */
  const ri_t last_block_idx = nbcols_C;

  for (j=0; j<nbrows_C; ++j) {
    /* const ri_t first_block_idx  = 0; */

    /*  set dense block entries to zero */
    for (k=0; k<D->bheight; ++k)
      memset(dense_block[k], 0, size);

    /*  copy sparse block data to dense representation */
    if (D->blocks[j][block_col_idx_D] != NULL)
      copy_sparse_to_dense_block(D->blocks[j][block_col_idx_D], dense_block,
          D->bheight, D->bwidth);

    for (k=0; k<last_block_idx; ++k) {
#if DDEBUG_C
      printf("j %d -- k %d -- lci %d\n", j, k, block_col_idx_D);
#endif
      if (C->blocks[j][k] != NULL) {
        red_with_rectangular_block(C->blocks[j][k], B->blocks[k][block_col_idx_D],
            dense_block, B->bheight, inv_scalars, modulus);
      }
      /*
         printf("RECTANGULAR DONE\n");
         for (int kk=0; kk<B->bheight; ++kk) {
         for (int ll=0; ll<B->bheight; ++ll) {
         printf("(%d,%d) %ld ",kk,ll,dense_block[kk][ll]);
         }
         printf("\n");
         }
         */
    }
    /* printf("OUT BEFORE %p\n",B->blocks[j][block_col_idx_B]); */
    copy_dense_block_to_sparse(
        dense_block, &D->blocks[j][block_col_idx_D],
        D->bheight, D->bwidth, modulus);
    /* printf("OUT AFTERWARDS %p\n",B->blocks[j][block_col_idx_B]); */
#if DDEBUG_C
    printf("after copying\n");
    if (D->blocks[j][block_col_idx_D] != NULL) {
      for (int kk=0; kk<D->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (D->blocks[j][block_col_idx_D][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<D->blocks[j][block_col_idx_D][kk].sz; ++ll) {
            printf("%d %d ",D->blocks[j][block_col_idx_D][kk].val[2*ll], D->blocks[j][block_col_idx_D][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }

  for (i=0; i<D->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#else
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C_in, sbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  sbm_fl_t *C = *C_in;

  ci_t i/*, rc*/;
  int tid ;
  ri_t j, k;
  /* const ci_t clD  = (ci_t) ceil((float) D->ncols / D->bwidth); */
  const ri_t rlC  = (ri_t) ceil((float) C->nrows / C->bheight);
  const ri_t clC  = (ci_t) ceil((float) C->ncols / C->bwidth);

  reduce_C_next_col_to_reduce = 0;
  omp_init_lock(&reduce_C_lock);

#pragma omp parallel shared(reduce_C_next_col_to_reduce) num_threads(nthrds)
  {
#pragma omp for nowait
    for (tid=0; tid<nthrds; ++tid) {
      /*rc  =*/ elim_fl_C_blocks_task(B, C,  D, tid, rlC, clC, inv_scalars, modulus);
    }
  }
  omp_destroy_lock(&reduce_C_lock);
  /*  free C */
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlC; ++j) {
      for (i=0; i<clC; ++i) {
        if (C->blocks[j][i] != NULL) {
          for (k=0; k<C->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(C->blocks[j][i][k].idx);
            C->blocks[j][i][k].idx  = NULL;
            free(C->blocks[j][i][k].val);
            C->blocks[j][i][k].val  = NULL;
          }
          free(C->blocks[j][i]);
          C->blocks[j][i] = NULL;
        }
      }
      free(C->blocks[j]);
      C->blocks[j]  = NULL;
    }
  }
  free(C->blocks);
  C->blocks = NULL;
  free(C);
  C = NULL;
  *C_in = C;
  return 0;
}

int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const int inv_scalars, const mod_t modulus) {

  const ci_t clD  = (ci_t) ceil((float) D->ncols / D->bwidth);
  int ret;
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[D->bheight] __attribute__((aligned(0x1000)));
  /* re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *)); */
  uint64_t size = D->bwidth * sizeof(re_l_t);
  for (i=0; i<D->bheight; ++i) {
    ret = posix_memalign((void **)&dense_block[i], ALIGNT, size);
  }

  uint32_t lci; /*  local column index */
  /*  take maximum number and check if blocks are NULL in loop */
  const ri_t last_block_idx = nbcols_C;

  while (1) {
    omp_set_lock(&reduce_C_lock);
    if (reduce_C_next_col_to_reduce < clD) {
      lci = reduce_C_next_col_to_reduce;
      ++reduce_C_next_col_to_reduce;
    } else {
      omp_unset_lock(&reduce_C_lock);
      break;
    }
    omp_unset_lock(&reduce_C_lock);

    for (j=0; j<nbrows_C; ++j) {
      /* const ri_t first_block_idx  = 0; */

      /*  set dense block entries to zero */
      for (k=0; k<D->bheight; ++k)
        memset(dense_block[k], 0, size);

      /*  copy sparse block data to dense representation */
      if (D->blocks[j][lci] != NULL) {
        copy_sparse_to_dense_block(D->blocks[j][lci], dense_block,
            D->bheight, D->bwidth);
      }
      for (k=0; k<last_block_idx; ++k) {
#if DDEBUG_C
        printf("j %d -- k %d -- lci %d\n", j, k, block_col_idx_D);
#endif
        if (C->blocks[j][k] != NULL) {
          red_with_rectangular_block(C->blocks[j][k], B->blocks[k][lci],
              dense_block, B->bheight, inv_scalars, modulus);
        }
      }
      copy_dense_block_to_sparse(
          dense_block, &D->blocks[j][lci],
          D->bheight, D->bwidth, modulus);
    }
  }

  for (i=0; i<D->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#endif

void red_with_triangular_block(mbl_t *block_A, re_l_t **dense_block,
    const ri_t bheight, const int inv_scalars, const mod_t modulus) {
  int  j;
  ri_t i, k ;

  int last_idx;

  for (i=0; i<bheight/2; ++i) {
    if (block_A[i].sz == 0)
      continue;

    last_idx  = -1;
    if (block_A[i].val[2*(block_A[i].sz-1)+1] == 0)
      last_idx  = block_A[i].sz-1;
    else
      last_idx  = block_A[i].sz-2;

    register re_m_t Av1_col1, Av2_col1;
    register re_m_t Av1_col2, Av2_col2;
    bi_t Ap1, Ap2;

#if DDDEBUG
    printf("lidx %d\n",last_idx);
#endif
    for (j=0; j<last_idx; ++j) {
      Ap1       = block_A[i].idx[j];
      Av1_col1  = block_A[i].val[2*j];
      Av2_col1  = block_A[i].val[2*j+1];

      if (inv_scalars == 1) {
        if (Av1_col1 != 0)
          Av1_col1  = (re_m_t)modulus - Av1_col1;
        if (Av2_col1 != 0)
          Av2_col1  = (re_m_t)modulus - Av2_col1;
      }

      if ((Ap1 % 2) == 0 && (j < last_idx-1)) {
        Ap2  = block_A[i].idx[j+1];
        if (Ap2 == Ap1+1) { /*  AXPY two rows */
          Av1_col2  = block_A[i].val[2*(j+1)];
          Av2_col2  = block_A[i].val[2*(j+1)+1];

          if (inv_scalars == 1) {
            if (Av1_col2 != 0)
              Av1_col2  = (re_m_t)modulus - Av1_col2;
            if (Av2_col2 != 0)
              Av2_col2  = (re_m_t)modulus - Av2_col2;
          }
          ++j;

          dense_scal_mul_sub_2_rows_array_array(
              Av1_col1, Av2_col1, Av1_col2, Av2_col2, bheight,
              dense_block[Ap1], dense_block[Ap1+1],
              dense_block[2*i], dense_block[2*i+1]);
        } else { /*  AXPY one row */
          dense_scal_mul_sub_1_row_array_array(
              Av1_col1, Av2_col1, bheight,
              dense_block[Ap1],
              dense_block[2*i], dense_block[2*i+1]);
        }
      } else { /*  AXPY one row */
        dense_scal_mul_sub_1_row_array_array(
            Av1_col1, Av2_col1, bheight,
            dense_block[Ap1],
            dense_block[2*i], dense_block[2*i+1]);
      }
    }

    /*  do modular reduction on dense array row */
    /* printf("outside %d\n", 2*i); */
    red_dense_array_modular(dense_block[2*i], bheight, modulus);

    /*  reduce lines within the same multiline */
    if (block_A[i].sz > 1) {
      j         = block_A[i].sz-2;
      Av1_col1  = block_A[i].val[2*j+1];
      Ap1       = block_A[i].idx[j];

      if(Av1_col1 != 0) {
        if (inv_scalars == 1)
          if (Av1_col1 != 0)
            Av1_col1  = (re_m_t)modulus - Av1_col1;

        const bi_t offset1  = 2*i+1;
        const bi_t offset2  = 2*i;

#ifdef GBLA_USE_AVX_XXX
        for (k=0; k<bheight; ++k)
          dense_block[offset1][k] +=  (re_m_t) Av1_col1 * dense_block[offset2][k];
#else
        for (k=0; k<bheight; ++k)
          dense_block[offset1][k] +=  (re_m_t) Av1_col1 * dense_block[offset2][k];
#endif

      }
    }

    /*  do modular reduction on dense array row */
    /* printf("outside2 %d\n", 2*i+1); */
    red_dense_array_modular(dense_block[2*i+1], bheight, modulus);
  }
}

void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_block,
    const ri_t bheight, const int inv_scalars, const mod_t modulus) {
  bi_t i, j;
  if (block_A == NULL || block_B == NULL)
    return;

  for (i=0; i<bheight/2; ++i) {
    const bi_t is_sparse  = block_A[i].dense == 0 ? 1 : 0;
    const bi_t N          = is_sparse == 1 ? block_A[i].sz : bheight;
    for (j=0; j<N; ++j) {
      /* printf("%d && %d\n",i,j); */
      /*
         for (int kk=0; kk<bheight; ++kk) {
         for (int ll=0; ll<bheight; ++ll) {
         printf("%ld ",dense_block[kk][ll]);
         }
         printf("\n");
         }
         */
      /* printf("%d // %d // %d\n",i,j, N); */
      const bi_t Ap1  = is_sparse == 1 ? block_A[i].idx[j] : j;
      register re_m_t Av1_col1  = block_A[i].val[2*j];
      register re_m_t Av2_col1  = block_A[i].val[2*j+1];

      if (inv_scalars == 1) {
        if (Av1_col1 != 0)
          Av1_col1  = (re_m_t)modulus - Av1_col1;
        if (Av2_col1 != 0)
          Av2_col1  = (re_m_t)modulus - Av2_col1;
      }
      if (((Ap1 % 2) == 0) && (j < (uint32_t)(N - 1))) {
        const bi_t Ap2  = (is_sparse == 1 ? block_A[i].idx[j+1] : j+1);
        if (Ap2 == Ap1+1) { /*  AXPY two rows */
          register re_m_t Av1_col2  = block_A[i].val[2*(j+1)];
          register re_m_t Av2_col2  = block_A[i].val[2*(j+1)+1];

          if (inv_scalars == 1) {
            if (Av1_col2 != 0)
              Av1_col2  = (re_m_t)modulus - Av1_col2;
            if (Av2_col2 != 0)
              Av2_col2  = (re_m_t)modulus - Av2_col2;
          }
          ++j;
          if (block_B[Ap1 / __GB_NROWS_MULTILINE].dense == 0) {
            /* printf("1S Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE); */
            sparse_scal_mul_sub_2_rows_vect_array(
                Av1_col1, Av2_col1, Av1_col2, Av2_col2,
                block_B[Ap1 / __GB_NROWS_MULTILINE],
                dense_block[2*i], dense_block[2*i+1]);
            /* printf("1 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]); */
          } else {
            /* printf("1D Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE); */
            dense_scal_mul_sub_2_rows_vect_array(
                Av1_col1, Av2_col1, Av1_col2, Av2_col2,
                block_B[Ap1 / __GB_NROWS_MULTILINE], bheight,
                dense_block[2*i], dense_block[2*i+1]);
            /* printf("2 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]); */
          }
        } else { /*  AXPY one row */
          if (block_B[Ap1 / __GB_NROWS_MULTILINE].dense == 0) {
            /* printf("2S Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE); */
            sparse_scal_mul_sub_1_row_vect_array(
                Av1_col1, Av2_col1,
                block_B[Ap1 / __GB_NROWS_MULTILINE],
                Ap1 % __GB_NROWS_MULTILINE,
                dense_block[2*i], dense_block[2*i+1]);
            /* printf("3 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]); */
          } else {
            /* printf("2D Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE); */
            dense_scal_mul_sub_1_row_vect_array(
                Av1_col1, Av2_col1,
                block_B[Ap1 / __GB_NROWS_MULTILINE],
                Ap1 % __GB_NROWS_MULTILINE, bheight,
                dense_block[2*i], dense_block[2*i+1]);
            /* printf("4 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]); */
          }
        }
      } else { /*  AXPY one row */
        if (block_B[Ap1 / __GB_NROWS_MULTILINE].dense == 0) {
          /* printf("3S Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE); */
          sparse_scal_mul_sub_1_row_vect_array(
              Av1_col1, Av2_col1,
              block_B[Ap1 / __GB_NROWS_MULTILINE],
              Ap1 % __GB_NROWS_MULTILINE,
              dense_block[2*i], dense_block[2*i+1]);
          /* printf("5 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]); */
        } else {
          /* printf("3D Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE); */
          dense_scal_mul_sub_1_row_vect_array(
              Av1_col1, Av2_col1,
              block_B[Ap1 / __GB_NROWS_MULTILINE],
              Ap1 % __GB_NROWS_MULTILINE, bheight,
              dense_block[2*i], dense_block[2*i+1]);
        }
      }
    }
    /*
    if (2*i == 83 || 2*i+1 == 83) {
      printf("---> %lu\n",dense_block[83][0]);
    }
    */
  }
}

void copy_dense_block_to_sparse(
    re_l_t **dense_block, mbl_t **sparse_block_in, int bheight, int bwidth, mod_t modulus) {

  mbl_t *sparse_block = *sparse_block_in;
  bi_t i,j, k, ctr, buffer, rows_empty;

  /*  if block was empty in the beginning, reallocate memory */
  if (sparse_block == NULL) {
    sparse_block  = (mbl_t *)malloc(bheight / __GB_NROWS_MULTILINE * sizeof(mbl_t));
    for (i=0; i<(bheight / __GB_NROWS_MULTILINE); ++i) {
      sparse_block[i].val = NULL;
      sparse_block[i].idx = NULL;
      sparse_block[i].sz  = sparse_block[i].dense =  0;
    }
  }
  rows_empty  = 0;
  for (i=0; i<bheight/2; ++i) {
    if (sparse_block[i].dense == 1) {
      sparse_block[i].idx = (bi_t *)malloc(bwidth * sizeof(bi_t));
    }
    ctr                 = 0;
    buffer              = sparse_block[i].sz;
    sparse_block[i].sz  = 0;
    sparse_block[i].dense = 0;
    for (j=0; j<bwidth; ++j) {
      if (dense_block[2*i][j] != 0) {
        dense_block[2*i][j] = (re_l_t)(MODP(dense_block[2*i][j], modulus));
      }
      if (dense_block[2*i+1][j] != 0) {
        dense_block[2*i+1][j] = (re_l_t)(MODP(dense_block[2*i+1][j] , modulus));
      }
      if (dense_block[2*i][j] != 0 || dense_block[2*i+1][j] != 0) {
        if (ctr >= buffer) {
          /*  if this happens just allocate memory for full block multiline row */
          sparse_block[i].idx = (bi_t*) realloc(
              sparse_block[i].idx, bwidth * sizeof(bi_t));
          sparse_block[i].val = (re_t*) realloc(
              sparse_block[i].val, 2 * bwidth * sizeof(re_t));
          buffer  = bwidth;
        }
        sparse_block[i].idx[ctr]      = j;
        sparse_block[i].val[2*ctr]    = (re_t) dense_block[2*i][j];
        sparse_block[i].val[2*ctr+1]  = (re_t) dense_block[2*i+1][j];
        /* printf("%d-vals %d .. %d   ",ctr, sparse_block[i].val[2*ctr], sparse_block[i].val[2*ctr+1]); */
        ctr++;
      }
    }
    sparse_block[i].sz  = ctr;
    /*  try to get hybrid representation, i.e. try to make multiline row dense */

    if ((float)sparse_block[i].sz / (float)bwidth < __GB_HYBRID_THRESHOLD) {
      /*  realloc memory, cut it down as much as possible */
      if (sparse_block[i].sz>0) {
        sparse_block[i].idx = (bi_t*) realloc(
            sparse_block[i].idx,
            sparse_block[i].sz * sizeof(bi_t));
        sparse_block[i].val = (re_t*) realloc(
            sparse_block[i].val,
            2 * sparse_block[i].sz  * sizeof(re_t));
      } else {
        free(sparse_block[i].idx);
        sparse_block[i].idx = NULL;
        free(sparse_block[i].val);
        sparse_block[i].val = NULL;
        rows_empty++;
      }
    } else { /*  dense multiline row */
      re_t *tmp  = (re_t *)malloc(2 * bwidth * sizeof(re_t));
      ctr  = 0;
      /* for (k=0; k<bwidth; ++k) { */
      k = 0;
      while (ctr<sparse_block[i].sz) {
        if (sparse_block[i].idx[ctr] == k) {
          tmp[2*k]   = sparse_block[i].val[2*ctr];
          tmp[2*k+1] = sparse_block[i].val[2*ctr+1];
          ctr++;
          k++;
        } else {
          tmp[2*k]   = 0;
          tmp[2*k+1] = 0;
          k++;
        }
      }
      for (; k<bwidth; ++k) {
        tmp[2*k]   = 0;
        tmp[2*k+1] = 0;
      }
      free(sparse_block[i].idx);
      sparse_block[i].idx   = NULL;
      free(sparse_block[i].val);
      sparse_block[i].val   = tmp;
      sparse_block[i].sz    = bwidth;
      sparse_block[i].dense = 1;
    }
    /* } */
  }
  /*  if block is completely empty remove memory */
  if (rows_empty == bheight/2) {
    free(sparse_block);
    sparse_block  = NULL;
  }
  *sparse_block_in  = sparse_block;
}


#if GBLA_WITH_FFLAS
ri_t elim_fl_D_fflas_ffpack(sbm_fl_t *D_old, mod_t modulus, int nthrds) {

  /*  generate DNS matrix D out of D_old */
  SAFE_MALLOC_DECL(D, 1, DNS);
  initDenseUnit(D);
  copyMetaData(D, D_old, modulus);

  /*  dense representation of D, alloc memory at once and set entries to zero */
  SAFE_CALLOC(D->ptr, (index_t)D->row * (index_t)D->ld, elemt_t);

  /*  copy multiline block matrix to DNS format */
  copy_block_ml_matrix_to_dns_matrix(&D_old, &D);

  /*  row reduce D with FFLAS-FFPACK */
  ri_t rank;
#if GBLA_WITH_FFLAS
  rank  = Mjoin(RowReduce,elemt_t)(D->mod,D->ptr,D->row,D->col,D->ld, nthrds);
#endif
  /*
     size_t i, j, k;
     ri_t ctr = 0;
     for (i=0; i<D->row; ++i) {
     ctr = 0;
     printf("\n%d || %d\n",i,ctr/256);
     for (j=0; j<D->col; ++j) {
     printf("%.1f ",D->ptr[i*D->ld+j]);
     ctr++;
     if (ctr % 256 == 0)
     printf("\n%d || %d\n",i,ctr/256);
     }
     }
     */

  return rank;
}
#endif

int elim_fl_dense_D_tasks(dm_t *D)
{
  ri_t curr_row_to_reduce;
  ri_t local_last_piv;
  ri_t local_last_row_fully_reduced;
  ri_t from_row;

  int ready_for_waiting_list  = 0;
  int curr_row_fully_reduced  = 0;
  int nred_consecutively      = 0;
  int ret;

  ri_t wl_idx = 0;
  ri_t wl_lp  = 0;

#if DEBUG_NEW_ELIM
  int tid = omp_get_thread_num();
#endif

  // computations done in while loop
  while (1) {
    if (global_last_row_fully_reduced >= D->nrows) {
      omp_unset_lock(&echelonize_lock);
      break;
    }
    local_last_piv                = global_last_piv;
    local_last_row_fully_reduced  = global_last_row_fully_reduced;
#if DEBUG_NEW_ELIM
    printf("setting llp\n");
    printf("llp %d\n",local_last_piv);
    printf("grr %d\n",global_next_row_to_reduce);
#endif
    omp_set_lock(&echelonize_lock);

    if ((ready_for_waiting_list == 0) && (global_next_row_to_reduce < D->nrows)) {
      curr_row_to_reduce  = global_next_row_to_reduce;
#if DEBUG_NEW_ELIM
      printf("%d -- inc grr\n", tid);
#endif
      global_next_row_to_reduce++;
      from_row  = 0;
    } else {
      ready_for_waiting_list  = 1;
    }

    if (ready_for_waiting_list == 1) {
      // no waiting rows in the list?
      if (get_smallest_waiting_row(&waiting_global, &wl_idx, &wl_lp) == 0) {
        // no more rows to be reduced?
        if (global_next_row_to_reduce >= D->nrows) {
          omp_unset_lock(&echelonize_lock);
          break;
        } else {
          // we do not have to check this since
          // global_last_piv == local_last_piv
          //if (local_last_piv >= D->nrows) {
          //  omp_unset_lock(&echelonize_lock);
          //  break;
          //} else {
            ready_for_waiting_list  = 0;
            omp_unset_lock(&echelonize_lock);
            continue;
          //}
        }
      }
      from_row            = wl_lp + 1;
      curr_row_to_reduce  = wl_idx;
#if DEBUG_NEW_ELIM
      printf("(%d) from row %d -- crr %d -- gllp %d\n", tid, from_row,
          curr_row_to_reduce, global_last_piv);
#endif
    }

    omp_unset_lock(&echelonize_lock);
    if (D->row[curr_row_to_reduce]->piv_val != NULL) {
#if DEBUG_NEW_ELIM
      printf("thread %d reduces %d with rows %d -- %d\n", tid,
          curr_row_to_reduce, from_row, local_last_piv);
#endif
      reduce_dense_row_task(D, curr_row_to_reduce, from_row, local_last_piv);
#if DEBUG_NEW_ELIM
      printf("thread %d done with rows %d -- %d\n", tid, from_row, local_last_piv);
      printf("row %u has lead index %u\n",curr_row_to_reduce, D->row[curr_row_to_reduce]->lead);
#endif
    }
    if (curr_row_to_reduce == local_last_row_fully_reduced + 1) {
    //if (curr_row_to_reduce == local_last_piv + 1) {
      curr_row_fully_reduced  = 1;
      ready_for_waiting_list  = 0;
    } else {
      curr_row_fully_reduced  = 0;
      nred_consecutively++;
      if (nred_consecutively > 4) {
        nred_consecutively      = 0;
        ready_for_waiting_list  = 1;
      } else {
        ready_for_waiting_list  = 0;
      }
    }

    // since we are using dense large D we do not have to store back rows
    // intermediately.
    // => we only save pivots, i.e. we normalize and reduce them w.r.t D->mod
    if (curr_row_fully_reduced == 1) { // new pivot found
      free(D->row[curr_row_to_reduce]->piv_val);
      D->row[curr_row_to_reduce]->piv_val = NULL;

      omp_set_lock(&echelonize_lock);
      global_last_row_fully_reduced++;
      if (D->row[curr_row_to_reduce]->val != NULL) {
        save_pivot(D, curr_row_to_reduce, global_last_piv+1);
      }
      if (D->row[global_last_piv+1]->piv_val == NULL)
        D->rank--;
      else
        ++global_last_piv;
#if DEBUG_NEW_ELIM
      printf("%d -- inc glp: %d\n", tid, global_last_piv);
#endif
      omp_unset_lock(&echelonize_lock);
    } else { // add row to waiting list
#if DEBUG_NEW_ELIM
      printf("%d -- pushes %d / %d\n", tid, curr_row_to_reduce, local_last_piv);
#endif
      omp_set_lock(&echelonize_lock);
      push_row_to_waiting_list(&waiting_global, curr_row_to_reduce, local_last_piv);
      omp_unset_lock(&echelonize_lock);
    }
  }
  return 0;
}

void reduce_upwards(dm_t *D, const ri_t last_row)
{
  ri_t i, j;
  re_t mult;
  
  i = last_row;
  while (i<(ri_t)(-1)) {
    j = 0;
    while (j<i) {
      mult  = MODP(D->row[j]->val[D->row[i]->lead], D->mod);
      if (mult != 0) {
        mult  = D->mod - mult;
        // also updates lead for row j
        reduce_dense_row(D, j, i, mult);
      }
      j++;
    }
    i--;
  }
  // modular computation of reduced rows
  for (i=0; i<=last_row; ++i)
    red_mod_dense_row(D, i);
}

void pre_elim_sequential(dm_t *D, const ri_t last_row)
{
  ri_t i, j, test_idx;
  re_t mult;

  copy_to_val(D, 0);
  save_pivot(D, 0, 0);

  i = 1;
  while (i<D->rank && i<=last_row) {
    j = 0;
    copy_to_val(D, i);
    while (j<i) {
      mult  = MODP(D->row[i]->val[D->row[j]->lead], D->mod);
      if (mult != 0) {
        //printf("mult in preelim for %u --> %u\n",i,mult);
        //inverse_val(&mult, D->mod);
        mult  = D->mod-mult;
        // also updates lead for row i
        reduce_dense_row(D, i, j, mult);
        // if reduced row i is zero row then swap row down and get a new
        // row from the bottom
        if (D->row[i]->val == NULL) {
          D->rank--;
          test_idx  = swap_zero_row(D, i);
          // if there is no nonzero row below we are done
          if (test_idx == i)
              return;
          // restart to reduce the new row i with j=0 in the next round of the
          // for loop
          j = 0;
        }
        if (D->row[i]->lead > i) {
          test_idx  = swap_row(D, i);
          // if there is no nonzero row below we are done
          if (test_idx == i)
              return;
          // restart to reduce the new row i with j=0 in the next round of the
          // for loop
          j = 0;
        }
      }
      j++;
    }
    //normalize_dense_row(D, i);
    save_pivot(D, i, i);
    i++;
  }
#if DEBUG_NEW_ELIM
  for (int ii=0; ii<=last_row; ++ii) {
    printf("ROW %d\n",ii);
    for (int jj=0; jj<D->ncols; ++jj)
      printf("%lu  ", D->row[ii]->val[jj]);
    printf("\n");
  }
#endif
  //reduce_upwards(D, last_row);
#if DEBUG_NEW_ELIM
  for (int ii=0; ii<=last_row; ++ii) {
    printf("ROW %d\n",ii);
    for (int jj=0; jj<D->ncols; ++jj)
      printf("%lu  ", D->row[ii]->val[jj]);
    printf("\n");
  }
#endif
}

ri_t elim_fl_dense_D(dm_t *D, int nthrds) {
  // setting global variables for open mp locks used later on
  global_next_row_to_reduce = nthrds * 2;
  global_last_piv           = global_next_row_to_reduce - 1;
  global_last_row_fully_reduced = global_last_piv;
  global_first_zero_row     = D->nrows;   

  /*  global waiting list */
  waiting_global.list = (wle_t *)malloc(D->nrows * sizeof(wle_t));
  waiting_global.sidx = 0;
  waiting_global.slp  = 0;
  waiting_global.sz   = 0;

#if COUNT_REDS
  nreductions = 0;
#endif

  int rc;

  // sort D w.r.t. first nonzero entry per row
  sort_dense_matrix_by_pivots(D);
#if DEBUG_NEW_ELIM
  printf("GLP %u\n",global_last_piv);
  for (int ii=0; ii<D->nrows; ++ii) {
    printf("%u - ", ii);
    if (D->row[ii] == NULL)
      printf("NULL\n");
    else
      printf("%u\n",D->row[ii]->lead);
  }
#endif
  // do sequential prereduction of first global_last_piv bunch of rows
  pre_elim_sequential(D, global_last_piv);
#if DEBUG_NEW_ELIM
  printf("GLP %u\n",global_last_piv);
  for (int ii=0; ii<D->nrows; ++ii) {
    printf("%u - ", ii);
    if (D->row[ii] == NULL)
      printf("NULL\n");
    else
      printf("%u\n",D->row[ii]->lead);
  }
#endif

  // if rank is smaller than the row until which we have been sequentially pre
  // eliminating then all other rows of D are already zero
  if (D->rank < global_last_piv)
    return D->rank;

  // do the parallel computation
  omp_init_lock(&echelonize_lock);
  int tid;

#pragma omp parallel shared(D, waiting_global, global_next_row_to_reduce, global_last_piv) num_threads(nthrds)
  {
#pragma omp for nowait
    for (tid=0; tid<nthrds; ++tid) {
      rc  = elim_fl_dense_D_tasks(D);
    }
  }
  omp_destroy_lock(&echelonize_lock);

#if COUNT_REDS
  printf("\n#REDUCTION STEPS: %lu\n",nreductions);
#endif
  return D->rank;
}

ri_t elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, mod_t modulus, int nthrds) {

  ri_t i;
  /* ci_t rc; */

  /*  row indices for subdividing echelonization parts in D_red */
  global_next_row_to_reduce     = nthrds * 2;
  global_last_piv               = global_next_row_to_reduce - 1;

  /*  meta data for the computation of the rank of D_red at the end */
  int head_line_1       = -1;
  int head_line_2       = -1;
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;
  const ci_t coldim     = D->ncols;
  ri_t wl_dim; /*  waiting list dimension */
  re_t h_a1;
#if COUNT_REDS
  nreductions = 0;
#endif

  wl_dim  = (D->nrows/__GB_NROWS_MULTILINE);
  if (D->nrows%__GB_NROWS_MULTILINE)
    wl_dim++;

  /*  global waiting list */
  waiting_global.list = (wle_t *)malloc(wl_dim * sizeof(wle_t));
  waiting_global.sidx = 0;
  waiting_global.slp  = 0;
  waiting_global.sz   = 0;

  /*  copy D to D_red and delete D */
  D_red = copy_block_matrix_to_multiline_matrix(&D, D_red, 1, nthrds);
#if DDEBUG_DD
  printf("BEFORE\n");
  const uint32_t rlD  = (uint32_t) ceil((float)D_red->nrows / __GB_NROWS_MULTILINE);
  int ii,jj,kk,ll;
  for (ii=0; ii<rlD; ++ii) {
    printf("%d .. \n",ii);
    /* printf("size %d\n", D_red->ml[ii].sz); */
    if (D_red->ml[ii].sz>0) {
      for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
        if (D_red->ml[ii].idx != NULL)
          printf("%d -- ", D_red->ml[ii].idx[ll]);
        else
          printf("%d -- ", ll);
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
      printf("\n");
    } else {
      printf("ml %d is zero! %p\n", ii, D_red->ml[ii]);
    }
  }
#endif
  echelonize_rows_sequential(D_red, 0, global_last_piv, modulus);
#if DDEBUG_DD
  printf("AFTER\n");
  for (int ii=0; ii<rlD; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", D_red->ml[ii].sz);
    if (D_red->ml[ii].sz>0) {
      for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
        /*
           if (D_red->ml[ii].idx != NULL)
           printf("%d -- ", D_red->ml[ii].idx[ll]);
           else
           printf("%d -- ", ll);
           */
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
    } else {
    }
    printf("\n");
  }
#endif

  const ri_t ml_nrows_D_red = (D_red->nrows % __GB_NROWS_MULTILINE == 0) ?
    D_red->nrows / __GB_NROWS_MULTILINE :
    D_red->nrows / __GB_NROWS_MULTILINE + 1;

  /*  define lock */
  /* omp_lock_t echelonize_lock; */
  omp_init_lock(&echelonize_lock);

  /*  if there are rows left do elimination in parallel */
  if (ml_nrows_D_red >= global_next_row_to_reduce) {
    int tid ;
    /*  TODO: parallel elimination with OpenMP */
#pragma omp parallel shared(D_red, waiting_global, global_next_row_to_reduce, global_last_piv) num_threads(nthrds)
    {
#pragma omp for nowait
      for (tid=0; tid<nthrds; ++tid) {
        /*rc  =*/ echelonize_rows_task(D_red, ml_nrows_D_red,
            /* global_next_row_to_reduce, global_last_piv, */
            /* &waiting_global, */
            modulus
            /* , echelonize_lock */
            );
      }
    }
  }
  omp_destroy_lock(&echelonize_lock);

  ri_t rank = 0;

  for (i=0; i<ml_nrows_D_red; ++i) {
    if (D_red->ml[i].sz == 0) {
      continue;
    }


    head_line_1 = get_head_multiline_hybrid(&(D_red->ml[i]), 0,
        &h_a1, &head_line_1_idx, coldim);
    head_line_2 = get_head_multiline_hybrid(&(D_red->ml[i]), 1,
        &h_a1, &head_line_2_idx, coldim);

    if (head_line_1 != -1)
      rank++;
    if (head_line_2 != -1)
      rank++;
  }
#if DDEBUG_DONE
  const uint32_t rlD  = (uint32_t) ceil((float)D_red->nrows / __GB_NROWS_MULTILINE);
  int ii,jj,kk,ll;
  printf("AFTER ALL\n");
  for (ii=0; ii<rlD; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", D_red->ml[ii].sz * 2);
    if (D_red->ml[ii].sz>0) {
      for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
        /*
           if (D_red->ml[ii].idx != NULL)
           printf("%d -- ", D_red->ml[ii].idx[ll]);
           else
           printf("%d -- ", ll);
           */
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
      printf("\n");
    }
  }
#endif
#if COUNT_REDS
  printf("#REDUCTION STEPS: %lu\n",nreductions);
#endif

  return rank;
}

ri_t echelonize_rows_sequential(sm_fl_ml_t *A, const ri_t from, const ri_t to,
    const mod_t modulus) {

  if (A->nrows == 0)
    return 0;

  int ret;
  ri_t npiv_real  = 0;
  ri_t N          = A->nrows / __GB_NROWS_MULTILINE +
    A->nrows % __GB_NROWS_MULTILINE;
  const ci_t coldim = A->ncols;

  ml_t *ml_row;
  re_l_t *dense_array_1, *dense_array_2;
  ret = posix_memalign((void **)&dense_array_1, ALIGNT, coldim * sizeof(re_l_t));
  ret = posix_memalign((void **)&dense_array_2, ALIGNT, coldim * sizeof(re_l_t));

  ri_t i;
  ci_t j, k;
  int head_line_1       = -1;
  int head_line_2       = -1;
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;

  re_t intermediate_val;
  re_m_t tmp_val;

  normalize_multiline(&A->ml[from], coldim, modulus);

  ri_t min_loop = to > N-1 ? N-1 : to;
  for (i=from; i<=min_loop; ++i) {
    if (A->ml[i].val != NULL) {
      memset(dense_array_1, 0, coldim * sizeof(re_l_t));
      memset(dense_array_2, 0, coldim * sizeof(re_l_t));
      copy_multiline_to_dense_array(A->ml[i], dense_array_1, dense_array_2, coldim);
#if DDEBUG_D
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("3-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("4-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
#endif
      re_t h_a1 = 1, h_a2 = 1;

      re_t v1_col1 = 0, v2_col1 = 0, v1_col2 = 0, v2_col2 = 0;
      /* ri_t tmp  = 0; */

      for (j=from; j<i; ++j) {
        ml_row  = &(A->ml[j]);
        if (ml_row->sz == 0)
          continue;

        head_line_1 = get_head_multiline_hybrid(ml_row, 0,
            &h_a1, &head_line_1_idx, coldim);
        head_line_2 = get_head_multiline_hybrid(ml_row, 1,
            &h_a2, &head_line_2_idx, coldim);

        if (head_line_1 != -1) {
          v1_col1 = MODP(dense_array_1[head_line_1] , modulus);
          v2_col1 = MODP(dense_array_2[head_line_1] , modulus);

          if (v1_col1 != 0)
            v1_col1 = modulus - v1_col1;
          if (v2_col1 != 0)
            v2_col1 = modulus - v2_col1;
        } else {
          v1_col1 = 0;
          v2_col1 = 0;
        }
#if DDEBUG_D_ONE
        printf("v11 %d\n",v1_col1);
        printf("v21 %d\n",v2_col1);
#endif
        if (head_line_2 != -1) {
          v1_col2 = MODP(dense_array_1[head_line_2] , modulus);
          v2_col2 = MODP(dense_array_2[head_line_2] , modulus);

          intermediate_val  = ml_row->val[2*head_line_2_idx];
          tmp_val = v1_col2 + (re_m_t)v1_col1 * intermediate_val;
          v1_col2 = MODP(tmp_val , modulus);
          tmp_val = v2_col2 + (re_m_t)v2_col1 * intermediate_val;
          v2_col2 = MODP(tmp_val , modulus);

          if (v1_col2 != 0)
            v1_col2 = modulus - v1_col2;
          if (v2_col2 != 0)
            v2_col2 = modulus - v2_col2;
        } else {
          v1_col2 = 0;
          v2_col2 = 0;
        }
#if DDEBUG_D_ONE
        printf("v12 %d\n",v1_col2);
        printf("v22 %d\n",v2_col2);
#endif
        if (ml_row->sz < coldim) {
          sparse_scal_mul_sub_2_rows_vect_array_multiline(
              v1_col1, v2_col1,
              v1_col2, v2_col2,
              *ml_row,
              dense_array_1,
              dense_array_2);
        } else {
          dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
              v1_col1, v2_col1,
              v1_col2, v2_col2,
              *ml_row,
              dense_array_1,
              dense_array_2,
              head_line_1,
              head_line_2);
        }
      }
#if DDEBUG_D
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("5-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("6-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
#endif
      /*  normalize dense arrays */
      head_line_1 = normalize_dense_array(dense_array_1, coldim, modulus);
      head_line_2 = normalize_dense_array(dense_array_2, coldim, modulus);
#if DDEBUG_D
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("7-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("8-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
#endif
      /*  reduce by same multiline */
      if (head_line_1 >= head_line_2 && head_line_1 != -1 && head_line_2 != -1) {
        dense_array_2[head_line_1] =MODP(dense_array_2[head_line_1], modulus);
        if (dense_array_2[head_line_1] != 0) {
          register re_m_t h = modulus - dense_array_2[head_line_1];
          register re_m_t v__;
#ifdef GBLA_USE_AVX_XXX
          for (k=head_line_1; k<coldim; ++k) {
            v__               =   CAST(dense_array_1[k]);
            dense_array_2[k]  +=  h * v__;
          }
#else
          for (k=head_line_1; k<coldim; ++k) {
            v__               =   CAST(dense_array_1[k]);
            dense_array_2[k]  +=  h * v__;
          }
#endif
        }
      }

      head_line_2 = get_head_dense_array(dense_array_2, &h_a2, coldim, modulus);

      /*  make A->ml[i] dense, i.e. free memory for idx and enlarge memory for val */
      /*  initialize all values in val with 0, set size of A->ml[i] to coldim */
      if (A->ml[i].dense == 0) {
        free (A->ml[i].idx);
        A->ml[i].idx  = NULL;
        A->ml[i].val  = (re_t *)realloc(A->ml[i].val, 2 * coldim * sizeof(re_t));
        A->ml[i].sz   = coldim;
      }
      memset(A->ml[i].val, 0, 2 * coldim * sizeof(re_t));

#if DDEBUG_Dd
      printf("DENSE1\n");
      for (int uu=0; uu<coldim; ++uu)
        printf("%lu :: ",dense_array_1[uu]);
      printf("DENSE2\n");
      for (int uu=0; uu<coldim; ++uu)
        printf("%lu :: ",dense_array_2[uu]);
#endif

      /*  save the line with the smallest column entry first */
      if (head_line_1 == -1) {
        copy_dense_array_to_zero_dense_multiline(
            dense_array_2, head_line_2, &A->ml[i], coldim, modulus);
      } else {
        if (head_line_2 == -1) {
          copy_dense_array_to_zero_dense_multiline(
              dense_array_1, head_line_1, &A->ml[i], coldim, modulus);
        } else { /*  both are not empty */
          if (head_line_1 > head_line_2) {
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_2, dense_array_1, head_line_2, &A->ml[i], coldim, modulus);
          } else {
            /* copy_dense_arrays_to_dense_multiline( */
            /*     dense_array_1, dense_array_2, ml, coldim, modulus); */
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_1, dense_array_2, head_line_1, &A->ml[i], coldim, modulus);
          }
        }
      }
#if DDEBUG_D
      printf("BEFORE NORMALIZE\n");
      const uint32_t rlD  = (uint32_t) ceil((float)A->nrows / __GB_NROWS_MULTILINE);
      int ii,jj,kk,ll;
      for (ii=0; ii<rlD; ++ii) {
        printf("%d .. \n",ii);
        printf("size %d\n", A->ml[ii].sz * 2);
        if (A->ml[ii].sz>0) {
          for (ll=0; ll<A->ml[ii].sz; ++ll) {
            /*
               if (D_red->ml[ii].idx != NULL)
               printf("%d -- ", D_red->ml[ii].idx[ll]);
               else
               printf("%d -- ", ll);
               */
            printf("%d %d ", A->ml[ii].val[2*ll], A->ml[ii].val[2*ll+1]);
          }
          printf("\n");
        }
      }
#endif
      /*  normalize multiline */
      normalize_multiline(&A->ml[i], coldim, modulus);
      if (head_line_1 != -1)
        npiv_real++;
      if (head_line_2 != -1)
        npiv_real++;
    }
  }
  free(dense_array_1);
  dense_array_1 = NULL;
  free(dense_array_2);
  dense_array_2 = NULL;

  return npiv_real;
}

int echelonize_rows_task(sm_fl_ml_t *A, const ri_t N,
    /* ri_t global_next_row_to_reduce, ri_t global_last_piv, */
    /* wl_t *waiting_global, */
    const mod_t modulus
    /* , omp_lock_t echelonize_lock */
    ) {

  const ci_t coldim = A->ncols;
  ri_t curr_row_to_reduce;
  ri_t local_last_piv;

  int ready_for_waiting_list = 0;
  int curr_row_fully_reduced = 0;
  int ret;

  ri_t from_row;
  int nreduced_consecutively  = 0;
  /*  local waiting list entries */
  ri_t wl_idx  = 0;
  ri_t wl_lp   = 0;

  re_l_t *dense_array_1, *dense_array_2;
  ret = posix_memalign((void **)&dense_array_1, ALIGNT, coldim * sizeof(re_l_t));
  ret = posix_memalign((void **)&dense_array_2, ALIGNT, coldim * sizeof(re_l_t));
#if DEBUG_ECHELONIZE
  int tid = omp_get_thread_num();
#endif

  /*  do the computations */
  while (1) {
    local_last_piv = global_last_piv;
#if DEBUG_ECHELONIZE
    printf("seeting llp\n");
    printf("llp %d\n",local_last_piv);
    printf("grr %d\n",global_next_row_to_reduce);
#endif
    if (global_last_piv >= N) {
      omp_unset_lock(&echelonize_lock);
      break;
    }
    omp_set_lock(&echelonize_lock);
    if ((ready_for_waiting_list == 0) && (global_next_row_to_reduce < N)) {
      curr_row_to_reduce = global_next_row_to_reduce;
#if DEBUG_ECHELONIZE
      printf("%d -- inc grr\n",tid);
#endif
      global_next_row_to_reduce++;
      from_row = 0;
    } else {
      ready_for_waiting_list = 1;
    }
    if (ready_for_waiting_list == 1) {
      if (get_smallest_waiting_row(&waiting_global, &wl_idx, &wl_lp) == 0) { /*  no waiting rows */
        if (global_next_row_to_reduce >= N) { /*  no more rows to reduce */
          omp_unset_lock(&echelonize_lock);
          break;
        } else {
          if (local_last_piv >= N) { /*  we are also done */
            omp_unset_lock(&echelonize_lock);
            break;
          } else {
            ready_for_waiting_list  = 0;
            omp_unset_lock(&echelonize_lock);
            continue;
          }
        }
      }
      from_row            = wl_lp + 1;
      curr_row_to_reduce  = wl_idx;
#if DEBUG_ECHELONIZE
      printf("(%d) from row %d -- crr %d -- gllp %d\n", tid,from_row, curr_row_to_reduce, global_last_piv);
#endif
    }

    omp_unset_lock(&echelonize_lock);
    if (A->ml[curr_row_to_reduce].val != NULL) {
      /*  set zero */
      memset(dense_array_1, 0, coldim * sizeof(re_l_t));
      memset(dense_array_2, 0, coldim * sizeof(re_l_t));

      copy_multiline_to_dense_array(A->ml[curr_row_to_reduce],
          dense_array_1, dense_array_2, coldim);
      /*  echelonize one row */
#if DEBUG_ECHELONIZE
      printf("thread %d reduces %d with rows %d -- %d\n",tid, curr_row_to_reduce, from_row, local_last_piv);
#endif
      echelonize_one_row(A, dense_array_1, dense_array_2, from_row,
          local_last_piv, modulus);
#if DEBUG_ECHELONIZE
      printf("thread %d done with rows %d -- %d\n",tid, from_row, local_last_piv);
#endif
    }
    if (curr_row_to_reduce == local_last_piv + 1) {
      curr_row_fully_reduced  = 1;
      ready_for_waiting_list  = 0;
    } else {
      curr_row_fully_reduced  = 0;
      ready_for_waiting_list  = 0;
      ++nreduced_consecutively;
      if (nreduced_consecutively > 4) {
        nreduced_consecutively  = 0;
        ready_for_waiting_list  = 1;
      }
    }

    /*  save back to multiline */
    /*
       printf("1 %p -- 2 %p\n",dense_array_1,dense_array_2);
       printf("BEFORE A.val %p\n",A->ml[curr_row_to_reduce].val);
       */
    if (A->ml[curr_row_to_reduce].val != NULL) {
      save_back_and_reduce(&(A->ml[curr_row_to_reduce]), dense_array_1,
          dense_array_2, coldim, modulus, curr_row_fully_reduced, curr_row_to_reduce);
    }
    /*
       printf("AFTER A.val %p\n",A->ml[curr_row_to_reduce].val);
       printf("1 %p -- 2 %p\n",dense_array_1,dense_array_2);
       printf("3-mlidx %p -- mlval %p\n",A->ml[curr_row_to_reduce].idx,A->ml[curr_row_to_reduce].val);
       printf("MULTILINE AFTER SAVE\n");
       for (int kk=0; kk<A->ml[curr_row_to_reduce].sz; ++kk) {
       printf("%d : %u -- %u ",kk,A->ml[curr_row_to_reduce].val[2*kk],A->ml[curr_row_to_reduce].val[2*kk+1]);
       }
       printf("\n");
       */

    if (curr_row_fully_reduced == 1) {
      omp_set_lock(&echelonize_lock);
      ++global_last_piv;
#if DEBUG_ECHELONIZE
      printf("%d -- inc glp: %d\n", tid, global_last_piv);
#endif
      omp_unset_lock(&echelonize_lock);
    } else {
      omp_set_lock(&echelonize_lock);
#if DEBUG_ECHELONIZE
      printf("%d -- pushes %d / %d\n",tid, curr_row_to_reduce, local_last_piv);
#endif
      push_row_to_waiting_list(&waiting_global, curr_row_to_reduce, local_last_piv);
      omp_unset_lock(&echelonize_lock);
    }
  }
  /*  free memory */
  free(dense_array_1);
  dense_array_1 = NULL;
  free(dense_array_2);
  dense_array_2 = NULL;

  return 0;
}

void echelonize_one_row(sm_fl_ml_t *A,
    re_l_t *dense_array_1, re_l_t *dense_array_2,
    const ri_t first_piv, const ri_t last_piv,
    const mod_t modulus) {

  ci_t j;
  int head_line_1       = -1;
  int head_line_2       = -1;
  int nzc               = 0; /*  tracks if there is a nonzero coeff at all */
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;

  re_t intermediate_val;
  re_m_t tmp_val;
  const ci_t coldim = A->ncols;

  ml_t *ml_row;

  re_t h_a1 = 1, h_a2 = 1;

  re_t v1_col1 = 0, v2_col1 = 0, v1_col2 = 0, v2_col2 = 0;
  /* ri_t tmp  = 0; */

#if DEBUG_ECHELONIZE
  printf("fp %d -- lp %d\n",first_piv, last_piv);
#endif
  for (j=first_piv; j<last_piv+1; ++j) {
    nzc = 0;
    ml_row  = &(A->ml[j]);
#if DEBUG_ECHELONIZE
    printf("j %d\n",j);
#endif
    if (ml_row->val == NULL || ml_row->sz == 0)
      continue;
#if DDEBUG_D
    if (ml_row != NULL) {
      printf("MULTILINE\n");
      for (int kk=0; kk<ml_row->sz; ++kk) {
        printf("%d : %u -- %u ",kk,ml_row->val[2*kk],ml_row->val[2*kk+1]);
      }
      printf("\n");
    }
#endif


    head_line_1 = get_head_multiline_hybrid(ml_row, 0,
        &h_a1, &head_line_1_idx, coldim);
    head_line_2 = get_head_multiline_hybrid(ml_row, 1,
        &h_a2, &head_line_2_idx, coldim);

#if DEBUG_ECHELONIZE
    printf("hl1 %d || h2 %d\n",head_line_1,head_line_2);
#endif

    if (head_line_1 != -1) {
#if DDEBUG_D_ONE
      printf("d11 %lu (%p)\n",dense_array_1[head_line_1],dense_array_1[head_line_1]);
      printf("d21 %lu\n",dense_array_2[head_line_1]);
      printf("d11 %lu\n",MODP(dense_array_1[head_line_1], modulus));
      printf("d21 %lu\n",MODP(dense_array_2[head_line_1], modulus));
#endif
      v1_col1 = MODP(dense_array_1[head_line_1], modulus);
      v2_col1 = MODP(dense_array_2[head_line_1], modulus);
#if DDEBUG_D_ONE
      printf("v11 %d\n",v1_col1);
      printf("v21 %d\n",v2_col1);
#endif

      if (v1_col1 != 0) {
        v1_col1 = modulus - v1_col1;
        nzc = 1;
      }
      if (v2_col1 != 0) {
        v2_col1 = modulus - v2_col1;
        nzc = 1;
      }
    } else {
      v1_col1 = 0;
      v2_col1 = 0;
    }
#if DDEBUG_D_ONE
    printf("v11 %d\n",v1_col1);
    printf("v21 %d\n",v2_col1);
#endif
    /* printf("hl2 %d => %d ?\n",head_line_2,head_line_2 != -1); */
    if (head_line_2 != -1) {
      v1_col2 = MODP(dense_array_1[head_line_2], modulus);
      v2_col2 = MODP(dense_array_2[head_line_2], modulus);
#if DDEBUG_D_ONE
      printf("v12 %d\n",v1_col2);
      printf("v22 %d\n",v2_col2);
#endif

      intermediate_val  = ml_row->val[2*head_line_2_idx];
      tmp_val = v1_col2 + (re_m_t)v1_col1 * intermediate_val;
      v1_col2 = MODP(tmp_val, modulus);
      tmp_val = v2_col2 + (re_m_t)v2_col1 * intermediate_val;
      v2_col2 = MODP(tmp_val, modulus);
#if DDEBUG_D_ONE
      printf("v12 %d\n",v1_col2);
      printf("v22 %d\n",v2_col2);
#endif

      if (v1_col2 != 0) {
        v1_col2 = modulus - v1_col2;
        nzc = 1;
      }
      if (v2_col2 != 0) {
        v2_col2 = modulus - v2_col2;
        nzc = 1;
      }
    } else {
      v1_col2 = 0;
      v2_col2 = 0;
    }
#if DDEBUG_D_ONE
    printf("v12 %d\n",v1_col2);
    printf("v22 %d\n",v2_col2);
#endif
    if (nzc == 1) {
      if (ml_row->sz < coldim) {
        sparse_scal_mul_sub_2_rows_vect_array_multiline(
            v1_col1, v2_col1,
            v1_col2, v2_col2,
            *ml_row,
            dense_array_1,
            dense_array_2);
      } else {
        dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
            v1_col1, v2_col1,
            v1_col2, v2_col2,
            *ml_row,
            dense_array_1,
            dense_array_2,
            head_line_1,
            head_line_2);
      }
    }
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("15-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("16-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
  }
}

void save_back_and_reduce(ml_t *ml, re_l_t *dense_array_1,
    re_l_t *dense_array_2, const ci_t coldim, const mod_t modulus,
    const int reduce, const ri_t curr_row_to_reduce) {


  if (reduce == 1) {
    int head_line_1 = -1;
    int head_line_2 = -1;
    re_t /* h_a1 = 1, */ h_a2 = 1;
    ci_t k;
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("17-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("18-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    /*  normalize dense arrays */
    head_line_1 = normalize_dense_array(dense_array_1, coldim, modulus);
    head_line_2 = normalize_dense_array(dense_array_2, coldim, modulus);
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("19-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("20-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    /*  reduce by same multiline */
    if (head_line_1 >= head_line_2 && head_line_1 != -1 && head_line_2 != -1) {
      dense_array_2[head_line_1] = MODP(dense_array_2[head_line_1], modulus);
      if (dense_array_2[head_line_1] != 0) {
        register re_m_t h = modulus - dense_array_2[head_line_1];
        register re_m_t v__;
#ifdef GBLA_USE_AVX_XXX
        for (k=head_line_1; k<coldim; ++k) {
          v__               =   CAST(dense_array_1[k]);
          dense_array_2[k]  +=  h * v__;
        }
#else
        for (k=head_line_1; k<coldim; ++k) {
          v__               =   CAST(dense_array_1[k]);
          dense_array_2[k]  +=  h * v__;
        }
#endif
      }
    }
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("21-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("22-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    head_line_2 = get_head_dense_array(dense_array_2, &h_a2, coldim, modulus);

#if DEBUG_ECHELONIZE
    printf("headl1 %d | headl2 %d\n",head_line_1,head_line_2);
    printf("mlidx %p -- mlval %p : sz %d\n",ml->idx,ml->val, ml->sz);
#endif
    /*  possibly we have a zero row */
    if (head_line_1 == -1 && head_line_2 == -1) {
      free(ml->val);
      ml->val  = NULL;
      ml->sz   = 0;
    } else {
      memset(ml->val, 0, 2 * coldim * sizeof(re_t));

      if (head_line_1 == -1) {
        copy_dense_array_to_zero_dense_multiline(
            dense_array_2, head_line_2, ml, coldim, modulus);
      } else {
        if (head_line_2 == -1) {
          copy_dense_array_to_zero_dense_multiline(
              dense_array_1, head_line_1, ml, coldim, modulus);
        } else { /*  both are not empty */
          if (head_line_1 > head_line_2) {
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_2, dense_array_1, head_line_2, ml, coldim, modulus);
          } else {
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_1, dense_array_2, head_line_1, ml, coldim, modulus);
          }
        }
      }
      /*  normalize multiline */
      normalize_multiline(ml, coldim, modulus);
    }
  } else { /*  do not reduce */
    /*  make A->ml[i] dense, i.e. free memory for idx and enlarge memory for val */
    /*  initialize all values in val with 0, set size of A->ml[i] to coldim */
    /* free(ml->idx); */
    /* ml->idx  = NULL; */
    /* ml->sz   = coldim; */
    memset(ml->val, 0, 2 * coldim * sizeof(re_t));
    /* copy_dense_arrays_to_dense_multiline( */
    /*  dense_array_1, dense_array_2, ml, coldim, modulus); */
    copy_dense_arrays_to_zero_dense_multiline(
        dense_array_1, dense_array_2, 0, ml, coldim, modulus);
  }
}


/******************************************************************************
 * MULTILINE VARIANT OF ELIMINATION PROCESS
 *****************************************************************************/

int elim_fl_C_ml(sm_fl_ml_t *C, sm_fl_ml_t *A, mod_t modulus, int nthrds) {
  ri_t i;
  /* ri_t |+ j, k, +| rc; */

  const ri_t rlC  = (ri_t) ceil((float) C->nrows / __GB_NROWS_MULTILINE);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    /*  each task takes one block column of B */
    for (i=0; i<rlC; ++i) {
#pragma omp task
      {
        /*rc  =*/ elim_fl_C_ml_task(C, A, i, modulus);
      }
    }
#pragma omp taskwait
  }

  return 0;
}

// choose the variant to be used in the following (only one of them should be 1,
// all others have to be 0
#define ONE 0
#define TWO 0
#define THREE 0
#define FOUR 0
#define FIVE 0
#define SIX 0
#define EIGHT 0
#define TEN 1

int elim_fl_C_sparse_dense_keep_A(sm_fl_t *C, sm_fl_t **A_in, const mod_t modulus,
    const int nthrds)
{
  sm_fl_t *A= *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ri_t rlC  = C->nrows;
  
  // TRY 10 LINES AT ONCE
  if (rlC > 5000) {
    int mod10 = 0, mod6 = 0, mod4 = 0, mod2 = 0;
    if (rlC > 9)
      mod2  = mod4  = mod6  = mod10  = rlC - rlC % 10;
    if (rlC-mod10 > 5)
      mod2  = mod4  = mod6  = rlC - ((mod10-rlC) % 6);
    if (rlC-mod6 > 3)
      mod2  = mod4  = rlC - ((mod6-rlC) % 4);
    if (rlC-mod4 > 1)
      mod2  = rlC - ((mod4-rlC) % 2);
#pragma omp parallel num_threads(nthrds) private(i)
    {
#pragma omp single
      {
        //#pragma omp task untied
        if (rlC > 9) {
          for (i=0; i<rlC-9; i=i+10) {
#pragma omp task untied
            rc  = elim_fl_C_sparse_dense_keep_A_tasks_ten(C, A, i, i+1, i+2, i+3, i+4, i+5, i+6, i+7, i+8, i+9, modulus);
          }
        }
        if (rlC-mod10 > 5) {
#pragma omp task untied
          rc  = elim_fl_C_sparse_dense_keep_A_tasks_six(C, A, mod10, mod10+1, mod10+2, mod10+3, mod10+4, mod10+5, modulus);
        }
        if (rlC-mod6 > 3) {
#pragma omp task untied
          rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, mod6, mod6+1, mod6+2, mod6+3, modulus);
        }
        if (rlC-mod4 > 1) {
#pragma omp task untied
          rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, mod4, mod4+1, modulus);
        }
        if (rlC-mod2 > 0) {
#pragma omp task untied
          rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
        }
      }
#pragma omp taskwait
    }
  } else {
    // TRY 6 LINES AT ONCE
    if (rlC > 3000) {
      int mod6 = 0, mod4 = 0, mod2 = 0;
      if (rlC > 5)
        mod2  = mod4  = mod6  = rlC - -rlC % 6;
      if (rlC-mod6 > 3)
        mod2  = mod4  = rlC - ((mod6-rlC) % 4);
      if (rlC-mod4 > 1)
        mod2  = rlC - ((mod4-rlC) % 2);
#pragma omp parallel num_threads(nthrds) private(i)
      {
#pragma omp single
        {
          //#pragma omp task untied
          if (rlC > 5) {
            for (i=0; i<rlC-5; i=i+6) {
#pragma omp task untied
              rc  = elim_fl_C_sparse_dense_keep_A_tasks_six(C, A, i, i+1, i+2, i+3, i+4, i+5, modulus);
            }
          }
          if (rlC-mod6 > 3) {
#pragma omp task untied
            rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, mod6, mod6+1, mod6+2, mod6+3, modulus);
          }
          if (rlC-mod4 > 1) {
#pragma omp task untied
            rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, mod4, mod4+1, modulus);
          }
          if (rlC-mod2 > 0) {
#pragma omp task untied
            rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
          }
        }
#pragma omp taskwait
      }
    } else {
    // TRY 4 LINES AT ONCE
      int mod4 = 0, mod2 = 0;
      if (rlC > 3)
        mod2  = mod4  = rlC - rlC % 4;
      if (rlC-mod4 > 1)
        mod2  = rlC - ((mod4-rlC) % 2);
#pragma omp parallel num_threads(nthrds) private(i)
      {
#pragma omp single
        {
          //#pragma omp task untied
          if (rlC > 3) {
            for (i=0; i<rlC-3; i=i+4) {
#pragma omp task untied
              rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, i, i+1, i+2, i+3, modulus);
            }
          }
          if (rlC-mod4 > 1) {
#pragma omp task untied
            rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, mod4, mod4+1, modulus);
          }
          if (rlC-mod2 > 0) {
#pragma omp task untied
            rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
          }
        }
#pragma omp taskwait
      }
    }
  }
/*
#if ONE
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; ++i) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, i, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if TWO
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; i=i+2) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, i, i+1, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if THREE
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; i=i+3) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, i, i+1, i+2, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if FOUR
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; i=i+4) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, i, i+1, i+2, i+3, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if FIVE
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; i=i+5) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_five(C, A, i, i+1, i+2, i+3, i+4, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if SIX
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; i=i+6) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_six(C, A, i, i+1, i+2, i+3, i+4, i+5, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if EIGHT
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; i=i+8) {
#pragma omp task
      {
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_eight(C, A, i, i+1, i+2, i+3, i+4, i+5, i+6, i+7, modulus);
      }
    }
#pragma omp taskwait
  }
#endif
#if TWO
    int mod2 = 0;
    if (rlC > 1)
      mod2  = rlC - rlC % 2;
#pragma omp parallel num_threads(nthrds) private(i)
{
#pragma omp single
  {
//#pragma omp task untied
    if (rlC > 1) {
      for (i=0; i<rlC-1; i=i+2) {
        #pragma omp task untied
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, i, i+1, modulus);
      }
    }
    if (rlC-mod2 > 0) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
    }
  }
#pragma omp taskwait
}
#endif
#if FOUR
    int mod4 = 0, mod2 = 0;
    if (rlC > 3)
      mod2  = mod4  = rlC - rlC % 4;
    if (rlC-mod4 > 1)
      mod2  = rlC - ((mod4-rlC) % 2);
#pragma omp parallel num_threads(nthrds) private(i)
{
#pragma omp single
  {
//#pragma omp task untied
    if (rlC > 3) {
      for (i=0; i<rlC-3; i=i+4) {
        #pragma omp task untied
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, i, i+1, i+2, i+3, modulus);
      }
    }
    if (rlC-mod4 > 1) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, mod4, mod4+1, modulus);
    }
    if (rlC-mod2 > 0) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
    }
  }
#pragma omp taskwait
}
#endif
#if SIX
    int mod6 = 0, mod4 = 0, mod2 = 0;
    if (rlC > 5)
      mod2  = mod4  = mod6  = rlC - -rlC % 6;
    if (rlC-mod6 > 3)
      mod2  = mod4  = rlC - ((mod6-rlC) % 4);
    if (rlC-mod4 > 1)
      mod2  = rlC - ((mod4-rlC) % 2);
#pragma omp parallel num_threads(nthrds) private(i)
{
#pragma omp single
  {
//#pragma omp task untied
    if (rlC > 5) {
      for (i=0; i<rlC-5; i=i+6) {
        #pragma omp task untied
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_six(C, A, i, i+1, i+2, i+3, i+4, i+5, modulus);
      }
    }
    if (rlC-mod6 > 3) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, mod6, mod6+1, mod6+2, mod6+3, modulus);
    }
    if (rlC-mod4 > 1) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, mod4, mod4+1, modulus);
    }
    if (rlC-mod2 > 0) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
    }
  }
#pragma omp taskwait
}
#endif
#if TEN
    int mod10 = 0, mod6 = 0, mod4 = 0, mod2 = 0;
    if (rlC > 9)
      mod2  = mod4  = mod6  = mod10  = rlC - rlC % 10;
    if (rlC-mod10 > 5)
      mod2  = mod4  = mod6  = rlC - ((mod10-rlC) % 6);
    if (rlC-mod6 > 3)
      mod2  = mod4  = rlC - ((mod6-rlC) % 4);
    if (rlC-mod4 > 1)
      mod2  = rlC - ((mod4-rlC) % 2);
#pragma omp parallel num_threads(nthrds) private(i)
{
#pragma omp single
  {
//#pragma omp task untied
    if (rlC > 9) {
      for (i=0; i<rlC-9; i=i+10) {
        #pragma omp task untied
        rc  = elim_fl_C_sparse_dense_keep_A_tasks_ten(C, A, i, i+1, i+2, i+3, i+4, i+5, i+6, i+7, i+8, i+9, modulus);
      }
    }
    if (rlC-mod10 > 5) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_six(C, A, mod10, mod10+1, mod10+2, mod10+3, mod10+4, mod10+5, modulus);
    }
    if (rlC-mod6 > 3) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, mod6, mod6+1, mod6+2, mod6+3, modulus);
    }
    if (rlC-mod4 > 1) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, mod4, mod4+1, modulus);
    }
    if (rlC-mod2 > 0) {
      #pragma omp task untied
      rc  = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, mod2, modulus);
    }
  }
#pragma omp taskwait
}
#endif
*/

  // free A
  free_sparse_matrix(&A, nthrds);
  *A_in = A;
  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_double(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const mod_t modulus)
{
  /*
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx1, modulus);
    return ret;
  }
  */

  // do two rows of C at once
  ci_t i;
  const ci_t c_ncols    = C->ncols;
  const ci_t start_idx  = C->pos[idx1][0] < C->pos[idx2][0] ? C->pos[idx1][0] : C->pos[idx2][0];
  ci_t nze_ctr1 = 0, nze_ctr2 = 0;
  re_l_t *wide_row1, *wide_row2;
  re_m_t multiplier1, multiplier2;
  init_wide_rows(&wide_row1, c_ncols);
  init_wide_rows(&wide_row2, c_ncols);
  memset(wide_row1, 0, c_ncols * sizeof(re_l_t));
  memset(wide_row2, 0, c_ncols * sizeof(re_l_t));
  
  copy_sparse_to_wide_row(wide_row1, C, idx1);
  copy_sparse_to_wide_row(wide_row2, C, idx2);
  for (i=start_idx; i<c_ncols-1; ++i) {
    multiplier1 = MODP(wide_row1[i], modulus);
    multiplier2 = MODP(wide_row2[i], modulus);
    if (multiplier1 != 0 && multiplier2 != 0) {
      multiplier1  = (re_m_t)(modulus - multiplier1);
      multiplier2  = (re_m_t)(modulus - multiplier2);
      update_wide_rows(wide_row1, wide_row2, A, multiplier1, multiplier2, i);
      nze_ctr1++;
      nze_ctr2++;
    } else {
      if (multiplier1 != 0) {
        multiplier1  = (re_m_t)(modulus - multiplier1);
        update_wide_row(wide_row1, A, multiplier1, i);
        nze_ctr1++;
      }
      if (multiplier2 != 0) {
        multiplier2  = (re_m_t)(modulus - multiplier2);
        update_wide_row(wide_row2, A, multiplier2, i);
        nze_ctr2++;
      }
    }
    wide_row1[i] = multiplier1;
    wide_row2[i] = multiplier2;
  }
  multiplier1          = MODP(wide_row1[c_ncols-1], modulus);
  multiplier2          = MODP(wide_row2[c_ncols-1], modulus);
  if (multiplier1 != 0) {
    multiplier1  = (re_m_t)(modulus - multiplier1);
    nze_ctr1++;
  }
  if (multiplier2 != 0) {
    multiplier2  = (re_m_t)(modulus - multiplier2);
    nze_ctr2++;
  }
  wide_row1[c_ncols-1] = multiplier1;
  wide_row2[c_ncols-1] = multiplier2;

  copy_wide_to_sparse_row(&C, wide_row2, idx2, nze_ctr2);
  copy_wide_to_sparse_row(&C, wide_row1, idx1, nze_ctr1);

  free(wide_row1);
  free(wide_row2);
  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_five(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx0, const ri_t idx1, const ri_t idx2, const ri_t idx3,
    const ri_t idx4, const mod_t modulus)
{
  if (idx1 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx0, modulus);
    return ret;
  }
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx0, idx1, modulus);
    return ret;
  }
  if (idx3 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx0, idx1, idx2, modulus);
    return ret;
  }
  if (idx4 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    return ret;
  }

  ci_t i;
  int j;
  const ci_t c_ncols    = C->ncols;
  ci_t min_01   = C->pos[idx0][0] < C->pos[idx1][0] ? C->pos[idx0][0] : C->pos[idx1][0];
  ci_t min_23   = C->pos[idx2][0] < C->pos[idx3][0] ? C->pos[idx2][0] : C->pos[idx3][0];
  ci_t min_0123  = min_01 < min_23 ? min_01 : min_23;
  const ci_t start_idx = min_0123 < C->pos[idx4][0] ? min_0123 : C->pos[idx4][0];
  ci_t nze_ctr[5] = {0};
  re_l_t **wide_row = (re_l_t **)malloc(5 * sizeof(re_l_t *));
  re_m_t multiplier[5];
  init_wide_rows(&wide_row[0], c_ncols);
  init_wide_rows(&wide_row[1], c_ncols);
  init_wide_rows(&wide_row[2], c_ncols);
  init_wide_rows(&wide_row[3], c_ncols);
  init_wide_rows(&wide_row[4], c_ncols);
  memset(wide_row[0], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[1], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[2], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[3], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[4], 0, c_ncols * sizeof(re_l_t));

  re_l_t **rows   = (re_l_t **)malloc(5 * sizeof(re_l_t *));
  re_m_t *mults   = (re_m_t *)malloc(5 * sizeof(re_m_t));
  int *row_to_id  = (int *)malloc(5 * sizeof(int));
  int ctr;

  copy_sparse_to_wide_row(wide_row[0], C, idx0);
  copy_sparse_to_wide_row(wide_row[1], C, idx1);
  copy_sparse_to_wide_row(wide_row[2], C, idx2);
  copy_sparse_to_wide_row(wide_row[3], C, idx3);
  copy_sparse_to_wide_row(wide_row[4], C, idx4);
  for (i=start_idx; i<c_ncols-1; ++i) {
    ctr = 0;
    multiplier[0] = MODP(wide_row[0][i], modulus);
    multiplier[1] = MODP(wide_row[1][i], modulus);
    multiplier[2] = MODP(wide_row[2][i], modulus);
    multiplier[3] = MODP(wide_row[3][i], modulus);
    multiplier[4] = MODP(wide_row[4][i], modulus);

    wide_row[0][i] = 0;
    wide_row[1][i] = 0;
    wide_row[2][i] = 0;
    wide_row[3][i] = 0;
    wide_row[4][i] = 0;

    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0 && multiplier[3] == 0 &&
        multiplier[4] == 0)
      continue;
    if (multiplier[0] != 0) {
      rows[ctr]       = wide_row[0];
      mults[ctr]      = (re_m_t)(modulus - multiplier[0]);
      row_to_id[ctr]  = 0;
      ++nze_ctr[0];
      ++ctr;
    }
    if (multiplier[1] != 0) {
      rows[ctr]       = wide_row[1];
      mults[ctr]      = (re_m_t)(modulus - multiplier[1]);
      row_to_id[ctr]  = 1;
      ++nze_ctr[1];
      ++ctr;
    }
    if (multiplier[2] != 0) {
      rows[ctr]       = wide_row[2];
      mults[ctr]      = (re_m_t)(modulus - multiplier[2]);
      row_to_id[ctr]  = 2;
      ++nze_ctr[2];
      ++ctr;
    }
    if (multiplier[3] != 0) {
      rows[ctr]       = wide_row[3];
      mults[ctr]      = (re_m_t)(modulus - multiplier[3]);
      row_to_id[ctr]  = 3;
      ++nze_ctr[3];
      ++ctr;
    }
    if (multiplier[4] != 0) {
      rows[ctr]       = wide_row[4];
      mults[ctr]      = (re_m_t)(modulus - multiplier[4]);
      row_to_id[ctr]  = 4;
      ++nze_ctr[4];
      ++ctr;
    }
    if (ctr == 5)
      update_wide_rows_five(rows[0], rows[1], rows[2], rows[3], rows[4], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], i);
    if (ctr == 4)
      update_wide_rows_four(rows[0], rows[1], rows[2], rows[3], A, mults[0],
          mults[1], mults[2], mults[3], i);
    if (ctr == 3)
      update_wide_rows_three(rows[0], rows[1], rows[2], A, mults[0],
          mults[1], mults[2], i);
    if (ctr == 2)
      update_wide_rows(rows[0], rows[1], A, mults[0], mults[1], i);
    if (ctr == 1)
      update_wide_row(rows[0], A, mults[0], i);

    for (j=0; j<ctr; ++j) {
      wide_row[row_to_id[j]]     = rows[j]; 
      wide_row[row_to_id[j]][i]  = mults[j]; 
    }
  }
  multiplier[0] = MODP(wide_row[0][c_ncols-1], modulus);
  multiplier[1] = MODP(wide_row[1][c_ncols-1], modulus);
  multiplier[2] = MODP(wide_row[2][c_ncols-1], modulus);
  multiplier[3] = MODP(wide_row[3][c_ncols-1], modulus);
  multiplier[4] = MODP(wide_row[4][c_ncols-1], modulus);
  if (multiplier[0] != 0) {
    multiplier[0]  = (re_m_t)(modulus - multiplier[0]);
    nze_ctr[0]++;
  }
  if (multiplier[1] != 0) {
    multiplier[1]  = (re_m_t)(modulus - multiplier[1]);
    nze_ctr[1]++;
  }
  if (multiplier[2] != 0) {
    multiplier[2]  = (re_m_t)(modulus - multiplier[2]);
    nze_ctr[2]++;
  }
  if (multiplier[3] != 0) {
    multiplier[3]  = (re_m_t)(modulus - multiplier[3]);
    nze_ctr[3]++;
  }
  if (multiplier[4] != 0) {
    multiplier[4]  = (re_m_t)(modulus - multiplier[4]);
    nze_ctr[4]++;
  }
  wide_row[0][c_ncols-1] = multiplier[0];
  wide_row[1][c_ncols-1] = multiplier[1];
  wide_row[2][c_ncols-1] = multiplier[2];
  wide_row[3][c_ncols-1] = multiplier[3];
  wide_row[4][c_ncols-1] = multiplier[4];

  //printf("idx %u -- nze[0] = %u\n",idx0, nze_ctr[0]);
  copy_wide_to_sparse_row(&C, wide_row[0], idx0, nze_ctr[0]);
  //printf("idx %u -- nze[1] = %u\n", idx1, nze_ctr[1]);
  copy_wide_to_sparse_row(&C, wide_row[1], idx1, nze_ctr[1]);
  //printf("idx %u -- nze[2] = %u\n", idx2, nze_ctr[2]);
  copy_wide_to_sparse_row(&C, wide_row[2], idx2, nze_ctr[2]);
  //printf("idx %u -- nze[3] = %u\n", idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[3], idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[4], idx4, nze_ctr[4]);

  free(wide_row[0]);
  free(wide_row[1]);
  free(wide_row[2]);
  free(wide_row[3]);
  free(wide_row[4]);
  free(wide_row);

  free(rows);
  free(row_to_id);

  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_six(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx0, const ri_t idx1, const ri_t idx2, const ri_t idx3,
    const ri_t idx4, const ri_t idx5, const mod_t modulus)
{
  /*
  if (idx1 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx0, modulus);
    return ret;
  }
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx0, idx1, modulus);
    return ret;
  }
  if (idx3 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx0, idx1, idx2, modulus);
    return ret;
  }
  if (idx4 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    return ret;
  }
  if (idx5 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx4, modulus);
    return ret;
  }
  */
  ci_t i;
  int j;
  const ci_t c_ncols    = C->ncols;
  ci_t min_01   = C->pos[idx0][0] < C->pos[idx1][0] ? C->pos[idx0][0] : C->pos[idx1][0];
  ci_t min_23   = C->pos[idx2][0] < C->pos[idx3][0] ? C->pos[idx2][0] : C->pos[idx3][0];
  ci_t min_45   = C->pos[idx4][0] < C->pos[idx5][0] ? C->pos[idx4][0] : C->pos[idx5][0];
  ci_t min_0123  = min_01 < min_23 ? min_01 : min_23;
  const ci_t start_idx = min_0123 < min_45 ? min_0123 : min_45;
  ci_t nze_ctr[6] = {0};
  re_l_t **wide_row = (re_l_t **)malloc(6 * sizeof(re_l_t *));
  re_m_t multiplier[6];
  init_wide_rows(&wide_row[0], c_ncols);
  init_wide_rows(&wide_row[1], c_ncols);
  init_wide_rows(&wide_row[2], c_ncols);
  init_wide_rows(&wide_row[3], c_ncols);
  init_wide_rows(&wide_row[4], c_ncols);
  init_wide_rows(&wide_row[5], c_ncols);
  memset(wide_row[0], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[1], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[2], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[3], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[4], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[5], 0, c_ncols * sizeof(re_l_t));

  re_l_t **rows   = (re_l_t **)malloc(6 * sizeof(re_l_t *));
  re_m_t *mults   = (re_m_t *)malloc(6 * sizeof(re_m_t));
  int *row_to_id  = (int *)malloc(6 * sizeof(int));
  int ctr;

  copy_sparse_to_wide_row(wide_row[0], C, idx0);
  copy_sparse_to_wide_row(wide_row[1], C, idx1);
  copy_sparse_to_wide_row(wide_row[2], C, idx2);
  copy_sparse_to_wide_row(wide_row[3], C, idx3);
  copy_sparse_to_wide_row(wide_row[4], C, idx4);
  copy_sparse_to_wide_row(wide_row[5], C, idx5);
  for (i=start_idx; i<c_ncols-1; ++i) {
    ctr = 0;
    multiplier[0] = MODP(wide_row[0][i], modulus);
    multiplier[1] = MODP(wide_row[1][i], modulus);
    multiplier[2] = MODP(wide_row[2][i], modulus);
    multiplier[3] = MODP(wide_row[3][i], modulus);
    multiplier[4] = MODP(wide_row[4][i], modulus);
    multiplier[5] = MODP(wide_row[5][i], modulus);

    wide_row[0][i] = 0;
    wide_row[1][i] = 0;
    wide_row[2][i] = 0;
    wide_row[3][i] = 0;
    wide_row[4][i] = 0;
    wide_row[5][i] = 0;

    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0 && multiplier[3] == 0 &&
        multiplier[4] == 0 && multiplier[5] == 0)
      continue;
    if (multiplier[0] != 0) {
      rows[ctr]       = wide_row[0];
      mults[ctr]      = (re_m_t)(modulus - multiplier[0]);
      row_to_id[ctr]  = 0;
      ++nze_ctr[0];
      ++ctr;
    }
    if (multiplier[1] != 0) {
      rows[ctr]       = wide_row[1];
      mults[ctr]      = (re_m_t)(modulus - multiplier[1]);
      row_to_id[ctr]  = 1;
      ++nze_ctr[1];
      ++ctr;
    }
    if (multiplier[2] != 0) {
      rows[ctr]       = wide_row[2];
      mults[ctr]      = (re_m_t)(modulus - multiplier[2]);
      row_to_id[ctr]  = 2;
      ++nze_ctr[2];
      ++ctr;
    }
    if (multiplier[3] != 0) {
      rows[ctr]       = wide_row[3];
      mults[ctr]      = (re_m_t)(modulus - multiplier[3]);
      row_to_id[ctr]  = 3;
      ++nze_ctr[3];
      ++ctr;
    }
    if (multiplier[4] != 0) {
      rows[ctr]       = wide_row[4];
      mults[ctr]      = (re_m_t)(modulus - multiplier[4]);
      row_to_id[ctr]  = 4;
      ++nze_ctr[4];
      ++ctr;
    }
    if (multiplier[5] != 0) {
      rows[ctr]       = wide_row[5];
      mults[ctr]      = (re_m_t)(modulus - multiplier[5]);
      row_to_id[ctr]  = 5;
      ++nze_ctr[5];
      ++ctr;
    }
    if (ctr == 6)
      update_wide_rows_six(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], mults[5], i);
    if (ctr == 5)
      update_wide_rows_five(rows[0], rows[1], rows[2], rows[3], rows[4], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], i);
    if (ctr == 4)
      update_wide_rows_four(rows[0], rows[1], rows[2], rows[3], A, mults[0],
          mults[1], mults[2], mults[3], i);
    if (ctr == 3)
      update_wide_rows_three(rows[0], rows[1], rows[2], A, mults[0],
          mults[1], mults[2], i);
    if (ctr == 2)
      update_wide_rows(rows[0], rows[1], A, mults[0], mults[1], i);
    if (ctr == 1)
      update_wide_row(rows[0], A, mults[0], i);

    for (j=0; j<ctr; ++j) {
      wide_row[row_to_id[j]]     = rows[j]; 
      wide_row[row_to_id[j]][i]  = mults[j]; 
    }
  }
  multiplier[0] = MODP(wide_row[0][c_ncols-1], modulus);
  multiplier[1] = MODP(wide_row[1][c_ncols-1], modulus);
  multiplier[2] = MODP(wide_row[2][c_ncols-1], modulus);
  multiplier[3] = MODP(wide_row[3][c_ncols-1], modulus);
  multiplier[4] = MODP(wide_row[4][c_ncols-1], modulus);
  multiplier[5] = MODP(wide_row[5][c_ncols-1], modulus);
  if (multiplier[0] != 0) {
    multiplier[0]  = (re_m_t)(modulus - multiplier[0]);
    nze_ctr[0]++;
  }
  if (multiplier[1] != 0) {
    multiplier[1]  = (re_m_t)(modulus - multiplier[1]);
    nze_ctr[1]++;
  }
  if (multiplier[2] != 0) {
    multiplier[2]  = (re_m_t)(modulus - multiplier[2]);
    nze_ctr[2]++;
  }
  if (multiplier[3] != 0) {
    multiplier[3]  = (re_m_t)(modulus - multiplier[3]);
    nze_ctr[3]++;
  }
  if (multiplier[4] != 0) {
    multiplier[4]  = (re_m_t)(modulus - multiplier[4]);
    nze_ctr[4]++;
  }
  if (multiplier[5] != 0) {
    multiplier[5]  = (re_m_t)(modulus - multiplier[5]);
    nze_ctr[5]++;
  }
  wide_row[0][c_ncols-1] = multiplier[0];
  wide_row[1][c_ncols-1] = multiplier[1];
  wide_row[2][c_ncols-1] = multiplier[2];
  wide_row[3][c_ncols-1] = multiplier[3];
  wide_row[4][c_ncols-1] = multiplier[4];
  wide_row[5][c_ncols-1] = multiplier[5];

  //printf("idx %u -- nze[0] = %u\n",idx0, nze_ctr[0]);
  copy_wide_to_sparse_row(&C, wide_row[0], idx0, nze_ctr[0]);
  //printf("idx %u -- nze[1] = %u\n", idx1, nze_ctr[1]);
  copy_wide_to_sparse_row(&C, wide_row[1], idx1, nze_ctr[1]);
  //printf("idx %u -- nze[2] = %u\n", idx2, nze_ctr[2]);
  copy_wide_to_sparse_row(&C, wide_row[2], idx2, nze_ctr[2]);
  //printf("idx %u -- nze[3] = %u\n", idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[3], idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[4], idx4, nze_ctr[4]);
  copy_wide_to_sparse_row(&C, wide_row[5], idx5, nze_ctr[5]);

  free(wide_row[0]);
  free(wide_row[1]);
  free(wide_row[2]);
  free(wide_row[3]);
  free(wide_row[4]);
  free(wide_row[5]);
  free(wide_row);

  free(rows);
  free(row_to_id);

  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_eight(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx0, const ri_t idx1, const ri_t idx2, const ri_t idx3,
    const ri_t idx4, const ri_t idx5, const ri_t idx6, const ri_t idx7,
    const mod_t modulus)
{
  if (idx1 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx0, modulus);
    return ret;
  }
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx0, idx1, modulus);
    return ret;
  }
  if (idx3 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx0, idx1, idx2, modulus);
    return ret;
  }
  if (idx4 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    return ret;
  }
  if (idx5 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx4, modulus);
    return ret;
  }
  if (idx6 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx4, idx5, modulus);
    return ret;
  }
  if (idx7 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx4, idx5, idx6, modulus);
    return ret;
  }

  ci_t i;
  int j;
  const ci_t c_ncols    = C->ncols;
  ci_t min_01   = C->pos[idx0][0] < C->pos[idx1][0] ? C->pos[idx0][0] : C->pos[idx1][0];
  ci_t min_23   = C->pos[idx2][0] < C->pos[idx3][0] ? C->pos[idx2][0] : C->pos[idx3][0];
  ci_t min_45   = C->pos[idx4][0] < C->pos[idx5][0] ? C->pos[idx4][0] : C->pos[idx5][0];
  ci_t min_67   = C->pos[idx6][0] < C->pos[idx7][0] ? C->pos[idx6][0] : C->pos[idx7][0];
  ci_t min_0123  = min_01 < min_23 ? min_01 : min_23;
  ci_t min_4567  = min_45 < min_67 ? min_45 : min_67;
  const ci_t start_idx = min_0123 < min_4567 ? min_0123 : min_4567;
  ci_t nze_ctr[8] = {0};
  re_l_t **wide_row = (re_l_t **)malloc(8 * sizeof(re_l_t *));
  re_m_t multiplier[8];
  init_wide_rows(&wide_row[0], c_ncols);
  init_wide_rows(&wide_row[1], c_ncols);
  init_wide_rows(&wide_row[2], c_ncols);
  init_wide_rows(&wide_row[3], c_ncols);
  init_wide_rows(&wide_row[4], c_ncols);
  init_wide_rows(&wide_row[5], c_ncols);
  init_wide_rows(&wide_row[6], c_ncols);
  init_wide_rows(&wide_row[7], c_ncols);
  memset(wide_row[0], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[1], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[2], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[3], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[4], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[5], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[6], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[7], 0, c_ncols * sizeof(re_l_t));

  re_l_t **rows   = (re_l_t **)malloc(8 * sizeof(re_l_t *));
  re_m_t *mults   = (re_m_t *)malloc(8 * sizeof(re_m_t));
  int *row_to_id  = (int *)malloc(8 * sizeof(int));
  int ctr;

  copy_sparse_to_wide_row(wide_row[0], C, idx0);
  copy_sparse_to_wide_row(wide_row[1], C, idx1);
  copy_sparse_to_wide_row(wide_row[2], C, idx2);
  copy_sparse_to_wide_row(wide_row[3], C, idx3);
  copy_sparse_to_wide_row(wide_row[4], C, idx4);
  copy_sparse_to_wide_row(wide_row[5], C, idx5);
  copy_sparse_to_wide_row(wide_row[6], C, idx6);
  copy_sparse_to_wide_row(wide_row[7], C, idx7);
  for (i=start_idx; i<c_ncols-1; ++i) {
    ctr = 0;
    multiplier[0] = MODP(wide_row[0][i], modulus);
    multiplier[1] = MODP(wide_row[1][i], modulus);
    multiplier[2] = MODP(wide_row[2][i], modulus);
    multiplier[3] = MODP(wide_row[3][i], modulus);
    multiplier[4] = MODP(wide_row[4][i], modulus);
    multiplier[5] = MODP(wide_row[5][i], modulus);
    multiplier[6] = MODP(wide_row[6][i], modulus);
    multiplier[7] = MODP(wide_row[7][i], modulus);

    wide_row[0][i] = 0;
    wide_row[1][i] = 0;
    wide_row[2][i] = 0;
    wide_row[3][i] = 0;
    wide_row[4][i] = 0;
    wide_row[5][i] = 0;
    wide_row[6][i] = 0;
    wide_row[7][i] = 0;

    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0 && multiplier[3] == 0 &&
        multiplier[4] == 0 && multiplier[5] == 0 && multiplier[6] == 0 && multiplier[7] == 0)
      continue;
    if (multiplier[0] != 0) {
      rows[ctr]       = wide_row[0];
      mults[ctr]      = (re_m_t)(modulus - multiplier[0]);
      row_to_id[ctr]  = 0;
      ++nze_ctr[0];
      ++ctr;
    }
    if (multiplier[1] != 0) {
      rows[ctr]       = wide_row[1];
      mults[ctr]      = (re_m_t)(modulus - multiplier[1]);
      row_to_id[ctr]  = 1;
      ++nze_ctr[1];
      ++ctr;
    }
    if (multiplier[2] != 0) {
      rows[ctr]       = wide_row[2];
      mults[ctr]      = (re_m_t)(modulus - multiplier[2]);
      row_to_id[ctr]  = 2;
      ++nze_ctr[2];
      ++ctr;
    }
    if (multiplier[3] != 0) {
      rows[ctr]       = wide_row[3];
      mults[ctr]      = (re_m_t)(modulus - multiplier[3]);
      row_to_id[ctr]  = 3;
      ++nze_ctr[3];
      ++ctr;
    }
    if (multiplier[4] != 0) {
      rows[ctr]       = wide_row[4];
      mults[ctr]      = (re_m_t)(modulus - multiplier[4]);
      row_to_id[ctr]  = 4;
      ++nze_ctr[4];
      ++ctr;
    }
    if (multiplier[5] != 0) {
      rows[ctr]       = wide_row[5];
      mults[ctr]      = (re_m_t)(modulus - multiplier[5]);
      row_to_id[ctr]  = 5;
      ++nze_ctr[5];
      ++ctr;
    }
    if (multiplier[6] != 0) {
      rows[ctr]       = wide_row[6];
      mults[ctr]      = (re_m_t)(modulus - multiplier[6]);
      row_to_id[ctr]  = 6;
      ++nze_ctr[6];
      ++ctr;
    }
    if (multiplier[7] != 0) {
      rows[ctr]       = wide_row[7];
      mults[ctr]      = (re_m_t)(modulus - multiplier[7]);
      row_to_id[ctr]  = 7;
      ++nze_ctr[7];
      ++ctr;
    }
    if (ctr == 8)
      update_wide_rows_eight(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5],
          rows[6], rows[7], A, mults[0], mults[1], mults[2], mults[3], mults[4],
          mults[5], mults[6], mults[7], i);
    if (ctr == 7)
      update_wide_rows_seven(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5],
          rows[6], A, mults[0], mults[1], mults[2], mults[3], mults[4], mults[5],
          mults[6], i);
    if (ctr == 6)
      update_wide_rows_six(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], mults[5], i);
    if (ctr == 5)
      update_wide_rows_five(rows[0], rows[1], rows[2], rows[3], rows[4], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], i);
    if (ctr == 4)
      update_wide_rows_four(rows[0], rows[1], rows[2], rows[3], A, mults[0],
          mults[1], mults[2], mults[3], i);
    if (ctr == 3)
      update_wide_rows_three(rows[0], rows[1], rows[2], A, mults[0],
          mults[1], mults[2], i);
    if (ctr == 2)
      update_wide_rows(rows[0], rows[1], A, mults[0], mults[1], i);
    if (ctr == 1)
      update_wide_row(rows[0], A, mults[0], i);

    for (j=0; j<ctr; ++j) {
      wide_row[row_to_id[j]]     = rows[j]; 
      wide_row[row_to_id[j]][i]  = mults[j]; 
    }
  }
  multiplier[0] = MODP(wide_row[0][c_ncols-1], modulus);
  multiplier[1] = MODP(wide_row[1][c_ncols-1], modulus);
  multiplier[2] = MODP(wide_row[2][c_ncols-1], modulus);
  multiplier[3] = MODP(wide_row[3][c_ncols-1], modulus);
  multiplier[4] = MODP(wide_row[4][c_ncols-1], modulus);
  multiplier[5] = MODP(wide_row[5][c_ncols-1], modulus);
  multiplier[6] = MODP(wide_row[6][c_ncols-1], modulus);
  multiplier[7] = MODP(wide_row[7][c_ncols-1], modulus);
  if (multiplier[0] != 0) {
    multiplier[0]  = (re_m_t)(modulus - multiplier[0]);
    nze_ctr[0]++;
  }
  if (multiplier[1] != 0) {
    multiplier[1]  = (re_m_t)(modulus - multiplier[1]);
    nze_ctr[1]++;
  }
  if (multiplier[2] != 0) {
    multiplier[2]  = (re_m_t)(modulus - multiplier[2]);
    nze_ctr[2]++;
  }
  if (multiplier[3] != 0) {
    multiplier[3]  = (re_m_t)(modulus - multiplier[3]);
    nze_ctr[3]++;
  }
  if (multiplier[4] != 0) {
    multiplier[4]  = (re_m_t)(modulus - multiplier[4]);
    nze_ctr[4]++;
  }
  if (multiplier[5] != 0) {
    multiplier[5]  = (re_m_t)(modulus - multiplier[5]);
    nze_ctr[5]++;
  }
  if (multiplier[6] != 0) {
    multiplier[6]  = (re_m_t)(modulus - multiplier[6]);
    nze_ctr[6]++;
  }
  if (multiplier[7] != 0) {
    multiplier[7]  = (re_m_t)(modulus - multiplier[7]);
    nze_ctr[7]++;
  }
  wide_row[0][c_ncols-1] = multiplier[0];
  wide_row[1][c_ncols-1] = multiplier[1];
  wide_row[2][c_ncols-1] = multiplier[2];
  wide_row[3][c_ncols-1] = multiplier[3];
  wide_row[4][c_ncols-1] = multiplier[4];
  wide_row[5][c_ncols-1] = multiplier[5];
  wide_row[6][c_ncols-1] = multiplier[6];
  wide_row[7][c_ncols-1] = multiplier[7];

  //printf("idx %u -- nze[0] = %u\n",idx0, nze_ctr[0]);
  copy_wide_to_sparse_row(&C, wide_row[0], idx0, nze_ctr[0]);
  //printf("idx %u -- nze[1] = %u\n", idx1, nze_ctr[1]);
  copy_wide_to_sparse_row(&C, wide_row[1], idx1, nze_ctr[1]);
  //printf("idx %u -- nze[2] = %u\n", idx2, nze_ctr[2]);
  copy_wide_to_sparse_row(&C, wide_row[2], idx2, nze_ctr[2]);
  //printf("idx %u -- nze[3] = %u\n", idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[3], idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[4], idx4, nze_ctr[4]);
  copy_wide_to_sparse_row(&C, wide_row[5], idx5, nze_ctr[5]);
  copy_wide_to_sparse_row(&C, wide_row[6], idx6, nze_ctr[6]);
  copy_wide_to_sparse_row(&C, wide_row[7], idx7, nze_ctr[7]);

  free(wide_row[0]);
  free(wide_row[1]);
  free(wide_row[2]);
  free(wide_row[3]);
  free(wide_row[4]);
  free(wide_row[5]);
  free(wide_row[6]);
  free(wide_row[7]);
  free(wide_row);

  free(rows);
  free(row_to_id);

  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_ten(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx0, const ri_t idx1, const ri_t idx2, const ri_t idx3,
    const ri_t idx4, const ri_t idx5, const ri_t idx6, const ri_t idx7,
    const ri_t idx8, const ri_t idx9, const mod_t modulus)
{
  /*
  if (idx1 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx0, modulus);
    return ret;
  }
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx0, idx1, modulus);
    return ret;
  }
  if (idx3 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx0, idx1, idx2, modulus);
    return ret;
  }
  if (idx4 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    return ret;
  }
  if (idx5 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx4, modulus);
    return ret;
  }
  if (idx6 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx4, idx5, modulus);
    return ret;
  }
  if (idx7 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx4, idx5, idx6, modulus);
    return ret;
  }
  if (idx8 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx4, idx5, idx6, idx7, modulus);
    return ret;
  }
  if (idx9 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx0, idx1, idx2, idx3, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_four(C, A, idx4, idx5, idx6, idx7, modulus);
    ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx8, modulus);
    return ret;
  }
  */
  ci_t i;
  int j;
  const ci_t c_ncols    = C->ncols;
  ci_t min_01   = C->pos[idx0][0] < C->pos[idx1][0] ? C->pos[idx0][0] : C->pos[idx1][0];
  ci_t min_23   = C->pos[idx2][0] < C->pos[idx3][0] ? C->pos[idx2][0] : C->pos[idx3][0];
  ci_t min_45   = C->pos[idx4][0] < C->pos[idx5][0] ? C->pos[idx4][0] : C->pos[idx5][0];
  ci_t min_67   = C->pos[idx6][0] < C->pos[idx7][0] ? C->pos[idx6][0] : C->pos[idx7][0];
  ci_t min_89   = C->pos[idx8][0] < C->pos[idx9][0] ? C->pos[idx8][0] : C->pos[idx9][0];
  ci_t min_0123  = min_01 < min_23 ? min_01 : min_23;
  ci_t min_4567  = min_45 < min_67 ? min_45 : min_67;
  ci_t min_01234567  = min_0123 < min_4567 ? min_0123 : min_4567;
  const ci_t start_idx  = min_01234567 < min_89 ? min_01234567 : min_89;
  ci_t nze_ctr[10] = {0};
  re_l_t **wide_row = (re_l_t **)malloc(10 * sizeof(re_l_t *));
  re_m_t multiplier[10];
  init_wide_rows(&wide_row[0], c_ncols);
  init_wide_rows(&wide_row[1], c_ncols);
  init_wide_rows(&wide_row[2], c_ncols);
  init_wide_rows(&wide_row[3], c_ncols);
  init_wide_rows(&wide_row[4], c_ncols);
  init_wide_rows(&wide_row[5], c_ncols);
  init_wide_rows(&wide_row[6], c_ncols);
  init_wide_rows(&wide_row[7], c_ncols);
  init_wide_rows(&wide_row[8], c_ncols);
  init_wide_rows(&wide_row[9], c_ncols);
  memset(wide_row[0], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[1], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[2], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[3], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[4], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[5], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[6], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[7], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[8], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[9], 0, c_ncols * sizeof(re_l_t));

  re_l_t **rows   = (re_l_t **)malloc(10 * sizeof(re_l_t *));
  re_m_t *mults   = (re_m_t *)malloc(10 * sizeof(re_m_t));
  int *row_to_id  = (int *)malloc(10 * sizeof(int));
  int ctr;

  copy_sparse_to_wide_row(wide_row[0], C, idx0);
  copy_sparse_to_wide_row(wide_row[1], C, idx1);
  copy_sparse_to_wide_row(wide_row[2], C, idx2);
  copy_sparse_to_wide_row(wide_row[3], C, idx3);
  copy_sparse_to_wide_row(wide_row[4], C, idx4);
  copy_sparse_to_wide_row(wide_row[5], C, idx5);
  copy_sparse_to_wide_row(wide_row[6], C, idx6);
  copy_sparse_to_wide_row(wide_row[7], C, idx7);
  copy_sparse_to_wide_row(wide_row[8], C, idx8);
  copy_sparse_to_wide_row(wide_row[9], C, idx9);
  for (i=start_idx; i<c_ncols-1; ++i) {
    ctr = 0;
    multiplier[0] = MODP(wide_row[0][i], modulus);
    multiplier[1] = MODP(wide_row[1][i], modulus);
    multiplier[2] = MODP(wide_row[2][i], modulus);
    multiplier[3] = MODP(wide_row[3][i], modulus);
    multiplier[4] = MODP(wide_row[4][i], modulus);
    multiplier[5] = MODP(wide_row[5][i], modulus);
    multiplier[6] = MODP(wide_row[6][i], modulus);
    multiplier[7] = MODP(wide_row[7][i], modulus);
    multiplier[8] = MODP(wide_row[8][i], modulus);
    multiplier[9] = MODP(wide_row[9][i], modulus);

    wide_row[0][i] = 0;
    wide_row[1][i] = 0;
    wide_row[2][i] = 0;
    wide_row[3][i] = 0;
    wide_row[4][i] = 0;
    wide_row[5][i] = 0;
    wide_row[6][i] = 0;
    wide_row[7][i] = 0;
    wide_row[8][i] = 0;
    wide_row[9][i] = 0;

    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0 && multiplier[3] == 0 &&
        multiplier[4] == 0 && multiplier[5] == 0 && multiplier[6] == 0 && multiplier[7] == 0 &&
        multiplier[8] == 0 && multiplier[9] == 0)
      continue;
    if (multiplier[0] != 0) {
      rows[ctr]       = wide_row[0];
      mults[ctr]      = (re_m_t)(modulus - multiplier[0]);
      row_to_id[ctr]  = 0;
      ++nze_ctr[0];
      ++ctr;
    }
    if (multiplier[1] != 0) {
      rows[ctr]       = wide_row[1];
      mults[ctr]      = (re_m_t)(modulus - multiplier[1]);
      row_to_id[ctr]  = 1;
      ++nze_ctr[1];
      ++ctr;
    }
    if (multiplier[2] != 0) {
      rows[ctr]       = wide_row[2];
      mults[ctr]      = (re_m_t)(modulus - multiplier[2]);
      row_to_id[ctr]  = 2;
      ++nze_ctr[2];
      ++ctr;
    }
    if (multiplier[3] != 0) {
      rows[ctr]       = wide_row[3];
      mults[ctr]      = (re_m_t)(modulus - multiplier[3]);
      row_to_id[ctr]  = 3;
      ++nze_ctr[3];
      ++ctr;
    }
    if (multiplier[4] != 0) {
      rows[ctr]       = wide_row[4];
      mults[ctr]      = (re_m_t)(modulus - multiplier[4]);
      row_to_id[ctr]  = 4;
      ++nze_ctr[4];
      ++ctr;
    }
    if (multiplier[5] != 0) {
      rows[ctr]       = wide_row[5];
      mults[ctr]      = (re_m_t)(modulus - multiplier[5]);
      row_to_id[ctr]  = 5;
      ++nze_ctr[5];
      ++ctr;
    }
    if (multiplier[6] != 0) {
      rows[ctr]       = wide_row[6];
      mults[ctr]      = (re_m_t)(modulus - multiplier[6]);
      row_to_id[ctr]  = 6;
      ++nze_ctr[6];
      ++ctr;
    }
    if (multiplier[7] != 0) {
      rows[ctr]       = wide_row[7];
      mults[ctr]      = (re_m_t)(modulus - multiplier[7]);
      row_to_id[ctr]  = 7;
      ++nze_ctr[7];
      ++ctr;
    }
    if (multiplier[8] != 0) {
      rows[ctr]       = wide_row[8];
      mults[ctr]      = (re_m_t)(modulus - multiplier[8]);
      row_to_id[ctr]  = 8;
      ++nze_ctr[8];
      ++ctr;
    }
    if (multiplier[9] != 0) {
      rows[ctr]       = wide_row[9];
      mults[ctr]      = (re_m_t)(modulus - multiplier[9]);
      row_to_id[ctr]  = 9;
      ++nze_ctr[9];
      ++ctr;
    }
    if (ctr == 10)
      update_wide_rows_ten(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5],
          rows[6], rows[7], rows[8], rows[9], A, mults[0], mults[1], mults[2], mults[3],
          mults[4], mults[5], mults[6], mults[7], mults[8], mults[9], i);
    if (ctr == 9)
      update_wide_rows_nine(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5],
          rows[6], rows[7], rows[8], A, mults[0], mults[1], mults[2], mults[3],
          mults[4], mults[5], mults[6], mults[7], mults[8], i);
    if (ctr == 8)
      update_wide_rows_eight(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5],
          rows[6], rows[7], A, mults[0], mults[1], mults[2], mults[3], mults[4],
          mults[5], mults[6], mults[7], i);
    if (ctr == 7)
      update_wide_rows_seven(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5],
          rows[6], A, mults[0], mults[1], mults[2], mults[3], mults[4], mults[5],
          mults[6], i);
    if (ctr == 6)
      update_wide_rows_six(rows[0], rows[1], rows[2], rows[3], rows[4], rows[5], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], mults[5], i);
    if (ctr == 5)
      update_wide_rows_five(rows[0], rows[1], rows[2], rows[3], rows[4], A,
          mults[0], mults[1], mults[2], mults[3], mults[4], i);
    if (ctr == 4)
      update_wide_rows_four(rows[0], rows[1], rows[2], rows[3], A, mults[0],
          mults[1], mults[2], mults[3], i);
    if (ctr == 3)
      update_wide_rows_three(rows[0], rows[1], rows[2], A, mults[0],
          mults[1], mults[2], i);
    if (ctr == 2)
      update_wide_rows(rows[0], rows[1], A, mults[0], mults[1], i);
    if (ctr == 1)
      update_wide_row(rows[0], A, mults[0], i);

    for (j=0; j<ctr; ++j) {
      wide_row[row_to_id[j]]     = rows[j]; 
      wide_row[row_to_id[j]][i]  = mults[j]; 
    }
  }
  multiplier[0] = MODP(wide_row[0][c_ncols-1], modulus);
  multiplier[1] = MODP(wide_row[1][c_ncols-1], modulus);
  multiplier[2] = MODP(wide_row[2][c_ncols-1], modulus);
  multiplier[3] = MODP(wide_row[3][c_ncols-1], modulus);
  multiplier[4] = MODP(wide_row[4][c_ncols-1], modulus);
  multiplier[5] = MODP(wide_row[5][c_ncols-1], modulus);
  multiplier[6] = MODP(wide_row[6][c_ncols-1], modulus);
  multiplier[7] = MODP(wide_row[7][c_ncols-1], modulus);
  multiplier[8] = MODP(wide_row[8][c_ncols-1], modulus);
  multiplier[9] = MODP(wide_row[9][c_ncols-1], modulus);
  if (multiplier[0] != 0) {
    multiplier[0]  = (re_m_t)(modulus - multiplier[0]);
    nze_ctr[0]++;
  }
  if (multiplier[1] != 0) {
    multiplier[1]  = (re_m_t)(modulus - multiplier[1]);
    nze_ctr[1]++;
  }
  if (multiplier[2] != 0) {
    multiplier[2]  = (re_m_t)(modulus - multiplier[2]);
    nze_ctr[2]++;
  }
  if (multiplier[3] != 0) {
    multiplier[3]  = (re_m_t)(modulus - multiplier[3]);
    nze_ctr[3]++;
  }
  if (multiplier[4] != 0) {
    multiplier[4]  = (re_m_t)(modulus - multiplier[4]);
    nze_ctr[4]++;
  }
  if (multiplier[5] != 0) {
    multiplier[5]  = (re_m_t)(modulus - multiplier[5]);
    nze_ctr[5]++;
  }
  if (multiplier[6] != 0) {
    multiplier[6]  = (re_m_t)(modulus - multiplier[6]);
    nze_ctr[6]++;
  }
  if (multiplier[7] != 0) {
    multiplier[7]  = (re_m_t)(modulus - multiplier[7]);
    nze_ctr[7]++;
  }
  if (multiplier[8] != 0) {
    multiplier[8]  = (re_m_t)(modulus - multiplier[8]);
    nze_ctr[8]++;
  }
  if (multiplier[9] != 0) {
    multiplier[9]  = (re_m_t)(modulus - multiplier[9]);
    nze_ctr[9]++;
  }
  wide_row[0][c_ncols-1] = multiplier[0];
  wide_row[1][c_ncols-1] = multiplier[1];
  wide_row[2][c_ncols-1] = multiplier[2];
  wide_row[3][c_ncols-1] = multiplier[3];
  wide_row[4][c_ncols-1] = multiplier[4];
  wide_row[5][c_ncols-1] = multiplier[5];
  wide_row[6][c_ncols-1] = multiplier[6];
  wide_row[7][c_ncols-1] = multiplier[7];
  wide_row[8][c_ncols-1] = multiplier[8];
  wide_row[9][c_ncols-1] = multiplier[9];

  //printf("idx %u -- nze[0] = %u\n",idx0, nze_ctr[0]);
  copy_wide_to_sparse_row(&C, wide_row[0], idx0, nze_ctr[0]);
  //printf("idx %u -- nze[1] = %u\n", idx1, nze_ctr[1]);
  copy_wide_to_sparse_row(&C, wide_row[1], idx1, nze_ctr[1]);
  //printf("idx %u -- nze[2] = %u\n", idx2, nze_ctr[2]);
  copy_wide_to_sparse_row(&C, wide_row[2], idx2, nze_ctr[2]);
  //printf("idx %u -- nze[3] = %u\n", idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[3], idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[4], idx4, nze_ctr[4]);
  copy_wide_to_sparse_row(&C, wide_row[5], idx5, nze_ctr[5]);
  copy_wide_to_sparse_row(&C, wide_row[6], idx6, nze_ctr[6]);
  copy_wide_to_sparse_row(&C, wide_row[7], idx7, nze_ctr[7]);
  copy_wide_to_sparse_row(&C, wide_row[8], idx8, nze_ctr[8]);
  copy_wide_to_sparse_row(&C, wide_row[9], idx9, nze_ctr[9]);

  free(wide_row[0]);
  free(wide_row[1]);
  free(wide_row[2]);
  free(wide_row[3]);
  free(wide_row[4]);
  free(wide_row[5]);
  free(wide_row[6]);
  free(wide_row[7]);
  free(wide_row[8]);
  free(wide_row[9]);
  free(wide_row);

  free(rows);
  free(row_to_id);

  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_four(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx0, const ri_t idx1, const ri_t idx2, const ri_t idx3,
    const mod_t modulus)
{
  /*
  if (idx1 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx0, modulus);
    return ret;
  }
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx0, idx1, modulus);
    return ret;
  }
  if (idx3 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_three(C, A, idx0, idx1, idx2, modulus);
    return ret;
  }
  */

  ci_t i;
  int j;
  const ci_t c_ncols    = C->ncols;
  ci_t min_01   = C->pos[idx0][0] < C->pos[idx1][0] ? C->pos[idx0][0] : C->pos[idx1][0];
  ci_t min_012  = C->pos[idx2][0] < min_01 ? C->pos[idx2][0] : min_01;
  const ci_t start_idx  = C->pos[idx3][0] < min_012 ? C->pos[idx3][0] : min_012;
  ci_t nze_ctr[4] = {0};
  re_l_t **wide_row = (re_l_t **)malloc(4 * sizeof(re_l_t *));
  re_m_t multiplier[4];
  init_wide_rows(&wide_row[0], c_ncols);
  init_wide_rows(&wide_row[1], c_ncols);
  init_wide_rows(&wide_row[2], c_ncols);
  init_wide_rows(&wide_row[3], c_ncols);
  memset(wide_row[0], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[1], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[2], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[3], 0, c_ncols * sizeof(re_l_t));

  re_l_t **rows   = (re_l_t **)malloc(4 * sizeof(re_l_t *));
  re_m_t *mults   = (re_m_t *)malloc(4 * sizeof(re_m_t));
  int *row_to_id  = (int *)malloc(4 * sizeof(int));
  int ctr;

  copy_sparse_to_wide_row(wide_row[0], C, idx0);
  copy_sparse_to_wide_row(wide_row[1], C, idx1);
  copy_sparse_to_wide_row(wide_row[2], C, idx2);
  copy_sparse_to_wide_row(wide_row[3], C, idx3);
  for (i=start_idx; i<c_ncols-1; ++i) {
    ctr = 0;
    multiplier[0] = MODP(wide_row[0][i], modulus);
    multiplier[1] = MODP(wide_row[1][i], modulus);
    multiplier[2] = MODP(wide_row[2][i], modulus);
    multiplier[3] = MODP(wide_row[3][i], modulus);

    wide_row[0][i] = 0;
    wide_row[1][i] = 0;
    wide_row[2][i] = 0;
    wide_row[3][i] = 0;
    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0 && multiplier[3] == 0)
      continue;
    if (multiplier[0] != 0) {
      rows[ctr]       = wide_row[0];
      mults[ctr]      = (re_m_t)(modulus - multiplier[0]);
      row_to_id[ctr]  = 0;
      ++nze_ctr[0];
      ++ctr;
    }
    if (multiplier[1] != 0) {
      rows[ctr]       = wide_row[1];
      mults[ctr]      = (re_m_t)(modulus - multiplier[1]);
      row_to_id[ctr]  = 1;
      ++nze_ctr[1];
      ++ctr;
    }
    if (multiplier[2] != 0) {
      rows[ctr]       = wide_row[2];
      mults[ctr]      = (re_m_t)(modulus - multiplier[2]);
      row_to_id[ctr]  = 2;
      ++nze_ctr[2];
      ++ctr;
    }
    if (multiplier[3] != 0) {
      rows[ctr]       = wide_row[3];
      mults[ctr]      = (re_m_t)(modulus - multiplier[3]);
      row_to_id[ctr]  = 3;
      ++nze_ctr[3];
      ++ctr;
    }
    if (ctr == 4)
      update_wide_rows_four(rows[0], rows[1], rows[2], rows[3], A, mults[0],
          mults[1], mults[2], mults[3], i);
    if (ctr == 3)
      update_wide_rows_three(rows[0], rows[1], rows[2], A, mults[0],
          mults[1], mults[2], i);
    if (ctr == 2)
      update_wide_rows(rows[0], rows[1], A, mults[0], mults[1], i);
    if (ctr == 1)
      update_wide_row(rows[0], A, mults[0], i);

    for (j=0; j<ctr; ++j) {
      wide_row[row_to_id[j]]     = rows[j]; 
      wide_row[row_to_id[j]][i]  = mults[j]; 
    }
  }
  multiplier[0] = MODP(wide_row[0][c_ncols-1], modulus);
  multiplier[1] = MODP(wide_row[1][c_ncols-1], modulus);
  multiplier[2] = MODP(wide_row[2][c_ncols-1], modulus);
  multiplier[3] = MODP(wide_row[3][c_ncols-1], modulus);
  if (multiplier[0] != 0) {
    multiplier[0]  = (re_m_t)(modulus - multiplier[0]);
    nze_ctr[0]++;
  }
  if (multiplier[1] != 0) {
    multiplier[1]  = (re_m_t)(modulus - multiplier[1]);
    nze_ctr[1]++;
  }
  if (multiplier[2] != 0) {
    multiplier[2]  = (re_m_t)(modulus - multiplier[2]);
    nze_ctr[2]++;
  }
  if (multiplier[3] != 0) {
    multiplier[3]  = (re_m_t)(modulus - multiplier[3]);
    nze_ctr[3]++;
  }
  wide_row[0][c_ncols-1] = multiplier[0];
  wide_row[1][c_ncols-1] = multiplier[1];
  wide_row[2][c_ncols-1] = multiplier[2];
  wide_row[3][c_ncols-1] = multiplier[3];

  //printf("idx %u -- nze[0] = %u\n",idx0, nze_ctr[0]);
  copy_wide_to_sparse_row(&C, wide_row[0], idx0, nze_ctr[0]);
  //printf("idx %u -- nze[1] = %u\n", idx1, nze_ctr[1]);
  copy_wide_to_sparse_row(&C, wide_row[1], idx1, nze_ctr[1]);
  //printf("idx %u -- nze[2] = %u\n", idx2, nze_ctr[2]);
  copy_wide_to_sparse_row(&C, wide_row[2], idx2, nze_ctr[2]);
  //printf("idx %u -- nze[3] = %u\n", idx3, nze_ctr[3]);
  copy_wide_to_sparse_row(&C, wide_row[3], idx3, nze_ctr[3]);

  free(wide_row[0]);
  free(wide_row[1]);
  free(wide_row[2]);
  free(wide_row[3]);
  free(wide_row);

  free(rows);
  free(row_to_id);

  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_three(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx0, const ri_t idx1, const ri_t idx2,
    const mod_t modulus)
{
  if (idx1 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_single(C, A, idx0, modulus);
    return ret;
  }
  if (idx2 >= C->nrows) {
    int ret = elim_fl_C_sparse_dense_keep_A_tasks_double(C, A, idx0, idx1, modulus);
    return ret;
  }

  // do two rows of C at once
  ci_t i;
  int j;
  const ci_t c_ncols    = C->ncols;
  ci_t min_01 = C->pos[idx0][0] < C->pos[idx1][0] ? C->pos[idx0][0] : C->pos[idx1][0];
  const ci_t start_idx  = C->pos[idx2][0] < min_01 ? C->pos[idx2][0] : min_01;
  ci_t nze_ctr0 = 0, nze_ctr1 = 0, nze_ctr2 = 0;
  re_l_t **wide_row = (re_l_t **)malloc(3 * sizeof(re_l_t *));
  re_m_t multiplier[3];
  init_wide_rows(&wide_row[0], c_ncols);
  init_wide_rows(&wide_row[1], c_ncols);
  init_wide_rows(&wide_row[2], c_ncols);
  memset(wide_row[0], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[1], 0, c_ncols * sizeof(re_l_t));
  memset(wide_row[2], 0, c_ncols * sizeof(re_l_t));

  re_l_t **rows   = (re_l_t **)malloc(3 * sizeof(re_l_t *));
  re_m_t *mults   = (re_m_t *)malloc(3 * sizeof(re_m_t));
  int *row_to_id  = (int *)malloc(3 * sizeof(int));
  int ctr;

  copy_sparse_to_wide_row(wide_row[0], C, idx0);
  copy_sparse_to_wide_row(wide_row[1], C, idx1);
  copy_sparse_to_wide_row(wide_row[2], C, idx2);
  for (i=start_idx; i<c_ncols-1; ++i) {
    rows[0] = rows[1] = rows[2] = NULL;
    ctr = 0;
    multiplier[0] = MODP(wide_row[0][i], modulus);
    multiplier[1] = MODP(wide_row[1][i], modulus);
    multiplier[2] = MODP(wide_row[2][i], modulus);
#if 1
    wide_row[0][i] = 0;
    wide_row[1][i] = 0;
    wide_row[2][i] = 0;
    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0)
      continue;
    if (multiplier[0] != 0) {
      rows[ctr]       = wide_row[0];
      mults[ctr]      = (re_m_t)(modulus - multiplier[0]);
      row_to_id[ctr]  = 0;
      ++nze_ctr0;
      ++ctr;
    }
    if (multiplier[1] != 0) {
      rows[ctr]       = wide_row[1];
      mults[ctr]      = (re_m_t)(modulus - multiplier[1]);
      row_to_id[ctr]  = 1;
      ++nze_ctr1;
      ++ctr;
    }
    if (multiplier[2] != 0) {
      rows[ctr]       = wide_row[2];
      mults[ctr]      = (re_m_t)(modulus - multiplier[2]);
      row_to_id[ctr]  = 2;
      ++nze_ctr2;
      ++ctr;
    }
    if (ctr == 3)
      update_wide_rows_three(rows[0], rows[1], rows[2], A, mults[0],
          mults[1], mults[2], i);
    if (ctr == 2)
      update_wide_rows(rows[0], rows[1], A, mults[0], mults[1], i);
    if (ctr == 1)
      update_wide_row(rows[0], A, mults[0], i);

    for (j=0; j<ctr; ++j) {
      wide_row[row_to_id[j]]     = rows[j]; 
      wide_row[row_to_id[j]][i]  = mults[j]; 
    }
  }
  multiplier[0] = MODP(wide_row[0][c_ncols-1], modulus);
  multiplier[1] = MODP(wide_row[1][c_ncols-1], modulus);
  multiplier[2] = MODP(wide_row[2][c_ncols-1], modulus);
  if (multiplier[0] != 0) {
    multiplier[0]  = (re_m_t)(modulus - multiplier[0]);
    nze_ctr0++;
  }
  if (multiplier[1] != 0) {
    multiplier[1]  = (re_m_t)(modulus - multiplier[1]);
    nze_ctr1++;
  }
  if (multiplier[2] != 0) {
    multiplier[2]  = (re_m_t)(modulus - multiplier[2]);
    nze_ctr2++;
  }
  wide_row[0][c_ncols-1] = multiplier[0];
  wide_row[1][c_ncols-1] = multiplier[1];
  wide_row[2][c_ncols-1] = multiplier[2];

  copy_wide_to_sparse_row(&C, wide_row[0], idx0, nze_ctr0);
  copy_wide_to_sparse_row(&C, wide_row[1], idx1, nze_ctr1);
  copy_wide_to_sparse_row(&C, wide_row[2], idx2, nze_ctr2);

  free(wide_row[0]);
  free(wide_row[1]);
  free(wide_row[2]);
  free(wide_row);

  free(rows);
  free(row_to_id);
#else
    if (multiplier[0] == 0 && multiplier[1] == 0 && multiplier[2] == 0) {
      wide_row[0][i] = multiplier[0];
      wide_row[1][i] = multiplier[1];
      wide_row[2][i] = multiplier[2];
      continue;
    }
    if (multiplier1 != 0 && multiplier2 != 0 && multiplier3 != 0) {
      multiplier1  = (re_m_t)(modulus - multiplier1);
      multiplier2  = (re_m_t)(modulus - multiplier2);
      multiplier3  = (re_m_t)(modulus - multiplier3);
      update_wide_rows_three(wide_row1, wide_row2, wide_row3, A,
          multiplier1, multiplier2, multiplier3, i);
      nze_ctr1++;
      nze_ctr2++;
      nze_ctr3++;
    } else {
      if (multiplier1 != 0) {
        if (multiplier2 != 0 && multiplier3 == 0) {
          multiplier1  = (re_m_t)(modulus - multiplier1);
          multiplier2  = (re_m_t)(modulus - multiplier2);
          update_wide_rows(wide_row1, wide_row2, A, multiplier1, multiplier2, i);
          nze_ctr1++;
          nze_ctr2++;
        }
        if (multiplier2 == 0 && multiplier3 != 0) {
          multiplier1  = (re_m_t)(modulus - multiplier1);
          multiplier3  = (re_m_t)(modulus - multiplier3);
          update_wide_rows(wide_row1, wide_row3, A, multiplier1, multiplier3, i);
          nze_ctr1++;
          nze_ctr3++;
        }
        if (multiplier2 == 0 && multiplier3 == 0) {
          multiplier1  = (re_m_t)(modulus - multiplier1);
          update_wide_row(wide_row1, A, multiplier1, i);
          nze_ctr1++;
        }
      } else {
        if (multiplier2 != 0 && multiplier3 != 0) {
          multiplier2  = (re_m_t)(modulus - multiplier2);
          multiplier3  = (re_m_t)(modulus - multiplier3);
          update_wide_rows(wide_row2, wide_row3, A, multiplier2, multiplier3, i);
          nze_ctr2++;
          nze_ctr3++;
        }
        if (multiplier2 != 0 && multiplier3 == 0) {
          multiplier2  = (re_m_t)(modulus - multiplier2);
          update_wide_row(wide_row2, A, multiplier2, i);
          nze_ctr2++;
        }
        if (multiplier2 == 0 && multiplier3 != 0) {
          multiplier3  = (re_m_t)(modulus - multiplier3);
          update_wide_row(wide_row3, A, multiplier3, i);
          nze_ctr3++;
        }
     }
    }
    wide_row1[i] = multiplier1;
    wide_row2[i] = multiplier2;
    wide_row3[i] = multiplier3;
  }
  multiplier1          = MODP(wide_row1[c_ncols-1], modulus);
  multiplier2          = MODP(wide_row2[c_ncols-1], modulus);
  multiplier3          = MODP(wide_row3[c_ncols-1], modulus);
  if (multiplier1 != 0) {
    multiplier1  = (re_m_t)(modulus - multiplier1);
    nze_ctr1++;
  }
  if (multiplier2 != 0) {
    multiplier2  = (re_m_t)(modulus - multiplier2);
    nze_ctr2++;
  }
  if (multiplier3 != 0) {
    multiplier3  = (re_m_t)(modulus - multiplier3);
    nze_ctr3++;
  }
  wide_row1[c_ncols-1] = multiplier1;
  wide_row2[c_ncols-1] = multiplier2;
  wide_row3[c_ncols-1] = multiplier3;
  copy_wide_to_sparse_row(&C, wide_row1, idx1, nze_ctr1);
  copy_wide_to_sparse_row(&C, wide_row2, idx2, nze_ctr2);
  copy_wide_to_sparse_row(&C, wide_row3, idx3, nze_ctr3);

  free(wide_row1);
  free(wide_row2);
  free(wide_row3);

  free(rows);
  free(row_to_id);
#endif

  return 0;
}

int elim_fl_C_sparse_dense_keep_A_tasks_single(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx, const mod_t modulus)
{
  ci_t i;
  const ci_t c_ncols    = C->ncols;
  const ci_t start_idx  = C->pos[idx][0];
  ci_t nze_ctr          = 0;
  re_l_t *wide_row;
  re_m_t multiplier;
  init_wide_rows(&wide_row, c_ncols);
  memset(wide_row, 0, c_ncols * sizeof(re_l_t));
  
  copy_sparse_to_wide_row(wide_row, C, idx);
  for (i=start_idx; i<c_ncols-1; ++i) {
    multiplier = MODP(wide_row[i], modulus);
    if (multiplier != 0) {
      multiplier  = (re_m_t)(modulus - multiplier);
      update_wide_row(wide_row, A, multiplier, i);
      nze_ctr++;
    }     
    wide_row[i] = multiplier;
  }
  multiplier          = MODP(wide_row[c_ncols-1], modulus);
  if (multiplier != 0) {
    multiplier  = (re_m_t)(modulus - multiplier);
    nze_ctr++;
  }
  wide_row[c_ncols-1] = multiplier;

  copy_wide_to_sparse_row(&C, wide_row, idx, nze_ctr);

  free(wide_row);

  return 0;
}

int elim_fl_C_ml_task(sm_fl_ml_t *C, sm_fl_ml_t *A, ri_t row_idx, mod_t modulus) {

  ci_t i;
  const ci_t coldim = C->ncols;

  int ret;

  ci_t start_idx;
  /*  Note that all our multilines are stored in a sparse fashion at the moment. */
  /*  For compatibility and later changes we keep this check in the code. */
  if (C->ml[row_idx].dense == 0)
    start_idx = C->ml[row_idx].idx[0];
  else
    start_idx = 0;

  re_m_t Cv1_col1, Cv2_col1;
  re_m_t Cv1_col2, Cv2_col2;
  uint32_t Cp1, Cp2;
  re_m_t tmp  = 0;
  ri_t row_in_A_mod, row_in_A;

  re_l_t *dense_array_C_1, *dense_array_C_2;
  ret = posix_memalign((void **)&dense_array_C_1, ALIGNT, coldim * sizeof(re_l_t));
  ret = posix_memalign((void **)&dense_array_C_2, ALIGNT, coldim * sizeof(re_l_t));

  memset(dense_array_C_1, 0, coldim * sizeof(re_l_t));
  memset(dense_array_C_2, 0, coldim * sizeof(re_l_t));

  copy_multiline_to_dense_array(C->ml[row_idx], dense_array_C_1,
      dense_array_C_2, coldim);
  for (i=start_idx; i<coldim; ++i) {
    Cp1 = i;
    Cv1_col1  = MODP(dense_array_C_1[i], modulus);
    Cv2_col1  = MODP(dense_array_C_2[i], modulus);

    if (Cv1_col1 == 0 && Cv2_col1 == 0)
      continue;

    if (Cv1_col1 != 0)
      Cv1_col1  = (re_m_t)modulus - Cv1_col1;
    if (Cv2_col1 != 0)
      Cv2_col1  = (re_m_t)modulus - Cv2_col1;

    row_in_A      = (A->nrows - 1 - Cp1) / __GB_NROWS_MULTILINE;
    row_in_A_mod  = (A->nrows - 1 - Cp1) % __GB_NROWS_MULTILINE;

    if (row_in_A_mod == 1) {
      Cp2 = i+1;
      Cv1_col2  = MODP(dense_array_C_1[i+1] , modulus);
      Cv2_col2  = MODP(dense_array_C_2[i+1] , modulus);

      register re_t v__ = A->ml[row_in_A].val[2*1+1];

      tmp = Cv1_col2 + (re_m_t)Cv1_col1 * v__;
      Cv1_col2 = MODP(tmp , modulus);
      tmp = Cv2_col2 + (re_m_t)Cv2_col1 * v__;
      Cv2_col2 = MODP(tmp , modulus);

      if (Cv1_col2 != 0)
        Cv1_col2  = (re_m_t)modulus - Cv1_col2;
      if (Cv2_col2 != 0)
        Cv2_col2  = (re_m_t)modulus - Cv2_col2;

      sparse_scal_mul_sub_2_rows_vect_array_multiline(
          Cv1_col2, Cv2_col2,
          Cv1_col1, Cv2_col1,
          A->ml[row_in_A],
          dense_array_C_1,
          dense_array_C_2);

      dense_array_C_1[Cp1]  = Cv1_col1;
      dense_array_C_1[Cp2]  = Cv1_col2;

      dense_array_C_2[Cp1]  = Cv2_col1;
      dense_array_C_2[Cp2]  = Cv2_col2;

      /*  since we have done two-in-one */
      ++i;
    } else {
      sparse_scal_mul_sub_1_row_vect_array_multiline(
          Cv1_col1, Cv2_col1,
          A->ml[row_in_A],
          row_in_A_mod,
          dense_array_C_1,
          dense_array_C_2);

      dense_array_C_1[Cp1]  = Cv1_col1;
      dense_array_C_2[Cp1]  = Cv2_col1;
    }
  }

  ml_t *ml = &C->ml[row_idx];
#if DENSE_MULTILINE_C_AFTER_ELIM
  /*  make multiline in C dense */
  free(ml->idx);
  ml->idx   = NULL;
  ml->dense = 1;
  ml->sz    = coldim;
  ml->val   = (re_t *)realloc(ml->val, 2 * coldim * sizeof(re_t));
  memset(ml->val, 0, 2 * coldim * sizeof(re_t));

  copy_dense_arrays_to_zero_dense_multiline(
      dense_array_C_1, dense_array_C_2, 0, ml, coldim, modulus);
#else
  copy_dense_arrays_to_multiline(
      dense_array_C_1, dense_array_C_2, 0, ml, coldim, modulus);
#endif

  free(dense_array_C_1);
  free(dense_array_C_2);

  return 0;
}


/* WAS IN HEADER */

void dense_scal_mul_sub_2_rows_array_array(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const re_m_t Av1_col2,
    const re_m_t Av2_col2,
    const bi_t bwidth,
    const re_l_t *dense_array_source1,
    const re_l_t *dense_array_source2,
    re_l_t *dense_array1,
    re_l_t *dense_array2) {

  /*  check cases where one pair of the elements is zero */
  if (Av1_col1 == 0 && Av2_col1 == 0) {
    dense_scal_mul_sub_1_row_array_array(
        Av1_col2, Av2_col2, bwidth,
        dense_array_source2,
        dense_array1, dense_array2);
    return;
  }
  if (Av1_col2 == 0 && Av2_col2 == 0) {
    dense_scal_mul_sub_1_row_array_array(
        Av1_col1, Av2_col1, bwidth,
        dense_array_source1,
        dense_array1, dense_array2);
    return;
  }

  bi_t i, j=0;
  register re_m_t v1__, v2__;

#ifdef GBLA_USE_AVX
  ELEM_4 av11 = SET1_4(Av1_col1);
  ELEM_4 av12 = SET1_4(Av1_col2);
  ELEM_4 av21 = SET1_4(Av2_col1);
  ELEM_4 av22 = SET1_4(Av2_col2);
  /* for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
     for (j=0; j<__GB_LOOP_UNROLL_BIG; j+=4) {       */
  for (i=0; i<bwidth; i+=4) {

    ELEM_4 ds1 = LOAD_4(dense_array_source1+i+j);
    ELEM_4 ds2 = LOAD_4(dense_array_source2+i+j);

    ELEM_4 da1=FMADD_4(LOAD_4(dense_array1+i+j),av11,ds1);
    STORE_4(dense_array1+i+j,FMADD_4(da1,av12,ds2));

    ELEM_4 da2=FMADD_4(LOAD_4(dense_array2+i+j),av21,ds1);
    STORE_4(dense_array2+i+j,FMADD_4(da2,av22,ds2));
  }
  /* }
     }      */
#else
  for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
    for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
      v1__  = CAST(dense_array_source1[i+j]);
      v2__  = CAST(dense_array_source2[i+j]);

      dense_array1[i+j] +=  v1__ * Av1_col1;
      dense_array1[i+j] +=  v2__ * Av1_col2;

      dense_array2[i+j] +=  v1__ * Av2_col1;
      dense_array2[i+j] +=  v2__ * Av2_col2;
    }
  }
#endif
}

void dense_scal_mul_sub_1_row_array_array(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const bi_t bwidth,
    const re_l_t *dense_array_source,
    re_l_t *dense_array1,
    re_l_t *dense_array2) {

  bi_t i, j=0;
#ifdef GBLA_USE_AVX
  if (Av1_col1 == 0) { /*  only one of them can be zero */
    ELEM_4 Av21 = SET1_4(Av2_col1);
    CHECK_ALIGN(dense_array2);
    /* for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
       for (j=0; j<__GB_LOOP_UNROLL_BIG; j+=4) {       */
    for (i=0; i<bwidth; i+=4) {
      STORE_4(dense_array2+i+j,FMADD_4(LOAD_4(dense_array2+i+j),Av21,LOAD_4(dense_array_source+i+j)));
    }
    /* }
       }    */
  } else {
    if (Av2_col1 == 0) { /*  only second one is zero */
    CHECK_ALIGN(dense_array1);
      ELEM_4 Av11 = SET1_4(Av1_col1);
      /* for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
         for (j=0; j<__GB_LOOP_UNROLL_BIG; j+=4) {       */
      for (i=0; i<bwidth; i+=4) {
        STORE_4(dense_array1+i+j,FMADD_4(LOAD_4(dense_array1+i),Av11,LOAD_4(dense_array_source+i)));
        /* }
           }      */
    }
    } else { /*  both are nonzero */
    CHECK_ALIGN(dense_array1);
    CHECK_ALIGN(dense_array2);
      ELEM_4 Av21 = SET1_4(Av2_col1);
      ELEM_4 Av11 = SET1_4(Av1_col1);
      /* for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
         for (j=0; j<__GB_LOOP_UNROLL_BIG; j+=4) {       */
      for (i=0; i<bwidth; i+=4) {
        STORE_4(dense_array2+i+j,FMADD_4(LOAD_4(dense_array2+i),Av21,LOAD_4(dense_array_source+i)));
        STORE_4(dense_array1+i+j,FMADD_4(LOAD_4(dense_array1+i),Av11,LOAD_4(dense_array_source+i)));
      }
      /* }
         }      */
    }
  }

#else
  assert(bwidth%__GB_LOOP_UNROLL_BIG == 0);
  register re_m_t v__;

  if (Av1_col1 == 0) { /*  only one of them can be zero */
    for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
      for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
        v__ = CAST(dense_array_source[i+j]) ;
        dense_array2[i+j] +=  v__ * Av2_col1;
      }
    }
  } else {
    if (Av2_col1 == 0) { /*  only second one is zero */
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
        for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
          v__ = CAST(dense_array_source[i+j]);
          dense_array1[i+j] +=  v__ * Av1_col1;
        }
      }
    } else { /*  both are nonzero */
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
        for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
          v__ = CAST(dense_array_source[i+j]);
          dense_array1[i+j] +=  v__ * Av1_col1;
          dense_array2[i+j] +=  v__ * Av2_col1;
        }
      }
    }
  }
#endif
}

/** TODO check align of dense_array[i0] */
void red_dense_array_modular(re_l_t *dense_array, const bi_t bwidth, const mod_t modulus) {
#ifdef GBLA_USE_AVX
  ELEM_4 P = SET1_4(modulus);
  ELEM_4 INVP = SET1_4(1./(double)modulus);
  ELEM_4 NEGP = SET1_4(-(double)modulus);
  ELEM_4 MIN = ZERO_4();
  ELEM_4 MAX = SET1_4(modulus-1);
  ELEM_4 Q,T ;

  bi_t i;
  CHECK_ALIGN(dense_array);
  for (i=0; i<bwidth/4*4; i+=4) {
    ELEM_4 C = LOAD_4(dense_array+i);
    FLOAT_MOD_4(C, P, INVP, Q);
    NORML_MOD_4(C, P, NEGP, MIN, MAX, Q, T);
    STORE_4(dense_array+i, C);
  }
  for (; i<bwidth; ++i) {
    dense_array[i]  = (re_l_t)(MODP(dense_array[i], modulus));
  }
#else
  bi_t i;
  for (i=0; i<bwidth; ++i) {
    dense_array[i]  = (re_l_t)(MODP(dense_array[i], modulus));
  }
#endif
}

/** TODO check align of dense_array[i0] */
int normalize_dense_array(re_l_t *dense_array,
    const ci_t coldim, const mod_t modulus) {
  re_t val;
  ci_t i;
  int h1  = get_head_dense_array(dense_array, &val, coldim, modulus);

  if (h1 == -1)
    return h1;

  inverse_val(&val, modulus);

#ifdef GBLA_USE_AVX_NO

  ELEM_4 P = SET1_4(modulus);
  ELEM_4 INVP = SET1_4(1./(double)modulus);
  ELEM_4 NEGP = SET1_4(-(double)modulus);
  ELEM_4 MIN = ZERO_4();
  ELEM_4 MAX = SET1_4(modulus-1);
  ELEM_4 Q,T ;
  ELEM_4 VAL = SET1_4(val);

  i = h1 ;
  while ( i < coldim && NOTALIGNED(dense_array+i) ) {
    dense_array[i]  =  MODP(dense_array[i]*val,modulus);
    ++i;
  }
  CHECK_ALIGN(dense_array+i);
  for (; i<(coldim)/4*4; i+=4) {
    ELEM_4 C = MUL_4(LOAD_4(dense_array+i),VAL);
    FLOAT_MOD_4(C, P, INVP, Q);
    NORML_MOD_4(C, P, NEGP, MIN, MAX, Q, T);
    STORE_4(dense_array+i, C);
  }


 for (; i<coldim; ++i) {
    dense_array[i]  =  MODP(dense_array[i]*val,modulus);
  }

#else
  for (i=h1; i<coldim; ++i) {
    dense_array[i]  *=  val;
    dense_array[i]  =  MODP(dense_array[i],modulus);
  }
#endif
  return h1;
}

int get_head_dense_array(re_l_t *dense_array,
    re_t *val, const ci_t coldim, const mod_t modulus) {
  ci_t i;
  for (i=0; i<coldim; ++i) {
    *val  = MODP(dense_array[i], modulus);
    if (*val != 0)
      return (int)i;
    else
      dense_array[i]  = 0;
  }
  return -1;
}

void sparse_scal_mul_sub_2_rows_vect_array(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const re_m_t Av1_col2,
    const re_m_t Av2_col2,
    const mbl_t multiline,
    re_l_t *dense_val1,
    re_l_t *dense_val2) {

  /*  check cases where one pair of the elements is zero */
  if (Av1_col1 == 0 && Av2_col1 == 0) {
    sparse_scal_mul_sub_1_row_vect_array(
        Av1_col2, Av2_col2, multiline, 1, dense_val1, dense_val2);
    return;
  }
  if (Av1_col2 == 0 && Av2_col2 == 0) {
    sparse_scal_mul_sub_1_row_vect_array(
        Av1_col1, Av2_col1, multiline, 0, dense_val1, dense_val2);
    return;
  }

  const uint32_t N    = multiline.sz;
  const bi_t *p_idx   = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val   = (N != 0) ? multiline.val : NULL;
  const re_t *p_val2  = p_val + 1;
  uint32_t i;
  bi_t j;

  register re_m_t v1__, v2__;
  register uint32_t idx;

#ifdef GBLA_USE_AVX_XXX
  for (i=0; i<__GB_ROUND_DOWN(N, __GB_LOOP_UNROLL_SMALL) ; i+=__GB_LOOP_UNROLL_SMALL) {
    SET2_2(Av1_d,Av1_col1,Av1_col2);
    SET2_2(Av2_d,Av2_col1,Av2_col2);
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; j+=2) {
      SET2_2(dv1,*(dense_val1+p_idx[i+j]),*(dense_val1+p_idx[i+j+1]));
      SET2_2(dv2,*(dense_val2+p_idx[i+j]),*(dense_val2+p_idx[i+j+1]));
      ELEM_2 v__ = LOAD_2(p_val+2*(i+j));
      AXPY_2_(dv1,v__,Av1_d);
      AXPY_2_(dv2,v__,Av2_d);
      STORE_L_2(dense_val1+p_idx[i+j],dv1);
      STORE_H_2(dense_val1+p_idx[i+j+1],dv1);
      STORE_L_2(dense_val2+p_idx[i+j],dv2);
      STORE_H_2(dense_val2+p_idx[i+j+1],dv2);

    }
  }
  for ( ; i<N; ++i) {
    idx   = p_idx[i];
    v1__  = p_val[2*i];
    v2__  = p_val2[2*i];

    dense_val1[idx] +=  Av1_col1 * v1__;
    dense_val1[idx] +=  Av1_col2 * v2__;

    dense_val2[idx] +=  Av2_col1 * v1__;
    dense_val2[idx] +=  Av2_col2 * v2__;
  }
#else
  for (i=0; i<__GB_ROUND_DOWN(N, __GB_LOOP_UNROLL_SMALL) ; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      idx   = p_idx[i+j];
      v1__  = p_val[2*(i+j)];
      v2__  = p_val2[2*(i+j)];

      dense_val1[idx] +=  Av1_col1 * v1__;
      dense_val1[idx] +=  Av1_col2 * v2__;

      dense_val2[idx] +=  Av2_col1 * v1__;
      dense_val2[idx] +=  Av2_col2 * v2__;
    }
  }
  for ( ; i<N; ++i) {
    idx   = p_idx[i];
    v1__  = p_val[2*i];
    v2__  = p_val2[2*i];

    dense_val1[idx] +=  Av1_col1 * v1__;
    dense_val1[idx] +=  Av1_col2 * v2__;

    dense_val2[idx] +=  Av2_col1 * v1__;
    dense_val2[idx] +=  Av2_col2 * v2__;
  }
#endif
}

void dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const re_m_t Av1_col2,
    const re_m_t Av2_col2,
    const ml_t multiline,
    re_l_t *dense_val1,
    re_l_t *dense_val2,
    const ci_t offset1,
    const ci_t offset2) {

  /*  if (offset1 < 0) return; */

  /*  check cases where one pair of the elements is zero */
  if (Av1_col1 == 0 && Av2_col1 == 0) {
    dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
        Av1_col2, Av2_col2, multiline, 1, dense_val1, dense_val2, offset2);
    return;
  }
  if (Av1_col2 == 0 && Av2_col2 == 0) {
    dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
        Av1_col1, Av2_col1, multiline, 0, dense_val1, dense_val2, offset1);
    return;
  }

  const re_t *p_val = multiline.val;
  uint32_t i;

#ifdef GBLA_USE_AVX_XXX
  const uint32_t outer_loop = multiline.sz - __GB_LOOP_UNROLL_SMALL;
  bi_t j;

  register re_m_t v1__, v2__;
  /*  register uint32_t idx; */

  for (i=offset1; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      v1__ = p_val[2*(i+j)];
      v2__ = p_val[2*(i+j)+1];

      dense_val1[i+j] +=  v1__ * Av1_col1;
      dense_val1[i+j] +=  v2__ * Av1_col2;


      dense_val2[i+j] +=  v1__ * Av2_col1;
      dense_val2[i+j] +=  v2__ * Av2_col2;
#if COUNT_REDS
      nreductions +=  4;
#endif
    }
  }
  for (; i<multiline.sz; ++i) {
    v1__ = p_val[2*i];
    v2__ = p_val[2*i+1];

    dense_val1[i] +=  v1__ * Av1_col1;
    dense_val1[i] +=  v2__ * Av1_col2;

    dense_val2[i] +=  v1__ * Av2_col1;
    dense_val2[i] +=  v2__ * Av2_col2;
#if COUNT_REDS
      nreductions +=  4;
#endif
  }
}

#else
const uint32_t outer_loop = multiline.sz - __GB_LOOP_UNROLL_SMALL;
bi_t j;

register re_m_t v1__, v2__;
/*  register uint32_t idx; */

for (i=offset1; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
  for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
    v1__ = p_val[2*(i+j)];
    v2__ = p_val[2*(i+j)+1];

    dense_val1[i+j] +=  v1__ * Av1_col1;
    dense_val1[i+j] +=  v2__ * Av1_col2;


    dense_val2[i+j] +=  v1__ * Av2_col1;
    dense_val2[i+j] +=  v2__ * Av2_col2;
#if COUNT_REDS
      nreductions +=  4;
#endif
  }
}
for (; i<multiline.sz; ++i) {
  v1__ = p_val[2*i];
  v2__ = p_val[2*i+1];

  dense_val1[i] +=  v1__ * Av1_col1;
  dense_val1[i] +=  v2__ * Av1_col2;

  dense_val2[i] +=  v1__ * Av2_col1;
  dense_val2[i] +=  v2__ * Av2_col2;
#if COUNT_REDS
      nreductions +=  4;
#endif
}
}
#endif

void sparse_scal_mul_sub_1_row_vect_array_multiline(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const ml_t multiline,
    const bi_t line_idx,
    re_l_t *dense_val1,
    re_l_t *dense_val2) {

  const uint32_t N    = multiline.sz;
  const mli_t *p_idx  = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val   = (N != 0) ? multiline.val : NULL;
  p_val +=  line_idx;
  uint32_t i  = 0;

  register re_m_t v__;
  register uint32_t idx;

#ifdef GBLA_USE_AVX_XXX
  /* printf("a11 %d -- a21 %d\n",Av1_col1,Av2_col1); */
  /*  both cannot be zero at the same time */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (; i<N; ++i) {
      idx = p_idx[i];
      v__ = p_val[2*i];

      dense_val1[idx] +=  v__ * Av1_col1;
      dense_val2[idx] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  2;
#endif
    }
  } else { /*  one of them is zero */
    if (Av1_col1 != 0) {
      /* printf("i %d -- N %d\n",i,N); */
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val1[idx] +=  v__ * Av1_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      }
    } else {
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val2[idx] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      }
    }
  }
#else
  /* printf("a11 %d -- a21 %d\n",Av1_col1,Av2_col1); */
  /*  both cannot be zero at the same time */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (; i<N; ++i) {
      idx = p_idx[i];
      v__ = p_val[2*i];

      dense_val1[idx] +=  v__ * Av1_col1;
      dense_val2[idx] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  2;
#endif
    }
  } else { /*  one of them is zero */
    if (Av1_col1 != 0) {
      /* printf("i %d -- N %d\n",i,N); */
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val1[idx] +=  v__ * Av1_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      }
    } else {
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val2[idx] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      }
    }
  }
#endif
}

void dense_scal_mul_sub_2_rows_vect_array(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const re_m_t Av1_col2,
    const re_m_t Av2_col2,
    const mbl_t multiline,
    const bi_t  bwidth,
    re_l_t *dense_val1,
    re_l_t *dense_val2) {

  /*  check cases where one pair of the elements is zero */
  if (Av1_col1 == 0 && Av2_col1 == 0) {
    dense_scal_mul_sub_1_row_vect_array(
        Av1_col2, Av2_col2, multiline, 1, bwidth, dense_val1, dense_val2);
    return;
  }
  if (Av1_col2 == 0 && Av2_col2 == 0) {
    dense_scal_mul_sub_1_row_vect_array(
        Av1_col1, Av2_col1, multiline, 0, bwidth, dense_val1, dense_val2);
    return;
  }

  const re_t *p_val = multiline.val;
  uint32_t i;
  bi_t j;

  register re_m_t v1__, v2__;
  /*  register uint32_t idx; */

#ifdef GBLA_USE_AVX_XXX
  for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      /* printf(";;%d..%d::",i,j); */
      v1__ = p_val[2*(i+j)];
      v2__ = p_val[2*(i+j)+1];
      /* printf("v11 %d v21 %d !! ",v1__,v2__); */

      dense_val1[i+j] +=  v1__ * Av1_col1;
      dense_val1[i+j] +=  v2__ * Av1_col2;

      v1__  *=  Av2_col1;
      v2__  *=  Av2_col2;
      /* printf("v12 %d v22 %d !! ",v1__,v2__); */

      dense_val2[i+j] +=  v1__;
      dense_val2[i+j] +=  v2__;
    }
  }
#else
  for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      /* printf(";;%d..%d::",i,j); */
      v1__ = p_val[2*(i+j)];
      v2__ = p_val[2*(i+j)+1];
      /* printf("v11 %d v21 %d !! ",v1__,v2__); */

      dense_val1[i+j] +=  v1__ * Av1_col1;
      dense_val1[i+j] +=  v2__ * Av1_col2;

      v1__  *=  Av2_col1;
      v2__  *=  Av2_col2;
      /* printf("v12 %d v22 %d !! ",v1__,v2__); */

      dense_val2[i+j] +=  v1__;
      dense_val2[i+j] +=  v2__;
    }
  }
#endif
  /* printf("\n"); */
}

void sparse_scal_mul_sub_1_row_vect_array(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const mbl_t multiline,
    const bi_t line_idx,
    re_l_t *dense_val1,
    re_l_t *dense_val2) {

  const uint32_t N  = multiline.sz;
  /* printf("mlsz %d\n",multiline.sz); */
  const bi_t *p_idx = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val = (N != 0) ? multiline.val : NULL;
  p_val +=  line_idx;
  uint32_t i  = 0;

  register re_m_t v__;
  register uint32_t idx;

#ifdef GBLA_USE_AVX_XXX
  /* printf("a11 %d -- a21 %d\n",Av1_col1,Av2_col1); */
  /*  both cannot be zero at the same time */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (; i<N; ++i) {
      idx = p_idx[i];
      v__ = p_val[2*i];

      dense_val1[idx] +=  v__ * Av1_col1;
      dense_val2[idx] +=  v__ * Av2_col1;
    }
  } else { /*  one of them is zero */
    if (Av1_col1 != 0) {
      /* printf("i %d -- N %d\n",i,N); */
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val1[idx] +=  v__ * Av1_col1;
      }
    } else {
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val2[idx] +=  v__ * Av2_col1;
      }
    }
  }
#else
  /* printf("a11 %d -- a21 %d\n",Av1_col1,Av2_col1); */
  /*  both cannot be zero at the same time */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (; i<N; ++i) {
      idx = p_idx[i];
      v__ = p_val[2*i];

      dense_val1[idx] +=  v__ * Av1_col1;
      dense_val2[idx] +=  v__ * Av2_col1;
    }
  } else { /*  one of them is zero */
    if (Av1_col1 != 0) {
      /* printf("i %d -- N %d\n",i,N); */
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val1[idx] +=  v__ * Av1_col1;
      }
    } else {
      for (; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val2[idx] +=  v__ * Av2_col1;
      }
    }
  }
#endif
}

void dense_scal_mul_sub_1_row_vect_array(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const mbl_t multiline,
    const bi_t line_idx,
    const bi_t  bwidth,
    re_l_t *dense_val1,
    re_l_t *dense_val2) {

  const re_t *p_val = multiline.val;
  p_val +=  line_idx;
  uint32_t i;
  bi_t j;

  register re_m_t v__;

#ifdef GBLA_USE_AVX_XXX
  /*  both cannot be zero at the same time */
  /* printf("mlsize %d\n",multiline.sz); */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
      for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
        v__ = p_val[2*(i+j)];

        dense_val1[i+j] +=  v__ * Av1_col1;
        dense_val2[i+j] +=  v__ * Av2_col1;
      }
    }
  } else { /*  one of them is zero */
    if (Av1_col1 == 0) {
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
        for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
          v__ = p_val[2*(i+j)];

          dense_val2[i+j] +=  v__ * Av2_col1;
        }
      }
    } else {
      if (Av2_col1 == 0) {
        for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
          for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
            v__ = p_val[2*(i+j)];

            dense_val1[i+j] +=  v__ * Av1_col1;
          }
        }
      }
    }
  }
#else
  /*  both cannot be zero at the same time */
  /* printf("mlsize %d\n",multiline.sz); */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
      for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
        v__ = p_val[2*(i+j)];

        dense_val1[i+j] +=  v__ * Av1_col1;
        dense_val2[i+j] +=  v__ * Av2_col1;
      }
    }
  } else { /*  one of them is zero */
    if (Av1_col1 == 0) {
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
        for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
          v__ = p_val[2*(i+j)];

          dense_val2[i+j] +=  v__ * Av2_col1;
        }
      }
    } else {
      if (Av2_col1 == 0) {
        for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
          for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
            v__ = p_val[2*(i+j)];

            dense_val1[i+j] +=  v__ * Av1_col1;
          }
        }
      }
    }
  }
#endif
}

void normalize_multiline(ml_t *m, const ci_t coldim, const mod_t modulus) {
  mli_t idx;
  re_t h1 = 0, h2 = 0;

  if (m->sz == 0)
    return;

  get_head_multiline_hybrid(m, 0, &h1, &idx, coldim);
  get_head_multiline_hybrid(m, 1, &h2, &idx, coldim);

  /*  invert values modulo modulus */
  inverse_val(&h1, modulus);
  inverse_val(&h2, modulus);


  /*  skip if both are 1 and/or 0 */
  if ((h1 == 0 || h1 == 1) && (h2 == 0 || h2 == 1))
    return;

  re_l_t tmp_val;
  /*  normalize h2 */
  if (h1 == 0 || h1 == 1) {
    for (idx=0; idx<m->sz; ++idx) {
      tmp_val         = (re_l_t)m->val[2*idx+1] * h2;
      m->val[2*idx+1] = MODP(tmp_val, modulus);
    }
  } else {
    /*  normalize h1 */
    if (h2 == 0 || h2 == 1) {
      for (idx=0; idx<m->sz; ++idx) {
        tmp_val       = (re_l_t)m->val[2*idx] * h1;
        m->val[2*idx] = MODP(tmp_val,  modulus);
      }
      /*  normalize h1 and h2 */
    } else {
#ifdef GBLA_USE_AVX
      ELEM_4 P = SET1_4(modulus);
      ELEM_4 INVP = SET1_4(1./(double)modulus);
      ELEM_4 NEGP = SET1_4(-(double)modulus);
      ELEM_4 MIN = ZERO_4();
      ELEM_4 MAX = SET1_4(modulus-1);
      ELEM_4 Q,T ;
      ELEM_4 H12= SET_4(h1,h2,h1,h2);

      idx = 0 ;
      while (idx < 2*m->sz &&  NOTALIGNED(m->val+idx)  ) {
        m->val[idx  ] = MODP(m->val[idx  ] * h1, modulus);
        m->val[idx+1] = MODP(m->val[idx+1] * h2, modulus);
        idx+=2;
      }
      CHECK_ALIGN(m->val+idx);
      for (; idx<2*(m->sz/4*4); idx+=4) {
        ELEM_4 C = MUL_4(LOAD_4(m->val+idx),H12);
        FLOAT_MOD_4(C, P, INVP, Q);
        NORML_MOD_4(C, P, NEGP, MIN, MAX, Q, T);
        STORE_4(m->val+idx, C);
      }

      for (; idx<2*m->sz; idx+=2) {
        m->val[idx  ] = MODP(m->val[idx  ] * h1, modulus);
        m->val[idx+1] = MODP(m->val[idx+1] * h2, modulus);
      }
#else
      for (idx=0; idx<m->sz; ++idx) {
        tmp_val         = (re_l_t)m->val[2*idx] * h1;
        m->val[2*idx]   = MODP(tmp_val, modulus);
        tmp_val         = (re_l_t)m->val[2*idx+1] * h2;
        m->val[2*idx+1] = MODP(tmp_val, modulus);
      }

#endif
    }
  }
}

void copy_multiline_to_dense_array(const ml_t m, re_l_t *dense_1,
    re_l_t *dense_2, const ci_t coldim) {
  if (m.sz == 0)
    return;

  register mli_t idx;
  ci_t i;

  if (m.sz < coldim) {
    for (i=0; i<m.sz; ++i) {
      idx           = m.idx[i];
      dense_1[idx]  = (re_l_t)m.val[2*i];
      dense_2[idx]  = (re_l_t)m.val[2*i+1];

    }
  } else {
    for (i=0; i<coldim; ++i) {
      dense_1[i]  = (re_l_t)m.val[2*i];
      dense_2[i]  = (re_l_t)m.val[2*i+1];
    }
  }
}

void sparse_scal_mul_sub_2_rows_vect_array_multiline(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const re_m_t Av1_col2,
    const re_m_t Av2_col2,
    const ml_t multiline,
    re_l_t *dense_val1,
    re_l_t *dense_val2) {

  /*  check cases where one pair of the elements is zero */
  if (Av1_col1 == 0 && Av2_col1 == 0) {
    sparse_scal_mul_sub_1_row_vect_array_multiline(
        Av1_col2, Av2_col2, multiline, 1, dense_val1, dense_val2);
    return;
  }
  if (Av1_col2 == 0 && Av2_col2 == 0) {
    sparse_scal_mul_sub_1_row_vect_array_multiline(
        Av1_col1, Av2_col1, multiline, 0, dense_val1, dense_val2);
    return;
  }

  const uint32_t N    = multiline.sz;
  const mli_t *p_idx  = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val   = (N != 0) ? multiline.val : NULL;
  const re_t *p_val2  = p_val + 1;
  uint32_t i;
  bi_t j;

  register re_m_t v1__, v2__;
  register uint32_t idx;

#ifdef GBLA_USE_AVX_XXX
  for (i=0; i<__GB_ROUND_DOWN(N, __GB_LOOP_UNROLL_SMALL) ; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      idx   = p_idx[i+j];
      v1__  = p_val[2*(i+j)];
      v2__  = p_val2[2*(i+j)];

      dense_val1[idx] +=  Av1_col1 * v1__;
      dense_val1[idx] +=  Av1_col2 * v2__;

      dense_val2[idx] +=  Av2_col1 * v1__;
      dense_val2[idx] +=  Av2_col2 * v2__;
#if COUNT_REDS
      nreductions +=  4;
#endif
    }
  }
  for ( ; i<N; ++i) {
    idx   = p_idx[i];
    v1__  = p_val[2*i];
    v2__  = p_val2[2*i];

    dense_val1[idx] +=  Av1_col1 * v1__;
    dense_val1[idx] +=  Av1_col2 * v2__;

    dense_val2[idx] +=  Av2_col1 * v1__;
    dense_val2[idx] +=  Av2_col2 * v2__;
#if COUNT_REDS
      nreductions +=  4;
#endif
  }
#else
  for (i=0; i<__GB_ROUND_DOWN(N, __GB_LOOP_UNROLL_SMALL) ; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      idx   = p_idx[i+j];
      v1__  = p_val[2*(i+j)];
      v2__  = p_val2[2*(i+j)];

      dense_val1[idx] +=  Av1_col1 * v1__;
      dense_val1[idx] +=  Av1_col2 * v2__;

      dense_val2[idx] +=  Av2_col1 * v1__;
      dense_val2[idx] +=  Av2_col2 * v2__;
#if COUNT_REDS
      nreductions +=  4;
#endif
    }
  }
  for ( ; i<N; ++i) {
    idx   = p_idx[i];
    v1__  = p_val[2*i];
    v2__  = p_val2[2*i];

    dense_val1[idx] +=  Av1_col1 * v1__;
    dense_val1[idx] +=  Av1_col2 * v2__;

    dense_val2[idx] +=  Av2_col1 * v1__;
    dense_val2[idx] +=  Av2_col2 * v2__;
#if COUNT_REDS
      nreductions +=  4;
#endif
  }
#endif
}

/* void check_inverse(re_t x, re_t x1, mod_t modulus)
{
  int64_t p = (int64_t)x * (int64_t)x1 ;
  assert (p % (int64_t) modulus == 1);
}                                                     */


void inverse_val(re_t *x, const mod_t modulus) {
  assert(*x);
  if ( *x == 1 ) return ;
  assert((int32_t)modulus > 0);
  int32_t u1 = 1, u2 = 0;
  int32_t v1 = 0, v3 = (int32_t)modulus;
  int32_t u3 = (int32_t)*x, v2 = 1;
  while (v3 != 0) {
    int32_t q  = u3 / v3;
    int32_t t1 = u1 - v1 * q;
    u1  = v1; v1  = t1;

    int32_t t3 = u3 - v3 * q;
    u3  = v3; v3  = t3;

    int32_t t2 = u2 - v2 * q;
    u2  = v2; v2  = t2;
  }
  if (u1 < 0) {
    u1  +=  modulus;
    /* check_inverse(*x,u1,modulus); */
    *x  =   (re_t)u1;
    return;
  }
  if (u1 > (int32_t)modulus) {
    u1  -=  modulus;
    /* check_inverse(*x,u1,modulus); */
    *x  =   (re_t) u1;
    return;
  }
  /* check_inverse(*x,u1,modulus); */
  *x  = (re_t)u1;
  return;
}

int get_smallest_waiting_row(wl_t *waiting_global,
    ri_t *wl_idx, ri_t *wl_lp) {
  if (waiting_global->sz == 0)
    return 0;

#if DEBUG_ECHELONIZE
  int tid = omp_get_thread_num();
  printf("BEFORE SORT\n");
  for (int ii=0; ii<waiting_global->sz; ++ii) {
    printf("T(%d) %d . %d\n",tid,waiting_global->list[ii].idx,waiting_global->list[ii].lp);
  }
#endif
  /*  sort the waiting list */
  qsort(waiting_global->list, waiting_global->sz, sizeof(wle_t), cmp_wle);
#if DEBUG_ECHELONIZE
  printf("AFTER SORT\n");
  for (int ii=0; ii<waiting_global->sz; ++ii) {
    if (ii<waiting_global->sz - 1) {
      if (waiting_global->list[ii].idx == waiting_global->list[ii+1].idx) {
        printf("SAME IN WAITING!\n");
      }
    }
    printf("ST(%d) %d/%d %d . %d\n",tid,ii,waiting_global->sz,waiting_global->list[ii].idx,waiting_global->list[ii].lp);
  }
#endif
  /*  store last (smallest index) element separately and */
  /*  remove it from the list */
  *wl_idx = waiting_global->list[waiting_global->sz-1].idx;
  *wl_lp  = waiting_global->list[waiting_global->sz-1].lp;
#if DEBUG_ECHELONIZE
  printf("wsidx %d -- wslp %d\n",*wl_idx, *wl_lp);
#endif
  waiting_global->sz--;

  return 1;
}

void copy_dense_arrays_to_multiline(const re_l_t *dense_1,
    const re_l_t *dense_2, const int start_pos, ml_t *m, const ci_t coldim,
    const mod_t modulus) {

#if DEBUG_ECHELONIZE
  printf("in copy mlv %p\n",m->val);
#endif

  ci_t buffer = m->sz;
  m->dense    = 0;
  m->sz       = 0;
  ci_t i;
  re_l_t tmp_1, tmp_2;
  for (i=start_pos; i<coldim; ++i) {
    tmp_1 = MODP( dense_1[i], modulus);
    tmp_2 = MODP( dense_2[i], modulus);
    /* printf("t1 %lu -- t2 %lu\n",tmp_1,tmp_2); */
    if (tmp_1 != 0 || tmp_2 != 0) {
      if (m->sz >= buffer) { /*  realloc memory for multiline ml */
        m->idx  =   (mli_t*)realloc(m->idx, 3*buffer*sizeof(mli_t));
        m->val  =   (re_t *)realloc(m->val, 6*buffer*sizeof(re_t));
        buffer  *=  3;
      }
      /* printf("B %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]); */
      m->idx[m->sz]     = i;
      m->val[2*m->sz]   = (re_t)tmp_1;
      /* printf("A1 %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]); */
      m->val[2*m->sz+1] = (re_t)tmp_2;
      m->sz++;
    }
    /* printf("A %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]); */
  }

  /*  fix sizes: we do not touch this anymore, so keep it as small as possible */
  m->idx  = (mli_t*)realloc(m->idx, m->sz * sizeof(mli_t));
  m->val  = (re_t *)realloc(m->val, 2 * m->sz * sizeof(re_t));
}

void copy_dense_array_to_zero_dense_multiline(const re_l_t *dense_1,
    const int start_pos, ml_t *m, const ci_t coldim, const mod_t modulus) {

#if DEBUG_ECHELONIZE
  printf("in copy mlv %p\n",m->val);
#endif
  ci_t i;
  re_l_t tmp_1;
  for (i=start_pos; i<coldim; ++i) {
    tmp_1 = MODP(dense_1[i], modulus);
    if (tmp_1 != 0) {
      m->val[2*i]   = (re_t)tmp_1;
    }
  }
}

void dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
    const re_m_t Av1_col1,
    const re_m_t Av2_col1,
    const ml_t multiline,
    const bi_t line_idx,
    re_l_t *dense_val1,
    re_l_t *dense_val2,
    const ci_t offset) {

  const re_t *p_val = multiline.val;
  p_val +=  line_idx;
  const uint32_t outer_loop = multiline.sz - __GB_LOOP_UNROLL_SMALL;
  uint32_t i;
  bi_t j;

  register re_m_t v__;

#ifdef GBLA_USE_AVX
  /*  both cannot be zero at the same time */
  /* printf("mlsize %d\n",multiline.sz); */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    j=0;
    ELEM_4 Av11 = SET1_4(Av1_col1);
    ELEM_4 Av21 = SET1_4(Av2_col1);
    /* for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
       for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {                */
    i = offset ;
    while ( i < multiline.sz && NOTALIGNED(dense_val1+i)) {
      v__ = p_val[2*(i+j)];

      dense_val1[i+j] +=  v__ * Av1_col1;
      dense_val2[i+j] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  2;
#endif
      ++i ;
    }
    CHECK_ALIGN(dense_val1);
    CHECK_ALIGN(dense_val2);
    for ( ; i<multiline.sz/4*4; i+=4) {
      ELEM_4 V = SET_4(p_val[2*(i+j)],p_val[2*(i+j)+2],p_val[2*(i+j)+4],p_val[2*(i+j)+6]);
      STORE_4(dense_val1+i+j,FMADD_4(LOAD_4(dense_val1+i+j),V,Av11));
      STORE_4(dense_val2+i+j,FMADD_4(LOAD_4(dense_val2+i+j),V,Av21));
    }
    for ( ; i < multiline.sz; ++i) {
      v__ = p_val[2*(i+j)];

      dense_val1[i+j] +=  v__ * Av1_col1;
      dense_val2[i+j] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  2;
#endif
    }

    /* }
       }      */
  } else { /*  one of them is zero */
    if (Av1_col1 == 0) {
    j=0;
    /* ELEM_4 Av11 = SET1_4(Av1_col1); */
    ELEM_4 Av21 = SET1_4(Av2_col1);
    /* for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
       for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {                */
    i = offset ;
    while ( i < multiline.sz && NOTALIGNED(dense_val2+i) ) {
      v__ = p_val[2*(i+j)];

      /* dense_val1[i+j] +=  v__ * Av1_col1; */
      dense_val2[i+j] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      ++i ;
    }
    CHECK_ALIGN(dense_val2);
    for ( ; i<multiline.sz/4*4; i+=4) {
      ELEM_4 V = SET_4(p_val[2*(i+j)],p_val[2*(i+j)+2],p_val[2*(i+j)+4],p_val[2*(i+j)+6]);
      /* STORE_4(dense_val1+i+j,FMADD_4(LOAD_4(dense_val1+i+j),V,Av11)); */
      STORE_4(dense_val2+i+j,FMADD_4(LOAD_4(dense_val2+i+j),V,Av21));
    }
    for ( ; i < multiline.sz; ++i) {
      v__ = p_val[2*(i+j)];

      /* dense_val1[i+j] +=  v__ * Av1_col1; */
      dense_val2[i+j] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
    }

    /* }
       }      */
    } else {
      if (Av2_col1 == 0) {
    j=0;
    ELEM_4 Av11 = SET1_4(Av1_col1);
    /* ELEM_4 Av21 = SET1_4(Av2_col1); */
    /* for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
       for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {                */
    i = offset ;
    while ( i < multiline.sz &&  NOTALIGNED(dense_val1+i)) {
      v__ = p_val[2*(i+j)];

      dense_val1[i+j] +=  v__ * Av1_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      /* dense_val2[i+j] +=  v__ * Av2_col1; */
      ++i ;
    }
    CHECK_ALIGN(dense_val1);
    for ( ; i<multiline.sz/4*4; i+=4) {
      ELEM_4 V = SET_4(p_val[2*(i+j)],p_val[2*(i+j)+2],p_val[2*(i+j)+4],p_val[2*(i+j)+6]);
      STORE_4(dense_val1+i+j,FMADD_4(LOAD_4(dense_val1+i+j),V,Av11));
      /* STORE_4(dense_val2+i+j,FMADD_4(LOAD_4(dense_val2+i+j),V,Av21)); */
    }
    for ( ; i < multiline.sz; ++i) {
      v__ = p_val[2*(i+j)];

      dense_val1[i+j] +=  v__ * Av1_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      /* dense_val2[i+j] +=  v__ * Av2_col1; */
    }

    /* }
       }      */
      }
    }
  }
#else
  /*  both cannot be zero at the same time */
  /* printf("mlsize %d\n",multiline.sz); */
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
      for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
        v__ = p_val[2*(i+j)];

        dense_val1[i+j] +=  v__ * Av1_col1;
        dense_val2[i+j] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  2;
#endif
      }
    }
    for (; i<multiline.sz; ++i) {
      v__ = p_val[2*i];

      dense_val1[i] +=  v__ * Av1_col1;
      dense_val2[i] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  2;
#endif
    }
  } else { /*  one of them is zero */
    if (Av1_col1 == 0) {
      for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
        for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
          v__ = p_val[2*(i+j)];

          dense_val2[i+j] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
        }
      }
      for (; i<multiline.sz; ++i) {
        v__ = p_val[2*i];

        /* dense_val1[i] +=  v__ * Av1_col1; */
        dense_val2[i] +=  v__ * Av2_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
      }
    } else {
      if (Av2_col1 == 0) {
        for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
          for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
            v__ = p_val[2*(i+j)];

            dense_val1[i+j] +=  v__ * Av1_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
          }
        }
        for (; i<multiline.sz; ++i) {
          v__ = p_val[2*i];

          dense_val1[i] +=  v__ * Av1_col1;
#if COUNT_REDS
      nreductions +=  1;
#endif
          /* dense_val2[i] +=  v__ * Av2_col1; */
        }
      }
    }
  }

#endif
}

void push_row_to_waiting_list(wl_t *waiting_global, const ri_t row_idx,
    const ri_t last_piv_reduced_by) {
#if DEBUG_ECHELONIZE
  int tid = omp_get_thread_num();
  printf("BEFORE PUSH\n");
  for (int ii=0; ii<waiting_global->sz; ++ii) {
    printf("T(%d) %d . %d\n",tid,waiting_global->list[ii].idx,waiting_global->list[ii].lp);
  }
#endif
  waiting_global->list[waiting_global->sz].idx  = row_idx;
  waiting_global->list[waiting_global->sz].lp   = last_piv_reduced_by;
#if DEBUG_ECHELONIZE
  printf("rowidx %u -- lprb %u\n", row_idx, last_piv_reduced_by);
#endif
  waiting_global->sz++;
#if DEBUG_ECHELONIZE
  printf("AFTER PUSH\n");
  for (int ii=0; ii<waiting_global->sz; ++ii) {
    printf("T(%d) %d . %d\n",tid,waiting_global->list[ii].idx,waiting_global->list[ii].lp);
  }
#endif
}

void copy_dense_arrays_to_zero_dense_multiline(const re_l_t *dense_1,
    const re_l_t *dense_2, const int start_pos, ml_t *m, const ci_t coldim,
    const mod_t modulus) {

#if DEBUG_ECHELONIZE
  printf("in copy mlv %p\n",m->val);
#endif
  ci_t i;
  re_l_t tmp_1, tmp_2;
  for (i=start_pos; i<coldim; ++i) {
    tmp_1 = MODP(dense_1[i], modulus);
    tmp_2 = MODP(dense_2[i], modulus);
    /* printf("t1 %lu -- t2 %lu\n",tmp_1,tmp_2); */
    if (tmp_1 != 0 || tmp_2 != 0) {
      /* printf("B %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]); */
      m->val[2*i]   = (re_t)tmp_1;
      /* printf("A1 %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]); */
      m->val[2*i+1] = (re_t)tmp_2;
    }
    /* printf("A %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]); */
  }
}

int cmp_wle(const void *a, const void *b) {
  return ((wle_t *)(b))->idx - ((wle_t *)(a))->idx;
}

void copy_sparse_to_dense_block(
    mbl_t *sparse_block, re_l_t **dense_block, int bheight, int bwidth) {

  bi_t i,j;
  for (i=0; i<bheight/2; ++i) {
    mli_t idx;
    /*  re_t val1, val2; */

    if (sparse_block[i].sz == 0)
      continue;

    if (sparse_block[i].dense == 0) {
      for (j=0; j<sparse_block[i].sz; ++j) {
        idx   = sparse_block[i].idx[j];
        dense_block[2*i  ][idx] = sparse_block[i].val[2*j  ];
        dense_block[2*i+1][idx] = sparse_block[i].val[2*j+1];
      }
    } else{
      for (j=0; j<bwidth; ++j) {
        dense_block[2*i  ][j] = sparse_block[i].val[2*j  ];
        dense_block[2*i+1][j] = sparse_block[i].val[2*j+1];
      }
    }
  }
}

mli_t get_head_multiline(const ml_t *m, const bi_t line_idx, re_t *h, mli_t *h_idx) {
  mli_t i;
  if (m->sz == 0)
    return -1;
  for (i=0; i<m->sz; ++i) {
    if (m->val[2*i+line_idx] != 0) {
      *h      = m->val[2*i+line_idx];
      *h_idx  = i;
      return m->idx[i];
    }
  }
  return -1;
}

mli_t get_head_multiline_hybrid(const ml_t *m,
    const bi_t line_idx, re_t *h, mli_t *h_idx, const ci_t coldim) {
  mli_t i;
  if (m->sz == 0)
    return -1;
  if (m->sz < coldim) {
    return get_head_multiline(m, line_idx, h, h_idx);
  } else {
    for (i=0; i<coldim; ++i) {
      if (m->val[2*i+line_idx] != 0) {
        *h      = m->val[2*i+line_idx];
        *h_idx  = i;
        return i;
      }
    }
  }
  return -1;
}

/* HEADER STUFF PUT HERE */

void init_fl_map_sizes(uint32_t col_size, uint32_t row_size, map_fl_t *map) {
  /*  initialize map arrays and */
  /*  set initial values to __GB_MINUS_ONE_8 */
  map->pc = (ci_t *)malloc(col_size * sizeof(ci_t));
  memset(map->pc, __GB_MINUS_ONE_8, col_size * sizeof(ci_t));

  map->npc  = (ci_t *)malloc(col_size * sizeof(ci_t));
  memset(map->npc, __GB_MINUS_ONE_8, col_size * sizeof(ci_t));

  map->pc_rev = (ci_t *)malloc(col_size * sizeof(ci_t));
  memset(map->pc_rev, __GB_MINUS_ONE_8, col_size * sizeof(ci_t));

  map->npc_rev  = (ci_t *)malloc(col_size * sizeof(ci_t));
  memset(map->npc_rev, __GB_MINUS_ONE_8, col_size * sizeof(ci_t));

  map->pri  = (ri_t *)malloc(col_size * sizeof(ri_t));
  memset(map->pri, __GB_MINUS_ONE_8, col_size * sizeof(ri_t));

  map->npri = (ri_t *)malloc(row_size * sizeof(ri_t));
  memset(map->npri, __GB_MINUS_ONE_8, row_size * sizeof(ri_t));
}

void init_fl_map(sm_t *M, map_fl_t *map) {
  /*  initialize map arrays and */
  /*  set initial values to __GB_MINUS_ONE_8 */

  ri_t max_length = M->ncols > M->nrows ? M->ncols : M->nrows;

  map->pc = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->pc, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));

  map->npc  = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->npc, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));

  map->pc_rev = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->pc_rev, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));

  map->npc_rev  = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->npc_rev, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));

  map->pri  = (ri_t *)malloc(M->ncols * sizeof(ri_t));
  memset(map->pri, __GB_MINUS_ONE_8, M->ncols * sizeof(ri_t));

  map->npri = (ri_t *)malloc(max_length * sizeof(ri_t));
  memset(map->npri, __GB_MINUS_ONE_8, M->nrows * sizeof(ri_t));
}

void realloc_rows_ml(sm_fl_ml_t *A, const mli_t mli,
    const bi_t init_buffer_A, mli_t *buffer_A) {
  *buffer_A +=  init_buffer_A;
  A->ml[mli].idx = (mli_t*) realloc(A->ml[mli].idx, (*buffer_A) * sizeof(mli_t));
  A->ml[mli].val = (re_t*)  realloc(A->ml[mli].val, 2 * (*buffer_A) * sizeof(re_t));
}

void realloc_block_rows(sbm_fl_t *A, const ri_t rbi, const ci_t bir,
    const bi_t lib, const bi_t init_buffer_A, bi_t *buffer_A) {
  *buffer_A +=  init_buffer_A;
  A->blocks[rbi][bir][lib].idx = realloc(
      A->blocks[rbi][bir][lib].idx,
      (*buffer_A) * sizeof(bi_t));
  A->blocks[rbi][bir][lib].val = realloc(
      A->blocks[rbi][bir][lib].val,
      2 * (*buffer_A) * sizeof(re_t));
}

void swap_block_data(sbm_fl_t *A, const ci_t clA, const bi_t rbi,
    const bi_t cvb) {
  uint32_t i ;
  int j, k, l;
  uint16_t ctr[clA];
  memset(ctr, 0, clA * sizeof(uint16_t));
  bi_t *old_idx_ptr = NULL;
  re_t *old_val_ptr = NULL;

  bi_t rounded_cvb_half = cvb/2;
  if (cvb % 2) /*  odd cvb */
    rounded_cvb_half  +=  1;
  for (i=0; i<rounded_cvb_half; ++i) {
    /*  allocate memory for swapping data */
    /* XXX LEAK HERE (when computing reduced form */
    bi_t *tmp_idx_ptr = (bi_t *)malloc(A->blocks[rbi][clA-1][i].sz * sizeof(bi_t));
    re_t *tmp_val_ptr = (re_t *)malloc(2 * A->blocks[rbi][clA-1][i].sz * sizeof(re_t));
    /*  swap data */
    l = 0;
    for (k=A->blocks[rbi][clA-1][i].sz-1; k>-1; --k) {
      tmp_idx_ptr[l]      = A->blocks[rbi][clA-1][i].idx[k];
      tmp_val_ptr[2*l]    = A->blocks[rbi][clA-1][i].val[2*k];
      tmp_val_ptr[2*l+1]  = A->blocks[rbi][clA-1][i].val[2*k+1];
      l++;
      ctr[clA-1]  = 1;
    }
    /*  keep track of old ptr */
    old_idx_ptr = A->blocks[rbi][clA-1][i].idx;
    old_val_ptr = A->blocks[rbi][clA-1][i].val;
    /*  swap starting ptrs */
    A->blocks[rbi][clA-1][i].idx  = tmp_idx_ptr;
    A->blocks[rbi][clA-1][i].val  = tmp_val_ptr;
    /*  try to reuse old, already allocated memory in the following */
    tmp_idx_ptr = old_idx_ptr;
    tmp_val_ptr = old_val_ptr;
    /*  now go through all remaining blocks in the corresponding row (bir) */
    /*  check if the already allocated memory is enough for executing the swap */
    for (j=clA-2; j>-1; --j) {
      if (A->blocks[rbi][j+1][i].sz == A->blocks[rbi][j][i].sz) {
        /*   Memory is not enough, so allocate new (reallocation is used */
        /*   at the moment: It might be prolematic since the old entries are */
        /*   useless and need not be copied, but often the already allocated memory */
        /*   can be extended in place, thus this operation is often faster than */
        /*   freeing memory and allocating new. */
        /*   Note that this reallocation is also used if the old memory allocated */
        /*   is bigger than what is needed. In this setting realloc just cuts the */
        /*   useless memory off. */
      } else {
        tmp_idx_ptr = realloc(tmp_idx_ptr,
            A->blocks[rbi][j][i].sz * sizeof(bi_t));
        tmp_val_ptr = realloc(tmp_val_ptr,
            2 * A->blocks[rbi][j][i].sz * sizeof(re_t));
      }
      l = 0;
      for (k=A->blocks[rbi][j][i].sz-1; k>-1; --k) {
        tmp_idx_ptr[l]      = A->blocks[rbi][j][i].idx[k];
        tmp_val_ptr[2*l]    = A->blocks[rbi][j][i].val[2*k];
        tmp_val_ptr[2*l+1]  = A->blocks[rbi][j][i].val[2*k+1];
        l++;
        ctr[j]  = 1;
      }
      /*  keep track of old ptr */
      old_idx_ptr = A->blocks[rbi][j][i].idx;
      old_val_ptr = A->blocks[rbi][j][i].val;
      /*  swap starting ptrs */
      A->blocks[rbi][j][i].idx  = tmp_idx_ptr;
      A->blocks[rbi][j][i].val  = tmp_val_ptr;
      /*  try to reuse old, already allocated memory in the following */
      tmp_idx_ptr = old_idx_ptr;
      tmp_val_ptr = old_val_ptr;
    }
    /*  finally remove temporary memory used for swapping */
    free(old_idx_ptr);
    free(old_val_ptr);
  }
  for (i=0; i<clA; ++i) {
    if (ctr[i]  ==  0) {
      free(A->blocks[rbi][i]);
      A->blocks[rbi][i] = NULL;
    }
  }
}

void insert_row_data_ml_1_1(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = M->rows[bi1][i1];
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = 0;
  A->ml[mli].sz++;
}

void insert_row_data_ml_1_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi2, const ci_t i2) {
  /* printf("sz %d -- %d\n", A->ml[mli].sz, eil); */
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = 0;
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = M->rows[bi2][i2];
  A->ml[mli].sz++;
}

void insert_row_data_ml_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1,
    const ci_t bi2, const ci_t i2) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = M->rows[bi1][i1];
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = M->rows[bi2][i2];
  A->ml[mli].sz++;
}


void insert_block_row_data_ml_1_1(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = M->rows[bi1][i1];
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = 0;
  A->blocks[rbi][bir][lib].sz++;
}

void insert_block_row_data_ml_1_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi2, const ci_t i2) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = 0;
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = M->rows[bi2][i2];
  A->blocks[rbi][bir][lib].sz++;
}

void insert_block_row_data_ml_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1, const ci_t bi2, const ci_t i2) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = M->rows[bi1][i1];
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = M->rows[bi2][i2];
  A->blocks[rbi][bir][lib].sz++;
}


/* vim:sts=2:sw=2:ts=2:et:
*/
