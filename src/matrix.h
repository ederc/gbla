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

/**
 * \file matrix.h
 * \brief Interfaces for matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_MATRIX_H
#define GB_MATRIX_H

#include <gbla_config.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <types.h>
#if GBLA_WITH_FFLAS
#include "../draft/matrix.h"
#endif

/**
 * \brief Sparse matrix structure for reading jcf matrices
 *
 */

typedef struct sm_t {
  ri_t nrows;     /*!<  number of rows */
  ci_t ncols;     /*!<  number of columns */
  nnz_t nnz;      /*!<  number of nonzero entries */
  float density;  /*!<  density used for adjusting memory allocation during
                        splicing and generation of ABCD blocks */
  mod_t mod;      /*!<  modulo/field characteristic */
  float fs;       /*!<  file size of input matrix */
  char *fsu;      /*!<  file size unit of input matrix, e.g. GB */
  re_t **rows;    /*!<  address of row: M->rows[i] gives first
                        address of ith row */
  ci_t **pos;     /*!<  position of entry in row: M->pos[i] gives first
                        address of first position of nonzero entry in row i */
  ci_t *rwidth;   /*!<  width of row: M->rwidth[i] gives number of nonzero
                        entries in row i */
  ci_t *buf;      /*!<  stores buffer of memory allocated for the given row*/
} sm_t;


/**
 * \brief A sparse matrix block
 */
typedef struct sbl_t {
  re_t **row; /*!< row entries */
  bi_t **pos; /*!< position in row */
  bi_t *sz;   /*!< size of row */
  bi_t *buf;  /*!< memory buffer already allocated */
} sbl_t;

/**
 * \brief Sparse block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small sparse blocks for the A part.
 */

typedef struct sb_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  sbl_t **blocks;   /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. */
} sb_fl_t;

/**
 * \brief Sparse matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small sparse blocks for the A part.
 */

typedef struct sm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  re_t **row;       /*!< row entries */
  ci_t **pos;       /*!< position in row */
  ci_t *sz;         /*!< size of row */
  ci_t *buf;        /*!< memory buffer already allocated */
} sm_fl_t;

/**
 * \brief A dense block is a block of size  __GBLA_SIMD_BLOCK_SIZE^2 of
 * matrix entries.
 */
typedef struct dbl_t {
  re_t *val;  /*!< rows */
} dbl_t;

/**
 * \brief A multiline block is a vector of an vector of multilines. It consists of
 * an index for the column of the value entries in rows. For each index rows stores
 * __GB_NROWS_MULTILINE elements, i.e. this many rows are taken care of at once.
 */
typedef struct mbl_t {
  bi_t *idx;  /*!< column index in the multiline vector */
  re_t *val;  /*!< multiline row, must be __GB_NROWS_MULTILINE * length(idx) */
  bi_t sz;    /*!< current length of the block row */
  bi_t dense; /*!< if 1 the multiline row is in dense representation */
} mbl_t;

/**
 * \brief A multiline is a vector of an vector of multilines. It consists of
 * an index for the column of the value entries in rows. For each index rows stores
 * __GB_NROWS_MULTILINE elements, i.e. this many rows are taken care of at once.
 *
 * \note In spite of type mbl_t the index idx as well as the size sz might
 * excced 2^16, thus bi_t is not enough and we need them to be of type ci_t
 * resp. uint32_t at least.
 *
 */
typedef struct ml_t {
  mli_t *idx; /*!< column index in the multiline vector */
  re_t  *val; /*!< multiline row, must be __GB_NROWS_MULTILINE * length(idx) */
  ci_t  sz;    /*!< current length of the multiline row */
  bi_t  dense; /*!< if 1 the multiline row is in dense representation */
} ml_t;


/**
 * \brief Enum of block alignments: Entries in the block submatrices are stored
 * w.r.t. different alignments, some from top to down, then from left to right,
 * other in combinations and inversions of these.
 */
enum ba_t {
  tdlr, /*!<  top-to-down, left-to-right */
  tdrl, /*!<  top-to-down, right-to-left */
  dtlr, /*!<  down-to-top, left-to-right */
  dtrl  /*!<  down-to-top, right-to-left */
};


/**
 * \brief Dense block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small dense blocks for exploiting SIMD
 * inctructions.
 */

typedef struct dbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  dbl_t **blocks;   /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. */
} dbm_fl_t;

/**
 * \brief Hybrid simd block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small inner dense blocks for exploiting
 * SIMD inctructions.
 *
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * 
 * The above picture represents one block of size __GBLA_SIMD_BLOCK_SIZE_RECT.
 * Inside we store small sub blocks of height 1 and width __GBLA_SIMD_INNER_SIZE.
 * Such an inner sub block is either NULL (if all entries are zero) or
 * represented in dense fashion.
 * All in all this leads to the fact that blocks is of type dbl_t ****blocks:
 * blocks[outer_row][outer_col][inner_row][inner_col].
 */
typedef struct hbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  dbl_t ****blocks; /*!<  address of blocks: M->blocks[i][j][k][l] as explained
                          above */
} hbm_fl_t;

/**
 * \brief Sparse block matrix structure for Faugère-Lachartre decompositions.
 * Can be used for usual line implementations and multi line implementations,
 * the corresponding multi line functions are labeled with an "_ml".
 *
 */

typedef struct sbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  ri_t bheight;     /*!<  number of rows per block */
  ci_t bwidth;      /*!<  number of columns per block */
  enum ba_t ba;     /*!<  memory alignment in the block: depending on
                          the block parts A,B,C and D the entries are
                          stored top-to-down, left-to-right, or any
                          combination and inversion of these two */
  int fe;           /*!<  if 1 then empty blocks are used to fill
                          block matrix
                          if 0 then no fill is done */
  int hr;           /*!<  if 0 then those rows are not accepted
                          if 1 then hybrid (sparse/dense) rows are accepted */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  mbl_t ***blocks;  /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. There are nrows/bheight * ncols/bwidth blocks. */
} sbm_fl_t;

/**
 * \brief Sparse matrix with multiline row structure for Faugère-Lachartre
 * decompositions. Can be used for multiline implementations, the corresponding
 * multiline functions are labeled with an "_ml".
 */

typedef struct sm_fl_ml_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  enum ba_t ba;     /*!<  memory alignment in the block: depending on
                          the block parts A,B,C and D the entries are
                          stored top-to-down, left-to-right, or any
                          combination and inversion of these two */
  int fe;           /*!<  if 1 then empty blocks are used to fill
                          block matrix
                          if 0 then no fill is done */
  int hr;           /*!<  if 1 then hybrid (sparse/dense) rows are accepted
                          if 0 then those rows are not accepted */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  ml_t *ml;         /*!<  address of multilines: M->ml[i] gives address of
                          multiline i. There are nrows/__GB_NROWS_MULTILINE
                          multilines. */
} sm_fl_ml_t;

/**
 * \brief returns number of row blocks for hybrid block submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \return number of row blocks for A
 */
inline ri_t get_number_hybrid_row_blocks(const hbm_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for hybrid block submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \return number of column blocks for A
 */
inline ci_t get_number_hybrid_col_blocks(const hbm_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of row blocks for dense block submatrix A
 *
 * \param dense block submatrix A
 *
 * \return number of row blocks for A
 */
inline ri_t get_number_dense_row_blocks(const dbm_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for dense block submatrix A
 *
 * \param dense block submatrix A
 *
 * \return number of column blocks for A
 */
inline ci_t get_number_dense_col_blocks(const dbm_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of row blocks for sparse block submatrix A
 *
 * \param sparse block submatrix A
 *
 * \return number of row blocks for A
 */
inline ri_t get_number_sparse_row_blocks(const sb_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for sparse block submatrix A
 *
 * \param sparse block submatrix A
 *
 * \return number of column blocks for A
 */
inline ci_t get_number_sparse_col_blocks(const sb_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief Initializes sparse submatrices
 *
 * \param sparse submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
inline void init_sm(sm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  int i;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for rows
  A->row  = (re_t **)malloc(nrows * sizeof(re_t *));
  A->pos  = (ci_t **)malloc(nrows * sizeof(ci_t *));
  A->sz   = (ci_t *)calloc(nrows, sizeof(ci_t));
  A->buf  = (ci_t *)calloc(nrows, sizeof(ci_t));
  for (i=0; i<nrows; ++i) {
    A->row[i] = NULL;
    A->pos[i] = NULL;
  }
}

/**
 * \brief Initializes sparse block submatrices
 *
 * \param sparse block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
inline void init_sb(sb_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  int i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_sparse_col_blocks(A);
  const ri_t rlA  = get_number_sparse_row_blocks(A);
  // we know that no block of A will be empty, thus we can already allocate
  // corresponding memory
  A->blocks = (sbl_t **)malloc(rlA * sizeof(sbl_t *));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (sbl_t *)malloc(clA * sizeof(sbl_t));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j].row = NULL;
      A->blocks[i][j].pos = NULL;
      A->blocks[i][j].sz  = NULL;
      A->blocks[i][j].buf = NULL;
      /*
      A->blocks[i][j].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
      A->blocks[i][j].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
      A->blocks[i][j].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
      A->blocks[i][j].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
        A->blocks[i][j].row[k]  = NULL;
        A->blocks[i][j].pos[k]  = NULL;
        A->blocks[i][j].sz[k]   = 0;
        A->blocks[i][j].buf[k]  = 0;
      }
      */
    }
  }
}

/**
 * \brief Initializes hybrid block submatrices
 *
 * \param hybrid block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
inline void init_hbm(hbm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  int i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_hybrid_col_blocks(A);
  const ri_t rlA  = get_number_hybrid_row_blocks(A);

  A->blocks = (dbl_t ****)malloc(rlA * sizeof(dbl_t ***));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (dbl_t ***)malloc(clA * sizeof(dbl_t **));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j] = NULL;
    }
  }
}

/**
 * \brief Initializes dense block submatrices
 *
 * \param dense block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
inline void init_dbm(dbm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  int i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_dense_col_blocks(A);
  const ri_t rlA  = get_number_dense_row_blocks(A);

  A->blocks = (dbl_t **)malloc(rlA * sizeof(dbl_t *));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (dbl_t *)malloc(clA * sizeof(dbl_t));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j].val = NULL;
    }
  }
}

/**
 * \brief Frees given sparse submatrix A
 *
 * \param sparse submatrix A
 *
 * \param number of threads for parallel computation
 */
inline void free_sparse_matrix(sm_fl_t **A_in, int nthrds)
{
  sm_fl_t *A      = *A_in;
  ri_t i, j, k, l;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i, j, k, l)
    for (i=0; i<A->nrows; ++i) {
      free(A->row[i]);
      free(A->pos[i]);
    }
  }
  free(A->row);
  free(A->pos);
  free(A->sz);
  free(A->buf);
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given sparse block submatrix A
 *
 * \param sparse block submatrix A
 *
 * \param number of threads for parallel computation
 */
inline void free_sparse_submatrix(sb_fl_t **A_in, int nthrds)
{
  sb_fl_t *A      = *A_in;
  const ci_t clA  = get_number_sparse_col_blocks(A);
  const ri_t rlA  = get_number_sparse_row_blocks(A);
  ri_t i, j, k, l;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i, j, k, l)
    for (i=0; i<rlA; ++i) {
      for (j=0; j<clA; ++j) {
        if (A->blocks[i][j].row != NULL) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            free(A->blocks[i][j].row[k]);
            free(A->blocks[i][j].pos[k]);
          }
          free(A->blocks[i][j].row);
          free(A->blocks[i][j].pos);
          free(A->blocks[i][j].sz);
          free(A->blocks[i][j].buf);
        }
      }
      free(A->blocks[i]);
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given hybrid block submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \param number of threads for parallel computation
 */
inline void free_hybrid_submatrix(hbm_fl_t **A_in, int nthrds)
{
  hbm_fl_t *A     = *A_in;
  const ci_t clA  = get_number_hybrid_col_blocks(A);
  const ri_t rlA  = get_number_hybrid_row_blocks(A);
  ri_t i, j, k, l;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i, j, k, l)
    for (i=0; i<rlA; ++i) {
      for (j=0; j<clA; ++j) {
        if (A->blocks[i][j] != NULL) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            if (A->blocks[i][j][k] != NULL) {
              for (l=0; l<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++l) {
                if (A->blocks[i][j][k][l].val != NULL) {
                  free(A->blocks[i][j][k][l].val);
                  A->blocks[i][j][k][l].val  = NULL;
                }
              }
              free(A->blocks[i][j][k]);
              A->blocks[i][j][k]  = NULL;
            }
          }
          free(A->blocks[i][j]);
          A->blocks[i][j] = NULL;
        }
      }
      free(A->blocks[i]);
      A->blocks[i]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given dense block submatrix A
 *
 * \param dense block submatrix A
 *
 * \param number of threads for parallel computation
 */
inline void free_dense_submatrix(dbm_fl_t **A_in, int nthrds)
{
  dbm_fl_t *A     = *A_in;
  const ci_t clA  = get_number_dense_col_blocks(A);
  const ri_t rlA  = get_number_dense_row_blocks(A);
  ri_t j;
  ci_t i;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i].val != NULL) {
          free(A->blocks[j][i].val);
          A->blocks[j][i].val  = NULL;
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
}

static inline sb_fl_t *copy_sparse_to_block_matrix(const sm_fl_t *A,
    const int nthrds)
{
  sb_fl_t *out  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  init_sb(out, A->nrows, A->ncols);

  const ri_t blocks_per_col = get_number_sparse_row_blocks(out);
  const ci_t blocks_per_row = get_number_sparse_col_blocks(out);

#pragma omp parallel num_threads(nthrds)
  {
    bi_t *block_length  = (bi_t *)calloc(blocks_per_row, sizeof(bi_t));
    bi_t col_idx;
    ri_t i, j;
    ci_t k, l, m, n;
    bi_t tmp1, tmp2;
    ci_t ctr;
    ri_t max;
#pragma omp for   
    for (i=0; i<blocks_per_col; ++i) {
      max = A->nrows < (i+1)*__GBLA_SIMD_BLOCK_SIZE ? A->nrows : (i+1)*__GBLA_SIMD_BLOCK_SIZE;
      for (j=i*__GBLA_SIMD_BLOCK_SIZE; j<max; ++j) {
        memset(block_length, 0, blocks_per_row * sizeof(bi_t));
        // get lengths for the sparse block rows
        for (k=0; k<A->sz[j]-1; ++k) {
          tmp1  = A->pos[j][k] / __GBLA_SIMD_BLOCK_SIZE;
          tmp2  = A->pos[j][k+1] / __GBLA_SIMD_BLOCK_SIZE;
          if (tmp1 == tmp2)
            block_length[tmp1]++;
        }
        if (tmp1 == tmp2)
          block_length[tmp1]++;
        else
          block_length[tmp2]++;
        
        col_idx = j % __GBLA_SIMD_BLOCK_SIZE;
        // allocate memory
        for (l=0; l<blocks_per_row; ++l) {
          if (block_length[l] != 0) {
            if (out->blocks[i][l].row == NULL) {
              out->blocks[i][l].row = (re_t **)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
              out->blocks[i][l].pos = (bi_t **)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
              out->blocks[i][l].sz = (bi_t *)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
              for (m=0; m<__GBLA_SIMD_BLOCK_SIZE; ++m) {
                out->blocks[i][l].row[m]  = NULL;
                out->blocks[i][l].pos[m]  = NULL;
                out->blocks[i][l].sz[m]   = 0;
              }
            }
          // memory for row j % __GBLA_SIMD_BLOCK_SIZE
          out->blocks[i][l].row[col_idx]  = (re_t *)malloc(
              block_length[l] * sizeof(re_t));
          out->blocks[i][l].pos[col_idx]  = (bi_t *)malloc(
              block_length[l] * sizeof(bi_t));
          }
        }

        // copy elements
        ctr = 0;
        for (l=0; l<blocks_per_row; ++l) {
          for (m=ctr; m<(ctr+block_length[l]); ++m) {
            out->blocks[i][l].row[col_idx][out->blocks[i][l].sz[col_idx]] =
              A->row[j][m];
            out->blocks[i][l].pos[col_idx][out->blocks[i][l].sz[col_idx]] =
              A->pos[j][m] % __GBLA_SIMD_BLOCK_SIZE;
            out->blocks[i][l].sz[col_idx]++;
          }
          ctr +=  block_length[l];
        }
      }
    }    
  }
  return out;
}
/**
 * \brief Copies data from block matrix in to input sparse matrix format. If deleteIn
 * is set the input block matrix in is deleted.
 *
 * \note This procedure is only used for complete reduced row echelon forms
 *
 * \param pointer to block matrix input_bl
 *
 * \param pointer to multiline matrix input_ml
 *
 * \param rank of multiline matrix rank_input_ml
 *
 * \param sparse matrix output
 *
 * \param deleteIn, if 1 the input matrix in is deleted, if 0 in is not deleted.
 *
 * \param number of threads to copy data in parallel nthrds
 */
void copy_block_ml_matrices_to_sparse_matrix(sbm_fl_t **input_bl,
    sm_fl_ml_t **input_ml, ri_t rank_input_ml, sm_t **output,
    int deleteIn, int nthrds);

#if GBLA_WITH_FFLAS
/**
 * \brief Copies data from block matrix into DNS matrix format.
 *
 * \param pointer to block matrix source
 *
 * \param pointer to DNS matrix destination
 */
void copy_block_ml_matrix_to_dns_matrix(sbm_fl_t **source, DNS **destination);
#endif

/**
 * \brief Copies data from block matrix in to multiline matrix out. If deleteIn
 * is set the input block matrix in is deleted.
 *
 * \param pointer to block matrix in, input matrix
 *
 * \param multiline matrix out, output matrix
 *
 * \param deleteIn, if 1 the input matrix in is deleted, if 0 in is not deleted.
 *
 * \param number of threads to copy data in parallel nthrds
 *
 * \return multiline matrix generated out of block input matrix
 */
sm_fl_ml_t *copy_block_matrix_to_multiline_matrix(sbm_fl_t **input,
     sm_fl_ml_t *out, int deleteIn, int nthrds);

/**
 * \brief Copies multiline format matrix A_in to block format matrix B
 *
 * \param pointer to pointer to multiline matrix A_in
 *
 * \param block height bheight
 *
 * \param block width bwidth
 -*
 * \param freeing memory? if 1 then memory is freed, otherwise not, free_memory
 *
 * \param number of threads for parallel computations nthrds
 *
 * \return matrix A in block format
 */
sbm_fl_t *copy_multiline_to_block_matrix_rl(sm_fl_ml_t **A_in,
    ri_t bheight, ci_t bwidth, int free_memory, int nthrds);
#endif

double compute_density(nnz_t nnz, ri_t nrows, ri_t ncols) ;


/* vim:sts=2:sw=2:ts=2:
 */
