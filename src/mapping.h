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
 * \file mapping.h
 * \brief Interfaces for sparse matrix index maps used for subdividing the matrix
 * in Faugère-Lachartre style
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_MAPPING_H
#define GB_MAPPING_H

#include <gbla_config.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <elimination.h>
#include <matrix.h>
#include <math.h>
#include <omp.h>

#define BUFFER  256

/**
 * \brief Indexer for subdividing sparse matrix into 4 parts as described by
 * Faugère and Lachartre in http://dx.doi.org/10.1145/1837210.1837225.
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 */

typedef struct map_fl_t {
  ri_t npiv;      /*!<  number of pivots found in input matrix */
  ci_t *pc;       /*!<  map of pivot columns: from input matrix M to
                        submatrix A; has length M->ncols, maps non-pivot
                        columns to __GB_MINUS_ONE_32 */
  ci_t *npc;      /*!<  map of non-pivot columns: from input matrix M to
                        submatrix B; has length M->ncols, maps pivot
                        columns to __GB_MINUS_ONE_32 */
  ci_t *pc_rev;   /*!<  reverse map of pivot columns: from submatrices
                        A and B to input matrix M */
  ci_t *npc_rev;  /*!<  reverse map of non-pivot columns: from submatrices
                        A and B to input matrix M */
  ri_t *pri;      /*!<  has length M->nrows, maps pivot columns to
                        their corresponding row index, maps non-pivot
                        columns to __GB_MINUS_ONE_32 */
  ri_t *npri;     /*!<  indices of non-pivot rows */

  int nthrds;     /*!<  number of threads to be used for the indexer */
} map_fl_t;


/**
 * \brief Initializes a map for the input matrix M with all entries set
 * to __GB_MINUS_ONE_8.
 *
 * \param size of columns col_size
 *
 * \param size of rows row_size
 *
 * \param map to be initialized
 *
 */
/* static */ /* inline */ void init_fl_map_sizes(uint32_t col_size, uint32_t row_size, map_fl_t *map);

/**
 * \brief Initializes a map for the input matrix M with all entries set
 * to __GB_MINUS_ONE_8.
 *
 * \param original matrix M
 *
 * \param map to be initialized
 *
 */
/* static */ /* inline */ void init_fl_map(sm_t *M, map_fl_t *map) ;

/**
 * \brief Process multiline matrix and applies mapping for it
 *
 * \param multiline matrix A
 *
 * \param map to be defined map
 *
 * \param block height of the corresponding block matrices bheight
 */
void process_matrix(sm_fl_ml_t *A, map_fl_t *map, const bi_t bheight);

/**
 * \brief Combines the mappings from the outer input matrix A and the
 * echelonized multiline matrix D
 *
 * \param mapping of the input matrix M outer_map
 *
 * \param pointer to mapping of the echelonized multiline matrix D. Will be
 * freed at the end. inner_map_in
 *
 * \param column dimension of input matrix M outer_coldim
 *
 * \param column dimension of multiline matrix D inner_coldim
 *
 * \param parameter that defines if only the new pivots are remap, i.e. no
 * columns are considered. This is done when the echelon form is wanted and not
 * the RREF
 */
void combine_maps(map_fl_t *outer_map, map_fl_t **inner_map_in,
    const ci_t outer_coldim, const ci_t inner_coldim,
    int only_rows);

/**
 * \brief Reconstructs matrix M after elimination process
 *
 * \param input and output matrix M
 *
 * \param block sub matrix A
 *
 * \param block sub matrix B
 *
 * \param multiline submatrix D
 *
 * \param outer mapping map
 *
 * \param column dimension of M coldim
 *
 * \param parameter to set/unset freeing of sub matrices free_matrices
 *
 * \param paramter indicating if M was already freed M_freed
 *
 * \param paramter indicating if A was already freed A_freed
 *
 * \param paramter indicating if D was already freed D_freed
 *
 * \param number of threads
 */
void reconstruct_matrix_block(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sm_fl_ml_t *D,
    map_fl_t *map, const ci_t coldim, int free_matrices,
    int M_freed, int A_freed, int D_freed, int nthrds);

/**
 * \brief Reconstructs matrix M after full reduced row echelon process
 *
 * \param input and output matrix M
 *
 * \param block sub matrix A
 *
 * \param block sub matrix B2
 *
 * \param block sub matrix D2
 *
 * \param outer mapping map
 *
 * \param column dimension of M coldim
 *
 * \param parameter to set/unset freeing of sub matrices free_matrices
 *
 * \param paramter indicating if M was already freed M_freed
 *
 * \param paramter indicating if A was already freed A_freed
 *
 * \param paramter indicating if D was already freed D_freed
 *
 * \param number of threads
 */
void reconstruct_matrix_block_reduced(sm_t *M, sbm_fl_t *A, sbm_fl_t *B2, sbm_fl_t *D2,
    map_fl_t *map, const ci_t coldim, int free_matrices,
    int M_freed, int A_freed, int D_freed, int nthrds);

/**
 * \brief Reconstructs matrix M after elimination process
 *
 * \param input and output matrix M
 *
 * \param multiline sub matrix A
 *
 * \param block sub matrix B
 *
 * \param multiline submatrix D
 *
 * \param outer mapping map
 *
 * \param column dimension of M coldim
 *
 * \param parameter to set/unset freeing of sub matrices free_matrices
 *
 * \param paramter indicating if M was already freed M_freed
 *
 * \param paramter indicating if A was already freed A_freed
 *
 * \param paramter indicating if D was already freed D_freed
 *
 * \param number of threads
 */
void reconstruct_matrix_ml(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, sm_fl_ml_t *D,
    map_fl_t *map, const ci_t coldim, int free_matrices,
    int M_freed, int A_freed, int D_freed, int nthrds);

/**
 * \brief Reallocates memory for the rows of the multiline during the splicing of
 * the input matrix. The buffer size buffer_A is doubled during this process
 *
 * \param multiline sub matrix A
 *
 * \param multiline index mli
 *
 * \param buffer size buffer_A
 */
/* static */ /* inline */ void realloc_rows_ml(sm_fl_ml_t *A, const mli_t mli,
    const bi_t init_buffer_A, mli_t *buffer_A);

/**
 * \brief Reallocates memory for the rows of the blocks during the splicing of
 * the input matrix. The buffer size buffer_A is doubled during this process
 *
 * \param block matrix A
 *
 * \param row block index rbi in A
 *
 * \param block index in row bir
 *
 * \param block index in row bir
 *
 * \param line index in block lib
 *
 * \param buffer size buffer_A
 *
 */
/* static */ /* inline */ void realloc_block_rows(sbm_fl_t *A, const ri_t rbi, const ci_t bir,
    const bi_t lib, const bi_t init_buffer_A, bi_t *buffer_A) ;

/**
 * \brief Swaps data arrangement in leftsided block matrices
 *
 * \param block matrix A
 *
 * \param number of blocks in the corresponding row clA
 *
 * \param current row index for blocks rbi
 *
 * \param current line in block lib
 *
 */
/* static */ /* inline */ void swap_block_data(sbm_fl_t *A, const ci_t clA, const bi_t rbi,
    const bi_t cvb) ;

/**
 * \brief Inserts elements from input matrix M in dense diagonal block of submatrix A
 *
 * \note For diagonal blocks memory is compressed to use only upper triangular
 * block space. We use the compressed format for upper resp. lowerd triangular
 * matrices also known as packed storage from LAPACK, see, for example,
 * http://www.netlib.org/lapack/lug/node123.html
 *
 * upper triangular, not inverted, aij is stored in AP(i+j(j-1)/2)
 * upper triangular, but inverted (for A in ABCD), aij is stored in AP( i+(2n-j)(j-1)/2)
 *
 * \param dense block submatrix A
 *
 * \param original matrix M
 *
 * \param current row block index rbi
 *
 * \param current block index in block row bir
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_dense_block_data_diagonalize(dbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const ci_t eil,
    const ci_t bi1, const ci_t i1)
{
  assert(rbi == bir);
  A->blocks[rbi][bir].val[lib + ((2*__GBLA_SIMD_BLOCK_SIZE-(eil+1))*eil)/2] =
    (re_t)((re_m_t)M->mod - M->rows[bi1][i1]);
}

/**
 * \brief Inserts elements from input matrix M in dense block of submatrix A and
 * inverts w.r.t. M->mod.
 *
 * \param dense block submatrix A
 *
 * \param original matrix M
 *
 * \param current row block index rbi
 *
 * \param current block index in block row bir
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_dense_block_data_inv(dbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const ci_t eil,
    const ci_t bi1, const ci_t i1)
{
  A->blocks[rbi][bir].val[(lib*__GBLA_SIMD_BLOCK_SIZE)+eil] =
    (re_t)((re_m_t)M->mod - M->rows[bi1][i1]);
}

/**
 * \brief Inserts elements from input matrix M in dense block of submatrix A
 *
 * \param dense block submatrix A
 *
 * \param original matrix M
 *
 * \param current row block index rbi
 *
 * \param current block index in block row bir
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_dense_block_data(dbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const ci_t eil,
    const ci_t bi1, const ci_t i1)
{
  A->blocks[rbi][bir].val[(lib*__GBLA_SIMD_BLOCK_SIZE)+eil] = M->rows[bi1][i1];
}

/**
 * \brief Inserts elements from input matrix M in hybrid submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_hbm(hbm_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  bi_t i, j;
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir] == NULL) {
    A->blocks[rbi][bir] = (dbl_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(dbl_t *));
    for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
      A->blocks[rbi][bir][i]  = NULL;
    }
    A->blocks[rbi][bir][lib]  = (dbl_t *)malloc(__GBLA_SIMD_INNER_BLOCKS_PER_ROW * sizeof(dbl_t));
    for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
      A->blocks[rbi][bir][lib][j].val = NULL;
    }
    A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val = (re_t *)calloc(
        __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
  } else {
    if (A->blocks[rbi][bir][lib] == NULL) {
      A->blocks[rbi][bir][lib]  = (dbl_t *)malloc(__GBLA_SIMD_INNER_BLOCKS_PER_ROW * sizeof(dbl_t));
      for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
        A->blocks[rbi][bir][lib][j].val = NULL;
      }
      A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val = (re_t *)calloc(
          __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
    } else {
      if (A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val  ==  NULL) {
        A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val = (re_t *)calloc(
            __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
      }
    }
  }
  // set values
  A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val[eil%__GBLA_SIMD_INNER_SIZE]  = 
    M->rows[bi][ri];
}

/**
 * \brief Swaps sorting of elements inside the blocks and cuts memory down for
 * sparse block matrix A.
 *
 * \param sparse block submatrix A
 */
static inline void swap_and_cut(sb_fl_t *A) {
  const ri_t rlA  = get_number_sparse_row_blocks(A);
  const ci_t clA  = get_number_sparse_col_blocks(A);
  bi_t sz, k, l;
  ri_t i;
  ci_t j;
  for (i=0; i<rlA; ++i) {
    if (A->blocks[i] != NULL) {
      for (j=0; j<clA; ++j) {
        if (A->blocks[i][j].row != NULL) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            sz  = A->blocks[i][j].sz[k];
            re_t *temp_row  = (re_t *)malloc(sz * sizeof(re_t));
            bi_t *temp_pos  = (bi_t *)malloc(sz * sizeof(bi_t));
            for (l=0; l<sz; ++l) {
              temp_row[l] = A->blocks[i][j].row[k][sz-l-1];
              temp_pos[l] = A->blocks[i][j].pos[k][sz-l-1];
            }
            free(A->blocks[i][j].row[k]);
            free(A->blocks[i][j].pos[k]);
            A->blocks[i][j].row[k]  = temp_row;
            A->blocks[i][j].pos[k]  = temp_pos;
          }
        }
      }
    }
  }
}

/**
 * \brief Cuts memory down for sparse block matrix A.
 *
 * \param sparse block submatrix A
 */
static inline void cut_blocks(sb_fl_t *A) {
  const ri_t rlA  = get_number_sparse_row_blocks(A);
  const ci_t clA  = get_number_sparse_col_blocks(A);
  bi_t sz, k, l;
  ri_t i;
  ci_t j;
  for (i=0; i<rlA; ++i) {
    if (A->blocks[i] != NULL) {
      for (j=0; j<clA; ++j) {
        if (A->blocks[i][j].row != NULL) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            sz  = A->blocks[i][j].sz[k];
            A->blocks[i][j].row[k] = realloc(A->blocks[i][j].row[k],sz * sizeof(re_t));
            A->blocks[i][j].pos[k] = realloc(A->blocks[i][j].pos[k],sz * sizeof(re_t));
          }
        }
      }
    }
  }
}

/**
 * \brief Cuts memory down in sparse submatrix A.
 *
 * \param sparse submatrix A
 */
static inline void cut(sm_fl_t *A) {
  ri_t i;
  for (i=0; i<A->nrows; ++i) {
    A->row[i] = realloc(A->row[i], A->sz[i] * sizeof(re_t));
    A->pos[i] = realloc(A->pos[i], A->sz[i] * sizeof(ci_t));
  }
}

/**
 * \brief Inserts elements from input matrix M in sparse submatrix A. Before
 * insertion the elements are inverted w.r.t. M->mod.
 *
 * \param sparse submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_sm_inv(sm_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const ri_t ridx = (rbi * __GBLA_SIMD_BLOCK_SIZE) + lib;

  bi_t i, j, k;
  // allocate memory if needed
  if (A->sz[ridx] == A->buf[ridx]) {
    A->buf[ridx]  *=  2;
    A->row[ridx]  =   realloc(A->row[ridx], A->buf[ridx] * sizeof(re_t));
    A->pos[ridx]  =   realloc(A->pos[ridx], A->buf[ridx] * sizeof(ci_t));
  }
  // set values
  A->row[ridx][A->sz[ridx]] = (re_t)((re_m_t)M->mod - M->rows[bi][ri]);
  A->pos[ridx][A->sz[ridx]] = shift;
  A->sz[ridx]++;
}

/**
 * \brief Inserts many elements from input matrix M in sparse submatrix A at
 * once.
 *
 * \param sparse submatrix A
 *
 * \param original matrix M
 *
 * \param vector of pivot resp. non-pivot column entries in map col
 *
 * \param array of column shift entries in map pos
 *
 * \param array of row entry positions from M ri
 *
 * \param length of arrays pos and ri length
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param block index bi
 */
static inline void insert_many_in_sm(sm_fl_t *A, const sm_t *M,
    const ci_t *pos, const ci_t *ri, const bi_t length, const ri_t rbi,
    const ri_t lib, const ri_t bi)
{
  const ri_t ridx = (rbi * __GBLA_SIMD_BLOCK_SIZE) + lib;

  bi_t i;
  // allocate memory
  if (A->sz[ridx]+length >= A->buf[ridx]) {
    A->buf[ridx]  *=  2;
    A->row[ridx]  =   realloc(A->row[ridx], A->buf[ridx] * sizeof(re_t));
    A->pos[ridx]  =   realloc(A->pos[ridx], A->buf[ridx] * sizeof(ci_t));
  }
  // set values
  for (i=0; i<length; ++i) {
    A->row[ridx][A->sz[ridx]] = (re_t)M->rows[bi][ri[i]];
    A->pos[ridx][A->sz[ridx]] = pos[i];
    A->sz[ridx]++;
  }
}

/**
 * \brief Inserts elements from input matrix M in sparse submatrix A.
 *
 * \param sparse submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_sm(sm_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const ri_t ridx = (rbi * __GBLA_SIMD_BLOCK_SIZE) + lib;

  bi_t i, j, k;
  // allocate memory if needed
  if (A->row[ridx] == NULL) {
    A->row[ridx]  = (re_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t));
    A->pos[ridx]  = (ci_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(ci_t));
    A->buf[ridx]  = __GBLA_SIMD_BLOCK_SIZE;
  } else {
    if (A->sz[ridx] == A->buf[ridx]) {
      A->buf[ridx]  *=  2;
      A->row[ridx]  =   realloc(A->row[ridx], A->buf[ridx] * sizeof(re_t));
      A->pos[ridx]  =   realloc(A->pos[ridx], A->buf[ridx] * sizeof(ci_t));
    }
  }
  // set values
  A->row[ridx][A->sz[ridx]] = (re_t)M->rows[bi][ri];
  A->pos[ridx][A->sz[ridx]] = shift;
  A->sz[ridx]++;
}

/**
 * \brief Inserts elements from input matrix M in sparse submatrix A in columns.
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_sb_by_column(sb_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir  = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil  = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  bi_t i, j, k;
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].row == NULL) {
    A->blocks[rbi][bir].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
    A->blocks[rbi][bir].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
    A->blocks[rbi][bir].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      A->blocks[rbi][bir].row[k]  = NULL;
      A->blocks[rbi][bir].pos[k]  = NULL;
      A->blocks[rbi][bir].sz[k]   = 0;
      A->blocks[rbi][bir].buf[k]  = 0;
    }
  }
  if (A->blocks[rbi][bir].row[eil] == NULL) {
    A->blocks[rbi][bir].row[eil]  = (re_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(re_t));
    A->blocks[rbi][bir].pos[eil]  = (bi_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf[eil]  = 2 * __GBLA_SIMD_INNER_SIZE;
  } else {
    if (A->blocks[rbi][bir].sz[eil] == A->blocks[rbi][bir].buf[eil]) {
      A->blocks[rbi][bir].buf[eil]  *=  2;
      A->blocks[rbi][bir].row[eil]  =   realloc(A->blocks[rbi][bir].row[eil],
          A->blocks[rbi][bir].buf[eil] * sizeof(re_t));
      A->blocks[rbi][bir].pos[eil]  =   realloc(A->blocks[rbi][bir].pos[eil],
          A->blocks[rbi][bir].buf[eil] * sizeof(bi_t));
    }
  }
  // set values
  A->blocks[rbi][bir].row[eil][A->blocks[rbi][bir].sz[eil]] = 
    (re_t)M->rows[bi][ri];
  A->blocks[rbi][bir].pos[eil][A->blocks[rbi][bir].sz[eil]] = lib;
  A->blocks[rbi][bir].sz[eil]++;
}

/**
 * \brief Inserts several elements from input matrix M in sparse submatrix A at
 * once.
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param array of column positions pos
 *
 * \param array of row indices ri
 *
 * \param length of pos resp. ri length
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param row index of corresponding element in M bi
 */
static inline void insert_many_in_sb_by_column(sb_fl_t *A, const sm_t *M, const ci_t *pos,
    const ci_t *ri, const bi_t length, const ri_t rbi, const ri_t lib, const ri_t bi)
{
  bi_t i, k;
  register bi_t eil;
  // bir is the same for all elements
  const register bi_t bir  = pos[i] / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].row == NULL) {
    A->blocks[rbi][bir].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
    A->blocks[rbi][bir].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
    A->blocks[rbi][bir].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      A->blocks[rbi][bir].row[k]  = (re_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t));
      A->blocks[rbi][bir].pos[k]  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
      A->blocks[rbi][bir].sz[k]   = 0;
      //A->blocks[rbi][bir].buf[k]  = 0;
    }
  }
  for (i=0; i<length; ++i) {
    eil  = pos[i] % __GBLA_SIMD_BLOCK_SIZE; // index in block line
    /*
    if (A->blocks[rbi][bir].row[eil] == NULL) {
      A->blocks[rbi][bir].row[eil]  = (re_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(re_t));
      A->blocks[rbi][bir].pos[eil]  = (bi_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(bi_t));
      A->blocks[rbi][bir].buf[eil]  = 2 * __GBLA_SIMD_INNER_SIZE;
    } else {
      if (A->blocks[rbi][bir].sz[eil] == A->blocks[rbi][bir].buf[eil]) {
        A->blocks[rbi][bir].buf[eil]  *=  2;
        A->blocks[rbi][bir].row[eil]  =   realloc(A->blocks[rbi][bir].row[eil],
            A->blocks[rbi][bir].buf[eil] * sizeof(re_t));
        A->blocks[rbi][bir].pos[eil]  =   realloc(A->blocks[rbi][bir].pos[eil],
            A->blocks[rbi][bir].buf[eil] * sizeof(bi_t));
      }
    }
    */
    // set values
    A->blocks[rbi][bir].row[eil][A->blocks[rbi][bir].sz[eil]] = 
      (re_t)M->rows[bi][ri[i]];
    A->blocks[rbi][bir].pos[eil][A->blocks[rbi][bir].sz[eil]] = lib;
    A->blocks[rbi][bir].sz[eil]++;
  }
}

/**
 * \brief Inserts several elements from input matrix M in sparse submatrix A at
 * once with inverted values w.r.t. M->mod.
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param array of column positions pos
 *
 * \param array of row indices ri
 *
 * \param length of pos resp. ri length
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param row index of corresponding element in M bi
 */
static inline void insert_many_in_sb_inv(sb_fl_t *A, const sm_t *M, const ci_t *pos,
    const ci_t *ri, const bi_t length, const ri_t rbi, const ri_t lib, const ri_t bi)
{
  bi_t i, k;
  
  register bi_t eil;
  // bir is the same for all elements
  register const bi_t bir  = pos[0] / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  // allocate memory
  if (A->blocks[rbi][bir].row == NULL) {
    A->blocks[rbi][bir].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
    A->blocks[rbi][bir].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
    A->blocks[rbi][bir].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      A->blocks[rbi][bir].row[k]  = NULL;
      A->blocks[rbi][bir].pos[k]  = NULL;
      A->blocks[rbi][bir].sz[k]   = 0;
      A->blocks[rbi][bir].buf[k]  = 0;
    }
  }
  A->blocks[rbi][bir].row[lib]  = (re_t *)malloc(length * sizeof(re_t));
  A->blocks[rbi][bir].pos[lib]  = (bi_t *)malloc(length * sizeof(bi_t));
  // insert entries in this block
  //printf("length %u | %u | %u | %u\n",length,bir,lib,A->blocks[rbi][bir].sz[lib]);
  for (i=0; i<length; ++i) {
    eil  = pos[i] % __GBLA_SIMD_BLOCK_SIZE; // index in block line
    // set values
    //if (A->blocks[rbi][bir].sz[lib]>=length)
      //printf("!!! %u\n",A->blocks[rbi][bir].sz[lib]);
    A->blocks[rbi][bir].row[lib][length - A->blocks[rbi][bir].sz[lib] - 1] = 
      (re_t)((re_m_t)M->mod - M->rows[bi][ri[i]]);
    A->blocks[rbi][bir].pos[lib][length - A->blocks[rbi][bir].sz[lib] - 1] = eil;
    A->blocks[rbi][bir].sz[lib]++;
  }
}

/**
 * \brief Inserts several elements from input matrix M in sparse submatrix A at
 * once.
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param array of column positions pos
 *
 * \param array of row indices ri
 *
 * \param length of pos resp. ri length
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param row index of corresponding element in M bi
 */
static inline void insert_many_in_sb(sb_fl_t *A, const sm_t *M, const ci_t *pos,
    const ci_t *ri, const bi_t length, const ri_t rbi, const ri_t lib, const ri_t bi)
{
  bi_t i, k;
  
  register bi_t eil;
  // bir is the same for all elements in col
  register const bi_t bir  = pos[0] / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].row == NULL) {
    A->blocks[rbi][bir].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
    A->blocks[rbi][bir].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
    A->blocks[rbi][bir].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      A->blocks[rbi][bir].row[k]  = NULL;
      A->blocks[rbi][bir].pos[k]  = NULL;
      A->blocks[rbi][bir].sz[k]   = 0;
      A->blocks[rbi][bir].buf[k]  = 0;
    }
  }
  A->blocks[rbi][bir].row[lib]  = (re_t *)malloc(length * sizeof(re_t));
  A->blocks[rbi][bir].pos[lib]  = (bi_t *)malloc(length * sizeof(bi_t));
  // insert entries in this block
  for (i=0; i<length; ++i) {
    eil  = pos[i] % __GBLA_SIMD_BLOCK_SIZE; // index in block line
    // set values
    A->blocks[rbi][bir].row[lib][A->blocks[rbi][bir].sz[lib]] = 
      (re_t)M->rows[bi][ri[i]];
    A->blocks[rbi][bir].pos[lib][A->blocks[rbi][bir].sz[lib]] = eil;
    A->blocks[rbi][bir].sz[lib]++;
  }
}

/**
 * \brief Inserts elements from input matrix M in sparse submatrix A.
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_sb(sb_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir  = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil  = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  bi_t i, j, k;
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].row == NULL) {
    A->blocks[rbi][bir].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
    A->blocks[rbi][bir].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
    A->blocks[rbi][bir].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      A->blocks[rbi][bir].row[k]  = NULL;
      A->blocks[rbi][bir].pos[k]  = NULL;
      A->blocks[rbi][bir].sz[k]   = 0;
      A->blocks[rbi][bir].buf[k]  = 0;
    }
  }
  if (A->blocks[rbi][bir].row[lib] == NULL) {
    A->blocks[rbi][bir].row[lib]  = (re_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(re_t));
    A->blocks[rbi][bir].pos[lib]  = (bi_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf[lib]  = 2 * __GBLA_SIMD_INNER_SIZE;
  } else {
    if (A->blocks[rbi][bir].sz[lib] == A->blocks[rbi][bir].buf[lib]) {
      A->blocks[rbi][bir].buf[lib]  *=  2;
      A->blocks[rbi][bir].row[lib]  =   realloc(A->blocks[rbi][bir].row[lib],
          A->blocks[rbi][bir].buf[lib] * sizeof(re_t));
      A->blocks[rbi][bir].pos[lib]  =   realloc(A->blocks[rbi][bir].pos[lib],
          A->blocks[rbi][bir].buf[lib] * sizeof(bi_t));
    }
  }
  // set values
  A->blocks[rbi][bir].row[lib][A->blocks[rbi][bir].sz[lib]] = 
    (re_t)M->rows[bi][ri];
  A->blocks[rbi][bir].pos[lib][A->blocks[rbi][bir].sz[lib]] = eil;
  A->blocks[rbi][bir].sz[lib]++;
}

/**
 * \brief Inserts elements from input matrix M in sparse submatrix A, invert
 * their values w.r.t. to addition and M->mod
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_sb_inv(sb_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir  = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil  = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  bi_t i, j, k;
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].row == NULL) {
    A->blocks[rbi][bir].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
    A->blocks[rbi][bir].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
    A->blocks[rbi][bir].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      A->blocks[rbi][bir].row[k]  = NULL;
      A->blocks[rbi][bir].pos[k]  = NULL;
      A->blocks[rbi][bir].sz[k]   = 0;
      A->blocks[rbi][bir].buf[k]  = 0;
    }
  }
  if (A->blocks[rbi][bir].row[lib] == NULL) {
    A->blocks[rbi][bir].row[lib]  = (re_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(re_t));
    A->blocks[rbi][bir].pos[lib]  = (bi_t *)malloc(2 * __GBLA_SIMD_INNER_SIZE * sizeof(bi_t));
    A->blocks[rbi][bir].buf[lib]  = 2 * __GBLA_SIMD_INNER_SIZE;
  } else {
    // allocate new memory if we are already full
    if (A->blocks[rbi][bir].sz[lib] == A->blocks[rbi][bir].buf[lib]) {
      A->blocks[rbi][bir].buf[lib]  *=  2;
      A->blocks[rbi][bir].row[lib]  =   realloc(A->blocks[rbi][bir].row[lib],
          A->blocks[rbi][bir].buf[lib] * sizeof(re_t));
      A->blocks[rbi][bir].pos[lib]  =   realloc(A->blocks[rbi][bir].pos[lib],
          A->blocks[rbi][bir].buf[lib] * sizeof(bi_t));
    }
  }
  // set values
  A->blocks[rbi][bir].row[lib][A->blocks[rbi][bir].sz[lib]] = 
    (re_t)((re_m_t)M->mod - M->rows[bi][ri]);
  A->blocks[rbi][bir].pos[lib][A->blocks[rbi][bir].sz[lib]] = eil;
  //printf("%u | %u , %u , %u , %u | %u\n",rbi, bir, lib, A->blocks[rbi][bir].row[lib][A->blocks[rbi][bir].sz[lib]], A->blocks[rbi][bir].sz[lib],A->blocks[rbi][bir].pos[lib][A->blocks[rbi][bir].sz[lib]]);
  A->blocks[rbi][bir].sz[lib]++;
}

/**
 * \brief Inserts elements from input matrix M in hybrid submatrix A, invert
 * their values w.r.t. to addition and M->mod
 *
 * \param hybrid block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_hbm_inv(hbm_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  bi_t i, j;
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir] == NULL) {
    A->blocks[rbi][bir] = (dbl_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(dbl_t *));
    for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
      A->blocks[rbi][bir][i]  = NULL;
    }
    A->blocks[rbi][bir][lib]  = (dbl_t *)malloc(__GBLA_SIMD_INNER_BLOCKS_PER_ROW * sizeof(dbl_t));
    for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
      A->blocks[rbi][bir][lib][j].val = NULL;
    }
    A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val = (re_t *)calloc(
        __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
  } else {
    if (A->blocks[rbi][bir][lib] == NULL) {
      A->blocks[rbi][bir][lib]  = (dbl_t *)malloc(__GBLA_SIMD_INNER_BLOCKS_PER_ROW * sizeof(dbl_t));
      for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
        A->blocks[rbi][bir][lib][j].val = NULL;
      }
      A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val = (re_t *)calloc(
          __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
    } else {
      if (A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val  ==  NULL) {
        A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val = (re_t *)calloc(
            __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
      }
    }
  }
  // set values
  A->blocks[rbi][bir][lib][eil/__GBLA_SIMD_INNER_SIZE].val[eil%__GBLA_SIMD_INNER_SIZE]  = 
    (re_t)((re_m_t)M->mod - M->rows[bi][ri]);
}

/**
 * \brief Inserts several elements from input matrix M in temporary dense blocks at
 * once.
 *
 * \param array of dense blocks
 *
 * \param original matrix M
 *
 * \param array of column positions pos
 *
 * \param array of row indices ri
 *
 * \param length of pos resp. ri length
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param row index of corresponding element in M bi
 */
static inline void insert_many_in_dense_blocks(re_t **db, const sm_t *M, const ci_t *pos,
    const ci_t *ri, const bi_t length, const ri_t rbi, const ri_t lib, const ri_t bi)
{
  bi_t i;
  register bi_t eil;
  // bir is the same for all entries added next
  register const bi_t bir = pos[0] / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  // set values
  for (i=0; i<length; ++i) {
    eil = pos[i] % __GBLA_SIMD_BLOCK_SIZE; // index in block line
    db[bir][(lib*__GBLA_SIMD_BLOCK_SIZE)+eil] = M->rows[bi][ri[i]];
  }
}

/**
 * \brief Inserts several elements from input matrix M in dense submatrix A at
 * once.
 *
 * \param sparse block submatrix A
 *
 * \param original matrix M
 *
 * \param array of column positions pos
 *
 * \param array of row indices ri
 *
 * \param length of pos resp. ri length
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param row index of corresponding element in M bi
 */
static inline void insert_many_in_dbm(dbm_fl_t *A, const sm_t *M, const ci_t *pos,
    const ci_t *ri, const bi_t length, const ri_t rbi, const ri_t lib, const ri_t bi)
{
  bi_t i;
  register bi_t eil;
  // bir is the same for all entries added next
  //printf("-pos[0] = %u | %u\n",pos[0], length);
  register const bi_t bir = pos[0] / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  //printf("-bir %u\n",bir);
  //printf("-rbi %u\n",rbi);
  //printf("-lib %u\n",lib);
  // allocate memory if needed, initialized to zero
  //printf("A->blocks[0][4].val = %p | %p\n", A->blocks[0][4].val, A->blocks[0][4]);
  //printf("A->blocks[0][5].val = %p | %p\n", A->blocks[0][5].val, A->blocks[0][5]);
  if (A->blocks[rbi][bir].val == NULL)
    A->blocks[rbi][bir].val = (re_t *)calloc(
        __GBLA_SIMD_BLOCK_SIZE_RECT, sizeof(re_t));
  for (i=0; i<length; ++i) {
    eil = pos[i] % __GBLA_SIMD_BLOCK_SIZE; // index in block line
    // set values
    insert_dense_block_data(A, M, rbi, bir, lib, eil, bi, ri[i]);
  }
}

/**
 * \brief Inserts elements from input matrix M in submatrix A
 *
 * \param dense block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_dbm(dbm_fl_t *A, const sm_t *M, const ci_t shift, const ri_t rbi,
    const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].val == NULL)
    A->blocks[rbi][bir].val = (re_t *)calloc(
        __GBLA_SIMD_BLOCK_SIZE * __GBLA_SIMD_BLOCK_SIZE, sizeof(re_t));
  // set values
  insert_dense_block_data(A, M, rbi, bir, lib, eil, bi, ri);
}

/**
 * \brief Inserts elements from input matrix M in submatrix A in inverse order.
 * The negative inverse w.r.t. M->mod is computed. Diagonalizes memory for
 * diagonal blocks in A.
 *
 * \note To be used only for A, not for C.
 *
 * \param dense block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_dbm_inv_diagonalize(dbm_fl_t *A, const sm_t *M,
    const ci_t shift, const ri_t rbi, const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].val == NULL) {
    // store only upper triangle for diagonal blocks
    if (rbi == bir) {
      A->blocks[rbi][bir].val = (re_t *)calloc(
          __GBLA_SIMD_BLOCK_SIZE_DIAG, sizeof(re_t));
      insert_dense_block_data_diagonalize(A, M, rbi, bir, lib,
          //__GBLA_SIMD_BLOCK_SIZE-1-eil, bi, ri);
          eil, bi, ri);
    } else {
      A->blocks[rbi][bir].val = (re_t *)calloc(
          __GBLA_SIMD_BLOCK_SIZE_RECT, sizeof(re_t));
      insert_dense_block_data_inv(A, M, rbi, bir, lib,
          //__GBLA_SIMD_BLOCK_SIZE-1-eil, bi, ri);
          eil, bi, ri);
    }
  } else {
    if (rbi == bir)
      insert_dense_block_data_diagonalize(A, M, rbi, bir, lib,
          //__GBLA_SIMD_BLOCK_SIZE-1-eil, bi, ri);
          eil, bi, ri);
    else
      insert_dense_block_data_inv(A, M, rbi, bir, lib,
          eil, bi, ri);
          //__GBLA_SIMD_BLOCK_SIZE-1-eil, bi, ri);
  }
}

/**
 * \brief Inserts elements from input matrix M in submatrix A in inverse order.
 * The negative inverse w.r.t. M->mod is computed.
 *
 * \param dense block submatrix A
 *
 * \param original matrix M
 *
 * \param shift to calculate correct coordinates of the corresponding block in A
 * and inside the block itself shift
 *
 * \param current row block index rbi
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
static inline void insert_in_dbm_inv(dbm_fl_t *A, const sm_t *M, const ci_t shift,
    const ri_t rbi, const ri_t lib, const ri_t bi, const ci_t ri)
{
  const bi_t bir = shift / __GBLA_SIMD_BLOCK_SIZE; // block index in block row
  const bi_t eil = shift % __GBLA_SIMD_BLOCK_SIZE; // index in block line
  // allocate memory if needed, initialized to zero
  if (A->blocks[rbi][bir].val == NULL)
    A->blocks[rbi][bir].val = (re_t *)calloc(
        __GBLA_SIMD_BLOCK_SIZE_RECT, sizeof(re_t));
  // write data
  insert_dense_block_data_inv(A, M, rbi, bir, lib,
      eil, bi, ri);
      //__GBLA_SIMD_BLOCK_SIZE-1-eil, bi, ri);
}


/**
 * \brief Inserts elements from input matrix M in multiline rows of A corresponding
 * to the given splicing and mapping precomputed. This is the version for multi
 * line rows inserting one entry in the first field, 0 in the second field.
 *
 * \param multiline sub matrix A
 *
 * \param original matrix M
 *
 * \param current multiline index mli
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
/* static */ /* inline */ void insert_row_data_ml_1_1(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1) ;

/**
 * \brief Inserts elements from input matrix M in multiline rows of A corresponding
 * to the given splicing and mapping precomputed. This is the version for multi
 * line rows inserting one entry in the second field, 0 in the first field.
 *
 * \param multiline sub matrix A
 *
 * \param original matrix M
 *
 * \param current multiline index mli
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi2
 *
 * \param index in row bi1 of corresponding element in M i2
 *
 */
/* static */ /* inline */ void insert_row_data_ml_1_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi2, const ci_t i2) ;

/**
 * \brief Inserts elements from input matrix M in multiline rows of A corresponding
 * to the given splicing and mapping precomputed. This is the version for multi
 * line ml inserting two entries.
 *
 * \param multiline sub matrix A
 *
 * \param original matrix M
 *
 * \param current multiline index mli
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 * \param row index of corresponding element in M bi2
 *
 * \param index in row bi1 of corresponding element in M i2
 *
 */
/* static */ /* inline */ void insert_row_data_ml_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1,
    const ci_t bi2, const ci_t i2) ;


/**
 * \brief Inserts elements from input matrix M in block rows of A corresponding
 * to the given splicing and mapping precomputed. This is the version for multi
 * line blocks inserting one entry in the first field, 0 in the second field.
 *
 * \param block matrix A
 *
 * \param original matrix M
 *
 * \param current row index for blocks rbi
 *
 * \param column index for the block in the current block row bir
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 */
/* static */ /* inline */ void insert_block_row_data_ml_1_1(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1) ;

/**
 * \brief Inserts elements from input matrix M in block rows of A corresponding
 * to the given splicing and mapping precomputed. This is the version for multi
 * line blocks inserting one entry in the second field, 0 in the first field.
 *
 * \param block matrix A
 *
 * \param original matrix M
 *
 * \param current row index for blocks rbi
 *
 * \param column index for the block in the current block row bir
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi2
 *
 * \param index in row bi1 of corresponding element in M i2
 *
 */

/* static */ /* inline */ void insert_block_row_data_ml_1_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi2, const ci_t i2) ;

/**
 * \brief Inserts elements from input matrix M in block rows of A corresponding
 * to the given splicing and mapping precomputed. This is the version for multi
 * line blocks inserting two entries.
 *
 * \param block matrix A
 *
 * \param original matrix M
 *
 * \param current row index for blocks rbi
 *
 * \param column index for the block in the current block row bir
 *
 * \param current line in block lib
 *
 * \param position of the element in line of the block eil
 *
 * \param row index of corresponding element in M bi1
 *
 * \param index in row bi1 of corresponding element in M i1
 *
 * \param row index of corresponding element in M bi2
 *
 * \param index in row bi1 of corresponding element in M i2
 *
 */
/* static */ /* inline */ void insert_block_row_data_ml_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1, const ci_t bi2, const ci_t i2) ;

/**
 * \brief Constructs a mapping for splicing the BD part for the second round of
 * FL reduction when performing a reduced row echelon form computations
 *
 * \param new indexer map map
 *
 * \param indexer map from D map_D
 *
 * \param number of rows of BD rows_BD
 *
 * \param rank of D rank_D
 *
 * \param column dimension of B and D coldim
 *
 * \param number of threads nthreads
 */
void construct_fl_map_reduced(map_fl_t *map, map_fl_t *map_D, ri_t rows_BD,
    ri_t rank_D, ci_t coldim, int nthreads);

/**
 * \brief Constructs an indexer map for a Faugère-Lachartre decomposition of the
 * original matrix M
 *
 * \param original matrix M
 *
 * \param indexer map for M
 */
void construct_fl_map(sm_t *M, map_fl_t *map);


/**
 * \brief Checks dimensions of sub matrices after constructing map and
 * initializing sub matrices
 *
 * \param initial input matrix M
 *
 * \param hybrid block submatrix A
 *
 * \param hybrid block submatrix B
 *
 * \param hybrid block submatrix C
 *
 * \param hybrid block submatrix D
 */
static inline void check_dimensions_hbm(sm_t *M, hbm_fl_t *A, hbm_fl_t *B, hbm_fl_t *C,
    hbm_fl_t *D)
{
  printf("M[%dx%d]\n",M->nrows,M->ncols);
  printf("A[%dx%d]\n",A->nrows,A->ncols);
  printf("B[%dx%d]\n",B->nrows,B->ncols);
  printf("C[%dx%d]\n",C->nrows,C->ncols);
  printf("D[%dx%d]\n",D->nrows,D->ncols);
  assert(A->nrows == A->ncols);
  assert(A->nrows == B->nrows);
  assert(A->ncols == C->ncols);
  assert(B->ncols == D->ncols);
  assert(C->nrows == D->nrows);
  assert(M->nrows == A->nrows + D->nrows);
  assert(M->ncols == B->ncols + C->ncols);
}

/**
 * \brief Checks dimensions of sub matrices after constructing map and
 * initializing sub matrices
 *
 * \param initial input matrix M
 *
 * \param dense block submatrix A
 *
 * \param dense block submatrix B
 *
 * \param dense block submatrix C
 *
 * \param dense block submatrix D
 */
static inline void check_dimensions(sm_t *M, dbm_fl_t *A, dbm_fl_t *B, dbm_fl_t *C,
    dbm_fl_t *D)
{
  printf("M[%dx%d]\n",M->nrows,M->ncols);
  printf("A[%dx%d]\n",A->nrows,A->ncols);
  printf("B[%dx%d]\n",B->nrows,B->ncols);
  printf("C[%dx%d]\n",C->nrows,C->ncols);
  printf("D[%dx%d]\n",D->nrows,D->ncols);
  assert(A->nrows == A->ncols);
  assert(A->nrows == B->nrows);
  assert(A->ncols == C->ncols);
  assert(B->ncols == D->ncols);
  assert(C->nrows == D->nrows);
  assert(M->nrows == A->nrows + D->nrows);
  assert(M->ncols == B->ncols + C->ncols);
}

/**
 * \brief Initializes pivot indices corresponding to the blocks
 * in A,B.
 *
 * \param array storing indices piv_start_idx
 *
 * \param number of pivots npiv
 *
 * \param splicer map privots or non-pivots map_piv
 *
 * \param number of columns in input matrix input_ncols
 *
 * \param number of rows needed for the submatrix splice nrows
 */
static inline void init_pivot_block_start_indices(ri_t **piv_start_idx,
    ri_t *npiv, const ri_t *map_piv, const ci_t range,
    const ri_t nrows)
{
  ri_t *psi = (ri_t *)malloc(
      ((nrows / __GBLA_SIMD_BLOCK_SIZE) + 2) * sizeof(ri_t));
  int i;
  *npiv  = 0;
  // find blocks for construction of A & B
  for (i = (int)range-1; i>-1; --i) {
    if (map_piv[i] != __GB_MINUS_ONE_32) {
      (*npiv)++;
    }
    if ((*npiv % __GBLA_SIMD_BLOCK_SIZE) == 0) {
      psi[*npiv/__GBLA_SIMD_BLOCK_SIZE]  = i;
    }
  }
  // loop might overwrite piv_start_idx[0] with a wrong index;
  // instead of checking "npiv > 0" in each if clause we just reset
  // piv_start_idx[0] after the for loop
  psi[0]  = range;

  // set leftout entries to zero
  for (i=*npiv/__GBLA_SIMD_BLOCK_SIZE+1;
      i < (nrows/__GBLA_SIMD_BLOCK_SIZE) + 2; ++i)
    psi[i] = 0;

  *piv_start_idx  = psi;
}

/**
 * \brief Frees given parts of the input matrix M.
 *
 * \param input matrix M
 *
 * \param rows in horizontal blocks that are already written to submatrices rihb
 *
 * \param current vector in block cvb
 */
static inline void free_input_matrix(sm_t **M_in, const uint32_t *rihb, const uint16_t cvb)
{
  sm_t *M = *M_in;
  int j;
  for (j=0; j<cvb; ++j) {
    free(M->pos[rihb[j]]);
    M->pos[rihb[j]] = NULL;
    free(M->rows[rihb[j]]);
    M->rows[rihb[j]] = NULL;
  }
  *M_in = M;
}

/**
 * \brief Writes corresponding entries of original matrix M into the sparse block
 * submatrix A and the dense block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note This procedure does not invert the entries in A. This is used for
 * constructing C when we keep A.
 *
 * \param original matrix M
 *
 * \param sparse block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_dense_blocks_matrix_no_inversion(const sm_t *M,
    sb_fl_t *A, dbm_fl_t *B, const map_fl_t *map, ri_t *rihb, const ri_t cvb,
    const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_sb(A, M, A->ncols-1-map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the sparse block
 * submatrix A and the dense block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note Data is buffered for a full block, and then written to a row in that
 * block at once.
 *
 * \param original matrix M
 *
 * \param sparse block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_dense_blocks_matrix_many(const sm_t *M, sb_fl_t *A,
    dbm_fl_t *B, const map_fl_t *map, ri_t *rihb, const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  ci_t itA[BUFFER];
  ci_t itB[BUFFER];
  ci_t riA[BUFFER];
  ci_t riB[BUFFER];
  bi_t ctrA, ctrB;
  ci_t old_col_posA;
  ci_t old_col_posB;
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    ctrA  = ctrB  = 0;
    old_col_posA  = old_col_posB  = UINT32_MAX;
    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32) {
        if (old_col_posA != ((A->ncols-1-map->pc[it]) / __GBLA_SIMD_BLOCK_SIZE)) {
          if (ctrA != 0) {
            insert_many_in_sb_inv(A, M, itA, riA, ctrA, rbi, i, bi); 
            ctrA  = 0;
          }
          //printf("ocpA %u | ctrA %u\n",old_col_posA, ctrA);
          itA[ctrA]     = A->ncols-1-map->pc[it];
          old_col_posA  = (A->ncols-1-map->pc[it]) / __GBLA_SIMD_BLOCK_SIZE;
          //printf("-> ocpA %u | itA %u\n",old_col_posA, itA[ctrA]);
          riA[ctrA]     = ri;
          ++ctrA;
          ++ri;
        } else {
          itA[ctrA]  = A->ncols-1-map->pc[it];
          riA[ctrA]  = ri;
          ++ctrA;
          ++ri;
        }
      } else {
        if (old_col_posB != (map->npc[it] / __GBLA_SIMD_BLOCK_SIZE)) {
          if (ctrB != 0) {
            insert_many_in_dbm(B, M, itB, riB, ctrB, rbi, i, bi); 
            ctrB  = 0;
          }
          itB[ctrB]     = map->npc[it];
          old_col_posB  = map->npc[it] / __GBLA_SIMD_BLOCK_SIZE;
          riB[ctrB]     = ri;
          ++ctrB;
          ++ri;
        } else {
          itB[ctrB]  = map->npc[it];
          riB[ctrB]  = ri;
          ++ctrB;
          ++ri;
        }
      }
    }
    // check possible lefovers
    if (ctrA > 0)
      insert_many_in_sb_inv(A, M, itA, riA, ctrA, rbi, i, bi); 
    if (ctrB > 0)
      insert_many_in_dbm(B, M, itB, riB, ctrB, rbi, i, bi); 
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the sparse block
 * submatrix A and the dense block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \param original matrix M
 *
 * \param sparse block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_dense_blocks_matrix(const sm_t *M, sb_fl_t *A,
    dbm_fl_t *B, const map_fl_t *map, ri_t *rihb, const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_sb_inv(A, M, A->ncols-1-map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the sparse block
 * submatrix A and the sparse block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note This procedure does not invert the entries in A. This is used for
 * constructing C when we keep A. Moreover, the blocks in A are also filled by
 * left to right order (instead of right to left order in block version).
 * It buffers entries for one full block row and writes then the whole block row
 * at once.
 *
 * \param original matrix M
 *
 * \param sparse submatrix A (left side)
 *
 * \param sparse block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_sparse_blocks_matrix_keep_A_many(
    const sm_t *M, sm_fl_t *A, sb_fl_t *B, const map_fl_t *map, ri_t *rihb,
    const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  ci_t itA[BUFFER];
  ci_t itB[BUFFER];
  ci_t riA[BUFFER];
  ci_t riB[BUFFER];
  bi_t ctrA, ctrB;
  ci_t old_col_posA;
  ci_t old_col_posB;

  const ri_t clB  = get_number_sparse_col_blocks(B);
  // store entries in B first dense and by row
  re_t **dense_blocks = (re_t **)malloc(clB *
      sizeof(re_t *));
  for (i=0; i<clB; ++i)
    dense_blocks[i] = (re_t *)calloc(__GBLA_SIMD_BLOCK_SIZE_RECT,
      sizeof(re_t));
  // allocate memory for row block splice of B in order to not reallocate too
  // often
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    ctrA  = ctrB  = 0;
    old_col_posA  = old_col_posB  = UINT32_MAX;
    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32) {
        if (old_col_posA != ((map->pc[it]) / __GBLA_SIMD_BLOCK_SIZE)) {
          if (ctrA != 0) {
            insert_many_in_sm(A, M, itA, riA, ctrA, rbi, i, bi); 
            ctrA  = 0;
          }
          itA[ctrA]     = map->pc[it];
          old_col_posA  = (map->pc[it]) / __GBLA_SIMD_BLOCK_SIZE;
          riA[ctrA]     = ri;
          ++ctrA;
          ++ri;
        } else {
          itA[ctrA]  = map->pc[it];
          riA[ctrA]  = ri;
          ++ctrA;
          ++ri;
        }
      } else {
        if (old_col_posB != (map->npc[it] / __GBLA_SIMD_BLOCK_SIZE)) {
          if (ctrB != 0) {
            insert_many_in_dense_blocks(dense_blocks, M, itB, riB, ctrB, rbi, i, bi);
            ctrB  = 0;
          }
          itB[ctrB]     = map->npc[it];
          old_col_posB  = map->npc[it] / __GBLA_SIMD_BLOCK_SIZE;
          riB[ctrB]     = ri;
          ++ctrB;
          ++ri;
        } else {
          itB[ctrB]  = map->npc[it];
          riB[ctrB]  = ri;
          ++ctrB;
          ++ri;
        }
      }
    }
    // check possible lefovers
    if (ctrA > 0)
      insert_many_in_sm(A, M, itA, riA, ctrA, rbi, i, bi); 
    if (ctrB > 0)
      insert_many_in_dense_blocks(dense_blocks, M, itB, riB, ctrB, rbi, i, bi);
  }

  // write data from dense blocks into B
#if __GBLA_COLUMN_B
  for (i=0; i<clB; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      ctrB  = 0;
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
        if (dense_blocks[i][k*__GBLA_SIMD_BLOCK_SIZE + j] != 0)
          ++ctrB;
      }
      if (ctrB > 0) {
        if (B->blocks[rbi][i].row == NULL) {
          B->blocks[rbi][i].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
          B->blocks[rbi][i].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
          B->blocks[rbi][i].sz = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            B->blocks[rbi][i].row[k]  = NULL;
            B->blocks[rbi][i].pos[k]  = NULL;
            B->blocks[rbi][i].sz[k]   = 0;
          }
        }
        B->blocks[rbi][i].row[j]  = (re_t *)malloc(ctrB * sizeof(re_t));
        B->blocks[rbi][i].pos[j]  = (bi_t *)malloc(ctrB * sizeof(bi_t));
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          if (dense_blocks[i][k*__GBLA_SIMD_BLOCK_SIZE + j] != 0) {
            B->blocks[rbi][i].row[j][B->blocks[rbi][i].sz[j]] =
              dense_blocks[i][k*__GBLA_SIMD_BLOCK_SIZE + j];
            B->blocks[rbi][i].pos[j][B->blocks[rbi][i].sz[j]] = k;
            B->blocks[rbi][i].sz[j]++;
          }
        }
      }
    }
  }

#else
  for (i=0; i<clB; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      ctrB  = 0;
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
        if (dense_blocks[i][j*__GBLA_SIMD_BLOCK_SIZE + k] != 0)
          ++ctrB;
      }
      if (ctrB > 0) {
        if (B->blocks[rbi][i].row == NULL) {
          B->blocks[rbi][i].row = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
          B->blocks[rbi][i].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
          B->blocks[rbi][i].sz = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            B->blocks[rbi][i].row[k]  = NULL;
            B->blocks[rbi][i].pos[k]  = NULL;
            B->blocks[rbi][i].sz[k]   = 0;
          }
        }
        B->blocks[rbi][i].row[j]  = (re_t *)malloc(ctrB * sizeof(re_t));
        B->blocks[rbi][i].pos[j]  = (bi_t *)malloc(ctrB * sizeof(bi_t));
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          if (dense_blocks[i][j*__GBLA_SIMD_BLOCK_SIZE + k] != 0) {
            B->blocks[rbi][i].row[j][B->blocks[rbi][i].sz[j]] =
              dense_blocks[i][j*__GBLA_SIMD_BLOCK_SIZE + k];
            B->blocks[rbi][i].pos[j][B->blocks[rbi][i].sz[j]] = k;
            B->blocks[rbi][i].sz[j]++;
          }
        }
      }
    }
  }
#endif

  // free dense_blocks
  for (i=0; i<clB; ++i)
    free(dense_blocks[i]);
  free(dense_blocks);
}


/**
 * \brief Writes corresponding entries of original matrix M into the sparse block
 * submatrix A and the sparse block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note This procedure does not invert the entries in A. This is used for
 * constructing C when we keep A. Moreover, the blocks in A are also filled by
 * left to right order (instead of right to left order in block version)
 *
 * \param original matrix M
 *
 * \param sparse submatrix A (left side)
 *
 * \param sparse block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_sparse_blocks_matrix_keep_A(
    const sm_t *M, sm_fl_t *A, sb_fl_t *B, const map_fl_t *map, ri_t *rihb,
    const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_sm(A, M, map->pc[it], rbi, i, bi, ri); 
      else
#if __GBLA_COLUMN_B
        insert_in_sb_by_column(B, M, map->npc[it], rbi, i, bi, ri); 
#else
      insert_in_sb(B, M, map->npc[it], rbi, i, bi, ri); 
#endif
      ri++;
    }
  }
}

/**
* \brief Writes corresponding entries of original matrix M into the sparse block
* submatrix A and the dense block submatrix B. The entries are defined by the
* mappings from M given by rihb, crb and rbi:
* parts of M --> A|B
*
* \note This procedure does not invert the entries in A. This is used for
* constructing C when we keep A. Moreover, the blocks in A are also filled by
 * left to right order (instead of right to left order in block version)
 *
 * \param original matrix M
 *
 * \param sparse block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_dense_blocks_matrix_keep_A_many(
    const sm_t *M, sm_fl_t *A, dbm_fl_t *B, const map_fl_t *map, ri_t *rihb,
    const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  ci_t itA[BUFFER];
  ci_t itB[BUFFER];
  ci_t riA[BUFFER];
  ci_t riB[BUFFER];
  bi_t ctrA, ctrB;
  ci_t old_col_posA;
  ci_t old_col_posB;
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    ctrA  = ctrB  = 0;
    old_col_posA  = old_col_posB  = UINT32_MAX;
    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32) {
        if (old_col_posA != ((map->pc[it]) / __GBLA_SIMD_BLOCK_SIZE)) {
          if (ctrA != 0) {
            insert_many_in_sm(A, M, itA, riA, ctrA, rbi, i, bi); 
            ctrA  = 0;
          }
          itA[ctrA]     = map->pc[it];
          old_col_posA  = (map->pc[it]) / __GBLA_SIMD_BLOCK_SIZE;
          riA[ctrA]     = ri;
          ++ctrA;
          ++ri;
        } else {
          itA[ctrA]  = map->pc[it];
          riA[ctrA]  = ri;
          ++ctrA;
          ++ri;
        }
      } else {
        if (old_col_posB != (map->npc[it] / __GBLA_SIMD_BLOCK_SIZE)) {
          if (ctrB != 0) {
            insert_many_in_dbm(B, M, itB, riB, ctrB, rbi, i, bi); 
            ctrB  = 0;
          }
          itB[ctrB]     = map->npc[it];
          old_col_posB  = map->npc[it] / __GBLA_SIMD_BLOCK_SIZE;
          riB[ctrB]     = ri;
          ++ctrB;
          ++ri;
        } else {
          itB[ctrB]  = map->npc[it];
          riB[ctrB]  = ri;
          ++ctrB;
          ++ri;
        }
      }
    }
    // check possible lefovers
    if (ctrA > 0)
      insert_many_in_sm(A, M, itA, riA, ctrA, rbi, i, bi); 
    if (ctrB > 0)
      insert_many_in_dbm(B, M, itB, riB, ctrB, rbi, i, bi); 
  }
}

/**
* \brief Writes corresponding entries of original matrix M into the sparse block
* submatrix A and the dense block submatrix B. The entries are defined by the
* mappings from M given by rihb, crb and rbi:
* parts of M --> A|B
*
* \note This procedure does not invert the entries in A. This is used for
* constructing C when we keep A. Moreover, the blocks in A are also filled by
 * left to right order (instead of right to left order in block version)
 *
 * \param original matrix M
 *
 * \param sparse block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_dense_blocks_matrix_keep_A(
    const sm_t *M, sm_fl_t *A, dbm_fl_t *B, const map_fl_t *map, ri_t *rihb,
    const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_sm(A, M, map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the sparse block
 * submatrix A and the dense block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note The blocks in A are filled by left to right order instead of the right
 * to left order used when not keeping A but updating B later on.
 *
 * \param original matrix M
 *
 * \param sparse block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_sparse_dense_blocks_matrix_inv_keep_A(const sm_t *M, sm_fl_t *A,
    dbm_fl_t *B, const map_fl_t *map, ri_t *rihb, const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_sm_inv(A, M, map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the hybrid block
 * submatrix A and the dense block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \param original matrix M
 *
 * \param hybrid block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_hybrid_dense_blocks_matrix(const sm_t *M, hbm_fl_t *A,
    dbm_fl_t *B, const map_fl_t *map, ri_t *rihb, const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  // ususally cvb is divisible by 2, but for the last row of blocks there might
  // be only an odd number of lines in the blocks
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_hbm_inv(A, M, A->ncols-1-map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the hybrid block
 * submatrices A and B. The entries are defined by the mappings from M given by
 * rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \param original matrix M
 *
 * \param hybrid block submatrix A (left side)
 *
 * \param hybrid block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_hybrid_blocks_matrix(const sm_t *M, hbm_fl_t *A, hbm_fl_t *B,
    const map_fl_t *map, ri_t *rihb, const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  // ususally cvb is divisible by 2, but for the last row of blocks there might
  // be only an odd number of lines in the blocks
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_hbm_inv(A, M, A->ncols-1-map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_hbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the dense block
 * submatrices A and B. The entries are defined by the mappings from M given by
 * rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \param original matrix M
 *
 * \param dense block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_dense_blocks_matrix(const sm_t *M, dbm_fl_t *A, dbm_fl_t *B,
    const map_fl_t *map, ri_t *rihb, const ri_t cvb, const ri_t rbi)
{
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  // ususally cvb is divisible by 2, but for the last row of blocks there might
  // be only an odd number of lines in the blocks
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_dbm_inv(A, M, A->ncols-1-map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Writes corresponding entries of original matrix M into the dense block
 * submatrices A and B. The entries are defined by the mappings from M given by
 * rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note Compresses memory for diagonal blocks for A, only useful for A|B
 * splicing, not for C|D splicing.
 *
 * \param original matrix M
 *
 * \param ense block submatrix A (left side)
 *
 * \param dense block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 */
static inline void write_dense_blocks_matrix_diagonalize(const sm_t *M,
    dbm_fl_t *A, dbm_fl_t *B, const map_fl_t *map, ri_t *rihb,
    const ri_t cvb, const ri_t rbi)
{

  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it, ri;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i])
  ri_t i, j, k, l, bi;

  // ususally cvb is divisible by 2, but for the last row of blocks there might
  // be only an odd number of lines in the blocks
  for (i=0; i<cvb; ++i) {
    bi  = rihb[i];
    ri  = 0;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (ri < M->rwidth[bi]) {
      it  = M->pos[bi][ri];
      if (map->pc[it] != __GB_MINUS_ONE_32)
        insert_in_dbm_inv_diagonalize(A, M, A->ncols-1-map->pc[it], rbi, i, bi, ri); 
      else
        insert_in_dbm(B, M, map->npc[it], rbi, i, bi, ri); 
      ri++;
    }
  }
}

/**
 * \brief Fills sparse submatrix A and dense submatrix B with values from M with
 * respect to the splicing stored in map. This version is for keeping A, i.e.
 * the entries inserted to C are not inverted w.r.t. M->mod.
 *
 * \param input matrix M
 *
 * \param left side sparse matrix A
 *
 * \param right side dense block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_sparse_dense_submatrices_keep_A(sm_t *M, sm_fl_t *A,
    dbm_fl_t *B, const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_sparse_dense_blocks_matrix_keep_A_many(
          //write_sparse_dense_blocks_matrix_keep_A(
              M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
  // cut memory down
  cut(A);
}

/**
 * \brief Fills sparse submatrix A and sparse block submatrix B with values
 * from M with respect to the splicing stored in map. This version is for
 * keeping A, i.e. the entries inserted to C are not inverted w.r.t. M->mod.
 *
 * \param input matrix M
 *
 * \param left side sparse matrix A
 *
 * \param right side sparse block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_sparse_sparse_submatrices_keep_A(sm_t *M, sm_fl_t *A,
    sb_fl_t *B, const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

    //#pragma omp for schedule(dynamic) nowait
#pragma omp for
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
#pragma omp task
      {
        // construct block submatrices A & B
        // Note: In the for loop we always construct block "block+1" and not block
        // "block".
        // TODO: Try to improve this rather strange looping.
        for (i = ((int)piv_start_idx[block_idx]-1);
            i > (int)piv_start_idx[block_idx+1]-1; --i) {
          if (range[i] != __GB_MINUS_ONE_32) {
            rihb[cvb] = range[i];
            cvb++;
          }
          if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
            write_sparse_sparse_blocks_matrix_keep_A_many(
            //write_sparse_sparse_blocks_matrix_keep_A(
                M, A, B, map, rihb, cvb, block_idx);

            // TODO: Destruct input matrix on the go
            if (destruct_input_matrix)
              free_input_matrix(&M, rihb, cvb);
            cvb = 0;
          }
        }
      }
    }
#pragma omp taskwait
  }
  // cut memory down
  cut(A);
  //cut_blocks(B);
}

/**
 * \brief Fills sparse submatrix A and dense submatrix B with values from M with
 * respect to the splicing stored in map.
 *
 * \param input matrix M
 *
 * \param left side sparse matrix A
 *
 * \param right side dense matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_sparse_dense_submatrices_inv_keep_A(sm_t *M, sm_fl_t *A, dbm_fl_t *B,
    const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_sparse_dense_blocks_matrix_inv_keep_A(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
  // cut memory down
  cut(A);
}

/**
 * \brief Fills sparse submatrix A and dense submatrix B with values from M with
 * respect to the splicing stored in map. This version is for keeping A, i.e.
 * the entries inserted to C are not inverted w.r.t. M->mod.
 *
 * \param input matrix M
 *
 * \param left side sparse block matrix A
 *
 * \param right side dense block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_sparse_dense_submatrices_no_inversion(sm_t *M, sb_fl_t *A,
    dbm_fl_t *B, const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_sparse_dense_blocks_matrix_no_inversion(
              M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
  // swap entries in A and cut memory at the same time
  swap_and_cut(A);
}

/**
 * \brief Fills sparse submatrix A and dense submatrix B with values from M with
 * respect to the splicing stored in map.
 *
 * \param input matrix M
 *
 * \param left side sparse block matrix A
 *
 * \param right side dense block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_sparse_dense_submatrices(sm_t *M, sb_fl_t *A, dbm_fl_t *B,
    const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_sparse_dense_blocks_matrix_many(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
  // no longer needed if write_sparse_dense_blocks_matrix_many is used
  //swap_and_cut(A);
}

/**
 * \brief Fills hybrid submatrix A and dense submatrix B with values from M with
 * respect to the splicing stored in map.
 *
 * \param input matrix M
 *
 * \param left side hybrid block matrix A
 *
 * \param right side dense block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_hybrid_dense_submatrices(sm_t *M, hbm_fl_t *A, dbm_fl_t *B,
    const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_hybrid_dense_blocks_matrix(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
}

/**
 * \brief Fills hybrid submatrices A and B with values from M with respect to the
 * splicing stored in map.
 *
 * \param input matrix M
 *
 * \param left side hybrid block matrix A
 *
 * \param right side hybrid block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_hybrid_submatrices(sm_t *M, hbm_fl_t *A, hbm_fl_t *B,
    const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_hybrid_blocks_matrix(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
}

/**
 * \brief Fills submatrices A and B with values from M with respect to the
 * splicing stored in map.
 *
 * \param input matrix M
 *
 * \param left side dense block matrix A
 *
 * \param right side dense block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_dense_submatrices(sm_t *M, dbm_fl_t *A, dbm_fl_t *B,
    const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_dense_blocks_matrix(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
}

/**
 * \brief Fills submatrices A and B with values from M with respect to the
 * splicing stored in map. Compresses diagonal blocks for A.
 *
 * \note Only useful for A|B splicing, not for C|D splicing
 *
 * \param input matrix M
 *
 * \param left side dense block matrix A
 *
 * \param right side dense block matrix B
 *
 * \param splicer map map
 *
 * \param range in map, either pivots or non-pivots range
 *
 * \param array storing indices piv_start_idx
 *
 * \param flag for destructing input matrix splices on the fly
 * destruct_input_matrix
 *
 * \param number of threads for parallel computations nthreads
 */
static inline void fill_dense_submatrices_diagonalize(sm_t *M, dbm_fl_t *A, dbm_fl_t *B,
    const map_fl_t *map, const ri_t *range, const ri_t *piv_start_idx,
    const int destruct_input_matrix, const int nthreads)
{
  int i;
  ri_t block_idx;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i) num_threads(nthreads)
  {
    ri_t rihb[__GBLA_SIMD_BLOCK_SIZE];  // rows indices horizontal block
    bi_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= A->nrows/__GBLA_SIMD_BLOCK_SIZE; ++block_idx) {
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping.
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (range[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = range[i];
          cvb++;
        }
        if (cvb == __GBLA_SIMD_BLOCK_SIZE || i == 0) {
          write_dense_blocks_matrix_diagonalize(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix)
            free_input_matrix(&M, rihb, cvb);
          cvb = 0;
        }
      }
    }
  }
}

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 *  \param original matrix M
 *
 *  \param block submatrix A
 *
 *  \param block submatrix B
 *
 *  \param block submatrix C
 *
 *  \param block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param dimension of blocks block_dim
 *
 *  \param number of rows per multiline rows_multiline
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param checks if map was already defined outside map_defined
 */
void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
                      map_fl_t *map, ri_t complete_nrows, ci_t complete_ncols,
                      int block_dim, int rows_multiline,
                      int nthreads, int destruct_input_matrix, int verbose,
                      int map_defined);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 * \note Sparse-dense version without multilines in order to exploit SIMD instructions.
 *
 *  \param original matrix M
 *
 *  \param sparse block submatrix A
 *
 *  \param sparse block submatrix B
 *
 *  \param sparse block submatrix C
 *
 *  \param dense block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param checks if map was already defined outside map_defined
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param level of verbosity
 *
 *  \param number of threads to be used nthreads
 */
void splice_fl_matrix_sparse_dense_keep_A(sm_t *M, sm_fl_t *A, sb_fl_t *B, sm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 * \note Sparse-dense version without multilines in order to exploit SIMD instructions.
 *
 *  \param original matrix M
 *
 *  \param sparse block submatrix A
 *
 *  \param dense block submatrix B
 *
 *  \param sparse block submatrix C
 *
 *  \param dense block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param checks if map was already defined outside map_defined
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param level of verbosity
 *
 *  \param number of threads to be used nthreads
 */
void splice_fl_matrix_sparse_dense_2(sm_t *M, sb_fl_t *A, dbm_fl_t *B, sb_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 * \note Sparse-dense version without multilines in order to exploit SIMD instructions.
 *
 *  \param original matrix M
 *
 *  \param sparse block submatrix A
 *
 *  \param dense block submatrix B
 *
 *  \param dense block submatrix C
 *
 *  \param dense block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param checks if map was already defined outside map_defined
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param level of verbosity
 *
 *  \param number of threads to be used nthreads
 */
void splice_fl_matrix_sparse_dense(sm_t *M, sb_fl_t *A, dbm_fl_t *B, dbm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 * \note Hybrid version without multilines in order to exploit SIMD instructions.
 *
 *  \param original matrix M
 *
 *  \param hybrid block submatrix A
 *
 *  \param dense block submatrix B
 *
 *  \param hybrid block submatrix C
 *
 *  \param dense block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param checks if map was already defined outside map_defined
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param level of verbosity
 *
 *  \param number of threads to be used nthreads
 */
void splice_fl_matrix_hybrid_dense(sm_t *M, hbm_fl_t *A, dbm_fl_t *B, hbm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 * \note Hybrid version without multilines in order to exploit SIMD instructions.
 *
 *  \param original matrix M
 *
 *  \param hybrid block submatrix A
 *
 *  \param hybrid block submatrix B
 *
 *  \param hybrid block submatrix C
 *
 *  \param hybrid block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param checks if map was already defined outside map_defined
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param level of verbosity
 *
 *  \param number of threads to be used nthreads
 */
void splice_fl_matrix_hybrid(sm_t *M, hbm_fl_t *A, hbm_fl_t *B, hbm_fl_t *C,
    hbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 * \note Dense version without multilines in order to exploit SIMD instructions.
 *
 *  \param original matrix M
 *
 *  \param dense block submatrix A
 *
 *  \param dense block submatrix B
 *
 *  \param dense block submatrix C
 *
 *  \param dense block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param checks if map was already defined outside map_defined
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param level of verbosity
 *
 *  \param number of threads to be used nthreads
 */
void splice_fl_matrix_dense(sm_t *M, dbm_fl_t *A, dbm_fl_t *B, dbm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 *  \note Submatrix A is in multiline format.
 *
 *  \param original matrix M
 *
 *  \param multiline submatrix A
 *
 *  \param block submatrix B
 *
 *  \param block submatrix C
 *
 *  \param block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param dimension of blocks block_dim
 *
 *  \param number of rows per multiline rows_multiline
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 */
void splice_fl_matrix_ml_A(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
                      map_fl_t *map, int block_dim, int rows_multiline,
                      int nthreads, int destruct_input_matrix, int verbose);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv
 *
 *  \note Submatrices A and C are in multiline format.
 *
 *  \param original matrix M
 *
 *  \param multiline submatrix A
 *
 *  \param block submatrix B
 *
 *  \param multiline submatrix C
 *
 *  \param block submatrix D
 *
 *  \param indexer mapping map
 *
 *  \param dimension of blocks block_dim
 *
 *  \param number of rows per multiline rows_multiline
 *
 *  \param destructing input matrix on the go? destruct_input_matrix
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 */
void splice_fl_matrix_ml_A_C(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, sm_fl_ml_t *C,
                      sbm_fl_t *D, map_fl_t *map, int block_dim, int rows_multiline,
                      int nthreads, int destruct_input_matrix, int verbose);


/**
 * \brief Writes corresponding entries of original matrix M into the block
 * submatrices A and B. The entries are defined by the mappings from M given by
 * rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \param original matrix M
 *
 * \param block submatrix A (left side)
 *
 * \param block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 *
 * \param level of density of initial buffer for submatrix splicing
 * density_level: If set to 0 then we assume the upper part of the splicing is
 * done, which means a very sparse left part and denser right part. If set to 1
 * then we assume that the lower part of the splicing is done, with a hybrid
 * left part and a very dense right part.
 */
void write_blocks_lr_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, map_fl_t *map,
        ri_t *rihb, const ri_t cvb, const ri_t rbi, const bi_t density_level);

/**
 * \brief Writes corresponding entries of original matrix M into the mutiline
 * submatrix A and the block submatrix B. The entries are defined by the
 * mappings from M given by rihb, crb and rbi:
 * parts of M --> A|B
 *
 * \note The main differene to the function write_blocks_lr_matrix() is that the
 * lefthand side sub matrix A is in multiline row format and not in block format.
 * The handling of the righthand side block sub matrix B is exactly the same.
 *
 * \param original matrix M
 *
 * \param multiline submatrix A (left side)
 *
 * \param block submatrix B (right side)
 *
 * \param splicer mapping map  that stores pivots and non pivots
 *
 * \param row indices in horizonal block rihb
 *
 * \param current row block crb
 *
 * \param row block index rbi
 *
 * \param level of density of initial buffer for submatrix splicing
 * density_level: If set to 0 then we assume the upper part of the splicing is
 * done, which means a very sparse left part and denser right part. If set to 1
 * then we assume that the lower part of the splicing is done, with a hybrid
 * left part and a very dense right part.
 */
void write_lr_matrix_ml(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, map_fl_t *map,
        ri_t *rihb, const ri_t cvb, const ri_t rbi, const bi_t density_level);
#endif

/* vim:sts=2:sw=2:ts=2:
 */
