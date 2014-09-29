/**
 * \file mapping.h
 * \brief Interfaces for sparse matrix index maps used for subdividing the matrix
 * in Faugère-Lachartre style
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_MAPPING_H
#define GB_MAPPING_H

#include <gb_config.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <matrix.h>
#include <math.h>
#include <omp.h>

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
 * \param original matrix M
 *
 * \param map to be initialized
 *
 */
static inline void init_fl_map(sm_t *M, map_fl_t *map) {
  // initialize map arrays and 
  // set initial values to __GB_MINUS_ONE_8
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
  
  map->npri = (ri_t *)malloc(M->ncols * sizeof(ri_t));
  memset(map->npri, __GB_MINUS_ONE_8, M->ncols * sizeof(ri_t));
}

/**
 * \brief Reallocates memory for the rows of the multiline during the splicing of
 * the input matrix. The buffer size bufferA is doubled during this process
 *
 * \param multiline sub matrix A
 *
 * \param multiline index mli
 *
 * \param buffer size bufferA
 *
 */
static inline void realloc_rows_ml(sm_fl_ml_t *A, const mli_t mli,
    const bi_t init_bufferA, mli_t *bufferA) {
  *bufferA +=  init_bufferA;
  A->ml[mli].idx = realloc(A->ml[mli].idx, (*bufferA) * sizeof(mli_t));
  A->ml[mli].val = realloc(A->ml[mli].val, 2 * (*bufferA) * sizeof(re_t));
}

/**
 * \brief Reallocates memory for the rows of the blocks during the splicing of
 * the input matrix. The buffer size bufferA is doubled during this process
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
 * \param buffer size bufferA
 *
 */
static inline void realloc_block_rows(sbm_fl_t *A, const ri_t rbi, const ci_t bir,
    const bi_t lib, const bi_t init_bufferA, bi_t *bufferA) {
  *bufferA +=  init_bufferA;
  A->blocks[rbi][bir][lib].idx = realloc(
      A->blocks[rbi][bir][lib].idx,
      (*bufferA) * sizeof(bi_t));
  A->blocks[rbi][bir][lib].val = realloc(
      A->blocks[rbi][bir][lib].val,
      2 * (*bufferA) * sizeof(re_t));
}

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
static inline void swap_block_data(sbm_fl_t *A, const ci_t clA, const bi_t rbi,
    const bi_t cvb) {
  int i, j, k, l;
  bi_t *old_idx_ptr = NULL;
  re_t *old_val_ptr = NULL;

  bi_t rounded_cvb_half = cvb/2;
  if (cvb % 2) // odd cvb
    rounded_cvb_half  +=  1;
  for (i=0; i<rounded_cvb_half; ++i) {
    // allocate memory for swapping data
    bi_t *tmp_idx_ptr = (bi_t *)malloc(A->blocks[rbi][clA-1][i].sz * sizeof(bi_t));
    re_t *tmp_val_ptr = (re_t *)malloc(2 * A->blocks[rbi][clA-1][i].sz * sizeof(re_t));
    // swap data
    l = 0;
    for (k=A->blocks[rbi][clA-1][i].sz-1; k>-1; --k) {
      tmp_idx_ptr[l]      = A->blocks[rbi][clA-1][i].idx[k];
      tmp_val_ptr[2*l]    = A->blocks[rbi][clA-1][i].val[2*k];
      tmp_val_ptr[2*l+1]  = A->blocks[rbi][clA-1][i].val[2*k+1];
      l++;
    }
    // keep track of old ptr
    old_idx_ptr = A->blocks[rbi][clA-1][i].idx;
    old_val_ptr = A->blocks[rbi][clA-1][i].val;
    // swap starting ptrs
    A->blocks[rbi][clA-1][i].idx  = tmp_idx_ptr;
    A->blocks[rbi][clA-1][i].val  = tmp_val_ptr;
    // try to reuse old, already allocated memory in the following
    tmp_idx_ptr = old_idx_ptr;
    tmp_val_ptr = old_val_ptr;
    // now go through all remaining blocks in the corresponding row (bir)
    // check if the already allocated memory is enough for executing the swap
    for (j=clA-2; j>-1; --j) {
      if (A->blocks[rbi][j+1][i].sz == A->blocks[rbi][j][i].sz) {
        //  Memory is not enough, so allocate new (reallocation is used
        //  at the moment: It might be prolematic since the old entries are
        //  useless and need not be copied, but often the already allocated memory
        //  can be extended in place, thus this operation is often faster than
        //  freeing memory and allocating new.
        //  Note that this reallocation is also used if the old memory allocated
        //  is bigger than what is needed. In this setting realloc just cuts the
        //  useless memory off.
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
      }
      // keep track of old ptr
      old_idx_ptr = A->blocks[rbi][j][i].idx;
      old_val_ptr = A->blocks[rbi][j][i].val;
      // swap starting ptrs
      A->blocks[rbi][j][i].idx  = tmp_idx_ptr;
      A->blocks[rbi][j][i].val  = tmp_val_ptr;
      // try to reuse old, already allocated memory in the following
      tmp_idx_ptr = old_idx_ptr;
      tmp_val_ptr = old_val_ptr;
    }
    // finally remove temporary memory used for swapping
    free(old_idx_ptr);
    free(old_val_ptr);
  }
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
static inline void insert_row_data_ml_1_1(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = M->rows[bi1][i1];
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = 0;
  A->ml[mli].sz++;
}

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
static inline void insert_row_data_ml_1_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi2, const ci_t i2) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = 0;
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = M->rows[bi2][i2];
  A->ml[mli].sz++;
}

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
static inline void insert_row_data_ml_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1,
    const ci_t bi2, const ci_t i2) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = M->rows[bi1][i1];
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = M->rows[bi2][i2];
  A->ml[mli].sz++;
}


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
static inline void insert_block_row_data_ml_1_1(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = M->rows[bi1][i1];
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = 0;
  A->blocks[rbi][bir][lib].sz++;
}

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

static inline void insert_block_row_data_ml_1_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi2, const ci_t i2) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = 0;
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = M->rows[bi2][i2];
  A->blocks[rbi][bir][lib].sz++;
}
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
static inline void insert_block_row_data_ml_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1, const ci_t bi2, const ci_t i2) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = M->rows[bi1][i1];
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = M->rows[bi2][i2];
  A->blocks[rbi][bir][lib].sz++;
}

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
 */
void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
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
 */
void write_blocks_lr_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, map_fl_t *map,
                            ri_t *rihb, const ri_t cvb, const ri_t rbi);


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
 */
void write_lr_matrix_ml(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, map_fl_t *map,
                            ri_t *rihb, const ri_t cvb, const ri_t rbi);
#endif
