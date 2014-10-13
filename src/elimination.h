/**
 * \file elimination.h
 * \brief Different Gaussian Elimination methods
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_ELIMINATION_H
#define GB_ELIMINATION_H

#include <mapping.h>

/**
 * \brief Copies entry from sparse block sparse_block to dense block dense_block
 * for elimination purposes.
 *
 * \param sparse block sparse_block
 *
 * \param dense block dense_block
 *
 * \param block height bheight
 *
 * \param block width bwidth
 */
static inline void copy_sparse_to_dense_block(
    mbl_t *sparse_block, re_l_t **dense_block, int bheight, int bwidth) {
  bi_t i,j;
  for (i=0; i<bheight/2; ++i) {
    mli_t idx;
    re_t val1, val2;

    if (sparse_block[i].sz == 0)
      continue;

    for (j=0; j<sparse_block[i].sz; ++j) {
      idx   = sparse_block[i].idx[j];
      dense_block[2*i][idx] = sparse_block[i].val[2*j];
      dense_block[2*i+1][idx] = sparse_block[i].val[2*j+1];
    }
  }
}

/**
 * \brief Reduces dense block block_B with rectangular sparse block block_A.
 *
 * \param sparse block block_A
 *
 * \param sparse block block_B
 *
 * \param dense block dense_B
 *
 * \param block height bheight
 */
void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight);

/**
 * \brief Reduces dense block block_B with triangular sparse block block_A.
 *
 * \param sparse block block_A
 *
 * \param sparse block block_B
 *
 * \param dense block dense_B
 *
 * \param block height bheight
 */
void red_with_triangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight);

/**
 * \brief Elimination procedure which reduces the block submatrix A to the unit
 * matrix. Corresponding changes in block submatrix B are carried out, too.
 *
 * \param block submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_block(sbm_fl_t *A, sbm_fl_t *B, int nthrds);

/**
 * \brief Different block tasks when reducing block submatrix A.
 *
 * \param block submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param column index of blocks in B block_col_idx_B
 *
 * \param number of block rows in A nblock_rows_A
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A);

/**
 * \brief Elimination procedure which reduces the multiline submatrix A
 * and the block submatrices B, C and D.
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_ml_block(sm_fl_ml_t *A, sbm_fl_t *B, int nthrds);

#endif
