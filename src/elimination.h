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

    if (sparse_block[i].dense == 0) {
      for (j=0; j<sparse_block[i].sz; ++j) {
        idx   = sparse_block[i].idx[j];
        dense_block[2*i][idx]   = sparse_block[i].val[2*j];
        dense_block[2*i+1][idx] = sparse_block[i].val[2*j+1];
      }
    } else{
      for (j=0; j<bwidth; ++j) {
        dense_block[2*i][j]   = sparse_block[i].val[2*j];
        dense_block[2*i+1][j] = sparse_block[i].val[2*j+1];
      }
    }
  }
}

/**
 * \brief Computes a dense AXPY for one row
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline block in B
 *
 * \param line index line_idx
 *
 * \param dense block width bwidth
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
static inline void dense_scal_mul_sub_1_row_vect_array(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const mbl_t multiline,
              const uint32_t line_idx,
              const bi_t  bwidth,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {
  const re_t *p_val = &(multiline.val[0]);
  p_val +=  line_idx;
  uint32_t i;
  bi_t j;

  register uint32_t v__;
  register uint32_t idx;

  // both cannot be zero at the same time
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
      for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
        v__ = p_val[2*(i+j)];

        dense_val1[i+j] +=  v__ * Av1_col1;
        dense_val2[i+j] +=  v__ * Av2_col1;
      }
    }
  } else { // one of them is zero
    if (Av1_col1 != 0) {
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
        for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
          v__ = p_val[2*(i+j)];

          dense_val1[i+j] +=  v__ * Av1_col1;
        }
      }
    } else {
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
        for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
          v__ = p_val[2*(i+j)];

          dense_val1[i+j] +=  v__ * Av1_col1;
        }
      }
    }
  }
}

/**
 * \brief Computes a sparse AXPY for one row
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline block in B
 *
 * \param line index line_idx
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
static inline void sparse_scal_mul_sub_1_row_vect_array(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const mbl_t multiline,
              const uint32_t line_idx,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {
  const uint32_t N  = multiline.sz;
  const bi_t *p_idx = (N != 0) ? &(multiline.idx[0]) : NULL;
  const re_t *p_val = (N != 0) ? &(multiline.val[0]) : NULL;
  p_val +=  line_idx;
  uint32_t i;

  register uint32_t v__;
  register uint32_t idx;

  // both cannot be zero at the same time
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=0; i<N; ++i) {
      idx = p_idx[i];
      v__ = p_val[2*i];

      dense_val1[idx] +=  v__ * Av1_col1;
      dense_val2[idx] +=  v__ * Av2_col1;
    }
  } else { // one of them is zero
    if (Av1_col1 != 0) {
      for (i=0; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val1[idx] +=  v__ * Av1_col1;
      }
    } else {
      for (i=0; i<N; ++i) {
        idx = p_idx[i];
        v__ = p_val[2*i];

        dense_val2[idx] +=  v__ * Av2_col1;
      }
    }
  }
}

/**
 * \brief Computes a dense AXPY for two rows
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline block in B
 *
 * \param dense block width bwidth
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
static inline void dense_scal_mul_sub_2_rows_vect_array(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const uint32_t Av1_col2,
              const uint32_t Av2_col2,
              const mbl_t multiline,
              const bi_t  bwidth,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {
  // check cases where one pair of the elements is zero
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
  
  const re_t *p_val = &(multiline.val[0]);
  uint32_t i;
  bi_t j;

  register uint32_t v1__, v2__;
  register uint32_t idx;

  for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      v1__ = p_val[2*(i+j)];
      v2__ = p_val[2*(i+j)+1];

      dense_val1[i+j] +=  v1__ * Av1_col1;
      dense_val1[i+j] +=  v2__ * Av1_col2;

      v1__  *=  Av2_col1;
      v2__  *=  Av2_col2;

      dense_val2[i+j] +=  v1__;
      dense_val2[i+j] +=  v2__;
    }
  }
}

/**
 * \brief Computes a sparse AXPY for two rows
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline block in B
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
static inline void sparse_scal_mul_sub_2_rows_vect_array(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const uint32_t Av1_col2,
              const uint32_t Av2_col2,
              const mbl_t multiline,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {
  // check cases where one pair of the elements is zero
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
  const bi_t *p_idx   = (N != 0) ? &(multiline.idx[0]) : NULL;
  const re_t *p_val   = (N != 0) ? &(multiline.val[0]) : NULL;
  const re_t *p_val2  = p_val + 1;
  uint32_t i;
  bi_t j;

  register uint32_t v1__, v2__;
  register uint32_t idx;

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
    idx   = p_idx[i+j];
    v1__  = p_val[2*(i+j)];
    v2__  = p_val2[2*(i+j)];

    dense_val1[idx] +=  Av1_col1 * v1__;
    dense_val1[idx] +=  Av1_col2 * v2__;

    dense_val2[idx] +=  Av2_col1 * v1__;
    dense_val2[idx] +=  Av2_col2 * v2__;
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
 *
 * \param invert scalars? inv_scalars
 *
 * \param characteristic of underlying field modulus
 */
void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight, mod_t modulus, int inv_scalars);

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
 *
 * \param characteristic of underlying field modulus
 */
void red_with_triangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight, mod_t modulus);

/**
 * \brief Elimination procedure which reduces the block submatrix A to the unit
 * matrix. Corresponding changes in block submatrix B are carried out, too.
 *
 * \param block submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_block(sbm_fl_t *A, sbm_fl_t *B, mod_t modulus, int nthrds);

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
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A, mod_t modulus);

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
