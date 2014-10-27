/**
 * \file elimination.h
 * \brief Different Gaussian Elimination methods
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_ELIMINATION_H
#define GB_ELIMINATION_H

#include <mapping.h>
#include <matrix.h>

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
 * \brief Copies entry from dense block to sparse rows
 * for elimination purposes.
 *
 * \param dense block dense_block
 *
 * \param sparse block sparse_block
 *
 * \param block height bheight
 *
 * \param block width bwidth
 *
 * \param modulus resp. field characteristic modulus
 */
mbl_t *copy_dense_block_to_sparse(
    re_l_t **dense_block, mbl_t *sparse_block, int bheight, int bwidth,
    mod_t modulus);

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
              const bi_t line_idx,
              const bi_t  bwidth,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {

  const re_t *p_val = multiline.val;
  p_val +=  line_idx;
  uint32_t i;
  bi_t j;

  register uint32_t v__;

  // both cannot be zero at the same time
  //printf("mlsize %d\n",multiline.sz);
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
      for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
        v__ = p_val[2*(i+j)];

        dense_val1[i+j] +=  v__ * Av1_col1;
        dense_val2[i+j] +=  v__ * Av2_col1;
      }
    }
  } else { // one of them is zero
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
              const bi_t line_idx,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {

  const uint32_t N  = multiline.sz;
  //printf("mlsz %d\n",multiline.sz);
  const bi_t *p_idx = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val = (N != 0) ? multiline.val : NULL;
  p_val +=  line_idx;
  uint32_t i  = 0;

  register uint32_t v__;
  register uint32_t idx;

  //printf("a11 %d -- a21 %d\n",Av1_col1,Av2_col1);
  // both cannot be zero at the same time
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (; i<N; ++i) {
      idx = p_idx[i];
      v__ = p_val[2*i];

      dense_val1[idx] +=  v__ * Av1_col1;
      dense_val2[idx] +=  v__ * Av2_col1;
    }
  } else { // one of them is zero
    if (Av1_col1 != 0) {
      //printf("i %d -- N %d\n",i,N);
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
  
  const re_t *p_val = multiline.val;
  uint32_t i;
  bi_t j;

  register uint32_t v1__, v2__;
  register uint32_t idx;

  for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      //printf(";;%d..%d::",i,j);
      v1__ = p_val[2*(i+j)];
      v2__ = p_val[2*(i+j)+1];
        //printf("v11 %d v21 %d !! ",v1__,v2__);

      dense_val1[i+j] +=  v1__ * Av1_col1;
      dense_val1[i+j] +=  v2__ * Av1_col2;

      v1__  *=  Av2_col1;
      v2__  *=  Av2_col2;
        //printf("v12 %d v22 %d !! ",v1__,v2__);

      dense_val2[i+j] +=  v1__;
      dense_val2[i+j] +=  v2__;
    }
  }
  //printf("\n");
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
  const bi_t *p_idx   = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val   = (N != 0) ? multiline.val : NULL;
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
    idx   = p_idx[i];
    v1__  = p_val[2*i];
    v2__  = p_val2[2*i];

    dense_val1[idx] +=  Av1_col1 * v1__;
    dense_val1[idx] +=  Av1_col2 * v2__;

    dense_val2[idx] +=  Av2_col1 * v1__;
    dense_val2[idx] +=  Av2_col2 * v2__;
  }
}

/**
 * \brief Computes a dense AXPY for one dense row (in triangular A^-1B situation)
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param block width bwidth
 *
 * \param dense source array dense_array_source
 *
 * \param dense array 1 holder for delayed modulus dense_array1
 *
 * \param dense array 2 holder for delayed modulus dense_array2
 */
static inline void dense_scal_mul_sub_1_row_array_array(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const bi_t bwidth,
              const re_l_t *dense_array_source,
              re_l_t *dense_array1,
              re_l_t *dense_array2) {
 
  bi_t i, j;
  register uint32_t v__;

  if (Av1_col1 == 0) { // only one of them can be zero
    for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
      for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
        v__ = dense_array_source[i+j] & 0x000000000000ffff;
        dense_array2[i+j] +=  v__ * Av2_col1;
      }
    }
  } else {
    if (Av2_col1 == 0) { // only second one is zero
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
        for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
          v__ = dense_array_source[i+j] & 0x000000000000ffff;
          dense_array1[i+j] +=  v__ * Av1_col1;
        }
      }
    } else { // both are nonzero
      for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
        for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
          v__ = dense_array_source[i+j] & 0x000000000000ffff;
          dense_array1[i+j] +=  v__ * Av1_col1;
          dense_array2[i+j] +=  v__ * Av2_col1;
        }
      }
    }
  }
}

/**
 * \brief Computes a dense AXPY for two dense rows (in triangular A^-1B situation)
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param block width bwidth
 *
 * \param dense source array dense_array_source1
 *
 * \param dense source array dense_array_source2
 *
 * \param dense array 1 holder for delayed modulus dense_array1
 *
 * \param dense array 2 holder for delayed modulus dense_array2
 */
static inline void dense_scal_mul_sub_2_rows_array_array(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const uint32_t Av1_col2,
              const uint32_t Av2_col2,
              const bi_t bwidth,
              const re_l_t *dense_array_source1,
              const re_l_t *dense_array_source2,
              re_l_t *dense_array1,
              re_l_t *dense_array2) {

  // check cases where one pair of the elements is zero
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

  bi_t i, j;
  register uint32_t v1__, v2__;

  for (i=0; i<bwidth; i+=__GB_LOOP_UNROLL_BIG) {
    for (j=0; j<__GB_LOOP_UNROLL_BIG; ++j) {
      v1__  = dense_array_source1[i+j] & 0x000000000000ffff;
      v2__  = dense_array_source2[i+j] & 0x000000000000ffff;
      
      dense_array1[i+j] +=  v1__ * Av1_col1;
      dense_array1[i+j] +=  v2__ * Av1_col2;
      
      dense_array2[i+j] +=  v1__ * Av2_col1;
      dense_array2[i+j] +=  v2__ * Av2_col2;
    }
  }
}

/**
 * \brief Modular reduction of dense row array
 *
 * \param dense row array to be reduced dense_array
 *
 * \param block width bwidth
 *
 * \param modulus resp. field characteristic modulus
 */
static inline void red_dense_array_modular(re_l_t *dense_array, bi_t bwidth, mod_t modulus) {
  bi_t i;
  for (i=0; i<bwidth; ++i) {
    dense_array[i]  = (re_l_t)(dense_array[i] % modulus);
  }
}

/**
 * \brief Gets first nonzero entry in multiline m at line index line_idx and
 * stores it in h resp. h_idx.
 *
 * \param multiline m
 *
 * \param line index line_idx
 *
 * \param storage for "head" of line h
 *
 * \param storage for "head index" of line h_idx
 *
 * \return index value corresponding to h from input matrix M
 */
mli_t get_head_multiline(const ml_t *m, const bi_t line_idx, re_t *h, mli_t *h_idx) {
  mli_t i;
  for (i=0; i<m.sz; ++i) {
    if (*h = m.val[2*i+line_idx] != 0) {
      *h_idx  = i;
      return m.idx[i];
    }
  }
  return -1;
}

/**
 * \brief Normalizes multiline vector.
 *
 * \param multiline m
 *
 * \param field characteristic modulus
 */
void normalize_multiline(ml_t *m, mod_t modulus) {
  mli_t idx;
  re_t h1 = 0, h2 = 0;

  if (m.sz == 0)
    return;

  get_head_multiline(m, 0, &h1, &idx);
  get_head_multiline(m, 0, &h2, &idx);

  // invert values modulo modulus
  inverse_val(&h1, modulus);
  inverse_val(&h2, modulus);

  // skip if both are 1 and/or 0
  if ((h1 == 0 || h1 == 1) && (h2 == 0 || h2 == 1))
    return;

  re_l_t tmp_val;
  // normalize h2
  if (h1 == 0 || h1 == 1) {
    for (idx=0; idx<m.sz; ++idx) {
      tmp_val         = (re_l_t)m.val[2*idx+1] * h2;
      m.val[2*idx+1]  = tmp_val % modulus;
    }
  } else {
    // normalize h1
    if (h2 == 0 || h2 == 1) {
      for (idx=0; idx<m.sz; ++idx) {
        tmp_val       = (re_l_t)m.val[2*idx] * h1;
        m.val[2*idx]  = tmp_val % modulus;
      }
    // normalize h1 and h2
    } else {
      for (idx=0; idx<m.sz; ++idx) {
        tmp_val       = (re_l_t)m.val[2*idx] * h1;
        m.val[2*idx]  = tmp_val % modulus;
        tmp_val         = (re_l_t)m.val[2*idx+1] * h2;
        m.val[2*idx+1]  = tmp_val % modulus;
      }
    }
  }
}

/**
 * \brief Computes inverse value of x modulo y:
 * We compute the inverse using the extended GCD. So we are only interested in x.
 * Note that internally we need signed types, but we return only unsigned type
 * re_t for x.
 *
 * \param x
 *
 * \param y
 */

void inverse_val(re_t *x, const re_t y) {
  int32_t u1 = 1, u2 = 0;
  int32_t v1 = 0, v3 = y;
  int32_t u3 = x, v2 = 1;
  while (v3 != 0) {
    int32_t q  = u3 / v3;
    int32_t t1 = u1 - v1 * q;
    u1  = v1; v1  = t1;

    int32_t t3 = u3 - v3 * q;
    u3  = v3; v3  = t3;

    int32_t t2 = u2 - v2 * q;
    u2  = v2; v2  = t2;
  }
  if (u1<0) {
    u1  +=  y;
    *x  =   u1;
    return;
  }
  if (u1 > y) {
    u1  -=  y;
    *x  =   u1;
    return;
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
 * \param characteristic of underlying field modulus
 *
 * \param invert scalars? inv_scalars
 */
void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight, mod_t modulus, int inv_scalars);

/**
 * \brief Reduces dense block block_B with triangular sparse block block_A.
 *
 * \param sparse block block_A
 *
 * \param dense block dense_B
 *
 * \param block height bheight
 *
 * \param characteristic of underlying field modulus
 *
 * \param invert scalars? inv_scalars
 */
void red_with_triangular_block(mbl_t *block_A, re_l_t **dense_B, ri_t bheight, mod_t modulus, int inv_scalars);

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
 * \brief Elimination procedure which reduces the block submatrix C to zero.
 * Corresponding changes in block submatrix D are carried out using B, too.
 *
 * \param block submatrix B (right upper side)
 *
 * \param block submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D, mod_t modulus, int nthrds);

/**
 * \brief Different block tasks when reducing block submatrix C.
 *
 * \param block submatrix B (right upper side)
 *
 * \param block submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \param column index of blocks in D block_col_idx_D
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
    ci_t block_col_idx_D, ri_t nbrows_C, ci_t nbcols_C, mod_t modulus);

/**
 * \brief Elimination procedure which reduces the block submatrix D to an
 * upper triangular matrix. Note that the input matrix D will be removed
 * later on and the output matrix D_red is in multiline format for further
 * reduction steps.
 *
 * \param block submatrix D (right lower side), input matrix
 *
 * \param multiline submatrix D_red, output matrix
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, mod_t modulus, int nthrds);

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
