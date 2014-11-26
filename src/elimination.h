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
 * \brief Structure of a waiting list element
 */
typedef struct wle_t {
  ri_t idx; /*!<  row index */
  ri_t lp;  /*!<  row index of last pivot that idx was reduced with */
} wle_t;

/**
 * \brief Structure defining the waiting list or queue for parallel dense
 * echelonization of the D part.
 *
 * \note We use a mutiline concept for the row index lists.
 */
typedef struct wl_t {
  wle_t *list;  /*!<  row indices of rows in waiting list */
  ri_t sidx;    /*!<  smallest row index in rows */
  ri_t slp;     /*!<  last pivot row index by which sidx is
                      reduced by already */
  ri_t sz;      /*!<  size of waiting list */
} wl_t;

/**
 * \brief Comparison function for qsort for waiting list
 *
 * \param waiting list element a
 *
 * \param waiting list element b
 *
 * \return *b.idx - *a.idx
 */
static inline int cmp_wle(const void *a, const void *b) {
  return ((wle_t *)(b))->idx - ((wle_t *)(a))->idx;
}

/**
 * \brief Adds a new row to end of waiting list
 *
 * \param waiting list wl
 *
 * \param row index to be added row_idx
 *
 * \param row index of last pivot row_idx was reduced by last_piv_reduced_by
 */
static inline void push_row_to_waiting_list(wl_t *wl, const ri_t row_idx,
    const ri_t last_piv_reduced_by) {
  wl->list[wl->sz].idx  = row_idx;
  wl->list[wl->sz].lp   = last_piv_reduced_by;
  wl->sz++;
}

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
void copy_dense_block_to_sparse(
    re_l_t **dense_block, mbl_t **sparse_block, int bheight, int bwidth,
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
 * \brief Computes a dense AXPY for one row using full multilines
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline
 *
 * \param line index line_idx
 *
 * \param dense block width bwidth
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 *
 * \param offset for first nonzero coefficient in multiline
 */
static inline void dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
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

  register uint32_t v__;

  // both cannot be zero at the same time
  //printf("mlsize %d\n",multiline.sz);
  if (Av1_col1 != 0 && Av2_col1 != 0) {
    for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
      for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
        v__ = p_val[2*(i+j)];

        dense_val1[i+j] +=  v__ * Av1_col1;
        dense_val2[i+j] +=  v__ * Av2_col1;
      }
    }
    for (; i<multiline.sz; ++i) {
      v__ = p_val[2*i];

      dense_val1[i] +=  v__ * Av1_col1;
      dense_val2[i] +=  v__ * Av2_col1;
    }
  } else { // one of them is zero
    if (Av1_col1 == 0) {
      for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
        for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
          v__ = p_val[2*(i+j)];

          dense_val2[i+j] +=  v__ * Av2_col1;
        }
      }
      for (; i<multiline.sz; ++i) {
        v__ = p_val[2*i];

        dense_val1[i] +=  v__ * Av1_col1;
        dense_val2[i] +=  v__ * Av2_col1;
      }
    } else {
      if (Av2_col1 == 0) {
        for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
          for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
            v__ = p_val[2*(i+j)];

            dense_val1[i+j] +=  v__ * Av1_col1;
          }
        }
        for (; i<multiline.sz; ++i) {
          v__ = p_val[2*i];

          dense_val1[i] +=  v__ * Av1_col1;
          dense_val2[i] +=  v__ * Av2_col1;
        }
      }
    }
  }
}

/**
 * \brief Computes a dense AXPY for two rows for multilines
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline 
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 *
 * \param offset for first nonzero coefficient in multiline
 */
static inline void dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const uint32_t Av1_col2,
              const uint32_t Av2_col2,
              const ml_t multiline,
              re_l_t *dense_val1,
              re_l_t *dense_val2,
              const ci_t offset) {

  if (offset < 0)
    return;

  // check cases where one pair of the elements is zero
  if (Av1_col1 == 0 && Av2_col1 == 0) {
    dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
        Av1_col2, Av2_col2, multiline, 1, dense_val1, dense_val2, offset);
    return;
  }
  if (Av1_col2 == 0 && Av2_col2 == 0) {
    dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
        Av1_col1, Av2_col1, multiline, 0, dense_val1, dense_val2, offset);
    return;
  }
  
  const re_t *p_val = multiline.val;
  const uint32_t outer_loop = multiline.sz - __GB_LOOP_UNROLL_SMALL;
  uint32_t i;
  bi_t j;

  register uint32_t v1__, v2__;
  register uint32_t idx;

  for (i=offset; i<outer_loop; i+=__GB_LOOP_UNROLL_SMALL) {
    for (j=0; j<__GB_LOOP_UNROLL_SMALL; ++j) {
      v1__ = p_val[2*(i+j)];
      v2__ = p_val[2*(i+j)+1];

      dense_val1[i+j] +=  v1__ * Av1_col1;
      dense_val1[i+j] +=  v2__ * Av1_col2;


      dense_val2[i+j] +=  v1__ * Av2_col1;
      dense_val2[i+j] +=  v2__ * Av2_col2;
    }
  }
  for (; i<multiline.sz; ++i) {
    v1__ = p_val[2*i];
    v2__ = p_val[2*i+1];

    dense_val1[i] +=  v1__ * Av1_col1;
    dense_val1[i] +=  v2__ * Av1_col2;

    dense_val2[i] +=  v1__ * Av2_col1;
    dense_val2[i] +=  v2__ * Av2_col2;
  }
}

/**
 * \brief Computes a sparse AXPY for one row for full multilines
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline
 *
 * \param line index line_idx
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
static inline void sparse_scal_mul_sub_1_row_vect_array_multiline(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const ml_t multiline,
              const bi_t line_idx,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {

  const uint32_t N    = multiline.sz;
  const mli_t *p_idx  = (N != 0) ? multiline.idx : NULL;
  const re_t *p_val   = (N != 0) ? multiline.val : NULL;
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
 * \brief Computes a sparse AXPY for two rows for multilines
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
static inline void sparse_scal_mul_sub_2_rows_vect_array_multiline(
              const uint32_t Av1_col1,
              const uint32_t Av2_col1,
              const uint32_t Av1_col2,
              const uint32_t Av2_col2,
              const ml_t multiline,
              re_l_t *dense_val1,
              re_l_t *dense_val2) {
 
  // check cases where one pair of the elements is zero
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
static inline void red_dense_array_modular(re_l_t *dense_array, const bi_t bwidth, const mod_t modulus) {
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
static inline mli_t get_head_multiline(const ml_t *m, const bi_t line_idx, re_t *h, mli_t *h_idx) {
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
}

/**
 * \brief Gets first nonzero entry in multiline m at line index line_idx and
 * stores it in h resp. h_idx. Use hybrid representation of multiline.
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
static inline mli_t get_head_multiline_hybrid(const ml_t *m,
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
}

/**
 * \brief Gets first nonzero entry in dense array. Reduces to zero on the go if
 * possible.
 *
 * \param dense array dense_array
 *
 * \param holder for nonzero head value
 *
 * \param column dimension coldim
 *
 * \param field characteristic modulus
 *
 * \return index value corresponding to head value index in dense array
 */
static inline int get_head_dense_array(re_l_t *dense_array,
    re_t *val, const ci_t coldim, const mod_t modulus) {
  ci_t i;
  for (i=0; i<coldim; ++i) {
    *val  = dense_array[i] % modulus;
    if (*val != 0)
      return (int)i;
    else
      dense_array[i]  = 0;
  }
  return -1;
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
static inline void inverse_val(re_t *x, const mod_t modulus) {
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
    *x  =   u1;
    return;
  }
  if (u1 > modulus) {
    u1  -=  modulus;
    *x  =   u1;
    return;
  }
  *x  = u1;
  return;
}

/**
 * \brief Normalizes dense array and returns index of head element in array
 *
 * \param dense array dense_array
 *
 * \param column dimension coldim
 *
 * \param field characteristic modulus
 *
 * \return index value corresponding to head value index in dense array
 */
static inline int normalize_dense_array(re_l_t *dense_array,
    const ci_t coldim, const mod_t modulus) {
  re_t val;
  ci_t i;
  int h1  = get_head_dense_array(dense_array, &val, coldim, modulus);

  if (h1 == -1)
    return h1;

  inverse_val(&val, modulus);
  for (i=h1; i<coldim; ++i) {
    dense_array[i]  *=  val;
    dense_array[i]  %=  modulus;
  }
  return h1;
}

/**
 * \brief Normalizes multiline vector.
 *
 * \param multiline m
 *
 * \param field characteristic modulus
 */
static inline void normalize_multiline(ml_t *m, const ci_t coldim, const mod_t modulus) {
  mli_t idx;
  re_t h1 = 0, h2 = 0;

  if (m->sz == 0)
    return;

  get_head_multiline_hybrid(m, 0, &h1, &idx, coldim);
  get_head_multiline_hybrid(m, 1, &h2, &idx, coldim);

  // invert values modulo modulus
  inverse_val(&h1, modulus);
  inverse_val(&h2, modulus);

  // skip if both are 1 and/or 0
  if ((h1 == 0 || h1 == 1) && (h2 == 0 || h2 == 1))
    return;

  re_l_t tmp_val;
  // normalize h2
  if (h1 == 0 || h1 == 1) {
    for (idx=0; idx<m->sz; ++idx) {
      tmp_val         = (re_l_t)m->val[2*idx+1] * h2;
      m->val[2*idx+1] = tmp_val % modulus;
    }
  } else {
    // normalize h1
    if (h2 == 0 || h2 == 1) {
      for (idx=0; idx<m->sz; ++idx) {
        tmp_val       = (re_l_t)m->val[2*idx] * h1;
        m->val[2*idx] = tmp_val % modulus;
      }
    // normalize h1 and h2
    } else {
      for (idx=0; idx<m->sz; ++idx) {
        tmp_val         = (re_l_t)m->val[2*idx] * h1;
        m->val[2*idx]   = tmp_val % modulus;
        tmp_val         = (re_l_t)m->val[2*idx+1] * h2;
        m->val[2*idx+1] = tmp_val % modulus;
      }
    }
  }
}

/**
 * \brief Copies two dense arrays to a dense multiline for further processing.
 * \note Multiline m must have already memory allocated correspondingly.
 * Moreover all entries of m->val must be initialized to zero beforehand.
 *
 * \param dense array dense_1
 *
 * \param dense array dense_2
 *
 * \param minimal first position of non-zero entries in dense_1 and dense_2
 * start_pos
 *
 * \param multiline m
 *
 * \param array size resp. column dimension coldim
 */
static inline void copy_dense_arrays_to_zero_dense_multiline(const re_l_t *dense_1,
    const re_l_t *dense_2, const int start_pos, ml_t *m, const ci_t coldim,
    const mod_t modulus) {

  printf("in copy mlv %p\n",m->val);
  ci_t i;
  re_l_t tmp_1, tmp_2;
  for (i=start_pos; i<coldim; ++i) {
    tmp_1 = dense_1[i] % modulus;
    tmp_2 = dense_2[i] % modulus;
    //printf("t1 %lu -- t2 %lu\n",tmp_1,tmp_2);
    if (tmp_1 != 0 || tmp_2 != 0) {
      //printf("B %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
      m->val[2*i]   = (re_t)tmp_1;
      //printf("A1 %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
      m->val[2*i+1] = (re_t)tmp_2;
    }
    //printf("A %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
  }
}

/**
 * \brief Copies one dense array to a dense multiline for further processing.
 * \note Multiline m must have already memory allocated correspondingly.
 * Moreover all entries of m->val must be initialized to zero beforehand.
 *
 * \param dense array dense_1
 *
 * \param position of first non-zero entry in array dense_1 start_pos
 *
 * \param multiline m
 *
 * \param array size resp. column dimension coldim
 */
static inline void copy_dense_array_to_zero_dense_multiline(const re_l_t *dense_1,
    const int start_pos, ml_t *m, const ci_t coldim, const mod_t modulus) {

  printf("in copy mlv %p\n",m->val);
  ci_t i;
  re_l_t tmp_1;
  for (i=start_pos; i<coldim; ++i) {
    tmp_1 = dense_1[i] % modulus;
    printf("t1 %lu\n",tmp_1);
    if (tmp_1 != 0) {
      printf("B %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
      m->val[2*i]   = (re_t)tmp_1;
      printf("A1 %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
    }
    printf("A %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
  }
}
static inline void copy_dense_arrays_to_dense_multiline(const re_l_t *dense_1,
    const re_l_t *dense_2, ml_t *m, const ci_t coldim, const mod_t modulus) {

  ci_t i;
  re_l_t tmp_1, tmp_2;
  for (i=0; i<coldim; ++i) {
    tmp_1 = dense_1[i] % modulus;
    tmp_2 = dense_2[i] % modulus;
    printf("t1 %lu -- t2 %lu\n",tmp_1,tmp_2);
      printf("B %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
      m->val[2*i]   = (re_t)tmp_1;
      m->val[2*i+1] = (re_t)tmp_2;
    printf("A %d v1 %u -- v2 %u\n",i,m->val[2*i],m->val[2*i+1]);
  }
}

/**
 * \brief Copies multiline to two dense arrays for further processing.
 *
 * \param multiline m
 *
 * \param dense array dense_1
 *
 * \param dense array dense_2
 *
 * \param array size resp. column dimension coldim
 */
static inline void copy_multiline_to_dense_array(const ml_t m, re_l_t *dense_1,
    re_l_t *dense_2, const ci_t coldim) {
  if (m.sz == 0)
    return;

  register mli_t idx;
  ci_t i;

  if (m.sz < coldim) {
    printf("m.sz %d\n",m.sz);
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

/**
 * \brief Returns smallest row index in waiting list for echelonization.
 *
 * \param waiting list waiting
 *
 * \return 1 if waiting list is not empty, 0 else
 */
static inline int get_smallest_waiting_row(wl_t *waiting) {
  if (waiting->sz == 0)
    return 0;

  // sort the waiting list
  qsort(waiting->list, waiting->sz, sizeof(wle_t), cmp_wle);
  // store last (smallest index) element separately and
  // remove it from the list
  waiting->sidx = waiting->list[waiting->sz-1].idx;
  waiting->slp  = waiting->list[waiting->sz-1].lp;
  waiting->sz--;

  return 1;
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
void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B,
    const ri_t bheight, const mod_t modulus, int inv_scalars);

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
void red_with_triangular_block(mbl_t *block_A, re_l_t **dense_B,
    const ri_t bheight, const mod_t modulus, int inv_scalars);

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
int elim_fl_A_block(sbm_fl_t **A, sbm_fl_t *B, const mod_t modulus, int nthrds);

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
int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, const ci_t block_col_idx_B,
    const ri_t nbrows_A, const mod_t modulus);

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
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C, sbm_fl_t *D, const mod_t modulus,
    const int nthrds);

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
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const mod_t modulus);

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
 * \return rank of D_red
 */
ri_t elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, const mod_t modulus, int nthrds);

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


/**
 * \brief Echelonizes the multiline rows of A from row index from up to row
 * index to sequential. This is done in order to prepare a basis for the
 * upcoming parallel structure echelonization
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param start row index from
 *
 * \param end row index to
 *
 * \param field characteristic modulus
 *
 * \return number of real pivots found
 */
ri_t echelonize_rows_sequential(sm_fl_ml_t *A, const ri_t from, const ri_t to,
    const mod_t modulus);

/**
 * \brief Echelonizes the multiline rows of A from row index from up to row
 * index to sequential. This is done in order to prepare a basis for the
 * upcoming parallel structure echelonization
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param global value of next row to be reduced (shared among the threads!)
 *
 * \param global last pivot known (shared among the threads!)
 *
 * \param list of rows waiting to be reduced waiting
 *
 * \param field characteristic modulus
 *
 * \return 0 if success, 1 if failure
 */
int echelonize_rows_task(sm_fl_ml_t *A, const ri_t ml_ncols,
    ri_t global_next_row_to_reduce, ri_t global_last_piv,
    wl_t *waiting, const mod_t modulus);

/**
 * \brief Echelonizes one multiline row (represented by dense_array_1 and
 * dense_array_2) by the already known pivot rows A->ml[first_piv] up to
 * A->ml[last_piv].
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param dense array for first part of multiline row dense_array_1
 *
 * \param dense array for second part of multiline row dense_array_2
 *
 * \param first known pivot first_piv
 *
 * \param last known pivot last_piv
 *
 * \param field characteristic modulus
 */
void echelonize_one_row(sm_fl_ml_t *A, re_l_t *dense_array_1,
    re_l_t *dense_array_2, const ri_t first_piv, const ri_t last_piv,
    const mod_t modulus);

/**
 * \brief Restores the two dense arrays in a multiline row in D after
 * echelonizing the corresponding multiline row
 *
 * \param multiline row in D ml
 *
 * \param dense array for first part of multiline row dense_array_1
 *
 * \param dense array for second part of multiline row dense_array_2
 *
 * \param columns dimension coldim
 *
 * \param field characteristic modulus
 *
 * \param reduce or just store in multiline row? reduce
 * if 1 then reduction is done, else we only store the multiline row
 */
void save_back_and_reduce(ml_t *ml, re_l_t *dense_array_1,
    re_l_t *dense_array_2, const ci_t coldim, const mod_t modulus,
    const int reduce);
#endif
