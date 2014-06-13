/**
 * \file sm.h
 * \brief Interfaces for sparse matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_SM_H
#define GB_SM_H

#include <stdint.h>

/// matrix row entry type
typedef uint16_t  re_t;
/// row and column index types
typedef uint32_t  ci_t;
typedef uint32_t  ri_t;
/// number of nonzero elements type
typedef uint64_t  nnz_t;
/// tyoe of field characteristic
typedef uint32_t  mod_t;

/**
 * \brief Sparse matrix structure for reading jcf matrices
 *
 */

typedef struct sm_t {
  ri_t nrows;     /*!<  number of rows */
  ci_t ncols;     /*!<  number of columns */
  nnz_t nnz;      /*!<  number of nonzero entries */
  mod_t mod;      /*!<  modulo/field characteristic */
  re_t **rows;    /*!<  address of row: M->row[i] gives first
                        address of ith row */
  ci_t **pos;     /*!<  position of entry in row: M->pos[i] gives first
                        address of first position of nonzero entry in row i */
  ci_t *rwidth;   /*!<  width of row: M->width[i] gives number of nonzero
                        entries in row i */
} sm_t;

/**
 * \brief Indexer for subdividing sparse matrix into 4 parts as described by
 * Faugère and Lachartre in http://dx.doi.org/10.1145/1837210.1837225.
 *
 *             A | B
 * M ---->     --+--
 *             C | D
 */
typedef struct sm_idx_t {
  ci_t *piv_col_map;            /*!<  map of pivot columns: from input matrix M to
                                      submatrix A; has length M->ncols, maps non-pivot
                                      columns to __GB_MINUS_ONE_32 */
  ci_t *non_piv_col_map;        /*!<  map of non-pivot columns: from input matrix M to
                                      submatrix B; has length M->ncols, maps pivot
                                      columns to __GB_MINUS_ONE_32 */
  ci_t *piv_col_rev_map;        /*!<  reverse map of pivot columns: from submatrices
                                      A and B to input matrix M */
  ci_t *non_piv_col_rev_map;    /*!<  reverse map of non-pivot columns: from submatrices
                                      A and B to input matrix M */
  ri_t *piv_rows_idxs_by_entry; /*!<  has length M->nrows, maps pivot columns to
                                      their corresponding row index, maps non-pivot
                                      columns to __GB_MINUS_ONE_32 */
  ri_t *non_piv_rows_idxs;      /*!<  indexes of non-pivot rows */
} sm_idx_t;

#endif
