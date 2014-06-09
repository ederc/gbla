/**
 * \file sparse_matrix.h
 * \brief Interface for sparse matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_SPARSE_MATRIX_H
#define GB_SPARSE_MATRIX_H

#include <typedefs.h>

// matrix typedefs
typedef uint16  entry_t;
typedef uint32  col_idx_t;
typedef uint32  row_idx_t;
typedef uint32  col_dim_t;
typedef uint32  row_dim_t;
typedef uint64  nnz_t;
typedef uint32  mod_t;

/**
 * \brief Sparse matrix structure for reading jcf matrices
 *
 */

typedef struct sparse_mat_t {
  row_dim_t nrows;    /*!< number of rows */
  col_dim_t ncols;    /*!< number of columns */
  nnz_t     nnz;      /*!< number of nonzero entries */
  mod_t     mod;      /*!< modulo/field characteristic */
  entry_t   **rows;   /*!< address of row: m->row[i] gives first
                           address of ith row */
  col_idx_t **pos;    /*!< position of entry in row: m->pos[i] gives first
                           address of first position of nonzero entry in row i */
} sparse_mat_t;

#endif
