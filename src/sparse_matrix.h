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

// rows
typedef struct row_t {
  entry_t   *first_entry;
  entry_t   *last_entry;
  col_idx_t *first_col;
  col_idx_t *last_col;
} row_t;

// sparse matrices
typedef struct sparse_mat_t {
  row_t **rows;
} sparse_mat_t;

#endif
