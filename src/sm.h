/**
 * \file sparse_matrix.h
 * \brief Interface for sparse matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_SPARSE_MATRIX_H
#define GB_SPARSE_MATRIX_H

#include <typedefs.h>

/// matrix row entry type
typedef uint16  re_t;
/// row and column index types
typedef uint32  ci_t;
typedef uint32  ri_t;
/// number of nonzero elements type
typedef uint64  nnz_t;
/// tyoe of field characteristic
typedef uint32  mod_t;

/**
 * \brief Sparse matrix structure for reading jcf matrices
 *
 */

typedef struct sm_t {
  ri_t nrows;     /*!<  number of rows */
  ci_t ncols;     /*!<  number of columns */
  nnz_t nnz;      /*!<  number of nonzero entries */
  mod_t mod;      /*!<  modulo/field characteristic */
  re_t **rows;    /*!<  address of row: m->row[i] gives first
                        address of ith row */
  ci_t **pos;     /*!<  position of entry in row: m->pos[i] gives first
                        address of first position of nonzero entry in row i */
  ci_t *rwidth;   /*!<  width of row: m->width[i] gives number of nonzero
                        entries in row i */
} sm_t;

#endif
