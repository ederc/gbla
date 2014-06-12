/**
 * \file sm.h
 * \brief Interface for sparse matrix structures
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
  re_t **rows;    /*!<  address of row: m->row[i] gives first
                        address of ith row */
  ci_t **pos;     /*!<  position of entry in row: m->pos[i] gives first
                        address of first position of nonzero entry in row i */
  ci_t *rwidth;   /*!<  width of row: m->width[i] gives number of nonzero
                        entries in row i */
} sm_t;

#endif
