/**
 * \file sm.h
 * \brief Interfaces for sparse matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_SM_H
#define GB_SM_H

#include <types.h>

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

#endif
