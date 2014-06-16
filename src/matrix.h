/**
 * \file matrix.h
 * \brief Interfaces for matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_MATRIX_H
#define GB_MATRIX_H

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


/**
 * \brief Enum of block alignments: Entries in the block submatrices are stored
 * w.r.t. different alignments, some from top to down, then from left to right,
 * other in combinations and inversions of these.
 */
enum ba_t {
  tdlr, /*!<  top-to-down, left-to-right */
  tdrl, /*!<  top-to-down, right-to-left */
  dtlr, /*!<  down-to-top, left-to-right */
  dtrl  /*!<  down-to-top, right-to-left */
};


/**
 * \brief Sparse block matrix structure for Faugère-Lachartre decompositions
 *
 */

typedef struct sbm_fl_t {
  ri_t nrows;                 /*!<  number of rows */
  ci_t ncols;                 /*!<  number of columns */
  ri_t bheight;               /*!<  number of rows per block */
  ci_t bwidth;                /*!<  number of columns per block */
  enum ba_t block_alignment;  /*!<  memory alignment in the block: depending on
                                    the block parts A,B,C and D the entries are
                                    stored top-to-down, left-to-right, or any
                                    combination and inversion of these two */
  int fill_empty;             /*!<  if 1 then empty blocks are used to fill
                                    block matrix
                                    if 0 then no fill is done */
  re_t **rows;                /*!<  address of row: M->row[i] gives first
                                    address of ith row */
  ci_t **pos;                 /*!<  position of entry in row: M->pos[i] gives
                                    first address of first position of nonzero
                                    entry in row i */
  ci_t *rwidth;               /*!<  width of row: M->width[i] gives number of
                                    nonzero entries in row i */
} sbm_fl_t;

#endif
