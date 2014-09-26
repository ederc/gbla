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
  float density;  /*!<  density used for adjusting memory allocation during
                        splicing and generation of ABCD blocks */
  mod_t mod;      /*!<  modulo/field characteristic */
  float fs;       /*!<  file size of input matrix */
  char *fsu;      /*!<  file size unit of input matrix, e.g. GB */
  re_t **rows;    /*!<  address of row: M->rows[i] gives first
                        address of ith row */
  ci_t **pos;     /*!<  position of entry in row: M->pos[i] gives first
                        address of first position of nonzero entry in row i */
  ci_t *rwidth;   /*!<  width of row: M->rwidth[i] gives number of nonzero
                        entries in row i */
} sm_t;


/**
 * \brief A multiline block is a vector of an vector of multilines. It consists of
 * an index for the column of the value entries in rows. For each index rows stores
 * __GB_NROWS_MULTILINE elements, i.e. this many rows are taken care of at once.
 */
typedef struct mbl_t {
  bi_t *idx;  /*!< column index in the multiline vector */
  re_t *val;  /*!< multiline row, must be __GB_NROWS_MULTILINE * length(idx) */
  bi_t sz;    /*!< current length of the block row */
} mbl_t;


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
 * \brief Sparse block matrix structure for Faugère-Lachartre decompositions.
 * Can be used for usual line implementations and multi line implementations,
 * the corresponding multi line functions are labeled with an "_ml".
 *
 */

typedef struct sbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  ri_t bheight;     /*!<  number of rows per block */
  ci_t bwidth;      /*!<  number of columns per block */
  enum ba_t ba;     /*!<  memory alignment in the block: depending on
                          the block parts A,B,C and D the entries are
                          stored top-to-down, left-to-right, or any
                          combination and inversion of these two */
  int fe;           /*!<  if 1 then empty blocks are used to fill
                          block matrix
                          if 0 then no fill is done */
  int hr;           /*!<  if 1 then hybrid (sparse/dense) rows are accepted
                          if 0 then those rows are not accepted */
  mbl_t ***blocks;  /*!<  address of blocks: M->blocks[i][j] gives address of 
                          block. There are nrows/bheight * ncols/bwidth blocks. */
} sbm_fl_t;

/**
 * \brief Sparse matrix with multiline row structure for Faugère-Lachartre
 * decompositions. Can be used for multiline implementations, the corresponding
 * multiline functions are labeled with an "_ml".
 */

typedef struct sm_fl_ml_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  enum ba_t ba;     /*!<  memory alignment in the block: depending on
                          the block parts A,B,C and D the entries are
                          stored top-to-down, left-to-right, or any
                          combination and inversion of these two */
  int fe;           /*!<  if 1 then empty blocks are used to fill
                          block matrix
                          if 0 then no fill is done */
  int hr;           /*!<  if 1 then hybrid (sparse/dense) rows are accepted
                          if 0 then those rows are not accepted */
  mbl_t *ml;        /*!<  address of multilines: M->ml[i] gives address of 
                          multiline i. There are nrows/__GB_NROWS_MULTILINE
                          multilines. */
} sm_fl_ml_t;

#endif
