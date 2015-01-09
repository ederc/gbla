/**
 * \file gb.h
 * \brief Main command line program
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_GB_H
#define GB_GB_H

#include "io.h"
#include <mapping.h>
#include <elimination.h>
#include <math.h>

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv.
 * Afterwards the row echelon form of M is computed.
 *
 *  \param original matrix M
 *
 *  \param dimension of blocks block_dimension
 *
 *  \param number of rows per multiline rows_multiline
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 */
int fl_block(sm_t *M, int block_dimension, int rows_multiline, int nthreads, int free_mem,
    int verbose, int reduce_completely);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv.
 * Afterwards the row echelon form of M is computed.
 *
 *  \note Submatrix A is in multiline format.
 *
 *  \param original matrix M
 *
 *  \param dimension of blocks block_dimension
 *
 *  \param number of rows per multiline rows_multiline
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 */
int fl_ml_A(sm_t *M, int block_dimension, int rows_multiline, int nthreads, int free_mem, int verbose);

/**
 * \brief Constructs the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 * In the subdivision the following dimensions hold:
 * A->nrows = B->nrows = map->npiv // number of pivots found
 * C->nrows = D->nrows = M->nrows - map->npiv // non-pivots
 * A->ncols = C->ncols = map->npiv
 * B->ncols = D->ncols = M->ncols - map->npiv.
 * Afterwards the row echelon form of M is computed.
 *
 *  \note Submatrices A and C are in multiline format.
 *
 *  \param original matrix M
 *
 *  \param dimension of blocks block_dimension
 *
 *  \param number of rows per multiline rows_multiline
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 */
int fl_ml_A_C(sm_t *M, int block_dimension, int rows_multiline, int nthreads, int free_mem,
    int verbose, int reduce_completely);

#endif
