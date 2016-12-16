/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */




/**
 * \file gbla.h
 * \brief Main command line program
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GBLA_GBLA_H
#define GBLA_GBLa_H

#include <math.h>
#include <unistd.h>
#include "io.h"
#include "mapping.h"
#include "elimination.h"

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
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block(sm_t *M, int block_dimension, int rows_multiline, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *  \note No multilines at all, but hybrid small block representations
 *  in order to exploit SIMD operations as much as possible.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block_hybrid(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *  \note No multilines at all, but hybrid small block representations
 *  in order to exploit SIMD operations as much as possible.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block_sparse_dense_keep_A(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *  \note No multilines at all, but hybrid small block representations
 *  in order to exploit SIMD operations as much as possible.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block_sparse_dense_old(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *  \note This implementation is a draft for the third version of gbla.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
int fl_v03(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);
*/

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
 *  \note No multilines at all, but hybrid small block representations
 *  in order to exploit SIMD operations as much as possible.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block_sparse_dense(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *  \note No multilines at all, but hybrid small block representations
 *  in order to exploit SIMD operations as much as possible.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block_hybrid_dense(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *  \note No multilines at all, but dense small block representations in order
 *  to exploit SIMD operations as much as possible.
 *
 *  \param original matrix M
 *
 *  \param destructing input matrix on the go? free_mem
 *
 *  \param number of threads to be used nthreads
 *
 *  \param level of verbosity
 *
 *  \param compute a complete reduced row echelon form? reduce_completely
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_block_dense(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

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
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_ml_A(sm_t *M, int block_dimension, int rows_multiline, int nthreads, int free_mem, int verbose, int dense_reducer);

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
 *
 *  \param use third party dense reducer for D dense_reducer
 */
int fl_ml_A_C(sm_t *M, int block_dimension, int rows_multiline, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer);

#endif

/* vim:sts=2:sw=2:ts=2:
 */
