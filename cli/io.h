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
 * \file io.h
 * \brief Input/output routines for matrices
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_IO_H
#define GB_IO_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <gbla_config.h>
#include <src/matrix.h>

// ========== TIMINGS ==========

/**
 * \brief Returns walltime since time t0.
 *
 * \param Start time stamp t0
 *
 * \return Difference in walltime between current time stamp and t0
 */

double walltime(struct timeval t_start);



// ========== READING ==========

/**
 * \brief Loads a matrix given in Jean-Charles Faugère's style to a sparse
 * matrix in (sm_t *) format.
 *
 * \param fn File name
 * \param vebose If 1: Printing of error messages
 *               If 2: Also printing of meta information
 * @param format 0 for old, 1 for new
 *
 * \return Corresponding sparse matrix in (sm_t *) format
 */
sm_t *load_jcf_matrix(const char *fn, int verbose, int new_format);



// ========== WRITING ==========

/**
 * \brief Writes a sparse matrix in (sm_t *) format to a file in Jean-Charles Faugère's style.
 *
 * \param M Matrix in (sm_t *) format
 * \param fn File name
 * \param vebose If 1: Printing of error messages
 *               If 2: Also printing meta information
 */
void write_jcf_matrix_to_file(sm_t *M, const char *fn, int verbose);

/**
 * \brief Writes a sparse matrix in (sm_t *) format to a file as portable bitmap
 * image.
 *
 * \param M Matrix in (sm_t *) format
 * \param fn File name
 * \param vebose If 1: Printing of error messages
 *               If 2: Also printing meta information
 */
void write_jcf_matrix_to_pbm(sm_t *M, const char *fn, int verbose);


// ========== PRINTING ==========

/**
 * \brief Computes number of nonzero elements in multiline submatrix A and
 * A's density. This is time consuming and thus only done for verbosity levels > 2.
 *
 * \param multiline submatrix A
 */

static inline void compute_density_ml_submatrix(sm_fl_ml_t *A) {
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GB_NROWS_MULTILINE);
  uint32_t i, j;
  for (i=0; i<rlA; ++i) {
    if (A->ml[i].sz>0) {
      for (j=0; j<A->ml[i].sz; ++j) {
        if (A->ml[i].val[2*j] != 0) {
          A->nnz++;
        }
        if (A->ml[i].val[2*j+1] != 0) {
          A->nnz++;
        }
      }
    }
  }
  A->density  = (double) (A->nnz * 100) / (A->nrows * (nnz_t)A->ncols);
}

/**
 * \brief Computes number of nonzero elements in block submatrix A and A's density.
 * This is time consuming and thus only done for verbosity levels > 2.
 *
 * \param block submatrix A
 */

static inline void compute_density_block_submatrix(sbm_fl_t *A) {
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / A->bheight);
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);
  uint32_t i, j, k, l;
  for (i=0; i<rlA; ++i) {
    for (j=0; j<clA; ++j) {
      if (A->blocks[i][j] != NULL) {
        for (k=0; k<A->bwidth/2; ++k) {
          for (l=0; l<A->blocks[i][j][k].sz; ++l) {
            if (A->blocks[i][j][k].val[2*l]) {
              A->nnz++;
            }
            if (A->blocks[i][j][k].val[2*l+1]) {
              A->nnz++;
            }
          }
        }
      }
    }
  }
  A->density  = (double) (A->nnz * 100) / (A->nrows * (nnz_t)A->ncols);
}

/**
 * \brief Prints memory usage by getting information from /proc/self/stat.
 *
 */
void print_mem_usage();
#endif
