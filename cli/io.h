/* gbla: Gröbner Basis Linear Algebra
 * This file is part of gbla.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA, or see  <http://www.gnu.org/licenses/>. */


/**
 * \file io.h
 * \brief Input/output routines for matrices
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GBLA_IO_H
#define GBLA_IO_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "../config.h"
#include "../src/matrix.h"
#include <omp.h>

/*  ========== TIMINGS ========== */

/**
 * \brief Returns walltime since time t0.
 *
 * \param Start time stamp t0
 *
 * \return Difference in walltime between current time stamp and t0
 */

double walltime(struct timeval t_start);



/*  ========== READING ========== */

/**
 * \brief Loads a matrix given in Jean-Charles Faugère's style to a sparse
 * matrix in (sm_t *) format.
 *
 * \param fn File name
 * \param vebose If 1: Printing of error messages
 *               If 2: Also printing of meta information
 * \param format 0 for old, 1 for new
 *
 * \param nthrds number of threads for construction of matrix (only possible in
 * new format
 *
 * \return Corresponding sparse matrix in (sm_t *) format
 */
sm_t *load_jcf_matrix(const char *fn, int verbose, int new_format, int nthrds);

/**
 * \brief Loads a matrix given in Schreyer's style to a sparse
 * matrix in (sm_t *) format. Matrices are arising from the Prym-Green
 * conjecture.
 *
 * \param fn File name
 * \param vebose If 1: Printing of error messages
 *               If 2: Also printing of meta information
 *
 * \return Corresponding sparse matrix in (sm_t *) format
 */
sm_t *load_schreyer_matrix(const char *fn, int verbose);



/*  ========== WRITING ========== */

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


/*  ========== PRINTING ========== */

/**
 * \brief Computes number of nonzero elements in multiline submatrix A and
 * A's density. This is time consuming and thus only done for verbosity levels > 2.
 *
 * \param multiline submatrix A
 */

 void compute_density_ml_submatrix(sm_fl_ml_t *A) ;

/**
 * \brief Computes number of nonzero elements in block submatrix A and A's density.
 * This is time consuming and thus only done for verbosity levels > 2.
 *
 * \param block submatrix A
 */

 void compute_density_block_submatrix(sbm_fl_t *A) ;



/**
 * \brief Prints memory usage by getting information from /proc/self/stat.
 *
 */
void print_mem_usage(void);
#endif

/* vim:sts=2:sw=2:ts=2:
 */
