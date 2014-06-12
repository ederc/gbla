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
#include <sm.h>

// ========== READING ==========

/**
 * \brief Loads a matrix given in Jean-Charles Faugère's style to a sparse
 * matrix in (sm_t *) format.
 *
 * \param fn File name
 * \param vebose If 1: Printing of error messages 
 *               If 2: Also printing of meta information
 *
 * \return Corresponding sparse matrix in (sm_t *) format
 */
sm_t *load_jcf_matrix(const char *fn, int verbose);



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
#endif
