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
#include <sparse_matrix.h>

// reading & writing
sparse_mat_t *load_matrix_jcf_format(const char *fn, int verbose);
void write_matrix_jcf_format(sparse_mat_t *mat, FILE *file);

#endif
