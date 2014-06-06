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
void load_matrix_jcf_format(FILE *file, sparse_mat_t *mat);
void write_matrix_jcf_format(sparse_mat_t *mat, FILE *file);

#endif
