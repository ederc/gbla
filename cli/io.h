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

// reading & writing
sm_t *load_jcf_matrix(const char *fn, int verbose);
void write_jcf_matrix_to_file(sm_t *mat, const char *fn);

void write_jcf_matrix_to_pbm(sm_t *mat, const char *fn);
#endif
