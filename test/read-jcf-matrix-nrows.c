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
 * \file read-jcf-matrix-nrows.c
 * \brief Unit test for gb
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include <assert.h>
#include "cli/io.h"

int main(int argc, char **argv) {
  const char *fn  = "examples/jcf-type-f4-kat12-mat1";
  int verbose     = 0;
  sm_t *M         = NULL;

  M = load_jcf_matrix(fn, verbose, 0);

  assert(M->nrows == 929);

  free(M);
  return 0;
}
