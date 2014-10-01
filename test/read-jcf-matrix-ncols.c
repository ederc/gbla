/**
 * \file read-jcf-matrix-ncols.c
 * \brief Unit test for gb
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include <assert.h>
#include "../cli/io.h"

int main(int argc, char **argv) {
  const char *fn  = "examples/jcf-type-f4-kat12-mat1";
  int verbose     = 0;
  sm_t *M         = NULL;

  M = load_jcf_matrix(fn, verbose);

  assert(M->ncols == 1657);

  free(M);
  return 0;
}
