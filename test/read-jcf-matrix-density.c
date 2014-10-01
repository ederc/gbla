/**
 * \file read-jcf-matrix-density.c
 * \brief Unit test for gb
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include <assert.h>
#include "../cli/io.h"

int main(int argc, char **argv) {
  return 0;
  const char *fn  = "examples/jcf-type-f4-kat12-mat1";
  int verbose     = 0;
  float density   = 6.62;
  sm_t *M         = NULL;

  M = load_jcf_matrix(fn, verbose);
  
  assert(M->density - density < 0.0000001);
  assert(density - M->density < 0.0000001);

  free(M);
  return 0;
}
