#include "gb.h"
void print_help() {
  printf("\n");
  printf("NAME\n");
  printf("    gb - computes a Groebner basis\n");
  printf("\n");
  printf("SYNOPSIS\n");
  printf("    gb [options] [file..]\n");
  printf("\n");
  printf("DESCRIPTION\n");
  printf("    Computes Gaussian Elimination of a structured hybrid sparse-\n");
  printf("    dense matrix coming from Groebner basis computations.\n");
  printf("\n");
  printf("OPTIONS\n");
  printf("    -b         Dimensions of a block in submatrices.\n");
  printf("               Default: 256.\n");
  printf("    -c         Result checked resp. validated with structured Gaussian\n");
  printf("               Elimination. By default there is no validation.\n");
  printf("    -h         Print help.\n");
  printf("    -f         Free memory on the go.\n");
  printf("               By default memory is freed only at the end.\n");
  printf("    -m         Number of rows per multiline.\n");
  printf("               Default: 1.\n");
  printf("    -p         Writes intermediate matrices in pbm format.\n");
  printf("    -t THRDS   Number of threads used.\n");
  printf("               Default: 1.\n");
  printf("    -v LEVEL   Level of verbosity:\n");
  printf("               1 -> only error messages printed\n");
  printf("               2 -> additionally meta information printed\n");
  printf("\n");

  return; 
}


int main(int argc, char *argv[]) {  
  const char *fn        = NULL;
  int free_mem          = 0;
  int validate_results  = 0;
  int verbose           = 0;
  int method            = 0;
  int nrows_multiline   = 1;
  int block_dimension   = 256;
  int write_pbm         = 0;
  int nthreads          = 1;

  int index;
  int opt;

  // timing structs
  struct timeval t_load_start;

  opterr  = 0;

  while ((opt = getopt(argc, argv, "b:cfhm:t:v:p")) != -1) {
    switch (opt) {
      case 'b': 
        block_dimension = atoi(optarg);
        if (block_dimension < 1)
          block_dimension = 256;
        break;
      case 'c': 
        validate_results  = 1;
        break;
      case 'f':
        free_mem  = 1;
        break;
      case 'h':
        print_help();
        return 0;
      case 'm': 
        nrows_multiline = atoi(optarg);
        if (nrows_multiline < 1)
          nrows_multiline = 1;
        break;
      case 'p': 
        write_pbm = 1;
        break;
      case 't': 
        nthreads  = atoi(optarg);
        break;
      case 'v': 
        verbose = atoi(optarg);
        if (verbose > 2)
          verbose = 2;
        break;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
              "Unknown option character `\\x%x'.\n",
              optopt);
        return 1;
      default:
        abort ();
    }
  }
  for (index = optind; index < argc; index++)
    fn = argv[index];

  if (fn == NULL) {
    fprintf(stderr, "File name is required.\nSee help using '-h' option.\n");
    return 1;
  }

  sm_t *M = NULL;

  if (verbose > 0) {
    printf("----------------------------------------------------------------\n");
    printf("----------- Computing a special Gaussian Elimination -----------\n");
    printf("--------------- with the following options set -----------------\n");
    printf("----------------------------------------------------------------\n");
    printf("number of threads   %4d\n", nthreads);
    printf("dimension of blocks %4d\n", block_dimension, block_dimension);
    printf("----------------------------------------------------------------\n");
  }

  M = load_jcf_matrix(fn, verbose);

  if (write_pbm) {
    const char *pbm_fn  = "input-matrix.pbm";
    write_jcf_matrix_to_pbm(M, pbm_fn, verbose);
  }

  // construct splicing of matrix M into A, B, C and D
  sbm_fl_t *A   = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  sbm_fl_t *B   = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  sbm_fl_t *C   = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  sbm_fl_t *D   = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  map_fl_t *map = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD
  /*
  sbm_fl_t A   = NULL;
  sbm_fl_t B   = NULL;
  sbm_fl_t C   = NULL;
  sbm_fl_t D   = NULL;
  map_fl_t map = NULL; // stores mappings from M <-> ABCD
  */
  if (verbose > 1)
    gettimeofday(&t_load_start, NULL);

  splice_fl_matrix(M, A, B, C, D, map, block_dimension, nrows_multiline, nthreads, verbose);

  if (verbose > 1) {
    printf("----------------------------------------------------------------\n");
    printf("Splicing and mapping of input matrix done.\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("----------------------------------------------------------------\n");
    printf("A [%d x %d]\n",A->nrows,A->ncols);
    printf("B [%d x %d]\n",B->nrows,B->ncols);
    printf("C [%d x %d]\n",C->nrows,C->ncols);
    printf("D [%d x %d]\n",D->nrows,D->ncols);
    printf("----------------------------------------------------------------\n");
    printf("||| Time for constructing submatrices A, B, C and D: %7.3f sec\n",
        walltime(t_load_start) / (1000000));
  }

  // computing Gaussian Elimination of A using methods of Faugère & Lachartre
  elim_fl(M);

  free(M);
  return 0;
}
