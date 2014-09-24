#include "gb.h"

#define __GB_CLI_DEBUG  0

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

#if __GB_CLI_DEBUG
  // column loops 
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / C->bwidth);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);
  // row loops 
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / A->bheight);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / C->bheight);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      for (kk=0; kk<block_dimension/2; ++kk) {
        printf("%d .. %d .. %d\n",ii,jj,kk);
        printf("size %d\n", A->blocks[ii][jj][kk].sz * 2);
        if (A->blocks[ii][jj][kk].sz>0) {
          for (ll=0; ll<A->blocks[ii][jj][kk].sz; ++ll) {
            printf("%d %d ", A->blocks[ii][jj][kk].val[2*ll], A->blocks[ii][jj][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
  }
  for (ii=0; ii<rlB; ++ii) {
    for (jj=0; jj<clB; ++jj) {
      for (kk=0; kk<block_dimension/2; ++kk) {
        printf("%d .. %d .. %d\n",ii,jj,kk);
        printf("size %d\n", B->blocks[ii][jj][kk].sz * 2);
        if (B->blocks[ii][jj][kk].sz>0) {
          for (ll=0; ll<B->blocks[ii][jj][kk].sz; ++ll) {
            printf("%d %d ", B->blocks[ii][jj][kk].val[2*ll], B->blocks[ii][jj][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
  }
  for (ii=0; ii<rlC; ++ii) {
    for (jj=0; jj<clC; ++jj) {
      for (kk=0; kk<block_dimension/2; ++kk) {
        printf("%d .. %d .. %d\n",ii,jj,kk);
        printf("size %d\n", C->blocks[ii][jj][kk].sz * 2);
        if (C->blocks[ii][jj][kk].sz>0) {
          for (ll=0; ll<C->blocks[ii][jj][kk].sz; ++ll) {
            printf("%d %d ", C->blocks[ii][jj][kk].val[2*ll], C->blocks[ii][jj][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
  }
  for (ii=0; ii<rlD; ++ii) {
    for (jj=0; jj<clD; ++jj) {
      for (kk=0; kk<block_dimension/2; ++kk) {
        printf("%d .. %d .. %d\n",ii,jj,kk);
        printf("size %d\n", D->blocks[ii][jj][kk].sz * 2);
        if (D->blocks[ii][jj][kk].sz>0) {
          for (ll=0; ll<D->blocks[ii][jj][kk].sz; ++ll) {
            printf("%d %d ", D->blocks[ii][jj][kk].val[2*ll], D->blocks[ii][jj][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
  }
#endif

  // computing Gaussian Elimination of A using methods of Faugère & Lachartre
  elim_fl(M);

  free(M);
  return 0;
}
