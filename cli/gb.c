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
  printf("    -h         Print help.\n");
  printf("    -m         Free memory on the go.\n");
  printf("               By default memory is freed only at the end.\n");
  printf("    -t THRDS   Number of threads used.\n");
  printf("               Default: 1.\n");
  printf("    -v         Result validated with structured Gaussian Elimination.\n");
  printf("               By default there is no validation.\n");
  printf("\n");

  return; 
}


int main(int argc, char *argv[]) {  
  const char *fn        = NULL;
  int free_mem          = 0;
  int validate_results  = 0;
  int verbose           = 0;
  int method            = 0;
  int write_pbm         = 0;
  int nthrds            = 1;

  int index;
  int opt;

  opterr  = 0;

  while ((opt = getopt(argc, argv, "f:hmt:cv:p")) != -1) {
    switch (opt) {
      case 'h':
        print_help();
        return 0;
      case 'f':
        printf("%s\n",optarg);
        fn = optarg;
        break;
      case 'm': 
        free_mem  = 1;
        break;
      case 't': 
        nthrds  = atoi(optarg);
        break;
      case 'v': 
        verbose = atoi(optarg);
        if (verbose > 2)
          verbose = 2;
        break;
      case 'c': 
        validate_results  = 1;
        break;
      case 'p': 
        write_pbm = 1;
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

  struct timeval t_load_start;
  gettimeofday(&t_load_start, NULL);

  sm_t *M = NULL;
  M = load_jcf_matrix(fn, verbose);

  // print walltime
  if (verbose > 1) {
    printf("Time for loading JCF matrix: %7.3f sec\n", walltime(t_load_start) / (1000000));
  }

  if (write_pbm) {
    const char *pbm_fn  = "input-matrix.pbm";
    write_jcf_matrix_to_pbm(M, pbm_fn, verbose);
  }

  // computing Gaussian Elimination of A using methods of Faugère & Lachartre
  elim_fl(M);

  free(M);
  return 0;
}
