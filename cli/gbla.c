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



#include "gbla.h"
#include <getopt.h>

#define __GBLA_CLI_DEBUG        0
#define __GBLA_CLI_DEBUG_1      0
#define __GBLA_CLI_DEBUG_A      0
#define __GBLA_CLI_DEBUG_D      0
#define __GBLA_CLI_DEBUG_D_TEST 0

// enable this only for some untested hybrid, etc. versions.
// NOTE: these implementations are not finished, not stable and you cannot rely
// on them.
#define __GBLA_UNTESTED_VERSIONS  0

#define SPARSE        0
#define INTERMEDIATE  0
#define DENSE         1

void print_help() {
  printf("\n");
  printf("NAME\n");
  printf("    gbla - Linear algebra for Groebner basis computations\n");
  printf("\n");
  printf("SYNOPSIS\n");
  printf("    gbla [options] [file..]\n");
  printf("\n");
  printf("DESCRIPTION\n");
  printf("    Computes Gaussian Elimination of a structured hybrid sparse-\n");
  printf("    dense matrix coming from Groebner basis computations.\n");
  printf("\n");
  printf("OPTIONS\n");
  printf("    -b          Dimensions of a block in submatrices.\n");
  printf("                Default: 256.\n");
  printf("                NOTE: The block dimension must be a multiple of __GBLA_LOOP_UNROLL_BIG\n");
  printf("                      which is by default 64. You can reset __GBLA_LOOP_UNROLL_BIG in \n");
  printf("                      src/gb_config.h.in (you need to rebuild the library afterwards).\n");
  printf("    -d          Reads in different matrix format (at the moment only for matrices arising\n");
  printf("                from the Prym-Green conjecture. \n");
  printf("    -f          Free memory on the go.\n");
  printf("                Default: 1, memory is freed on the go.\n");
  printf("                Use 0 for keeping complete memory until the end.\n");
  printf("    -g GIT      Outputs git commit hash if verbosity level is >0.\n");
  printf("    -h          Print help.\n");
  /*
  printf("    -m          Number of rows per multiline.\n");
  printf("                Default: 1.\n");
  */
  printf("    -n          Use the new format for matrices.\n");
  printf("    -o          Use third party dense reducer for D. By default this is FFLAS-FFPACK\n");
  printf("    -p          Writes intermediate matrices in pbm format.\n");
  printf("    -r          Compute a REDUCED row echelon form.\n");
  printf("    -s          Splicing options:\n");
  printf("                0: standard Faugère-Lachartre block method.\n");
  printf("                1: not reducing A, for rank profiling.\n");
  printf("                2: standard Faugère-Lachartre block method, completely multiline.   (OLD VERSION)\n");
  printf("                3: A and C are multiline, B and D are blocks, completely multiline. (OLD VERSION)\n");
  printf("                Default: 1.\n");
  printf("    -t THRDS    Number of threads used.\n");
  printf("                Default: 1.\n");
  printf("    -v LEVEL    Level of verbosity:\n");
  printf("                1 -> only error messages printed\n");
  printf("                2 -> some meta information printed\n");
  printf("                3 -> additionally meta information printed\n");
  printf("                Note: Everything >2 is time consuming and slows down\n");
  printf("                      the overall computations.\n");
  printf("\n");

  return;
}


int main(int argc, char *argv[]) {
  const char *fn        = NULL;
  int free_mem          = 1;
  int reduce_completely = 0;
  int verbose           = 0;
	/* int method            = 0; */
  int nrows_multiline   = 1;
  int block_dimension   = 256;
  int write_pbm         = 0;
  int nthreads          = 1;
  int splicing          = 1;
  int dense_reducer     = 0;
  int new_format        = 0;
  int schreyer_matrix   = 0;
  int git_hash          = 0;

  int index;
  int opt;

  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;

	opterr  = 0;

  while ((opt = getopt(argc, argv, "db:f:ghm:t:v:oprs:n")) != -1) {
    switch (opt) {
      case 'b':
        block_dimension = atoi(optarg);
        if (block_dimension < 1)
          block_dimension = 256;
        break;
      case 'd':
        schreyer_matrix = 1;
        break;
      case 'f':
        free_mem = atoi(optarg);
        if (free_mem < 0)
          free_mem  = 1;
        break;
      case 'g':
        git_hash  = 1;
        break;
      case 'h':
        print_help();
        return 0;
      case 'm':
        nrows_multiline = atoi(optarg);
        if (nrows_multiline < 1)
          nrows_multiline = 1;
        break;
      case 'o':
        dense_reducer = 1;
        break;
      case 'p':
        write_pbm = 1;
        break;
      case 't':
        nthreads  = atoi(optarg);
        break;
      case 'r':
        reduce_completely = 1;
        break;
      case 's':
        splicing = atoi(optarg);
        if (splicing < 0)
          splicing  = 1;
        if (splicing > 7)
          splicing  = 1;
        break;
      case 'v':
        verbose = atoi(optarg);
        if (verbose > 3)
          verbose = 2;
        break;
      case 'n':
	      new_format = 1;
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

  if (reduce_completely == 1)
    splicing  = 2;

  sm_t *M = NULL;

  if (verbose > 1) {
    printf("---------------------------------------------------------------------\n");
    printf("------------- Computing an FL-style Gaussian Elimination ------------\n");
    printf("------------------ with the following options set -------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("number of threads           %4d\n", nthreads);
    printf("splicing option             %4d\n", splicing);
    printf("third party dense reducer   %4d\n", dense_reducer);
    printf("dimension of blocks         %4d\n", block_dimension);
    printf("free memory on the go       %4d\n", free_mem);
    printf("reduced row echelon form    %4d\n", reduce_completely);
    printf("write PBM file              %4d\n", write_pbm);
    printf("---------------------------------------------------------------------\n");
  }


  if (verbose > 0) {
    if (git_hash == 1) {
      printf("---------------------------------------------------------------------\n");
      FILE *fp;
      char hash_string[200];
      fp  = popen("git rev-parse HEAD", "r");
      if (fp == NULL) {
        printf("Failed to run command\n");
        exit(1);
      }
      printf("git commit hash             ");
      while (fgets(hash_string, sizeof(hash_string)-1, fp) != NULL)
        printf("%s", hash_string);
      pclose(fp);
    }
    printf("---------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Loading matrix ...");
    fflush(stdout);
  }
  if (schreyer_matrix) {
    /*  load Schreyer matrix */
    M = load_schreyer_matrix(fn, verbose);
  } else {
    /*  load JCF matrix */
    M = load_jcf_matrix(fn, verbose, new_format, nthreads);
  }
  if (verbose > 0) {
    printf("%9.3f sec (%.3f %s/sec)\n",
        walltime(t_load_start) / (1000000),
        M->fs / (walltime(t_load_start) / (1000000)), M->fsu);
  }  
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("Data for %s\n", fn);
    printf("---------------------------------------------------------------------\n");
    printf("modulus                     %14d\n", (uint32_t) M->mod);
    printf("number of rows              %14d\n", (uint32_t)M->nrows);
    printf("number of columns           %14d\n", (uint32_t)M->ncols);
    printf("number of nonzero elements  %14ld\n", M->nnz);
    printf("density                     %14.2f %%\n", M->density);
    printf("size                        %14.2f %s\n", M->fs, M->fsu);
    printf("---------------------------------------------------------------------\n");
  }

  /*  write JCF matrix to PBM file */
  if (write_pbm) {
    if (verbose > 0) {
      gettimeofday(&t_load_start, NULL);
      printf("%-38s","Wrintg matrix to PBM file ...");
      fflush(stdout);
    }
    const char *pbm_fn  = "input-matrix.pbm";
    write_jcf_matrix_to_pbm(M, pbm_fn, verbose);
    if (verbose > 0) {
      printf("t%9.3f sec (%.3f %s/sec)\n",
          walltime(t_load_start) / (1000000),
          M->fs / (walltime(t_load_start) / (1000000)), M->fsu);
    }
  }
  if (schreyer_matrix == 1) {
    M = sort_schreyer_matrix(M);
    normalize_schreyer_input_rows(M);
    if (write_pbm) {
      if (verbose > 1) {
        gettimeofday(&t_load_start, NULL);
        printf("---------------------------------------------------------------------\n");
        printf(">>>>\tSTART writing input matrix to PBM file ...\n");
      }
      const char *pbm_fn  = "input-2-matrix.pbm";
      write_jcf_matrix_to_pbm(M, pbm_fn, verbose);
      if (verbose > 1) {
        printf("<<<<\tDONE writing input matrix to PBM file.\n");
        // print walltime
        printf("TIME\t%.3f sec (%.3f %s/sec)\n",
            walltime(t_load_start) / (1000000),
            M->fs / (walltime(t_load_start) / (1000000)), M->fsu);
        print_mem_usage();
        printf("---------------------------------------------------------------------\n");
        printf("\n");
      }
    }
  }
  /*
  printf("\n");
  for (int ii=0; ii<M->nrows; ++ii) {
    printf("ROW %d\n",ii);
      printf("%u || ", M->pos[ii][0]);
      for (int jj=0; jj<M->rwidth[ii]; ++jj)
        printf("%u  ", M->rows[ii][M->pos[ii][jj]]);
    printf("\n");
  }
  */
  /*  track time for the complete reduction process (excluding load) */
  if (verbose > 0)
    gettimeofday(&t_complete, NULL);


  switch (splicing) {
    // A, C sparse, B & D dense matrices, no multilines
    case 0:
      if (fl_block_sparse_dense(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
    // A, C sparse, B & D dense matrices, no multilines, directly reducing C
    case 1:
      if (fl_block_sparse_dense_keep_A(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
    /*  all submatrices are of block type */
    case 2:
      if (fl_block(M, block_dimension, nrows_multiline, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
    /*  submatrices A and C of multiline type, submatrices B and D are of block type */
    case 3:
      if (fl_ml_A_C(M, block_dimension, nrows_multiline, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in A and C multiline,\n",fn);
        printf("B and D block type mode.\n");
      }
      break;
    // A, C sparse, B & D dense matrices, no multilines, new draft for v0.3
      /*
    case 4:
      if (fl_v03(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
      */
    /*
     * ------------------------------------------------------------------------------
     * untested versions - start
     * ------------------------------------------------------------------------------
     */
#if __GBLA_UNTESTED_VERSIONS
    // all submatrices are of small dense block type without multilines
    case 4:
      if (fl_block_dense_old(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
    // all submatrices are of small hybrid block type without multilines and
    // inner sub blocks of that are either dense or NULL
    case 5:
      if (fl_block_hybrid(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
    // A & C hybrid, B & D dense matrices, no multilines
    case 6:
      if (fl_block_hybrid_dense(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
    // A sparse, C, B & D dense matrices, no multilines
    case 7:
      if (fl_block_sparse_dense(M, nthreads, free_mem, verbose,
            reduce_completely, dense_reducer)) {
        printf("Error while trying to eliminate matrix from file '%s' in all block type mode.\n",fn);
      }
      break;
#endif
     /*
     * ------------------------------------------------------------------------------
     * untested versions - end
     * ------------------------------------------------------------------------------
     */
    default:
      abort ();
  }

	/* freeing memory */
	ri_t	ii = 0 ;
	for ( ; ii < M->nrows ; ++ii) {
		if (M->rows[ii] != NULL)
		free(M->rows[ii]);
		if (M->pos[ii] != NULL)
		free(M->pos[ii]);
	}

  free(M->rows);
  free(M->pos);
  free(M->rwidth);
  free(M);
  M = NULL;

  if (verbose > 0) {
    printf("-------------------------------------------------------------------\n");
    printf("%-38s","Reduction completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    if (verbose > 1) 
      print_mem_usage();
  }
    printf("-------------------------------------------------------------------\n");
  return 0;
}
#if 0
int fl_v03(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer)
{

#if __GBLA_HAVE_OPENCL
  const cl_long n   = 1000;
  const int ntrips  = 100;

  cl_context context;
  cl_command_queue queue;
  create_context_on(CHOOSE_INTERACTIVELY, CHOOSE_INTERACTIVELY, 0, &context, &queue, 0);

  if (verbose > 1)
    print_device_info_from_queue(queue);

  // --------------------------------------------------------------------------
  // load kernels
  // --------------------------------------------------------------------------
  // kernels are loaded relative to the binary: even if we are in cli/ the
  // binary lives in ../. for building in other directories we have to take care
  // of where to put the kernels.
  char *knl_text = read_file("src/vec-add-soln.cl");
  cl_kernel knl = kernel_from_string(context, knl_text, "sum", NULL);
  free(knl_text);

  printf("%l | %d\n", n, ntrips);
  // --------------------------------------------------------------------------
  // allocate and initialize CPU memory
  // --------------------------------------------------------------------------
  float *a = (float *) malloc(sizeof(float) * n);
  if (!a) { perror("alloc x"); abort(); }
  float *b = (float *) malloc(sizeof(float) * n);
  if (!b) { perror("alloc y"); abort(); }
  float *c = (float *) malloc(sizeof(float) * n);
  if (!c) { perror("alloc z"); abort(); }

  for (size_t i = 0; i < n; ++i)
  {
    a[i] = i;
    b[i] = 2*i;
  }

  // --------------------------------------------------------------------------
  // allocate device memory
  // --------------------------------------------------------------------------
  cl_int status;
  cl_mem buf_a = clCreateBuffer(context, CL_MEM_READ_WRITE,
      sizeof(float) * n, 0, &status);
  CHECK_CL_ERROR(status, "clCreateBuffer");

  cl_mem buf_b = clCreateBuffer(context, CL_MEM_READ_WRITE,
      sizeof(float) * n, 0, &status);
  CHECK_CL_ERROR(status, "clCreateBuffer");

  cl_mem buf_c = clCreateBuffer(context, CL_MEM_READ_WRITE,
      sizeof(float) * n, 0, &status);
  CHECK_CL_ERROR(status, "clCreateBuffer");

  // --------------------------------------------------------------------------
  // transfer to device
  // --------------------------------------------------------------------------
  CALL_CL_GUARDED(clEnqueueWriteBuffer, (
        queue, buf_a, /*blocking*/ CL_TRUE, /*offset*/ 0,
        n * sizeof(float), a,
        0, NULL, NULL));

  CALL_CL_GUARDED(clEnqueueWriteBuffer, (
        queue, buf_b, /*blocking*/ CL_TRUE, /*offset*/ 0,
        n * sizeof(float), b,
        0, NULL, NULL));

  // --------------------------------------------------------------------------
  // run code on device
  // --------------------------------------------------------------------------

  CALL_CL_GUARDED(clFinish, (queue));

  timestamp_type time1, time2;
  get_timestamp(&time1);

  for (int trip = 0; trip < ntrips; ++trip)
  {
    SET_4_KERNEL_ARGS(knl, buf_a, buf_b, buf_c, n);
    size_t ldim[] = { 32 };
    size_t gdim[] = { ((n + ldim[0] - 1)/ldim[0])*ldim[0] };
    CALL_CL_GUARDED(clEnqueueNDRangeKernel,
        (queue, knl,
         /*dimensions*/ 1, NULL, gdim, ldim,
         0, NULL, NULL));
  }

  CALL_CL_GUARDED(clFinish, (queue));

  get_timestamp(&time2);
  double elapsed = timestamp_diff_in_seconds(time1,time2)/ntrips;
  printf("%f s\n", elapsed);
  printf("%f GB/s\n",
      3*n*sizeof(float)/1e9/elapsed);

  // --------------------------------------------------------------------------
  // transfer back & check
  // --------------------------------------------------------------------------
  CALL_CL_GUARDED(clEnqueueReadBuffer, (
        queue, buf_c, /*blocking*/ CL_TRUE, /*offset*/ 0,
        n * sizeof(float), c,
        0, NULL, NULL));

  for (size_t i = 0; i < n; ++i)
    if (c[i] != 3*i)
    {
      printf("BAD %ld %f %f!\n", i, c[i], c[i] - 3*i);
      abort();
    }
  puts("GOOD");

  // --------------------------------------------------------------------------
  // clean up
  // --------------------------------------------------------------------------
  CALL_CL_GUARDED(clReleaseMemObject, (buf_a));
  CALL_CL_GUARDED(clReleaseMemObject, (buf_b));
  CALL_CL_GUARDED(clReleaseMemObject, (buf_c));
  CALL_CL_GUARDED(clReleaseKernel, (knl));
  CALL_CL_GUARDED(clReleaseCommandQueue, (queue));
  CALL_CL_GUARDED(clReleaseContext, (context));
#endif
  return 0;
}
#endif

int fl_block_sparse_dense_keep_A(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  // all submatrices of block type
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s", "Splicing of input matrix ...");
    fflush(stdout);
  }
  // construct splicing of matrix M into A, B, C and D
  sm_fl_t *A      = (sm_fl_t *)malloc(sizeof(sm_fl_t));
  sb_fl_t *B      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  sm_fl_t *C      = (sm_fl_t *)malloc(sizeof(sm_fl_t));
  dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD

  /*
  if (nthreads > 8)
    splice_fl_matrix_sparse_dense_keep_A(M, A, B, C, D, map, 0, free_mem, verbose, 8);
  else
  */
    splice_fl_matrix_sparse_dense_keep_A(M, A, B, C, D, map, 0, free_mem, verbose, nthreads);

  if (verbose > 0) {
    //printf("<<<<\tDONE  splicing and mapping of input matrix.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("---------------------------------------------------------------------\n");
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG_A
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  int ii,jj,kk,ll;
  for (ii=0; ii<A->nrows; ++ii) {
    if (A->row[ii] != NULL) {
      printf("%d|||\n", ii);
      for (jj=0; jj<A->sz[ii]; ++jj) {
        printf("%u | %u || ", A->row[ii][jj], A->pos[ii][jj]);
      }
      printf("\n");
    }
  }
  printf("==================================\n");
  for (ii=0; ii<C->nrows; ++ii) {
    if (C->row[ii] != NULL) {
      printf("%d|||\n", ii);
      for (jj=0; jj<C->sz[ii]; ++jj) {
        printf("%u | %u || ", C->row[ii][jj], C->pos[ii][jj]);
      }
      printf("\n");
    }
  }
#endif
  // reducing submatrix A using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s","Storing A in C ...");
    fflush(stdout);
  }
  if (elim_fl_C_sparse_dense_keep_A(C, &A, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 0) {
    //printf("<<<<\tDONE  reducing A.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    //printf("---------------------------------------------------------------------\n");
    //printf("\n");
  }
  /*  copying sparse matrix C to a sparse block matrix C_block */
#if SPARSE
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s", "Copying C to sparse block representation ...");
    fflush(stdout);
  }
  sb_fl_t *C_block = copy_sparse_to_block_matrix(C, nthreads);
  free_sparse_matrix(&C,nthreads);
  if (verbose > 0) {
    //printf("<<<<\tDONE  copying sparse C to sparse block representation.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s","Reducing C to zero ...");
    fflush(stdout);
  }
  if (elim_fl_C_sparse_sparse_block(B, &C_block, D, 0, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 0) {
    //printf("<<<<\tDONE  reducing C to zero.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#endif
#if INTERMEDIATE
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s", "Copying C to intermediate block representation ...");
    fflush(stdout);
  }
  ibm_fl_t *C_block = copy_sparse_to_intermediate_block_matrix(C, nthreads);
  free_sparse_matrix(&C,nthreads);
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s","Reducing C to zero ...");
    fflush(stdout);
  }
  if (elim_fl_C_intermediate_block(B, &C_block, D, 1, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 0) {
    //printf("<<<<\tDONE  reducing C to zero.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#endif
#if DENSE
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s", "Transforming C ...");
    fflush(stdout);
  }
  dbm_fl_t *C_block = copy_sparse_to_dense_block_matrix(C, nthreads);
  free_sparse_matrix(&C,nthreads);
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s","Reducing C to zero ...");
    fflush(stdout);
  }
  if (elim_fl_C_dense_sparse_block(B, &C_block, D, 1, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 0) {
    //printf("<<<<\tDONE  reducing C to zero.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#endif
#if __GBLA_CLI_DEBUG_D_TEST
  printf("DDDD\n");
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlD; ++ii) {
    for (int jj=0; jj<clD; ++jj) {
      if (D->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", D->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif

  // copy block D to dense wide (re_l_t) representation
  dm_t *D_red = copy_block_to_dense_matrix(&D, nthreads, 1);
  D_red->mod  = M->mod;

  // eliminate D_red using a structured Gaussian Elimination process on the rows
  ri_t rank_D = 0;
  // echelonizing D to zero using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (D_red->nrows > 0)
    rank_D = elim_fl_dense_D(D_red, nthreads);
  if (verbose > 0) {
    printf("%9.3f sec (rank D: %u)\n",
        walltime(t_load_start) / (1000000), rank_D);
  }
  if (verbose > 1) {
    print_mem_usage();
  }
#if __GBLA_CLI_DEBUG_D_TEST
  for (int ii=0; ii<D_red->nrows; ++ii) {
    printf("ROW %d\n",ii);
    if (D_red->row[ii]->piv_val == NULL)
      printf("NULL!");
    else {
      for (int jj=0; jj<D_red->ncols; ++jj)
        printf("%lu  ", D_red->row[ii]->piv_val[jj]);
    }
    printf("\n");
  }
#endif
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reconstructing M ...");
    fflush(stdout);
  }
  reconstruct_matrix_no_multiline_keep_A(M, A, B, D_red, map, nthreads);
  if (verbose > 0) {
    printf("%9.3f sec (rank M: %u)\n",
        walltime(t_load_start) / (1000000), M->nrows);
  }
  if (verbose > 1) {
    print_mem_usage();
  }
  return 0;
}

int fl_block_sparse_dense(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer)
{
  struct timeval t_load_start;
  // all submatrices of block type
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s", "Splicing of input matrix ...");
    fflush(stdout);
  }
  // construct splicing of matrix M into A, B, C and D
  sb_fl_t *A      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *B     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  sb_fl_t *C      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD

  splice_fl_matrix_sparse_dense_2(M, A, B, C, D, map, 0, free_mem, verbose, nthreads);

  // free elements in M, but keep general meta data for later references
	ri_t ii;
	for (ii=0; ii < M->nrows; ++ii) {
		if (M->rows[ii] != NULL)
		free(M->rows[ii]);
		if (M->pos[ii] != NULL)
		free(M->pos[ii]);
	}
  free(M->rows);
  free(M->pos);

  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG_A
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  printf("============= AAAAAAAAAAAAAAAA ===================\n");
  for (int ii=0; ii<rlA; ++ii) {
    for (int jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj].val != NULL) {
        printf("%d .. %d\n", ii, jj);
        for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
          for (int ll=0; ll<A->blocks[ii][jj].sz[kk]; ++ll) {
            printf("%d | %d || ", A->blocks[ii][jj].val[kk][ll], A->blocks[ii][jj].pos[kk][ll]);
          }
          printf("\n");
        }
      }
    }
  }
  printf("==================================\n");
#endif
#if __GBLA_CLI_DEBUG
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj].val != NULL) {
        printf("%d .. %d\n", ii, jj);
        if (ii != jj) {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", A->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
        } else {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE_DIAG; ++kk) {
            printf("%d | ", A->blocks[ii][jj].val[kk]);
          }
          printf("\n");
        }
        printf("\n");
      }
    }
  }
#endif
  // reducing submatrix A using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing A ...");
    fflush(stdout);
  }
  if (elim_fl_A_sparse_dense_block(&A, B, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
  }
#if __GBLA_CLI_DEBUG_11
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", B->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif

  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing C ...");
    fflush(stdout);
  }
  if (elim_fl_C_sparse_dense_block(B, &C, D, 1, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
  }

#if __GBLA_CLI_DEBUG_D_TEST
  printf("DDDD\n");
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlD; ++ii) {
    for (int jj=0; jj<clD; ++jj) {
      if (D->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", D->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif

  // copy block D to dense wide (re_l_t) representation
  dm_t *D_red = copy_block_to_dense_matrix(&D, nthreads, 1);
  D_red->mod  = M->mod;

  // eliminate D_red using a structured Gaussian Elimination process on the rows
  ri_t rank_D = 0;
  // echelonizing D to zero using methods of Faugère & Lachartre
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (D_red->nrows > 0)
    rank_D = elim_fl_dense_D(D_red, nthreads);
  if (verbose > 0) {
    printf("%9.3f sec (rank D: %u)\n",
        walltime(t_load_start) / (1000000), rank_D);
  }
  if (verbose > 1) {
    print_mem_usage();
  }
#if __GBLA_CLI_DEBUG_D_TEST
  for (int ii=0; ii<D_red->rank; ++ii) {
    printf("%u | ",D_red->row[ii]->piv_lead);
  }
  printf("\n");
  for (int ii=0; ii<D_red->rank; ++ii) {
    printf("ROW %d\n",ii);
    if (D_red->row[ii]->piv_val == NULL)
      printf("NULL!");
    else {
      printf("%u || ", D_red->row[ii]->piv_lead);
      for (int jj=0; jj<D_red->ncols; ++jj)
#if defined(GBLA_USE_UINT16) || defined(GBLA_USE_UINT32)
        printf("%u  ", D_red->row[ii]->piv_val[jj]);
#else
        printf("%.0f  ", D_red->row[ii]->piv_val[jj]);
#endif
    }
    printf("\n");
  }
#endif
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reconstructing M ...");
    fflush(stdout);
  }
  reconstruct_matrix_block_no_multiline(M, A, B, D_red, map, nthreads);
  if (verbose > 0) {
    printf("%9.3f sec (rank M: %u)\n",
        walltime(t_load_start) / (1000000), M->nrows);
  }
  if (verbose > 1) {
    print_mem_usage();
  }
  return 0;
}

int fl_block(sm_t *M, int block_dimension, int nrows_multiline, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  /*  all submatrices of block type */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s", "Splicing of input matrix ...");
    fflush(stdout);
  }
  /*  construct splicing of matrix M into A, B, C and D */
  sbm_fl_t *A     = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  sbm_fl_t *B     = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  sbm_fl_t *C     = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  sbm_fl_t *D     = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); /*  stores mappings from M <-> ABCD */
  map_fl_t *map_D = (map_fl_t *)malloc(sizeof(map_fl_t)); /*  stores mappings for reduced D */

  splice_fl_matrix(M, A, B, C, D, map,
      M->nrows, M->ncols, block_dimension,
      nrows_multiline, nthreads, free_mem, verbose, 0);

  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG
  /*  column loops */
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / C->bwidth);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);
  /*  row loops */
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / A->bheight);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / C->bheight);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj] != NULL) {
        for (kk=0; kk<block_dimension/2; ++kk) {
          printf("%d .. %d .. %d\n",ii,jj,kk);
          if (A->blocks[ii][jj][kk].dense == 0) {
            for (ll=0; ll<A->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", A->blocks[ii][jj][kk].idx[ll]);
              printf("%d %d ", A->blocks[ii][jj][kk].val[2*ll], A->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          } else {
            for (ll=0; ll<A->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", ll);
              printf("%d %d ", A->blocks[ii][jj][kk].val[2*ll], A->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          }
        }
      }
    }
  }
#endif

  /*  reducing submatrix A using methods of Faugère & Lachartre */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing A ...");
    fflush(stdout);
  }
  if (elim_fl_A_block(&A, B, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
  }
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
 printf("++++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj] != NULL) {
        for (int kk=0; kk<block_dimension/2; ++kk) {
          printf("%d .. %d .. %d\n",ii,jj,kk);
          if (B->blocks[ii][jj][kk].dense == 0) {
            for (int ll=0; ll<B->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", B->blocks[ii][jj][kk].idx[ll]);
              printf("%d %d ", B->blocks[ii][jj][kk].val[2*ll], B->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          } else {
            for (int ll=0; ll<B->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", ll);
              printf("%d %d ", B->blocks[ii][jj][kk].val[2*ll], B->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          }
        }
      }
    }
  }
#endif

  /*  reducing submatrix C to zero using methods of Faugère & Lachartre */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing C ...");
    fflush(stdout);
  }
  if (elim_fl_C_block(B, &C, D, 1, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
  }

#if __GBLA_CLI_DEBUG_D_TEST
  printf("DDDD\n");
  // column loops
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);
  // row loops
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);
 printf("++++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlD; ++ii) {
    for (int jj=0; jj<clD; ++jj) {
      if (D->blocks[ii][jj] != NULL) {
        for (int kk=0; kk<block_dimension/2; ++kk) {
          printf("%d .. %d .. %d\n",ii,jj,kk);
          if (D->blocks[ii][jj][kk].dense == 0) {
            for (int ll=0; ll<D->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", D->blocks[ii][jj][kk].idx[ll]);
              printf("%d %d ", D->blocks[ii][jj][kk].val[2*ll], D->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          } else {
            for (int ll=0; ll<D->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", ll);
              printf("%d %d ", D->blocks[ii][jj][kk].val[2*ll], D->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          }
        }
      }
    }
  }
#endif

  /*  echelonize D using methods of Faugère & Lachartre */

  /*  We need to do at least two rounds on D, possibly even reducing B further */
  /*  with the reduced D later on. For this purpose we use a multiline version of */
  /*  D as output of the Gaussian elimination. */
  sm_fl_ml_t *D_red = (sm_fl_ml_t *)malloc(sizeof(sm_fl_ml_t));
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  ri_t rank_D;
#if GBLA_WITH_FFLAS
  if (dense_reducer == 1)
    rank_D = elim_fl_D_fflas_ffpack(D, M->mod, nthreads);
  else
#endif
    rank_D = elim_fl_D_block(D, D_red, M->mod, nthreads);
  if (rank_D == (ri_t)-1) {
    printf("Error while reducing D to upper triangular matrix.\n");
    return 1;
  }
  if (verbose > 0) {
    printf("%9.3f sec (rank D: %u)\n",
        walltime(t_load_start) / (1000000), rank_D);
  }
  if (verbose > 1) {
    print_mem_usage();
  }
#if __GBLA_CLI_DEBUG_D_TEST
  printf("DDDD11111111111111\n");
  for (int ii=0; ii<D_red->nrows / 2; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", D_red->ml[ii].sz * 2);
    if (D_red->ml[ii].sz>0) {
      if (D_red->ml[ii].dense == 0) {
        for (int ll=0; ll<D_red->ml[ii].sz; ++ll) {
          printf("%d -- ", D_red->ml[ii].idx[ll]);
          printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
        }
        printf("\n");
      } else {
        for (int ll=0; ll<D_red->ml[ii].sz; ++ll) {
          printf("%d -- ", ll);
          printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
        }
      }
    }
  }
#endif
  ri_t rank_M = map->npiv + rank_D;
  M->nrows  = rank_M;
  if (reduce_completely == 0) {
    if (verbose > 0) {
      gettimeofday(&t_load_start, NULL);
      printf("%-38s","Reconstructing M ...");
      fflush(stdout);
    }
    process_matrix(D_red, map_D, block_dimension);
    combine_maps(map, &map_D, M->ncols, D_red->ncols, 1);
    reconstruct_matrix_block(M, A, B, D_red, map, M->ncols, 1, 1, 1, 0, nthreads);
    if (verbose > 0) {
      printf("%9.3f sec (rank M: %u)\n",
          walltime(t_load_start) / (1000000), M->nrows);
    }
    if (verbose > 1) {
      print_mem_usage();
    }
  } else { /*  compute reduced row echelon form of input matrix */
    if (verbose > 0) {
      printf("Starting second iteration of complete reduction\n");
    }
    process_matrix(D_red, map_D, block_dimension);
    sm_t *BD        = (sm_t *)malloc(sizeof(sm_t));

    /*  copy block matrix B back to a sparse matrix in order to use */
    /*  splice_fl_matrix() procedure again */
    copy_block_ml_matrices_to_sparse_matrix(&B, &D_red, rank_D, &BD, 1, nthreads);

    map_fl_t *map2  = (map_fl_t *)malloc(sizeof(map_fl_t));
    construct_fl_map_reduced(map2, map_D, BD->nrows, rank_D, BD->ncols, nthreads);
    sbm_fl_t *B1    = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
    sbm_fl_t *B2    = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
    sbm_fl_t *D1    = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
    sbm_fl_t *D2    = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));

    splice_fl_matrix(BD, D1, D2, B1, B2, map2, BD->nrows, BD->ncols,
        block_dimension, nrows_multiline, nthreads, free_mem, verbose, 1);
    if (verbose > 0) {
      gettimeofday(&t_load_start, NULL);
      printf("%-38s","Reducing D1 ...");
      fflush(stdout);
    }
    if (elim_fl_A_block(&D1, D2, M->mod, nthreads)) {
      printf("Error while reducing D1.\n");
      return 1;
    }
    if (verbose > 0) {
      printf("%9.3f sec\n",
          walltime(t_load_start) / (1000000));
    }
    if (verbose > 1) {
      print_mem_usage();
    }

    /*  reducing submatrix C to zero using methods of Faugère & Lachartre */
    if (verbose > 0) {
      gettimeofday(&t_load_start, NULL);
      printf("%-38s","Reducing D2 ...");
      fflush(stdout);
    }
    if (elim_fl_C_block(D2, &B1, B2, 1, M->mod, nthreads)) {
      printf("Error while reducing D1.\n");
      return 1;
    }
    if (verbose > 0) {
      printf("%9.3f sec\n",
          walltime(t_load_start) / (1000000));
    }
    if (verbose > 1) {
      print_mem_usage();
    }

    /*  reconstruct matrix */
    if (verbose > 0) {
      gettimeofday(&t_load_start, NULL);
      printf("%-38s","Reconstructing M ...");
      fflush(stdout);
    }
    combine_maps(map, &map_D, M->ncols, BD->ncols, 0);
    reconstruct_matrix_block_reduced(M, A, B2, D2, map, M->ncols, 1, 1, 1, 0, nthreads);
    if (verbose > 0) {
      printf("%9.3f sec (rank M: %u)\n",
          walltime(t_load_start) / (1000000), M->nrows);
    }
    if (verbose > 1) {
      print_mem_usage();
    }
  }

  return 0;
}

int fl_ml_A_C(sm_t *M, int block_dimension, int nrows_multiline, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  /*  submatrices A and C of multiline type, B and D of block type */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s", "Splicing of input matrix ...");
    fflush(stdout);
  }
  /*  construct splicing of matrix M into A, B, C and D */
  sm_fl_ml_t *A   = (sm_fl_ml_t *)malloc(sizeof(sm_fl_ml_t));
  sbm_fl_t *B     = (sbm_fl_t   *)malloc(sizeof(sbm_fl_t));
  sm_fl_ml_t *C   = (sm_fl_ml_t *)malloc(sizeof(sm_fl_ml_t));
  sbm_fl_t *D     = (sbm_fl_t   *)malloc(sizeof(sbm_fl_t));
  map_fl_t *map   = (map_fl_t   *)malloc(sizeof(map_fl_t)); /*  stores mappings from M <-> ABCD */
  map_fl_t *map_D = (map_fl_t   *)malloc(sizeof(map_fl_t)); /*  stores mappings for reduced D */

  splice_fl_matrix_ml_A_C(M, A, B, C, D, map, block_dimension, nrows_multiline, nthreads,
      free_mem, verbose);

  if (verbose > 0) {
    //printf("<<<<\tDONE  splicing and mapping of input matrix.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("---------------------------------------------------------------------\n");
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG_A

  /*  column loops */
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);
  /*  row loops */
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_NROWS_MULTILINE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_NROWS_MULTILINE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", A->ml[ii].sz * 2);
    if (A->ml[ii].sz>0) {
      for (ll=0; ll<A->ml[ii].sz; ++ll) {
        printf("%d -- ", A->ml[ii].idx[ll]);
        printf("%d %d ", A->ml[ii].val[2*ll], A->ml[ii].val[2*ll+1]);
      }
      printf("\n");
    }
  }
  for (ii=0; ii<rlB; ++ii) {
    for (jj=0; jj<clB; ++jj) {
      for (kk=0; kk<block_dimension/2; ++kk) {
        printf("%d .. %d .. %d\n",ii,jj,kk);
        if (B->blocks[ii][jj] != NULL) {
          for (ll=0; ll<B->blocks[ii][jj][kk].sz; ++ll) {
            printf("%d -- ", B->blocks[ii][jj][kk].idx[ll]);
            printf("%d %d ", B->blocks[ii][jj][kk].val[2*ll], B->blocks[ii][jj][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
  }
  for (ii=0; ii<rlC; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", C->ml[ii].sz * 2);
    if (C->ml[ii].sz>0) {
      for (ll=0; ll<C->ml[ii].sz; ++ll) {
        printf("%d -- ", C->ml[ii].idx[ll]);
        printf("%d %d ", C->ml[ii].val[2*ll], C->ml[ii].val[2*ll+1]);
      }
      printf("\n");
    }
  }
  for (ii=0; ii<rlD; ++ii) {
    for (jj=0; jj<clD; ++jj) {
      for (kk=0; kk<block_dimension/2; ++kk) {
        printf("%d .. %d .. %d\n",ii,jj,kk);
        if (D->blocks[ii][jj] != NULL) {
          for (ll=0; ll<D->blocks[ii][jj][kk].sz; ++ll) {
            printf("%d -- ", D->blocks[ii][jj][kk].idx[ll]);
            printf("%d %d ", D->blocks[ii][jj][kk].val[2*ll], D->blocks[ii][jj][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
  }
#endif



  /*  reducing submatrix A using methods of Faugère & Lachartre */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s","Storing A in C ...");
    fflush(stdout);
  }
  if (elim_fl_C_ml(C, A, M->mod, nthreads)) {
    printf("Error storing data from A in C.\n");
    return 1;
  }
  if (verbose > 0) {
    //printf("<<<<\tDONE  reducing A.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    //printf("---------------------------------------------------------------------\n");
    //printf("\n");
  }
  /*  copying multiline matrix C to a block matrix C_block */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s", "Transforming C ...");
    fflush(stdout);
  }
  sbm_fl_t *C_block = copy_multiline_to_block_matrix_rl(&C, block_dimension, block_dimension, 1, nthreads);
  if (verbose > 0) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
  /*  reducing submatrix C to zero using methods of Faugère & Lachartre */
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    //printf("---------------------------------------------------------------------\n");
    printf("%-38s","Reducing C to zero ...");
    fflush(stdout);
  }
  if (elim_fl_C_block(B, &C_block, D, 0, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 0) {
    //printf("<<<<\tDONE  reducing C to zero.\n");
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#if __GBLA_CLI_DEBUG_D_TEST
  printf("DDDD\n");
  // column loops
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);
  // row loops
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);
 printf("++++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlD; ++ii) {
    for (int jj=0; jj<clD; ++jj) {
      if (D->blocks[ii][jj] != NULL) {
        for (int kk=0; kk<block_dimension/2; ++kk) {
          printf("%d .. %d .. %d\n",ii,jj,kk);
          if (D->blocks[ii][jj][kk].dense == 0) {
            for (int ll=0; ll<D->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", D->blocks[ii][jj][kk].idx[ll]);
              printf("%d %d ", D->blocks[ii][jj][kk].val[2*ll], D->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          } else {
            for (int ll=0; ll<D->blocks[ii][jj][kk].sz; ++ll) {
              printf("%d -- ", ll);
              printf("%d %d ", D->blocks[ii][jj][kk].val[2*ll], D->blocks[ii][jj][kk].val[2*ll+1]);
            }
            printf("\n");
          }
        }
      }
    }
  }
#endif
  /*  echelonize D using methods of Faugère & Lachartre */

  /*  We need to do at least two rounds on D, possibly even reducing B further */
  /*  with the reduced D later on. For this purpose we use a multiline version of */
  /*  D as output of the Gaussian elimination. */
  sm_fl_ml_t *D_red = (sm_fl_ml_t *)malloc(sizeof(sm_fl_ml_t));
  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  ri_t rank_D;
#if GBLA_WITH_FFLAS
  if (dense_reducer == 1)
    rank_D = elim_fl_D_fflas_ffpack(D, M->mod, nthreads);
  else
#endif
    rank_D = elim_fl_D_block(D, D_red, M->mod, nthreads);
  if (rank_D == (ri_t)-1) {
    printf("Error while reducing D to upper triangular matrix.\n");
    return 1;
  }
  if (verbose > 0) {
    printf("%9.3f sec (rank D: %u)\n",
        walltime(t_load_start) / (1000000), rank_D);
  }
  if (verbose > 1) {
    print_mem_usage();
  }

  ri_t rank_M = map->npiv + rank_D;
  M->nrows  = rank_M;

  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reconstructing M ...");
    fflush(stdout);
  }
  process_matrix(D_red, map_D, block_dimension);
  combine_maps(map, &map_D, M->ncols, D_red->ncols, 1);
  reconstruct_matrix_ml(M, A, B, D_red, map, M->ncols, 1, 1, 0, 0, nthreads);

  if (verbose > 0) {
    printf("%9.3f sec (rank M: %u)\n",
        walltime(t_load_start) / (1000000), M->nrows);
  }
  if (verbose > 1) {
    print_mem_usage();
  }

  return 0;
}

#if __GBLA_UNTESTED_VERSIONS
int fl_block_sparse_dense_old(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  // all submatrices of block type
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART splicing and mapping of input matrix ...\n");
  }
  // construct splicing of matrix M into A, B, C and D
  sb_fl_t *A      = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *B     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  dbm_fl_t *C     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD
  map_fl_t *map_D = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings for reduced D

  splice_fl_matrix_sparse_dense(M, A, B, C, D, map, 0, free_mem, verbose, nthreads);

  if (verbose > 1) {
    printf("<<<<\tDONE  splicing and mapping of input matrix.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj].val != NULL) {
        printf("%d .. %d\n", ii, jj);
        if (ii != jj) {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", A->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
        } else {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE_DIAG; ++kk) {
            printf("%d | ", A->blocks[ii][jj].val[kk]);
          }
          printf("\n");
        }
        printf("\n");
      }
    }
  }
#endif
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj] != NULL) {
        //printf("%d .. %d\n", ii, jj);
        for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
          if (B->blocks[ii][jj][kk] != NULL) {
            for (int ll=0; ll<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++ll) {
              printf("\n++ %d | %d || %d | %d ++\n", ii, jj, kk, ll);
              if (B->blocks[ii][jj][kk][ll].val != NULL) {
                for (int mm=0; mm<__GBLA_SIMD_INNER_SIZE; ++mm) {
                  printf("%d | ", B->blocks[ii][jj][kk][ll].val[mm]);
                }
              }
            }
            printf("\n");
          }
          printf("\n");
        }
      }
    }
  }
#endif
  // reducing submatrix A using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART reducing A ...\n");
  }
  if (elim_fl_A_sparse_dense_block(&A, B, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 1) {
    printf("<<<<\tDONE  reducing A.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", B->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif

  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART reducing C to zero ...\n");
  }
  if (elim_fl_C_dense_block(B, &C, D, 1, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 1) {
    printf("<<<<\tDONE  reducing C to zero.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlD; ++ii) {
    for (int jj=0; jj<clD; ++jj) {
      if (D->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", D->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif


  return 0;
}

int fl_block_hybrid_dense(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  // all submatrices of block type
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART splicing and mapping of input matrix ...\n");
  }
  // construct splicing of matrix M into A, B, C and D
  hbm_fl_t *A     = (hbm_fl_t *)malloc(sizeof(hbm_fl_t));
  dbm_fl_t *B     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  hbm_fl_t *C     = (hbm_fl_t *)malloc(sizeof(hbm_fl_t));
  dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD
  map_fl_t *map_D = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings for reduced D

  splice_fl_matrix_hybrid_dense(M, A, B, C, D, map, 0, free_mem, verbose, nthreads);

  if (verbose > 1) {
    printf("<<<<\tDONE  splicing and mapping of input matrix.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj].val != NULL) {
        printf("%d .. %d\n", ii, jj);
        if (ii != jj) {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", A->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
        } else {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE_DIAG; ++kk) {
            printf("%d | ", A->blocks[ii][jj].val[kk]);
          }
          printf("\n");
        }
        printf("\n");
      }
    }
  }
#endif
#if __GBLA_CLI_DEBUG_2
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj] != NULL) {
        //printf("%d .. %d\n", ii, jj);
        for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
          if (B->blocks[ii][jj][kk] != NULL) {
            for (int ll=0; ll<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++ll) {
              printf("\n++ %d | %d || %d | %d ++\n", ii, jj, kk, ll);
              if (B->blocks[ii][jj][kk][ll].val != NULL) {
                for (int mm=0; mm<__GBLA_SIMD_INNER_SIZE; ++mm) {
                  printf("%d | ", B->blocks[ii][jj][kk][ll].val[mm]);
                }
              }
            }
            printf("\n");
          }
          printf("\n");
        }
      }
    }
  }
#endif
  // reducing submatrix A using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART reducing A ...\n");
  }
  if (elim_fl_A_hybrid_dense_block(&A, B, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 1) {
    printf("<<<<\tDONE  reducing A.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }

  return 0;
}

int fl_block_hybrid(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  // all submatrices of block type
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART splicing and mapping of input matrix ...\n");
  }
  // construct splicing of matrix M into A, B, C and D
  hbm_fl_t *A     = (hbm_fl_t *)malloc(sizeof(hbm_fl_t));
  hbm_fl_t *B     = (hbm_fl_t *)malloc(sizeof(hbm_fl_t));
  hbm_fl_t *C     = (hbm_fl_t *)malloc(sizeof(hbm_fl_t));
  hbm_fl_t *D     = (hbm_fl_t *)malloc(sizeof(hbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD
  map_fl_t *map_D = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings for reduced D

  splice_fl_matrix_hybrid(M, A, B, C, D, map, 0, free_mem, verbose, nthreads);

  if (verbose > 1) {
    printf("<<<<\tDONE  splicing and mapping of input matrix.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj].val != NULL) {
        printf("%d .. %d\n", ii, jj);
        if (ii != jj) {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", A->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
        } else {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE_DIAG; ++kk) {
            printf("%d | ", A->blocks[ii][jj].val[kk]);
          }
          printf("\n");
        }
        printf("\n");
      }
    }
  }
#endif
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj] != NULL) {
        //printf("%d .. %d\n", ii, jj);
        for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
          if (B->blocks[ii][jj][kk] != NULL) {
            for (int ll=0; ll<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++ll) {
              printf("\n++ %d | %d || %d | %d ++\n", ii, jj, kk, ll);
              if (B->blocks[ii][jj][kk][ll].val != NULL) {
                for (int mm=0; mm<__GBLA_SIMD_INNER_SIZE; ++mm) {
                  printf("%d | ", B->blocks[ii][jj][kk][ll].val[mm]);
                }
              }
            }
            printf("\n");
          }
          printf("\n");
        }
      }
    }
  }
#endif
  // reducing submatrix A using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART reducing A ...\n");
  }
  if (elim_fl_A_hybrid_block(&A, B, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 1) {
    printf("<<<<\tDONE  reducing A.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }

  return 0;
}
int fl_block_dense(sm_t *M, int nthreads, int free_mem,
    int verbose, int reduce_completely, int dense_reducer) {
  struct timeval t_load_start;
  // all submatrices of block type
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART splicing and mapping of input matrix ...\n");
  }
  // construct splicing of matrix M into A, B, C and D
  dbm_fl_t *A     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  dbm_fl_t *B     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  dbm_fl_t *C     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  dbm_fl_t *D     = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  map_fl_t *map   = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings from M <-> ABCD
  map_fl_t *map_D = (map_fl_t *)malloc(sizeof(map_fl_t)); // stores mappings for reduced D

  splice_fl_matrix_dense(M, A, B, C, D, map, 0, free_mem, verbose, nthreads);

  if (verbose > 1) {
    printf("<<<<\tDONE  splicing and mapping of input matrix.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%9d x %9d]\n",
        A->nrows, A->ncols);
    printf("B [%9d x %9d]\n",
        B->nrows, B->ncols);
    printf("C [%9d x %9d]\n",
        C->nrows, C->ncols);
    printf("D [%9d x %9d]\n",
        D->nrows, D->ncols);
    printf("---------------------------------------------------------------------\n");
  }
#if __GBLA_CLI_DEBUG
  // column loops
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);

  int ii,jj,kk,ll;
  for (ii=0; ii<rlA; ++ii) {
    for (jj=0; jj<clA; ++jj) {
      if (A->blocks[ii][jj].val != NULL) {
        printf("%d .. %d\n", ii, jj);
        if (ii != jj) {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", A->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
        } else {
          for (kk=0; kk<__GBLA_SIMD_BLOCK_SIZE_DIAG; ++kk) {
            printf("%d | ", A->blocks[ii][jj].val[kk]);
          }
          printf("\n");
        }
        printf("\n");
      }
    }
  }
#endif
  // reducing submatrix A using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART reducing A ...\n");
  }
  if (elim_fl_A_dense_block(&A, B, M->mod, nthreads)) {
    printf("Error while reducing A.\n");
    return 1;
  }
  if (verbose > 1) {
    printf("<<<<\tDONE  reducing A.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlB; ++ii) {
    for (int jj=0; jj<clB; ++jj) {
      if (B->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", B->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif

  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("---------------------------------------------------------------------\n");
    printf(">>>>\tSTART reducing C to zero ...\n");
  }
  if (elim_fl_C_dense_block(B, &C, D, 1, M->mod, nthreads)) {
    printf("Error while reducing C.\n");
    return 1;
  }
  if (verbose > 1) {
    printf("<<<<\tDONE  reducing C to zero.\n");
    printf("TIME\t%.3f sec\n",
        walltime(t_load_start) / (1000000));
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("\n");
  }
#if __GBLA_CLI_DEBUG_1
  // column loops
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / __GBLA_SIMD_BLOCK_SIZE);
  // row loops
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / __GBLA_SIMD_BLOCK_SIZE);
  printf("+++$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
  for (int ii=0; ii<rlD; ++ii) {
    for (int jj=0; jj<clD; ++jj) {
      if (D->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
            for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
              printf("%d | ", D->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll]);
            }
            printf("\n");
          }
          printf("\n");
      }
    }
  }
#endif


  return 0;
}
#endif


/* vim:sts=2:sw=2:ts=2:
 */
