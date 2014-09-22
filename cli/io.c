#include "io.h"


// ========== TIMINGS ==========
double walltime(struct timeval t_start) {
  struct timeval t_end;
  gettimeofday(&t_end, NULL);
  return (double)((t_end.tv_sec - t_start.tv_sec) * 1000000 + t_end.tv_usec - t_start.tv_usec);
}


// ========== READING ==========
sm_t *load_jcf_matrix(const char *fn, int verbose) {

  // meta information of matrix
  double density;
  double fs;
  char *fsu;
  // timing structs
  struct timeval t_load_start;

  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
  }

  // start loading the matrix
  ri_t m;
  ci_t n;
  mod_t     mod;
  nnz_t     nnz;
  uint64_t  fl;

  // open in binary mode first to get file size with fseek
  FILE *fh        = fopen(fn,"rb");
  if (fh == NULL) {
    if (verbose > 0)
      printf("File not found!\n");
    return NULL;
  } else {
    fseek(fh, 0L, SEEK_END);
    fl  = ftell(fh);
    fclose(fh);
  }

  // now read data from file
  fh  = fopen(fn,"r");
  // get columns
  if (fread(&m, sizeof(uint32_t), 1, fh) != 1) {
    if (verbose > 0)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }
  // get rows
  if (fread(&n, sizeof(uint32_t), 1, fh) != 1) {
    if (verbose > 0)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }
  // get modulus
  if ((fread(&mod, sizeof(uint32_t), 1, fh) != 1) || (mod == 1)) {
    if (verbose > 0)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }
  // get number of nonzero elements
  if (fread(&nnz, sizeof(uint64_t), 1, fh) != 1) {
    if (verbose > 0)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }

  if (verbose > 1) {
    // density of matrix
    density =   (double) n * (double) m;
    density =   (double) (nnz) / density;
    density *=  100.0;
    // file size of matrix
    fs  = (double) fl / 1024 / 1024;
    fsu = "MB";
    if (fs > 1000) {
      fs  = fs / 1024;
      fsu = "GB";
    }
    printf("--------------------------------------------------------------\n");
    printf("Data for %s\n", fn);
    printf("--------------------------------------------------------------\n");
    printf("modulus                     %14d\n", mod);
    printf("number of rows              %14d\n", m);
    printf("number of columns           %14d\n", n);
    printf("number of nonzero elements  %14ld\n", nnz);
    printf("density                     %14.2f %\n", density);
    printf("size                        %14.2f %s\n", fs, fsu);
    printf("--------------------------------------------------------------\n");
  }

  // read entries from file
  sm_t *M   = (sm_t *)malloc(sizeof(sm_t));
  M->rows   = (re_t **)malloc(m*sizeof(re_t *));
  M->pos    = (ci_t **)malloc(m*sizeof(ci_t *));
  M->rwidth = (ci_t *)malloc(m*sizeof(ci_t));
  
  // get meta data
  M->nrows   = m;
  M->ncols   = n;
  M->nnz     = nnz;
  M->mod     = mod;
  M->density = (float)density;
  
  // maximal possible nonzero entries per row is n*sizeof(entry_t)
  re_t *nze = (re_t *)malloc(n * sizeof(re_t));
  ci_t *pos = (ci_t *)malloc(n * sizeof(ci_t));

  // store header size of file
  // size of m, n, mod and nb
  uint32_t hs  = 3 * sizeof(uint32_t) + sizeof(uint64_t);

  // offsets for file handling
  uint64_t row_size_offset  = nnz * sizeof(re_t) + nnz * sizeof(uint32_t) + hs;
  uint64_t row_val_offset   = hs;
  uint64_t row_pos_offset   = nnz * sizeof(re_t) + hs;

  ri_t i;
  ci_t j;
  ci_t sz;

  for (i = 0; i < m; ++i) {
    fseek(fh, row_size_offset, SEEK_SET);
    if (fread(&sz, sizeof(ci_t), 1, fh) != 1) {
      if (verbose > 0)
        printf("Error while reading file '%s'\n",fn);
      free(M);
      fclose(fh);
      return NULL;
    }

    row_size_offset +=  sizeof(ci_t);

    fseek(fh, row_val_offset, SEEK_SET);
    if (fread(nze, sizeof(re_t), sz, fh) != sz) {
      if (verbose > 0)
        printf("Error while reading file '%s'\n",fn);
      free(M);
      fclose(fh);
      return NULL;
    }

    row_val_offset +=  sz * sizeof(re_t);

    fseek(fh, row_pos_offset, SEEK_SET);
    if (fread(pos, sizeof(ci_t), sz, fh) != sz) {
      if (verbose > 0)
        printf("Error while reading file '%s'\n",fn);
      free(M);
      fclose(fh);
      return NULL;
    }

    row_pos_offset +=  sz * sizeof(ci_t);

    // reserve memory in matrix M for rows[i]
    M->rows[i] = (re_t *)malloc(sz * sizeof(re_t));
    M->pos[i]     = (ci_t *)malloc(sz * sizeof(ci_t));
    for (j = 0; j < sz; ++j) {
      M->rows[i][j] = nze[j];
      M->pos[i][j]  = pos[j];
    }
    M->rwidth[i]  = sz;
  }

  fclose(fh); 
  // print walltime
  if (verbose > 1) {
    printf("||| Time for loading JCF matrix: %7.3f sec (%7.3f %s/sec)\n",
        walltime(t_load_start) / (1000000),
        fs / (walltime(t_load_start) / (1000000)), fsu);
  }
  return M;
}


// ========== WRITING ==========

void write_jcf_matrix_to_file(sm_t *M, const char *fn, int verbose) {
}

void write_jcf_matrix_to_pbm(sm_t *M, const char *fn, int verbose) {
  char buffer[512];
  unsigned char out_byte  = 0;

  ri_t m = M->nrows;
  ci_t n = M->ncols;

  FILE *fh  = fopen(fn, "wb");

  // magic PBM header
#ifdef __LP64__ // 64bit machine
  sprintf(buffer, "P4\n# matrix size(%lu, %lu)\n%lu %lu\n", m, n, n, m);
#else // 32bit machine
  sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#endif

  fwrite(buffer, sizeof(char), strlen(buffer), fh);

  ri_t i;
  ci_t j, k;
  // row width: number of nonzero elements in current row
  ci_t sz;

  for (i = 0; i < m; ++i) {
    k   = 0;
    sz  = M->rwidth[i];
    for (j = 0; j < n; ++j) {
      if (k < sz && M->pos[i][k] == j) {
        out_byte  |=  (1 << (7 - (j % 8)));
        k++;
      } else {
        out_byte  &=  ~(1 << (7 - (j % 8)));
      }
      if (j % 8 == 7) {
        fwrite(&out_byte, sizeof(unsigned char), 1, fh);
        out_byte  = 0;
      }
    }
    if (j % 8 != 0)
      fwrite(&out_byte, sizeof(unsigned char), 1, fh);

    fflush(fh);
  }
  fclose(fh);
  
}
