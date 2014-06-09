#include "io.h"

// reading & writing
sparse_mat_t *load_matrix_jcf_format(const char *fn, int verbose) {

  // start loading the matrix
  row_dim_t m;
  col_dim_t n;
  mod_t     mod;
  nnz_t     nnz;
  uint64    fl;

  // open in binary mode first to get file size with fseek
  FILE *fh        = fopen(fn,"rb");
  if (fh == NULL) {
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
  if (fread(&m, sizeof(uint32), 1, fh) != 1) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }
  // get rows
  if (fread(&n, sizeof(uint32), 1, fh) != 1) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }
  // get modulus
  if ((fread(&mod, sizeof(uint32), 1, fh) != 1) || (mod == 1)) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }
  // get number of nonzero elements
  if (fread(&nnz, sizeof(uint64), 1, fh) != 1) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    fclose(fh);
    return NULL;
  }

  if (verbose) {
    // density of matrix
    double density  = (double) n * (double) m;
    density   =   (double) (nnz) / density;
    density   *=  100.0;
    // file size of matrix
    double fs =   (double) fl / 1024 / 1024;
    char *fsu =   "MB";
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
    printf("density                     %14.2f%\n", density);
    printf("size                        %14.2f%s\n", fs, fsu);
    printf("--------------------------------------------------------------\n");
  }

  // read entries from file
  sparse_mat_t *A = (sparse_mat_t *)malloc(sizeof(sparse_mat_t));
  A->rows         = (entry_t **)malloc(sizeof(entry_t **));
  A->pos          = (col_idx_t **)malloc(sizeof(col_idx_t **));
  // get meta data
  A->nrows   = m;
  A->ncols   = n;
  A->nnz     = nnz;
  A->mod     = mod;
  
  // maximal possible nonzero entries per row is n*sizeof(entry_t)
  entry_t   *nze  = (entry_t *)malloc(n * sizeof(entry_t));
  col_idx_t *pos  = (col_idx_t *)malloc(n * sizeof(col_idx_t));

  // store header size of file
  // size of m, n, mod and nb
  uint32 hs  = 3 * sizeof(uint32) + sizeof(uint64);

  // offsets for file handling
  uint64 row_size_offset = nnz * sizeof(entry_t) + nnz * sizeof(uint32) + hs;
  uint64 row_val_offset = hs;
  uint64 row_pos_offset = nnz * sizeof(entry_t) + hs;

  uint32 sz;
  uint32 i, j;

  for (i = 0; i < m; ++i) {
    fseek(fh, row_size_offset, SEEK_SET);
    if (fread(&sz, sizeof(col_idx_t), 1, fh) != 1) {
      if (verbose)
        printf("Error while reading file '%s'\n",fn);
      free(A);
      fclose(fh);
      return NULL;
    }

    row_size_offset +=  sizeof(uint32);

    fseek(fh, row_val_offset, SEEK_SET);
    if (fread(nze, sizeof(entry_t), sz, fh) != sz) {
      if (verbose)
        printf("Error while reading file '%s'\n",fn);
      free(A);
      fclose(fh);
      return NULL;
    }

    row_val_offset +=  sz * sizeof(entry_t);

    fseek(fh, row_pos_offset, SEEK_SET);
    if (fread(pos, sizeof(col_idx_t), sz, fh) != sz) {
      if (verbose)
        printf("Error while reading file '%s'\n",fn);
      free(A);
      fclose(fh);
      return NULL;
    }

    row_pos_offset +=  sz * sizeof(uint32);

    // reserve memory in matrix A for row[i]
    A->rows[i] = (entry_t *)malloc(sz * sizeof(entry_t));
    A->pos[i]  = (col_idx_t *)malloc(sz * sizeof(col_idx_t));
    for (j = 0; j < sz; ++j) {
      A->rows[i][j]  = nze[j];
      A->pos[i][j]   = pos[j];
    }
  }

  fclose(fh); 
  return A;
}

void write_matrix_jcf_format(sparse_mat_t *mat, FILE *file) {
}
