#include "io.h"

// reading & writing
sparse_mat_t *load_matrix_jcf_format(const char *fn, int verbose) {

  // start loading the matrix
  uint16 *nz;
  uint32 *pos;
  uint32 sz;
  uint32 n;
  uint32 m;
  uint32 mod;
  uint32 nb;
  uint64 fl;

  sparse_mat_t *A = NULL;

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
  if (fread(&n, sizeof(uint32), 1, fh) != 1) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    return NULL;
  }
  // get rows
  if (fread(&m, sizeof(uint32), 1, fh) != 1) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    return NULL;
  }
  // get modulus
  if ((fread(&mod, sizeof(uint32), 1, fh) != 1) || (mod == 1)) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    return NULL;
  }
  // get number of nonzero elements
  if (fread(&nb, sizeof(uint64), 1, fh) != 1) {
    if (verbose)
      printf("Error while reading file '%s'\n",fn);
    return NULL;
  }

  if (verbose) {
    // density of matrix
    double density  = (double) n * (double) m;
    density   =   (double) (nb) / density;
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
    printf("number of rows              %14d\n", n);
    printf("number of columns           %14d\n", m);
    printf("number of nonzero elements  %14ld\n", nb);
    printf("density                     %14.2f%\n", density);
    printf("size                        %14.2f%s\n", fs, fsu);
    printf("--------------------------------------------------------------\n");
  }
  return A;
}

void write_matrix_jcf_format(sparse_mat_t *mat, FILE *file) {
}
