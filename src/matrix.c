#include "matrix.h"

void copyBlockMatrixToMultilineMatrix(sbm_fl_t *in, sm_fl_ml_t *out,
    int deleteIn, int nthrds) {

  ri_t i;
  ci_t j;
  bi_t k;
  bi_t l;

  // initialize meta data for multiline matrix out
  out->nrows    = in->nrows;  // row dimension
  out->ncols    = in->ncols;  // col dimension
  out->ba       = dtrl;       // block alignment
  out->fe       = 0;          // fill empty blocks?
  out->hr       = 0;          // allow hybrid rows?
  out->nnz      = 0;          // number nonzero elements

  // allocate memory for blocks
  ri_t rl = out->nrows / 2;
  if (out->nrows % 2)
    rl++;

  out->ml = (ml_t *)malloc(rl * sizeof(ml_t));
  for (i=0; i<rl; ++i) {
    out->ml[i].val  = NULL;
    out->ml[i].idx  = NULL;
    out->ml[i].sz   = 0;
  }

  const ri_t rlin = (ri_t) ceil((float) in->nrows / in->bheight);
  const ci_t clin = (ci_t) ceil((float) in->ncols / in->bwidth);
  mli_t init_buffer = 4 * in->bwidth;
  mli_t buffer      = 0;

  mli_t crb;
  #pragma omp parallel num_threads(nthrds)
  {
    #pragma omp for private(j,k) nowait
    for (i=0; i<rlin; ++i) {
      // curr_row_base in LELA
      ci_t crb  = i * in->bheight / __GB_NROWS_MULTILINE;

      for (j=0; j<clin; ++j) {
        if (in->blocks[i][j] == NULL)
          continue;

        ci_t block_idx  = in->bwidth * j;
        re_t v1, v2;
        for (k=0; in->bheight; ++k) {
          if (in->blocks[i][j][k].sz == 0)
            continue;
          if (in->blocks[i][j][k].dense == 0) {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              // possibly realloc more memory
              if (out->ml[crb+k].sz == buffer) {
                buffer  +=  init_buffer;
                out->ml[crb+k].idx = realloc(out->ml[crb+k].idx,
                    buffer * sizeof(mli_t));
                out->ml[crb+k].val = realloc(out->ml[crb+k].val,
                    2 * buffer * sizeof(re_t));
              }
              // fill in data
              out->ml[crb+k].idx[out->ml[crb+k].sz] =
                block_idx + in->blocks[i][j][k].idx[l];
              out->ml[crb+k].val[2*out->ml[crb+k].sz] =
                in->blocks[i][j][k].val[2*l];
              out->ml[crb+k].val[2*out->ml[crb+k].sz+1] =
                in->blocks[i][j][k].val[2*l+1];
              out->ml[crb+k].sz++;
            }
          } else { // block submatrix in is dense
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              // possibly realloc more memory
              if (out->ml[crb+k].sz == buffer) {
                buffer  +=  init_buffer;
                out->ml[crb+k].idx = realloc(out->ml[crb+k].idx,
                    buffer * sizeof(mli_t));
                out->ml[crb+k].val = realloc(out->ml[crb+k].val,
                    2 * buffer * sizeof(re_t));
              }
              // fill in data
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].idx[out->ml[crb+k].sz] =
                  block_idx + l;
                out->ml[crb+k].val[2*out->ml[crb+k].sz] =
                  in->blocks[i][j][k].val[2*l];
                out->ml[crb+k].val[2*out->ml[crb+k].sz+1] =
                  in->blocks[i][j][k].val[2*l+1];
                out->ml[crb+k].sz++;
              }
            }
          }
          // realloc memory
          out->ml[crb+k].idx  = realloc(out->ml[crb+k].idx,
              out->ml[crb+k].sz * sizeof(bi_t));
          out->ml[crb+k].val  = realloc(out->ml[crb+k].val,
              2 * out->ml[crb+k].sz * sizeof(re_t));

          // destruct input matrix?
          if (deleteIn) {
            free(in->blocks[i][j][k].idx);
            free(in->blocks[i][j][k].val);
          }
        }
        if (deleteIn) {
          free(in->blocks[i][j]);
        }
      }
      if (deleteIn) {
        free(in->blocks[i]);
      }
    }
  }
  free(in);
  in  = NULL;
}
