#include "matrix.h"

void copy_block_matrix_to_multiline_matrix(sbm_fl_t *in, sm_fl_ml_t *out,
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
  mli_t init_buffer = 2 * in->bwidth;
  // we need buffers for all multiline entries since the copying process
  // revisits already filled up multilines.if we reset the buffer to zero, we
  // might realloc only init_buffer memory and lose what we have already in the
  // multiline stored.
  mli_t *buffer = (mli_t *)malloc(rl * sizeof(mli_t));
  memset(buffer, 0, rl * sizeof(mli_t));

  mli_t crb;
  ci_t block_idx;
  re_t v1, v2;
  #pragma omp parallel num_threads(nthrds)
  {
    #pragma omp for private(j,k,l,crb) nowait
    for (i=0; i<rlin; ++i) {
      // curr_row_base in LELA
      ci_t crb  = i * in->bheight / __GB_NROWS_MULTILINE;

      for (j=0; j<clin; ++j) {
        if (in->blocks[i][j] == NULL) {
          continue;
        }

        block_idx  = in->bwidth * j;
        for (k=0; k<in->bheight/__GB_NROWS_MULTILINE; ++k) {
          if (in->blocks[i][j][k].sz == 0) {
            continue;
          }
          if (in->blocks[i][j][k].dense == 0) {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              // possibly realloc more memory
              if (out->ml[crb+k].sz == buffer[crb+k]) {
                buffer[crb+k] +=  init_buffer;
                out->ml[crb+k].idx = realloc(out->ml[crb+k].idx,
                    buffer[crb+k] * sizeof(mli_t));
                out->ml[crb+k].val = realloc(out->ml[crb+k].val,
                    2 * buffer[crb+k] * sizeof(re_t));
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
              if (out->ml[crb+k].sz == buffer[crb+k]) {
                buffer[crb+k]  +=  init_buffer;
                out->ml[crb+k].idx = realloc(out->ml[crb+k].idx,
                    buffer[crb+k] * sizeof(mli_t));
                out->ml[crb+k].val = realloc(out->ml[crb+k].val,
                    2 * buffer[crb+k] * sizeof(re_t));
              }
              // fill in data
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].idx[out->ml[crb+k].sz] =
                  block_idx + l;
                out->ml[crb+k].val[2*out->ml[crb+k].sz]   = v1;
                out->ml[crb+k].val[2*out->ml[crb+k].sz+1] = v2;
                out->ml[crb+k].sz++;
              }
            }
          }

          // destruct input matrix?
          if (deleteIn) {
            if (in->blocks[i][j][k].idx != NULL) {
              free(in->blocks[i][j][k].idx);
              in->blocks[i][j][k].idx = NULL;
            }
            if (in->blocks[i][j][k].val != NULL) {
              free(in->blocks[i][j][k].val);
              in->blocks[i][j][k].val = NULL;
            }
          }
        }
        if (deleteIn) {
          if (in->blocks[i][j] != NULL) {
            free(in->blocks[i][j]);
            in->blocks[i][j] = NULL;
          }
        }
      }
      if (deleteIn) {
        if (in->blocks[i] != NULL) {
          free(in->blocks[i]);
          in->blocks[i] = NULL;
        }
      }
    }
    // realloc memory
    if (out->ml[crb+k].sz > 0) {
      out->ml[crb+k].idx  = realloc(out->ml[crb+k].idx,
          out->ml[crb+k].sz * sizeof(mli_t));
      out->ml[crb+k].val  = realloc(out->ml[crb+k].val,
          2 * out->ml[crb+k].sz * sizeof(re_t));
    }
  }
  free(in);
  in  = NULL;
}
