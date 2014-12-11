#include "matrix.h"

#define NOT_DENSE_COPYING 0

void copy_block_matrix_to_sparse_matrix(sbm_fl_t **input,
    sm_t **output, int deleteIn, int nthrds) {

  sbm_fl_t *in  = *input;
  sm_t *out     = *output;

  ri_t i, ii;
  ci_t j;
  bi_t k;
  bi_t l;

  // initialize meta data for multiline matrix out
  out->nrows    = in->nrows;  // row dimension
  out->ncols    = in->ncols;  // col dimension

  // allocate memory for blocks
  ri_t rl = out->nrows;

  out->rows   = (re_t **)malloc(rl * sizeof(re_t *));
  out->pos    = (ci_t **)malloc(rl * sizeof(ci_t *));
  out->rwidth = (ci_t *)malloc(rl * sizeof(ci_t));
  memset(out->rwidth, 0, rl * sizeof(ci_t));
  #pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i) nowait
    for (i=0; i<rl; ++i) {
      out->rows[i]  = (re_t *)malloc(out->ncols * sizeof(re_t));
      memset(out->rows[i], 0, out->ncols * sizeof(re_t));
      out->pos[i]   = (ci_t *)malloc(out->ncols * sizeof(ci_t));
      memset(out->pos[i], 0, out->ncols * sizeof(ci_t));
    }
  }

  const ri_t rlin = (ri_t) ceil((float) in->nrows / in->bheight);
  const ci_t clin = (ci_t) ceil((float) in->ncols / in->bwidth);
  // we need buffers for all multiline entries since the copying process
  // revisits already filled up multilines.if we reset the buffer to zero, we
  // might realloc only init_buffer memory and lose what we have already in the
  // multiline stored.

  ci_t block_idx;
  ri_t block_vert;
  ri_t sparse_idx;
  re_t v1, v2;
  mli_t crb = 0;
  const bi_t ml_bheight = in->bheight / __GB_NROWS_MULTILINE;
  #pragma omp parallel num_threads(nthrds)
  {
    #pragma omp for private(sparse_idx,i,ii,j,k,l,block_idx,block_vert,v1,v2) nowait
    for (i=0; i<rlin; ++i) { // loop over multilines
      block_vert = i * in->bheight;
      for (j=0; j<clin; ++j) {
        block_idx = j * in->bwidth;
        if (in->blocks[i][j] == NULL) {
          continue;
        }
        for (k=0;k<ml_bheight; ++k) {
          if (in->blocks[i][j][k].sz == 0) {
            continue;
          }
          sparse_idx  = block_vert+2*k;
          if (in->blocks[i][j][k].dense == 0) {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              // fill in data
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0) {
                out->rows[sparse_idx][out->rwidth[sparse_idx]]  = v1;
                out->pos[sparse_idx][out->rwidth[sparse_idx]]   = block_idx + in->blocks[i][j][k].idx[l];
                out->rwidth[sparse_idx]++;
              }
              if (v2 != 0) {
                out->rows[sparse_idx+1][out->rwidth[sparse_idx+1]]  = v2;
                out->pos[sparse_idx+1][out->rwidth[sparse_idx+1]]   = block_idx + in->blocks[i][j][k].idx[l];
                out->rwidth[sparse_idx+1]++;
              }
            }
          } else {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              // fill in data
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0) {
                out->rows[sparse_idx][out->rwidth[sparse_idx]]  = v1;
                out->pos[sparse_idx][out->rwidth[sparse_idx]]   = block_idx + l;
                out->rwidth[sparse_idx]++;
              }
              if (v2 != 0) {
                out->rows[sparse_idx+1][out->rwidth[sparse_idx+1]]  = v2;
                out->pos[sparse_idx+1][out->rwidth[sparse_idx+1]]   = block_idx + l;
                out->rwidth[sparse_idx+1]++;
              }
            }
          }
        }
      }
    }
  }
  // realloc memory
  #pragma omp parallel num_threads(nthrds)
    {
#pragma omp for private(i) nowait
      for (i=0; i<out->nrows; ++i) {
        out->rows[i]  = realloc(out->rows[i],out->rwidth[i]*sizeof(re_t));
        out->pos[i]   = realloc(out->pos[i],out->rwidth[i]*sizeof(ci_t));
      }
    }

  if (deleteIn) {
  // free memory for input matrix
  #pragma omp parallel num_threads(nthrds)
    {
#pragma omp for private(i,j,k) nowait
      for (i=0; i<rlin; ++i) {
        for (j=0; j<clin; ++j) {
          if (in->blocks[i][j] != NULL) {
            for (k=0; k<ml_bheight; ++k) {
              free(in->blocks[i][j][k].idx);
              in->blocks[i][j][k].idx = NULL;

              free(in->blocks[i][j][k].val);
              in->blocks[i][j][k].val = NULL;
            }
            free(in->blocks[i][j]);
            in->blocks[i][j] = NULL;
          }
        }
        free(in->blocks[i]);
        in->blocks[i] = NULL;
      }
    }
    free(in->blocks);
    in->blocks  = NULL;
    free(in);
    in  = NULL;
  }
  *input  = in;
  *output = out;
}

sm_fl_ml_t *copy_block_matrix_to_multiline_matrix(sbm_fl_t **input,
    sm_fl_ml_t *out, int deleteIn, int nthrds) {

  sbm_fl_t *in  = *input;

  ri_t i, ii;
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
  out->density  = (double)0;

  // allocate memory for blocks
  ri_t rl = out->nrows / 2;
  if (out->nrows % 2)
    rl++;

  out->ml = (ml_t *)malloc(rl * sizeof(ml_t));
  for (i=0; i<rl; ++i) {
#if NOT_DENSE_COPYING
    out->ml[i].val    = NULL;
    out->ml[i].idx    = NULL;
    out->ml[i].sz     = 0;
    out->ml[i].dense  = 0;
#else
    out->ml[i].val    = (re_t *)malloc(2 * out->ncols * sizeof(re_t));
    out->ml[i].idx    = NULL;
    out->ml[i].sz     = out->ncols;
    out->ml[i].dense  = 1;
    //out->ml[i].val  = realloc(out->ml[i].val,2 * out->ncols * sizeof(re_t));
#endif
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

  ci_t block_idx;
  re_t v1, v2;
  mli_t crb = 0;
  const bi_t ml_bheight = in->bheight / __GB_NROWS_MULTILINE;
  #pragma omp parallel shared(buffer) num_threads(nthrds)
  {
    #pragma omp for private(i,ii,j,k,l,block_idx,v1,v2) nowait
    for (i=0; i<rl; ++i) { // loop over multilines
      memset(out->ml[i].val, 0, 2 * out->ncols * sizeof(re_t));
      ii  = i / ml_bheight;
      k   = i % ml_bheight;
      for (j=0; j<clin; ++j) {
        block_idx = j * in->bwidth;
        if (in->blocks[ii][j] == NULL) {
          continue;
        }
        if (in->blocks[ii][j][k].sz == 0) {
          continue;
        }
        if (in->blocks[ii][j][k].dense == 0) {
          for (l=0; l<in->blocks[ii][j][k].sz; ++l) {
            // fill in data
            v1  = in->blocks[ii][j][k].val[2*l];
            v2  = in->blocks[ii][j][k].val[2*l+1];
            if (v1 != 0 || v2 != 0) {
              out->ml[i].val[2*(block_idx + in->blocks[ii][j][k].idx[l])]   = v1;
              out->ml[i].val[2*(block_idx + in->blocks[ii][j][k].idx[l])+1] = v2;
            }
            buffer[i]++;
          }
        } else {
          for (l=0; l<in->blocks[ii][j][k].sz; ++l) {
            // fill in data
            v1  = in->blocks[ii][j][k].val[2*l];
            v2  = in->blocks[ii][j][k].val[2*l+1];
            if (v1 != 0 || v2 != 0) {
              out->ml[i].val[2*(block_idx+l)]   = v1;
              out->ml[i].val[2*(block_idx+l)+1] = v2;
            }
            buffer[i]++;
          }
        }
      }
      if (buffer[i] == 0) {
        free (out->ml[i].val);
        out->ml[i].val  = NULL;
        out->ml[i].sz   = 0;
      }
    }
  }
  if (deleteIn) {
  // free memory for input matrix
  #pragma omp parallel num_threads(nthrds)
    {
#pragma omp for private(i,j,k) nowait
      for (i=0; i<rlin; ++i) {
        for (j=0; j<clin; ++j) {
          if (in->blocks[i][j] != NULL) {
            for (k=0; k<ml_bheight; ++k) {
              free(in->blocks[i][j][k].idx);
              in->blocks[i][j][k].idx = NULL;

              free(in->blocks[i][j][k].val);
              in->blocks[i][j][k].val = NULL;
            }
            free(in->blocks[i][j]);
            in->blocks[i][j] = NULL;
          }
        }
        free(in->blocks[i]);
        in->blocks[i] = NULL;
      }
    }
    free(in->blocks);
    in->blocks  = NULL;
    free(in);
    in  = NULL;
  }
/*
 *  NOTE:
 *  old and deprecated code: keep for possible tests with sparse-hybrid-dense
 *  representations for D in multiline structures later on
 *
  #pragma omp parallel shared(buffer) num_threads(nthrds)
  {
    mli_t crb = 0;
    #pragma omp for private(i,j,k,l,block_idx,crb) nowait
    for (i=0; i<rlin; ++i) {
      // curr_row_base in LELA
      crb  = i * in->bheight / __GB_NROWS_MULTILINE;

      for (j=0; j<clin; ++j) {
        if (in->blocks[i][j] == NULL) {
          continue;
        }

        block_idx  = in->bwidth * j;
        for (k=0; k<in->bheight/__GB_NROWS_MULTILINE; ++k) {
#if NOT_DENSE_COPYING
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
#else
          if (in->blocks[i][j][k].sz == 0) {
            continue;
          }
          if (in->blocks[i][j][k].dense == 0) {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              printf("%d - %d - %d - %d - %d - %d\n",i,j,k,l,crb,block_idx);
              // fill in data
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].val[2*(block_idx + in->blocks[i][j][k].idx[l])]   = v1;
                out->ml[crb+k].val[2*(block_idx + in->blocks[i][j][k].idx[l])+1] = v2;
              }
              buffer[crb+k]++;
            }
          } else {
            for (l=0; l<in->blocks[i][j][k].sz; ++l) {
              printf("%d - %d - %d - %d - %d - %d\n",i,j,k,l,crb,block_idx);
              // fill in data
              v1  = in->blocks[i][j][k].val[2*l];
              v2  = in->blocks[i][j][k].val[2*l+1];
              if (v1 != 0 || v2 != 0) {
                out->ml[crb+k].val[2*(block_idx+l)]   = v1;
                out->ml[crb+k].val[2*(block_idx+l)+1] = v2;
              }
              buffer[crb+k]++;
            }
          }
#endif
          // destruct input matrix?
          if (deleteIn) {
            free(in->blocks[i][j][k].idx);
            in->blocks[i][j][k].idx = NULL;

            free(in->blocks[i][j][k].val);
            in->blocks[i][j][k].val = NULL;
          }
        }
        if (deleteIn) {
          free(in->blocks[i][j]);
          in->blocks[i][j] = NULL;
        }
      }
      if (deleteIn) {
        free(in->blocks[i]);
        in->blocks[i] = NULL;
      }
    }
    for (i=0; i<rlin; ++i) {
#if NOT_DENSE_COPYING
      // realloc memory, only needed if the multilines are not dense copied
      if (out->ml[i].sz > 0) {
        if (out->ml[i].sz < out->ncols) {
          out->ml[i].idx  = realloc(out->ml[i].idx,
              out->ml[i].sz * sizeof(mli_t));
        } else {
          free (out->ml[i].idx);
          out->ml[i].idx = NULL;
        }
        out->ml[i].val  = realloc(out->ml[i].val,
            2 * out->ml[i].sz * sizeof(re_t));
      }
#else
      if (buffer[i] == 0) {
        free (out->ml[i].val);
        out->ml[i].val  = NULL;
        out->ml[i].sz   = 0;
      }
#endif
    }
  }
  if (deleteIn) {
    free(in->blocks);
    in->blocks  = NULL;
    printf("addr in %p\n",in);
    free(in);
    in  = NULL;
    printf("addr in %p\n",in);
  }
  */
  free(buffer);
  *input  = in;
  return out;
}

sbm_fl_t *copy_multiline_to_block_matrix_rl(sm_fl_ml_t **A_in,
    ri_t bheight, ci_t bwidth, int free_memory, int nthrds) {

  sm_fl_ml_t *A = *A_in;
  sbm_fl_t *B = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));

  ci_t i;
  ri_t j;
  ri_t k;

  const ri_t rlB            = (ri_t) ceil((float) A->nrows / bheight);
  const ci_t clB            = (ci_t) ceil((float) A->ncols / bwidth);
  const ri_t ml_rlB         = (A->nrows % __GB_NROWS_MULTILINE == 0) ?
    A->nrows / __GB_NROWS_MULTILINE :
    A->nrows / __GB_NROWS_MULTILINE + 1;
  const bi_t ml_bheight     = bheight / __GB_NROWS_MULTILINE;

  // initialize meta data for block submatrices
  B->nrows    = A->nrows; // row dimension
  B->ncols    = A->ncols; // col dimension
  B->bheight  = bheight;  // block height
  B->bwidth   = bwidth;   // block width
  B->ba       = dtlr;     // block alignment
  B->fe       = 1;        // fill empty blocks?
  B->hr       = 1;        // allow hybrid rows?
  B->nnz      = 0;        // number nonzero elements
  // initialize B
  B->blocks = (mbl_t ***)malloc(rlB * sizeof(mbl_t **));
  for (i=0; i<rlB; ++i) {
    B->blocks[i]  = (mbl_t **)malloc(clB * sizeof(mbl_t *));
    for (j=0; j<clB; ++j) {
      B->blocks[i][j] = (mbl_t *)malloc(
          ml_bheight * sizeof(mbl_t));
      for (k=0; k<ml_bheight; ++k) {
        B->blocks[i][j][k].val  = (re_t *)malloc(2*bwidth*sizeof(re_t));
        B->blocks[i][j][k].idx  = (bi_t *)malloc(bwidth*sizeof(bi_t));
        B->blocks[i][j][k].sz   = B->blocks[i][j][k].dense  = 0;
      }
    }
  }
  const bi_t init_buffer_B  = (bi_t)(B->bwidth/2);
#pragma omp parallel private(j,k) num_threads(nthrds)
  {
    ri_t rbi;           // row block index
    mli_t idx;          // multiline index of values
    bi_t eil;           // element index in line
    bi_t bir;           // block index in row
    bi_t lib;           // line index in block
    bi_t buffer_B[clB]; // buffer status for blocks in B
    memset(buffer_B, 0, clB * sizeof(bi_t));

#pragma omp for nowait ordered
    for (i=0; i<ml_rlB; ++i) {
      if (A->ml[i].sz == 0)
        continue;

      rbi = i / ml_bheight;
      lib = i % ml_bheight;

      for (j=A->ml[i].sz; j>0; --j) {
        idx = A->ml[i].idx[j-1];
        bir = (A->ncols - 1 - idx) / bwidth;
        eil = (A->ncols - 1 - idx) % bwidth;
        // realloc memory if needed
        /*
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir]) {
          printf("init_buffer_B %d || buffer_B[%d] %d\n",init_buffer_B,bir,buffer_B[bir]);
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        }
        */
        // set values
        B->blocks[rbi][bir][lib].idx[B->blocks[rbi][bir][lib].sz]   = eil;
        B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz]   =
          A->ml[i].val[2*(j-1)];
        B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz+1] =
          A->ml[i].val[2*(j-1)+1];
        B->blocks[rbi][bir][lib].sz++;
      }
      if (free_memory == 1) {
        free(A->ml[i].idx);
        //A->ml[i].idx  = NULL;
        free(A->ml[i].val);
        //A->ml[i].val  = NULL;
      }
    }
  }
  if (free_memory == 1) {
    free(A->ml);
    A->ml = NULL;
    free(A);
    A = NULL;
  }
  A_in  = &A;
  // hybrid multirows for the righthand side block matrices?
  if (B->hr) {
#pragma omp parallel num_threads(nthrds)
    {
      ri_t idx;
      bi_t l;
#pragma omp for private(i,j,k,l) schedule(dynamic) nowait ordered
      // TODO: Implement hybrid stuff
      for (i=0; i<rlB; ++i) {
        for (j=0; j<clB; ++j) {
          for (k=0; k<ml_bheight; ++k) {
            if (B->blocks[i][j][k].sz == 0) {
              free(B->blocks[i][j][k].idx);
              B->blocks[i][j][k].idx  = NULL;
              free(B->blocks[i][j][k].val);
              B->blocks[i][j][k].val  = NULL;
            }
            if ((float)B->blocks[i][j][k].sz / (float)B->bwidth
                < __GB_HYBRID_THRESHOLD) {
              B->blocks[i][j][k].idx =  realloc(
                  B->blocks[i][j][k].idx,
                  B->blocks[i][j][k].sz * sizeof(bi_t));
              B->blocks[i][j][k].val =  realloc(
                  B->blocks[i][j][k].val,
                  2 * B->blocks[i][j][k].sz * sizeof(bi_t));
              continue;
            }
            re_t *tmp_val_ptr = (re_t *)malloc(2 * B->bwidth * sizeof(re_t));
            idx  = 0;
            for (l=0; l<B->bwidth; ++l) {
              if (idx < B->blocks[i][j][k].sz && B->blocks[i][j][k].idx[idx] == l) {
                tmp_val_ptr[2*l]    = B->blocks[i][j][k].val[2*idx];
                tmp_val_ptr[2*l+1]  = B->blocks[i][j][k].val[2*idx+1];
                idx++;
              } else {
                tmp_val_ptr[2*l]    = 0;
                tmp_val_ptr[2*l+1]  = 0;
              }
            }
            free(B->blocks[i][j][k].idx);
            B->blocks[i][j][k].idx    = NULL;
            free(B->blocks[i][j][k].val);
            B->blocks[i][j][k].val    = tmp_val_ptr;
            B->blocks[i][j][k].sz     = B->bwidth;
            B->blocks[i][j][k].dense  = 1;
          }
        }
      }
    }
  } else { // cut down memory usage
    // Realloc memory usage:
    // Note that A is reallocated during the swapping of the data, so we
    // reallocate only B here.
#pragma omp parallel num_threads(nthrds)
    {
#pragma omp for schedule(dynamic) nowait ordered
      for (i=0; i<rlB; ++i) {
        int ctr;
        for (j=0; j<clB; ++j) {
          ctr = 0;
          for (k=0; k<ml_bheight; ++k) {
            if (B->blocks[i][j][k].sz > 0) {
              ctr = 1;
              B->blocks[i][j][k].idx = realloc(
                  B->blocks[i][j][k].idx,
                  B->blocks[i][j][k].sz * sizeof(bi_t));
              B->blocks[i][j][k].val = realloc(
                  B->blocks[i][j][k].val,
                  2 * B->blocks[i][j][k].sz  * sizeof(re_t));
            }
          }
          // if full block is empty, remove it!
          if (ctr == 0) {
            free(B->blocks[i][j]);
            B->blocks[i][j] = NULL;
          }
        }
      }
    }
  }
  return B;
}
