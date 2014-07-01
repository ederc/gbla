#include <mapping.h>

map_fl_t *construct_fl_map(sm_t *M) {
  // initialize all map entries to __GB_MINUS_ONE_8
  map_fl_t *map = init_fl_map(M);

  uint32_t npiv = 0;  // number of pivots
  ri_t i = 0;         // current row index
  re_t entry;         // possible pivot entry to be checked

  // sweeps the rows to identify the row pivots and column pivots
  for (i=0; i<M->nrows; ++i) {
    if (M->rwidth[i] != 0) {
      entry = M->rows[i][M->pos[i][0]];
      if (map->pri[entry] == __GB_MINUS_ONE_32) {
        map->pri[entry] = i;
        npiv++;
      } else { // check for a sparser pivot row (see ELAGB talk from Lachartre)
        if (M->rwidth[map->pri[entry]] > M->rwidth[i]) {
          map->npri[map->pri[entry]]  = map->pri[entry];
          map->pri[entry] = i;
        } else {
          map->npri[i]  = i;
        }
      }
    } else {
      map->npri[i]  = i;
    }
  }
  map->npiv = npiv;

  ci_t pc_idx = 0, npc_idx = 0, j;

  // construct pivot columns and non-pivot columns maps and the corresponding
  // reverse maps
  for (j=0; j<M->ncols; ++j) {
    if (map->pri[j] !=  __GB_MINUS_ONE_32) {
      map->pc[j]           = pc_idx;
      map->pc_rev[pc_idx]  = j;
      pc_idx++;
    } else {
      map->npc[j]           = npc_idx;
      map->npc_rev[npc_idx] = j;
      npc_idx++;
    }
  }
  return map;
}

void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
                      map_fl_t *map, int block_dim, int rows_multiline,
                      int nthreads) {
  
  A = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  B = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  C = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  D = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));

  int bdim  = block_dim / rows_multiline;
  // construct index map for Faugère-Lachartre decomposition of matrix M
  map  = construct_fl_map(M);

  // initialize meta data for block submatrices
  A->nrows    = map->npiv;            // row dimension
  A->ncols    = map->npiv;            // col dimension
  A->bheight  = bdim;                 // block height
  A->bwidth   = bdim;                 // block width
  A->ba       = dtrl;                 // block alignment
  A->fe       = 0;                    // fill empty blocks?
  A->hr       = 0;                    // allow hybrid rows?

  B->nrows    = map->npiv;            // row dimension
  B->ncols    = M->nrows - map->npiv; // col dimension
  B->bheight  = bdim;                 // block height
  B->bwidth   = bdim;                 // block width
  B->ba       = dtlr;                 // block alignment
  B->fe       = 1;                    // fill empty blocks?
  B->hr       = 1;                    // allow hybrid rows?

  C->nrows    = M->nrows - map->npiv; // row dimension
  C->ncols    = map->npiv;            // col dimension
  C->bheight  = bdim;                 // block height
  C->bwidth   = bdim;                 // block width
  C->ba       = dtrl;                 // block alignment
  C->fe       = 0;                    // fill empty blocks?
  C->hr       = 1;                    // allow hybrid rows?

  D->nrows    = M->nrows - map->npiv; // row dimension
  D->ncols    = M->nrows - map->npiv; // col dimension
  D->bheight  = bdim;                 // block height
  D->bwidth   = bdim;                 // block width
  D->ba       = dtlr;                 // block alignment
  D->fe       = 1;                    // fill empty blocks?
  D->hr       = 1;                    // allow hybrid rows?

  assert(A->nrows == A->ncols);
  assert(A->nrows == B->nrows);
  assert(A->ncols == C->ncols);
  assert(B->ncols == D->ncols);
  assert(C->nrows == D->nrows);
  assert(M->nrows == A->nrows + D->nrows);
  assert(M->ncols == B->ncols + C->ncols);

  uint32_t npiv = 0; // number pivots handled
  // Take maximal number of rows in blocks to be able to use array
  // piv_start_idx for construction of A & B as well as for the construction of
  // C & D.
  uint32_t max_nrows =  (A->nrows > C->nrows) ? A->nrows : C->nrows; 
  uint32_t piv_start_idx[(max_nrows / B->bheight) + 2];
  piv_start_idx[0]  = M->nrows;
  uint32_t block;

  int i;

  // find blocks for construction of A & B
  for (i = (int)M->nrows-1; i > -1; --i) {
    if (map->pri[i] != __GB_MINUS_ONE_32)
      npiv++;

    if (npiv % B->bheight == 0)
      piv_start_idx[npiv / B->bheight]  = i;
  }

  piv_start_idx[(map->npiv / B->bheight) + 1] = 0;

  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  // rows indices horizontal block
    uint16_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block=0; block<=(map->npiv/B->bheight); ++block) {
      // construct block submatrices A & B
      for (i = ((int)piv_start_idx[block]-1); i >= (int)piv_start_idx[block+1];
          --i) {
        if (map->pri[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        if (cvb == B->bheight || i == 0) {
          write_left_multiline_right_block(M, A, B, rihb, cvb, block);

          // TODO: Destruct input matrix on the go
          /*
          if (destruct_input_matrix) {
          // free memory
          }
          */
          cvb = 0;
        }
      }
    }
  } 

  // find blocks for construction of C & D
  npiv  = 0;
  piv_start_idx[0]  = M->nrows;
  for (i = (int)M->nrows-1; i > -1; --i) {
    if (map->npri[i] != __GB_MINUS_ONE_32)
      npiv++;

    if (npiv % B->bheight == 0)
      piv_start_idx[npiv / B->bheight]  = i;
  }

  piv_start_idx[(C->nrows / B->bheight) + 1] = 0;

  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
#pragma omp for schedule(dynamic) nowait
    for (block=0; block<=(C->nrows/B->bheight); ++block) {
      // construct block submatrices A & B
      for (i = ((int)piv_start_idx[block]-1); i >= (int)piv_start_idx[block+1];
          --i) {
        if (map->npri[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          write_left_multiline_right_block(M, C, D, rihb, cvb, block);

          // TODO: Destruct input matrix on the go
          /*
          if (destruct_input_matrix) {
          // free memory
          }
          */
          cvb = 0;
        }
      }
    }
  } 
}


void write_blocks_lr_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, ri_t *rihb,
                            ri_t crb, ri_t rbi) {
 
  uint32_t npivs  = 0;
  uint32_t piv_start_idxs[
}
