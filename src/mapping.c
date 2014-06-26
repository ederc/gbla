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
                      map_fl_t *map, int block_dim, int rows_multiline) {
  
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

  ri_t rbi  = 0; // row block index
  ri_t crb  = 0; // current row in block
  ri_t rihb[A->bheight]; // row indices in horizontal block

  /*
   * TODO:  1.  The following should be done in parallel, use M read
   *            only and construct A,B and C,D in parallel.
   *        2.  Destroy M on the fly or at the end of the construction of
   *            A,B,C,D.
   */

  ri_t i;
  // construct matrices A and B
  for (i=M->nrows-1; i>-1; --i) {
    if (map->pri[i] == __GB_MINUS_ONE_32) // no pivot row
      continue;

    rihb[crb] = map->pri[i];
    crb++;

    if (crb == A->bheight)
      write_blocks_lr_matrix(M, A, B, rihb, crb, rbi);
    // TODO: destruct original matrix on the fly here

    rbi++;
    crb = 0;
  }

  if (crb != 0)
    write_blocks_lr_matrix(M, A, B, rihb, crb, rbi);
  // TODO: destruct original matrix on the fly here
  

  // construct matrices C and D
  rbi = 0;
  crb = 0;

  for (i=M->ncols-1; i>-1; --i) {
    if (map->npri[i] == __GB_MINUS_ONE_32)
      continue;

    rihb[crb] = map->npri[i];
    crb++;

    if (crb == C->bheight)
      write_blocks_lr_matrix(M, C, D, rihb, crb, rbi);
    // TODO: destruct original matrix on the fly here

    rbi++;
    crb = 0;
  }

  if (crb != 0)
    write_blocks_lr_matrix(M, C, D, rihb, crb, rbi);
  // TODO: destruct original matrix on the fly here
  
}
