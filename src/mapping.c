#include <mapping.h>

map_fl_t *construct_fl_map(sm_t *M) {
  // initialize all map entries to __GB_MINUS_ONE_8
  map_fl_t *map = init_fl_map(M);

  uint32_t npiv = 0;  // number of pivots
  ri_t i = 0;         // current row index
  ri_t idx;          // possible pivot entry index

  // sweeps the rows to identify the row pivots and column pivots
  for (i=0; i<M->nrows; ++i) {
    if (M->rwidth[i] != 0) {
      idx = M->pos[i][0];
      if (map->pri[idx] == __GB_MINUS_ONE_32) {
        map->pri[idx] = i;
        npiv++;
      } else { // check for a sparser pivot row (see ELAGB talk from Lachartre)
        if (M->rwidth[map->pri[idx]] > M->rwidth[i]) {
          map->npri[map->pri[idx]]  = map->pri[idx];
          map->pri[idx] = i;
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
                      int nthreads, int verbose) {
  
  A = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  B = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  C = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));
  D = (sbm_fl_t *)malloc(sizeof(sbm_fl_t));

  int bdim  = block_dim / rows_multiline;
  // construct index map for Faugère-Lachartre decomposition of matrix M
  map  = construct_fl_map(M);
  
  printf("H1\n");
  int i, j, k, l, m;

  // initialize meta data for block submatrices
  A->nrows    = map->npiv;            // row dimension
  A->ncols    = map->npiv;            // col dimension
  A->bheight  = bdim;                 // block height
  A->bwidth   = bdim;                 // block width
  A->ba       = dtrl;                 // block alignment
  A->fe       = 0;                    // fill empty blocks?
  A->hr       = 0;                    // allow hybrid rows?

  // allocate memory for blocks
  
  // row and column loops 
  const uint32_t rlA  = (uint32_t) floor((float)A->nrows / A->bheight);
  const uint32_t clA  = (uint32_t) floor((float)A->ncols / A->bwidth);

  A->blocks = (mbl_t **)malloc(rlA * sizeof(mbl_t *));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (mbl_t *)malloc(clA * sizeof(mbl_t));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j].rows  = (re_t ***)malloc(
          (A->bheight / __GB_NROWS_MULTILINE) * sizeof(re_t **));
      /*
      for (k=0; k<(A->bheight / __GB_NROWS_MULTILINE); ++k) {
        A->blocks[i][j].rows[k] = (re_t **)malloc(
          A->bwidth * sizeof(re_t *));
        for (l=0; l<A->bwidth; ++l) {
          A->blocks[i][j].rows[k][l] = (re_t *)malloc(
            __GB_NROWS_MULTILINE * sizeof(re_t));
        }
      }
      */
    }
  }

  B->nrows    = map->npiv;            // row dimension
  B->ncols    = M->ncols - map->npiv; // col dimension
  B->bheight  = bdim;                 // block height
  B->bwidth   = bdim;                 // block width
  B->ba       = dtlr;                 // block alignment
  B->fe       = 1;                    // fill empty blocks?
  B->hr       = 1;                    // allow hybrid rows?

  // allocate memory for blocks

  // row and column loops 
  const uint32_t rlB  = (uint32_t) floor((float)B->nrows / B->bheight);
  const uint32_t clB  = (uint32_t) floor((float)B->ncols / B->bwidth);

  B->blocks = (mbl_t **)malloc(rlB * sizeof(mbl_t *));
  for (i=0; i<rlB; ++i) {
    B->blocks[i]  = (mbl_t *)malloc(clB * sizeof(mbl_t));
    for (j=0; j<clB; ++j) {
      B->blocks[i][j].rows  = (re_t ***)malloc(
          (B->bheight / __GB_NROWS_MULTILINE) * sizeof(re_t **));
      /*
      for (k=0; k<(B->bheight / __GB_NROWS_MULTILINE); ++k) {
        B->blocks[i][j].rows[k] = (re_t **)malloc(
          B->bwidth * sizeof(re_t *));
        for (l=0; l<B->bwidth; ++l) {
          B->blocks[i][j].rows[k][l] = (re_t *)malloc(
            __GB_NROWS_MULTILINE * sizeof(re_t));
        }
      }
      */
    }
  }

  C->nrows    = M->nrows - map->npiv; // row dimension
  C->ncols    = map->npiv;            // col dimension
  C->bheight  = bdim;                 // block height
  C->bwidth   = bdim;                 // block width
  C->ba       = dtrl;                 // block alignment
  C->fe       = 0;                    // fill empty blocks?
  C->hr       = 1;                    // allow hybrid rows?

  // allocate memory for blocks

  // row and column loops 
  const uint32_t rlC  = (uint32_t) floor((float)C->nrows / C->bheight);
  const uint32_t clC  = (uint32_t) floor((float)C->ncols / C->bwidth);

  C->blocks = (mbl_t **)malloc(rlC * sizeof(mbl_t *));
  for (i=0; i<rlC; ++i) {
    C->blocks[i]  = (mbl_t *)malloc(clC * sizeof(mbl_t));
    for (j=0; j<clC; ++j) {
      C->blocks[i][j].rows  = (re_t ***)malloc(
          (C->bheight / __GB_NROWS_MULTILINE) * sizeof(re_t **));
      /*
      for (k=0; k<(C->bheight / __GB_NROWS_MULTILINE); ++k) {
        C->blocks[i][j].rows[k] = (re_t **)malloc(
          C->bwidth * sizeof(re_t *));
        for (l=0; l<C->bwidth; ++l) {
          C->blocks[i][j].rows[k][l] = (re_t *)malloc(
            __GB_NROWS_MULTILINE * sizeof(re_t));
        }
      }
      */
    }
  }

  D->nrows    = M->nrows - map->npiv; // row dimension
  D->ncols    = M->ncols - map->npiv; // col dimension
  D->bheight  = bdim;                 // block height
  D->bwidth   = bdim;                 // block width
  D->ba       = dtlr;                 // block alignment
  D->fe       = 1;                    // fill empty blocks?
  D->hr       = 1;                    // allow hybrid rows?

  // allocate memory for blocks

  // row and column loops 
  const uint32_t rlD  = (uint32_t) floor((float)D->nrows / D->bheight);
  const uint32_t clD  = (uint32_t) floor((float)D->ncols / D->bwidth);

  D->blocks = (mbl_t **)malloc(rlD * sizeof(mbl_t *));
  for (i=0; i<rlD; ++i) {
    D->blocks[i]  = (mbl_t *)malloc(clD * sizeof(mbl_t));
    for (j=0; j<clD; ++j) {
      D->blocks[i][j].rows  = (re_t ***)malloc(
          (D->bheight / __GB_NROWS_MULTILINE) * sizeof(re_t **));
      /*
      for (k=0; k<(D->bheight / __GB_NROWS_MULTILINE); ++k) {
        D->blocks[i][j].rows[k] = (re_t **)malloc(
          D->bwidth * sizeof(re_t *));
        for (l=0; l<D->bwidth; ++l) {
          D->blocks[i][j].rows[k][l] = (re_t *)malloc(
            __GB_NROWS_MULTILINE * sizeof(re_t));
        }
      }
      */
    }
  }

  printf("M[%dx%d]\n",M->nrows,M->ncols);
  printf("A[%dx%d]\n",A->nrows,A->ncols);
  printf("B[%dx%d]\n",B->nrows,B->ncols);
  printf("C[%dx%d]\n",C->nrows,C->ncols);
  printf("D[%dx%d]\n",D->nrows,D->ncols);
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
  printf("size psi %d\n",(max_nrows / B->bheight) + 2);
  piv_start_idx[0]  = M->nrows;
  uint32_t block_idx;
  printf("npiv %d\n",map->npiv);

  // find blocks for construction of A & B
  printf("map %p\n",map);
  printf("npiv %d\n", npiv);
  for (i = (int)M->nrows-1; i > -1; --i) {
    if (map->pri[i] != __GB_MINUS_ONE_32)
      npiv++;

  printf("npiv %d\n", npiv);
    if (npiv % B->bheight == 0) {
      //printf("idx %d\n", idx);
      piv_start_idx[npiv/B->bheight]  = i;
    }
    //printf("(%d / %d) npiv2 %d\n",i, M->nrows, map->npiv);
  }
  printf("map2 %p\n",map);
    printf("npiv22 %d\n",map->npiv);

  // set leftout entries to zero
  for (i=npiv/B->bheight; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

    printf("npiv\t%d\n", map->npiv);
    printf("div\t%d\n", map->npiv/B->bheight);
  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  // rows indices horizontal block
    uint16_t cvb  = 0;          // current vector in block

    printf("npiv\t%d\n", map->npiv);
    printf("div\t%d\n", map->npiv/B->bheight);
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= (map->npiv/B->bheight); ++block_idx) {
      printf("bi %d\n", block_idx);
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping. 
      for (i = ((int)piv_start_idx[block_idx]-1);
           i >= (int)piv_start_idx[block_idx+1]; --i) {
        if (map->pri[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        if (cvb == B->bheight || i == 0) {
          write_blocks_lr_matrix(M, A, B, map, rihb, cvb, block_idx);

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
  piv_start_idx[0]  = M->nrows;
  npiv  = 0;
  for (i = (int)M->nrows-1; i > -1; --i) {
    if (map->npri[i] != __GB_MINUS_ONE_32)
      npiv++;

    if (npiv % B->bheight == 0) {
      piv_start_idx[npiv/B->bheight]  = i;
    }
  }

  // set leftout entries to zero
  for (i=npiv/B->bheight; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  // rows indices horizontal block
    uint16_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= (C->nrows/B->bheight); ++block_idx) {
      // construct block submatrices C & D
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping. 
      for (i = ((int)piv_start_idx[block_idx]-1);
           i >= (int)piv_start_idx[block_idx+1]; --i) {
        if (map->npri[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          write_blocks_lr_matrix(M, C, D, map, rihb, cvb, block_idx);

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

  // print inforamtion if verbose level >= 2
  if (verbose > 1) {
    printf("Number of pivots found: %d\n", map->npiv);
    printf("A [%d x %d]\n",A->nrows,A->ncols);
    printf("B [%d x %d]\n",B->nrows,B->ncols);
    printf("C [%d x %d]\n",C->nrows,C->ncols);
    printf("D [%d x %d]\n",D->nrows,D->ncols);
  }
}


void write_blocks_lr_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, map_fl_t *map,
                            ri_t *rihb, const ri_t cvb, const ri_t rbi) {

  uint32_t  bir;  // block index in row
  uint32_t  eil;  // element index in line
  uint32_t  lib;  // line index in block
  ci_t      it1, it2, i1, i2;

  // memory for block entries is already allocated in splice_fl_matrix()
  uint32_t i;
  const uint32_t loop_size  = (uint32_t) floor(cvb / __GB_NROWS_MULTILINE);

  // NOTE: In LELA Martani uses a stack vector to store the values for sub
  // block matrix A and first writes to B and this stack vector. Afterwards a
  // looping is done over the stack vector and all values are written into A in
  // one run.
  // Here we do not use this feature, but write directly into A.

  // NOTE: The following implementation (copied from Martani's LELA
  // implementation) works only for __GB_NROWS_MULTILINE = 2 !
  i1  = 0;
  i2  = 0;
  for (i=0; i<loop_size; i=i+2) {
    lib  = i/2;

    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (i1 != (M->rwidth[i]-1) && i2 != (M->rwidth[i+1]-1)) {
      it1 = M->pos[i][i1];
      it2 = M->pos[i+1][i2];
      if (it1 < it2) {
        if (map->pc[it1] != __GB_MINUS_ONE_32) {
          bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
          eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
          A->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
          A->blocks[rbi][bir].rows[lib][eil][1] = 0;
        } else {
          bir = map->npc[it1] / B->bwidth;
          eil = map->npc[it1] % B->bwidth;
          B->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
          B->blocks[rbi][bir].rows[lib][eil][1] = 0;
        }
        i1++;
      } else {
        if (it1 > it2) {
          if (map->pc[it1] != __GB_MINUS_ONE_32) {
            bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
            eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
            A->blocks[rbi][bir].rows[lib][eil][0] = 0;
            A->blocks[rbi][bir].rows[lib][eil][1] = M->rows[i+1][i2];
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            B->blocks[rbi][bir].rows[lib][eil][0] = 0;
            B->blocks[rbi][bir].rows[lib][eil][1] = M->rows[i+1][i2];
          }
          i2++;
        } else { // it1 == it2
          if (map->pc[it1] != __GB_MINUS_ONE_32) { // holds also for it2
            bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
            eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
            A->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
            A->blocks[rbi][bir].rows[lib][eil][1] = M->rows[i+1][i2];
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            B->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
            B->blocks[rbi][bir].rows[lib][eil][1] = M->rows[i+1][i2];
          }
          i1++;
          i2++;
        }
      }
    }
    // Now we splice rows i and i+1 separately (one of the rows we have already
    // completely spliced since we are out of the above while loop. Thus maximal
    // one of the following two while loops are executed (if rows i and i+1 have
    // the same number of elements none of the following while loops is
    // executed).
    while (i1 != (M->rwidth[i]-1)) {
      it1 = M->pos[i][i1];
      if (map->pc[it1] != __GB_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
        A->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
        A->blocks[rbi][bir].rows[lib][eil][1] = 0;
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        B->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
        B->blocks[rbi][bir].rows[lib][eil][1] = 0;
      }
      i1++;
    }
    while (i2 != (M->rwidth[i+1]-1)) {
      it2 = M->pos[i][i2];
      if (map->pc[it2] != __GB_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
        A->blocks[rbi][bir].rows[lib][eil][0] = 0;
        A->blocks[rbi][bir].rows[lib][eil][1] = M->rows[i+1][i2];
      } else {
        bir = map->npc[it2] / B->bwidth;
        eil = map->npc[it2] % B->bwidth;
        B->blocks[rbi][bir].rows[lib][eil][0] = 0;
        B->blocks[rbi][bir].rows[lib][eil][1] = M->rows[i+1][i2];
      }
      i2++;
    }
  }

  // Backup check if block size is not a multiple of 2:
  // There is maximal one line left thus we only work with the current i.
  if (i < cvb) {
    lib = i/2;
    while (i1 != (M->rwidth[i]-1)) {
      it1 = M->pos[i][i1];
      if (map->pc[it1] != __GB_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
        A->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
        A->blocks[rbi][bir].rows[lib][eil][1] = 0;
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        B->blocks[rbi][bir].rows[lib][eil][0] = M->rows[i][i1];
        B->blocks[rbi][bir].rows[lib][eil][1] = 0;
      }
      i1++;
    }
  }

  // hybrid multirows?
  if (B->hr) {
    // TODO: Implement hybrid stuff
  }
}
