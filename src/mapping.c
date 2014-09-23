#include <mapping.h>

#define __GB_DEBUG  0

void construct_fl_map(sm_t *M, map_fl_t *map) {
  // initialize all map entries to __GB_MINUS_ONE_8
  init_fl_map(M, map);

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
}

void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
                    map_fl_t *map, int block_dim, int rows_multiline,
                    int nthreads, int verbose) {

  int bdim  = block_dim / rows_multiline;
  // construct index map for Faugère-Lachartre decomposition of matrix M
  construct_fl_map(M, map);

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
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / A->bheight);
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);

  A->blocks = (mbl_t ***)malloc(rlA * sizeof(mbl_t **));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (mbl_t **)malloc(clA * sizeof(mbl_t *));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j] = (mbl_t *)malloc(
          (A->bheight / __GB_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(A->bheight / __GB_NROWS_MULTILINE); ++k) {
        A->blocks[i][j][k].val  = NULL;
        A->blocks[i][j][k].idx  = NULL;
        A->blocks[i][j][k].sz   = 0;
      }
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
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  B->blocks = (mbl_t ***)malloc(rlB * sizeof(mbl_t **));
  for (i=0; i<rlB; ++i) {
    B->blocks[i]  = (mbl_t **)malloc(clB * sizeof(mbl_t *));
    for (j=0; j<clB; ++j) {
      B->blocks[i][j] = (mbl_t *)malloc(
          (B->bheight / __GB_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(B->bheight / __GB_NROWS_MULTILINE); ++k) {
        B->blocks[i][j][k].val  = NULL;
        B->blocks[i][j][k].idx  = NULL;
        B->blocks[i][j][k].sz   = 0;
      }
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
  C->hr       = 0;                    // allow hybrid rows?

  // allocate memory for blocks

  // row and column loops 
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / C->bheight);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / C->bwidth);

  C->blocks = (mbl_t ***)malloc(rlC * sizeof(mbl_t **));
  for (i=0; i<rlC; ++i) {
    C->blocks[i]  = (mbl_t **)malloc(clC * sizeof(mbl_t *));
    for (j=0; j<clC; ++j) {
      C->blocks[i][j] = (mbl_t *)malloc(
          (C->bheight / __GB_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(C->bheight / __GB_NROWS_MULTILINE); ++k) {
        C->blocks[i][j][k].val  = NULL;
        C->blocks[i][j][k].idx  = NULL;
        C->blocks[i][j][k].sz   = 0;
      }
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
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);

  D->blocks = (mbl_t ***)malloc(rlD * sizeof(mbl_t **));
  for (i=0; i<rlD; ++i) {
    D->blocks[i]  = (mbl_t **)malloc(clD * sizeof(mbl_t *));
    for (j=0; j<clD; ++j) {
      D->blocks[i][j] = (mbl_t *)malloc(
          (D->bheight / __GB_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(D->bheight / __GB_NROWS_MULTILINE); ++k) {
        D->blocks[i][j][k].val  = NULL;
        D->blocks[i][j][k].idx  = NULL;
        D->blocks[i][j][k].sz   = 0;
      }
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

#if __GB_DEBUG
  printf("M[%dx%d]\n",M->nrows,M->ncols);
  printf("A[%dx%d]\n",A->nrows,A->ncols);
  printf("B[%dx%d]\n",B->nrows,B->ncols);
  printf("C[%dx%d]\n",C->nrows,C->ncols);
  printf("D[%dx%d]\n",D->nrows,D->ncols);
#endif
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
  uint32_t block_idx;

  uint16_t destruct_input_matrix  = 1;

  // find blocks for construction of A & B
  for (i = (int)M->ncols-1; i > -1; --i) {
    if (map->pri[i] != __GB_MINUS_ONE_32) {
      npiv++;
    }
    if (npiv % B->bheight == 0) {
      //printf("idx %d\n", idx);
      piv_start_idx[npiv/B->bheight]  = i;
    }
    //printf("(%d / %d) npiv2 %d\n",i, M->nrows, map->npiv);
  }

  // set leftout entries to zero
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  // rows indices horizontal block
    uint16_t cvb  = 0;          // current vector in block

#if __GB_DEBUG
    printf("nthreads %d\n",nthreads);
    printf("npiv %d -- bheight %d\n",map->npiv,B->bheight);
    printf("div %d\n", map->npiv/B->bheight);
#endif
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= npiv/A->bheight; ++block_idx) {
#if __GB_DEBUG
      printf("bi %d\n", block_idx);
      printf("piv_idx[block] %d\n",piv_start_idx[block_idx]);
      printf("piv_idx[block+1] %d\n",piv_start_idx[block_idx+1]);
#endif
      // construct block submatrices A & B
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping. 
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (map->pri[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        //printf("block_idx %d -- i %d\n",block_idx, i);
        //printf(" cvb %d / %d\n",cvb, B->bheight);
        if (cvb == B->bheight || i == 0) {
          write_blocks_lr_matrix_ml(M, A, B, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              free(M->rows[rihb[j]]);
            }
          // free memory
          }
          cvb = 0;
        }
      }
    }
  } 

  // find blocks for construction of C & D
  piv_start_idx[0]  = M->nrows;
  npiv  = 0;
  for (i = (int)M->ncols-1; i > -1; --i) {
    if (map->npri[i] != __GB_MINUS_ONE_32)
      npiv++;

    if (npiv % B->bheight == 0) {
      piv_start_idx[npiv/B->bheight]  = i;
    }
  }

  // set leftout entries to zero
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  // rows indices horizontal block
    uint16_t cvb  = 0;          // current vector in block

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= npiv/C->bheight; ++block_idx) {
      // construct block submatrices C & D
      // Note: In the for loop we always construct block "block+1" and not block
      // "block".
      // TODO: Try to improve this rather strange looping. 
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (map->npri[i] != __GB_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          write_blocks_lr_matrix_ml(M, C, D, map, rihb, cvb, block_idx);

          // TODO: Destruct input matrix on the go
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              free(M->rows[rihb[j]]);
            }
          }
          cvb = 0;
        }
      }
    }
  } 
}


void write_blocks_lr_matrix_ml(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, map_fl_t *map,
    ri_t *rihb, const ri_t cvb, const ri_t rbi) {

  bi_t  bir;    // block index in row
  bi_t  eil;    // element index in line
  bi_t  lib;    // line index in block
  bi_t  length; // local helper for block line length arithmetic
  ci_t  it1, it2, i1, i2;

  // memory for block entries is already allocated in splice_fl_matrix()

  // current loop variable i, block indices 1 (rihb[i]) and 2 (rihb[i+1])
  uint32_t i, j, k, l, bi1, bi2;
  //const uint32_t loop_size  = (uint32_t) ceil(cvb / __GB_NROWS_MULTILINE);


  // column loops 
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  // usually blocks in B tend to be denser than blocks in A, thus we allocate
  // already at the beginning more memory for those lines.
  // we use a density threshold (which is globally defined in gb_config.h.in) to
  // specify how much memory we allocate initially.
  bi_t init_bufferA, init_bufferB;
  if (M->density > __GB_DENSITY_THRESHOLD) {
    init_bufferA = A->bwidth / (A->bwidth/16);
    init_bufferB = B->bwidth / (B->bwidth/64);
  } else {
    init_bufferA = A->bwidth / (A->bwidth/8);
    init_bufferB = B->bwidth / (B->bwidth/16);
  }

  bi_t bufferA[clA];
  bi_t bufferB[clB];

  // NOTE: In LELA Martani uses a stack vector to store the values for sub
  // block matrix A and first writes to B and this stack vector. Afterwards a
  // looping is done over the stack vector and all values are written into A in
  // one run.
  // Here we do not use this feature, but write directly into A.

  // ususally cvb is divisible by 2, but for the last row of blocks there might
  // be only an odd number of lines in the blocks
  uint16_t rounded_cvb  = cvb;
  if (cvb % 2)
    rounded_cvb = cvb-1;
  for (i=0; i<rounded_cvb; i=i+2) {
    bi1 = rihb[i];
    bi2 = rihb[i+1];
    i1  = 0;
    i2  = 0;
    lib = i/2;
    memset(bufferA, 0, clA * sizeof(bi_t));
    memset(bufferB, 0, clB * sizeof(bi_t));

    //printf("ridxs %d -- %d -- %d\n",bi1,bi2,lib);

    // Note: The blocks of A are on or above the diagonal. Thus in heigth
    /*
    for (j=0; j<clA-rbi; ++j) {
      bufferA[j]  = init_bufferA;
      //printf("A_idx %p\n",A->blocks[rbi][j][lib].idx);
      //printf("A_val %p\n",A->blocks[rbi][j][lib].val);
      //printf("%d - %d - %d\n",rbi,j,lib);
      A->blocks[rbi][j][lib].val = (re_t *)malloc(2 * bufferA[j] * sizeof(re_t));
      A->blocks[rbi][j][lib].idx = (bi_t *)malloc(bufferA[j] * sizeof(bi_t));
    }
    for (j=0; j<clB; ++j) {
      bufferB[j]  = init_bufferB;
      //printf("B_idx %p\n",B->blocks[rbi][j][lib].idx);
      //printf("B_val %p\n",B->blocks[rbi][j][lib].val);
      printf("-- %d - %d - %d\n",rbi,j,lib);
      B->blocks[rbi][j][lib].idx = (bi_t *)malloc(bufferB[j] * sizeof(bi_t));
      B->blocks[rbi][j][lib].val = (re_t *)malloc(2 * bufferB[j] * sizeof(re_t));
    }
    */
    // loop over rows i and i+1 of M and splice correspondingly into A & B
    while (i1 < M->rwidth[bi1] && i2 < M->rwidth[bi2]) {
      it1 = M->pos[bi1][i1];
      it2 = M->pos[bi2][i2];
      if (it1 < it2) {
        if (map->pc[it1] != __GB_MINUS_ONE_32) {
          bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
          eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
          // realloc memory if needed
          if (A->blocks[rbi][bir][lib].sz == bufferA[bir]) {
            //printf("bufferA[bir] before %d\n",bufferA[bir]);
            realloc_block_rows(A, rbi, bir, lib, init_bufferA, &bufferA[bir]);
            //printf("bufferA[bir] after  %d\n",bufferA[bir]);
          }
          // set values
          insert_block_row_data_ml(A, M, rbi, bir, lib, eil, bi1, i1, 0, 0);
#if __GB_DEBUG
          printf("1 %d -- %d -- %d\n", rbi,bir,lib);
          printf("1 eil %d\n",eil);
#endif
        } else {
          bir = map->npc[it1] / B->bwidth;
          eil = map->npc[it1] % B->bwidth;
          // realloc memory if needed
          if (B->blocks[rbi][bir][lib].sz == bufferB[bir])
            realloc_block_rows(B, rbi, bir, lib, init_bufferB, &bufferB[bir]);
          // set values
          insert_block_row_data_ml(B, M, rbi, bir, lib, eil, bi1, i1, 0, 0);
#if __GB_DEBUG
          printf("2 %d -- %d -- %d\n", rbi,bir,lib);
          printf("2 eil %d\n",eil);
#endif
        }
        i1++;
      } else {
        if (it1 > it2) {
          if (map->pc[it2] != __GB_MINUS_ONE_32) {
            bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
            eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
            // realloc memory if needed
            if (A->blocks[rbi][bir][lib].sz == bufferA[bir])
              realloc_block_rows(A, rbi, bir, lib, init_bufferA, &bufferA[bir]);
            // set values
            insert_block_row_data_ml(A, M, rbi, bir, lib, eil, 0, 0, bi2, i2);
#if __GB_DEBUG
            printf("3 %d -- %d -- %d\n", rbi,bir,lib);
            printf("3 eil %d\n",eil);
#endif
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            // realloc memory if needed
            if (B->blocks[rbi][bir][lib].sz == bufferB[bir])
              realloc_block_rows(B, rbi, bir, lib, init_bufferB, &bufferB[bir]);
            // set values
            insert_block_row_data_ml(B, M, rbi, bir, lib, eil, 0, 0, bi2, i2);
#if __GB_DEBUG
            printf("4 %d -- %d -- %d\n", rbi,bir,lib);
            printf("4 eil %d\n",eil);
#endif
          }
          i2++;
        } else { // it1 == it2
          if (map->pc[it1] != __GB_MINUS_ONE_32) { // holds also for it2
            bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
            eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
            // realloc memory if needed
            if (A->blocks[rbi][bir][lib].sz == bufferA[bir])
              realloc_block_rows(A, rbi, bir, lib, init_bufferA, &bufferA[bir]);
            // set values
            insert_block_row_data_ml(A, M, rbi, bir, lib, eil, bi1, i1, bi2, i2);
#if __GB_DEBUG
            printf("5 %d -- %d -- %d\n", rbi,bir,lib);
            printf("5 eil %d\n",eil);
#endif
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            // realloc memory if needed
            if (B->blocks[rbi][bir][lib].sz == bufferB[bir])
              realloc_block_rows(B, rbi, bir, lib, init_bufferB, &bufferB[bir]);
            // set values
            insert_block_row_data_ml(B, M, rbi, bir, lib, eil, bi1, i1, bi2, i2);
#if __GB_DEBUG
            printf("6 %d -- %d -- %d\n", rbi,bir,lib);
            printf("6 eil %d\n",eil);
            printf("1 %d\n",M->rows[bi1][i1]);
            printf("2 %d\n",M->rows[bi2][i2]);
#endif
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
    while (i1 < (M->rwidth[bi1])) {
      it1 = M->pos[bi1][i1];
      if (map->pc[it1] != __GB_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
        // realloc memory if needed
        if (A->blocks[rbi][bir][lib].sz == bufferA[bir])
          realloc_block_rows(A, rbi, bir, lib, init_bufferA, &bufferA[bir]);
        // set values
        insert_block_row_data_ml(A, M, rbi, bir, lib, eil, bi1, i1, 0, 0);
#if __GB_DEBUG
        printf("7 %d -- %d -- %d\n", rbi,bir,lib);
        printf("7 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        // realloc memory if needed
        if (B->blocks[rbi][bir][lib].sz == bufferB[bir])
          realloc_block_rows(B, rbi, bir, lib, init_bufferB, &bufferB[bir]);
        // set values
        insert_block_row_data_ml(B, M, rbi, bir, lib, eil, bi1, i1, 0, 0);
#if __GB_DEBUG
        printf("8 %d -- %d -- %d\n", rbi,bir,lib);
        printf("8 eil %d\n",eil);
#endif
      }
      i1++;
    }
    while (i2 < (M->rwidth[bi2])) {
      it2 = M->pos[bi2][i2];
      //printf("it2 %d\n",it2);
      if (map->pc[it2] != __GB_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
        // realloc memory if needed
        if (A->blocks[rbi][bir][lib].sz == bufferA[bir])
          realloc_block_rows(A, rbi, bir, lib, init_bufferA, &bufferA[bir]);
        // set values
        insert_block_row_data_ml(A, M, rbi, bir, lib, eil, 0, 0, bi2, i2);
#if __GB_DEBUG
        printf("9 %d -- %d -- %d\n", rbi,bir,lib);
        printf("9 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it2] / B->bwidth;
        eil = map->npc[it2] % B->bwidth;
        // realloc memory if needed
        if (B->blocks[rbi][bir][lib].sz == bufferB[bir])
          realloc_block_rows(B, rbi, bir, lib, init_bufferB, &bufferB[bir]);
        // set values
        insert_block_row_data_ml(B, M, rbi, bir, lib, eil, 0, 0, bi2, i2);
#if __GB_DEBUG
        printf("10 %d -- %d -- %d\n", rbi,bir,lib);
        printf("10 eil %d\n",eil);
#endif
      }
      i2++;
    }

  }

#if __GB_DEBUG
  printf("i %d -- cvb %d\n",i,cvb);
#endif
  if (i<cvb) {
    bi1 = rihb[i];
    i1  = 0;
    lib = i/2;
    memset(bufferA, 0, clA * sizeof(bi_t));
    memset(bufferB, 0, clB * sizeof(bi_t));
    /*
    for (j=0; j<clA; ++j) {
      bufferA[j]  = init_bufferA;
      //printf("A_idx %p\n",A->blocks[rbi][j][lib].idx);
      //printf("A_val %p\n",A->blocks[rbi][j][lib].val);
      //printf("%d - %d - %d\n",rbi,j,lib);
      A->blocks[rbi][j][lib].val = (re_t *)malloc(2 * bufferA[j] * sizeof(re_t));
      A->blocks[rbi][j][lib].idx = (bi_t *)malloc(bufferA[j] * sizeof(bi_t));
    }
    for (j=0; j<clB; ++j) {
      bufferB[j]  = init_bufferB;
      //printf("B_idx %p\n",B->blocks[rbi][j][lib].idx);
      //printf("B_val %p\n",B->blocks[rbi][j][lib].val);
      //printf("%d - %d - %d\n",rbi,j,lib);
      B->blocks[rbi][j][lib].idx = (bi_t *)malloc(bufferB[j] * sizeof(bi_t));
      B->blocks[rbi][j][lib].val = (re_t *)malloc(2 * bufferB[j] * sizeof(re_t));
    }
    */
    while (i1 < (M->rwidth[bi1])) {
      it1 = M->pos[bi1][i1];
      if (map->pc[it1] != __GB_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
        // realloc memory if needed
        if (A->blocks[rbi][bir][lib].sz == bufferA[bir])
          realloc_block_rows(A, rbi, bir, lib, init_bufferA, &bufferA[bir]);
        // set values
        insert_block_row_data_ml(A, M, rbi, bir, lib, eil, bi1, i1, 0, 0);
#if __GB_DEBUG
        printf("11 %d -- %d -- %d\n", rbi,bir,lib);
        printf("11 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        // realloc memory if needed
        if (B->blocks[rbi][bir][lib].sz == bufferB[bir])
          realloc_block_rows(B, rbi, bir, lib, init_bufferB, &bufferB[bir]);
        // set values
        insert_block_row_data_ml(B, M, rbi, bir, lib, eil, bi1, i1, 0, 0);
#if __GB_DEBUG
        printf("12 %d -- %d -- %d\n", rbi,bir,lib);
        printf("12 eil %d\n",eil);
#endif
      }
      i1++;
    }
  }
  // write data to A (from right to left)
  swap_block_data(A, clA, rbi, cvb);
  // Realloc memory usage:
  // We do not need to realloc memory for A since we have already fixed memory
  // for A in the swapping of the data above.
  for (j=0; j<clB; ++j) {
    B->blocks[rbi][j][lib].idx = realloc(
        B->blocks[rbi][j][lib].idx,
        B->blocks[rbi][j][lib].sz * sizeof(bi_t));
    B->blocks[rbi][j][lib].val = realloc(
        B->blocks[rbi][j][lib].val,
        2 * B->blocks[rbi][j][lib].sz  * sizeof(re_t));
  }

  // hybrid multirows for the righthand side block matrices?
  if (B->hr) {
    uint32_t idx  = 0;
    // TODO: Implement hybrid stuff
    for (i=0; i<clB; ++i) {
      for (j=0; j<B->bheight/__GB_NROWS_MULTILINE; ++j) {
        if ((float)B->blocks[rbi][i][j].sz / (float)B->bwidth
            < __GB_HYBRID_THRESHOLD)
          continue;
        re_t *tmp_val_ptr = (re_t *)malloc(2 * B->bwidth * sizeof(re_t));
        idx  = 0;
        for (k=0; k<B->bwidth; ++k) {
          if (idx < B->blocks[rbi][i][j].sz && B->blocks[rbi][i][j].idx[idx] == k) {
            tmp_val_ptr[2*k]    = B->blocks[rbi][i][j].val[2*idx];
            tmp_val_ptr[2*k+1]  = B->blocks[rbi][i][j].val[2*idx+1];
            idx++;
          } else {
            tmp_val_ptr[2*k]    = 0;
            tmp_val_ptr[2*k+1]  = 0;
          }
        }
        free(B->blocks[rbi][i][j].idx);
        free(B->blocks[rbi][i][j].val);
        B->blocks[rbi][i][j].val  = tmp_val_ptr;
      }
    }
  }
}
