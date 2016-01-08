/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */


#include <mapping.h>

#define __GBLA_DEBUG_LL     0
#define __GBLA_DEBUG        0
#define __GBLA_NEW_SPLICER  0
#define __GBLA_BENCHMARKS  1

void construct_fl_map(sm_t *M, map_fl_t *map) {
  /*  initialize all map entries to __GBLA_MINUS_ONE_8 */
  init_fl_map(M, map);

  /* ri_t max_length = M->ncols > M->nrows ? M->ncols : M->n */
  uint32_t npiv   = 0;  /*  number of pivots */
  ri_t i = 0;         /*  current row index */
  ri_t idx;          /*  possible pivot entry index */

  /*  sweeps the rows to identify the row pivots and column pivots */
  for (i=0; i<M->nrows; ++i) {
    if (M->rwidth[i] != 0) {
      idx = M->pos[i][0];
      if (map->pri[idx] == __GBLA_MINUS_ONE_32) {
        map->pri[idx] = i;
        npiv++;
      } else { /*  check for a sparser pivot row (see ELAGB talk from Lachartre) */
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

  /*  construct pivot columns and non-pivot columns maps and the corresponding */
  /*  reverse maps */
  for (j=0; j<M->ncols; ++j) {
    if (map->pri[j] !=  __GBLA_MINUS_ONE_32) {
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

void construct_fl_map_reduced(map_fl_t *map, map_fl_t *map_D, ri_t rows_BD,
    ri_t rank_D, ci_t coldim, int nthreads) {
  ri_t i;
	/* ci_t j; */

	/* ri_t rows_B = rows_BD - rank_D; */
  /*  pivot and non-pivot columns are exactly the ones from map_D, so we just */
  /*  point to them */
  map->pc       = map_D->pc;
  map->npc      = map_D->npc;
  map->pc_rev   = map_D->pc_rev;
  map->npc_rev  = map_D->npc_rev;

  /*  for the rows we have to take the values from map_D, but we have to adjust */
  /*  them by rows_B down. So we have to create new mappings */
  map->npri = (ri_t *)malloc((rows_BD) * sizeof(ri_t));
  /* memset(map->npri, __GBLA_MINUS_ONE_8, (rows_BD) * sizeof(ri_t)); */

  /*  number of pivots stays the same, i.e. rank_D */
  map->npiv = map_D->npiv;

  /*  pivot rows stay the same, too */
  map->pri  = map_D->pri;

  /*  all rows in B are non-pivot rows */
#pragma omp parallel num_threads(nthreads)
  {
#pragma omp for private(i) schedule(dynamic)
  for (i=0; i<rows_BD; ++i)
    map->npri[i]  = i;
  }
#pragma omp parallel num_threads(nthreads)
  {
#pragma omp for private(i) schedule(dynamic)
  for (i=0; i<coldim; ++i)
    if (map->pri[i] != __GBLA_MINUS_ONE_32)
      map->npri[map->pri[i]]  = __GBLA_MINUS_ONE_32;
  }
}

void process_matrix(sm_fl_ml_t *A, map_fl_t *map, const bi_t bheight) {
  uint32_t i;
  uint32_t curr_row_idx = 0;
  int h1, h2;
  uint32_t h1_idx, h2_idx;
  re_t h1_val, h2_val;
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GBLA_NROWS_MULTILINE);
	/* const uint32_t clA  = (uint32_t) ceil((float)A->ncols / __GBLA_NROWS_MULTILINE); */
  const ci_t coldim = A->ncols;

  init_fl_map_sizes(A->ncols + (A->ncols % bheight), A->nrows + (A->nrows % bheight), map);

  map->npiv = 0;

  for (i=0; i<rlA; ++i) {
    if (A->ml[i].sz > 0) {
      h1 = get_head_multiline_hybrid(&(A->ml[i]), 0,
          &h1_val, &h1_idx, coldim);
      h2 = get_head_multiline_hybrid(&(A->ml[i]), 1,
          &h2_val, &h2_idx, coldim);
      /*  Line 1 */
      if (h1 == -1) {
        map->npri[curr_row_idx] = curr_row_idx;
      } else {
        if (map->pri[h1] == __GBLA_MINUS_ONE_32) {
          map->pri[h1]  = curr_row_idx;
          map->npiv++;
        } else {
          map->npri[curr_row_idx] = curr_row_idx;
        }
      }
      ++curr_row_idx;
      /*  Line 2 */
      if (h2 == -1) {
        map->npri[curr_row_idx] = curr_row_idx;
      } else {
        if (map->pri[h2] == __GBLA_MINUS_ONE_32) {
          map->pri[h2]  = curr_row_idx;
          map->npiv++;
        } else {
          map->npri[curr_row_idx] = curr_row_idx;
        }
      }
      ++curr_row_idx;
    } else {
      map->npri[curr_row_idx]   = curr_row_idx;
      map->npri[curr_row_idx+1] = curr_row_idx + 1;
      curr_row_idx  +=  2;
    }
  }

  /*  construct column pivots and their reverses */
  uint32_t piv_col_idx  = 0, non_piv_col_idx  = 0;
  for (i=0; i<coldim; ++i) {
    if (map->pri[i] != __GBLA_MINUS_ONE_32) {
      map->pc[i]                = piv_col_idx;
      map->pc_rev[piv_col_idx]  = i;
      piv_col_idx++;
    } else {
      map->npc[i]                   = non_piv_col_idx;
      map->npc_rev[non_piv_col_idx] = i;
      non_piv_col_idx++;
    }
  }
}

void combine_maps(map_fl_t *outer_map, map_fl_t **inner_map_in,
     const ci_t outer_coldim, const ci_t inner_coldim, int only_rows) {

#if DEBUG_RECONSTRUCT
  printf("inner coldim %d\n",inner_coldim);
#endif
  map_fl_t *inner_map = *inner_map_in;

  uint32_t i, idx_B, idx_B2, rev_idx_A;

  /*  remap pivot rows indices (lines of matrix A) */
  uint32_t next_piv = 0;
  for (i=0; i<outer_coldim; ++i) {
    if (outer_map->pri[i] == __GBLA_MINUS_ONE_32)
      continue;
    outer_map->pri[i] = next_piv++;
  }

  if (only_rows == 1) {
    for (i=0; i<inner_coldim; ++i) {
#if DEBUG_RECONSTRUCT
      printf("i %d, %d\n",i,inner_map->pri[i]);
#endif
      if (inner_map->pri[i] == __GBLA_MINUS_ONE_32)
        continue;
      idx_B     = i;
      rev_idx_A = outer_map->npc_rev[idx_B];
#if DEBUG_RECONSTRUCT
      printf("outer_map->pri[%d] = %d\n",rev_idx_A,outer_map->npiv + inner_map->pri[i]);
#endif
      outer_map->pri[rev_idx_A] = outer_map->npiv + inner_map->pri[i];
    }
		/* return; */
    goto cleanup ;
  }

  /*  remap the pivot columns of D */
  for (i=0; i<inner_map->npiv; ++i) {
    idx_B     = inner_map->pc_rev[i];
    rev_idx_A = outer_map->npc_rev[idx_B];

    /*  the pivot column at index rev_idx_A points to a row in D, */
    /*  i. e. outer_map->npiv rows + relative index in D */
    outer_map->pri[rev_idx_A] = i + outer_map->npiv;
  }

  /*  remap the nonpivot columns from inner_map to this map */
  for (i=0; i<outer_coldim; ++i) {
    if (outer_map->npc[i] == __GBLA_MINUS_ONE_32)
      continue;

    /*  point index of the nonpivot column in matrix B = B1|B2 to the index in */
    /*  the matrix B2 */
    idx_B = outer_map->npc[i];

    /*  column is a pivot column in B, skip it */
    if (inner_map->pc[idx_B] != __GBLA_MINUS_ONE_32)
      continue;

    idx_B2  = inner_map->npc[idx_B];

    outer_map->npc[i]           = idx_B2;
    outer_map->npc_rev[idx_B2]  = i;
  }

cleanup:
  /*  free memory for inner map */
  free(inner_map->pc);
  free(inner_map->npc);
  free(inner_map->pc_rev);
  free(inner_map->npc_rev);
  free(inner_map->pri);
  free(inner_map->npri);
  free(inner_map);
  inner_map     = NULL;
  *inner_map_in = inner_map;
}

void reconstruct_matrix_block_reduced(sm_t *M, sbm_fl_t *A, sbm_fl_t *B2, sbm_fl_t *D2,
    map_fl_t *map, const ci_t coldim, int free_matrices,
    int M_freed, int A_freed, int D_freed, int nthrds) {

  int ret;
  uint8_t *vec_free_B2, *vec_free_D2;
  /*  1 pivot from A and otherwise at most D->ncols = B->ncols elements */
  /*  ==> no reallocations! */
  /* ci_t buffer = (ci_t) ceil((float)M->ncols / 10); */

  if (free_matrices == 1) {
    ret = posix_memalign((void *)&vec_free_B2, ALIGNT, (B2->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    memset(vec_free_B2, 0, (B2->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    ret = posix_memalign((void *)&vec_free_D2, ALIGNT, (D2->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    memset(vec_free_D2, 0, (D2->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
  }

  uint32_t i, j, k;
	/* uint32_t new_piv_to_row[coldim]; */
  uint32_t new_piv  = 0;
  for (i=0; i<coldim; ++i) {
    if (map->pri[i] == __GBLA_MINUS_ONE_32)
      continue;
		/* new_piv_to_row[i] = new_piv; */
    new_piv++;
  }

  /*  maximal row length in M is 1 + B2->ncols = 1 + D2->ncols */
  ci_t max_row_length = 1 + B2->ncols;
  const uint32_t rlB  = (uint32_t) ceil((float)B2->nrows / B2->bheight);
  const uint32_t clB  = (uint32_t) ceil((float)B2->ncols / B2->bwidth);
  const uint32_t rlD  = (uint32_t) ceil((float)D2->nrows / D2->bheight);
  const uint32_t clD  = (uint32_t) ceil((float)D2->ncols / D2->bwidth);
  M->nnz  = 0;

  int ctr = 0;
	/* ci_t row_M_width; */
	/* ci_t tmp_width; */
  uint32_t local_piv = 0;
  uint16_t line;
  for (i=0; i<coldim; ++i) {
    if (map->pri[i] == __GBLA_MINUS_ONE_32)
      continue;
    ctr++;

    M->rows[local_piv]  = realloc(M->rows[local_piv], max_row_length  * sizeof(re_t));
    M->pos[local_piv]   = realloc(M->pos[local_piv], max_row_length * sizeof(ci_t));
    M->rwidth[local_piv]  = 0;
    ci_t block_row_idx, row_idx_block, block_start_idx, idx;
    re_t val;
    if (map->pri[i] >= map->npiv) { /*  we are in D */
      block_row_idx = (D2->nrows - 1 - (map->pri[i] - map->npiv)) / D2->bheight;
      row_idx_block = ((D2->nrows - 1 - (map->pri[i] - map->npiv)) % D2->bheight)
        / __GBLA_NROWS_MULTILINE;

      line  = (D2->nrows - 1 - (map->pri[i] - map->npiv)) % __GBLA_NROWS_MULTILINE;

      M->rows[local_piv][M->rwidth[local_piv]] = (re_t)1;
      M->pos[local_piv][M->rwidth[local_piv]]  = i;
      M->rwidth[local_piv]++;

      for (j=0; j<clD; ++j) {
        if (D2->blocks[block_row_idx][j] == NULL)
          continue;
        if ((D2->blocks[block_row_idx][j][row_idx_block].val == NULL) ||
            (D2->blocks[block_row_idx][j][row_idx_block].sz == 0))
          continue;

        block_start_idx = D2->bwidth * j;

        if (D2->blocks[block_row_idx][j][row_idx_block].dense == 0) {
          for (k=0; k<D2->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
            idx = D2->blocks[block_row_idx][j][row_idx_block].idx[k];
            val = D2->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
            if (val != 0) {
              M->rows[local_piv][M->rwidth[local_piv]] = val;
              M->pos[local_piv][M->rwidth[local_piv]]  = block_start_idx + idx;
              M->rwidth[local_piv]++;
            }
          }
        } else { /*  D2 block multiline is dense */
          for (k=0; k<D2->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
            val = D2->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
            if (val != 0) {
              M->rows[local_piv][M->rwidth[local_piv]] = val;
              M->pos[local_piv][M->rwidth[local_piv]]  = block_start_idx + k;
              M->rwidth[local_piv]++;
            }
          }
        }
        if (free_matrices == 1) {
          if (vec_free_D2[(D2->nrows - 1 - (map->pri[i] - map->npiv))/__GBLA_NROWS_MULTILINE] == 1) {
            if (D2->blocks[block_row_idx][j][row_idx_block].dense == 0) {
              free (D2->blocks[block_row_idx][j][row_idx_block].idx);
              D2->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
            }
            free (D2->blocks[block_row_idx][j][row_idx_block].val);
            D2->blocks[block_row_idx][j][row_idx_block].val  = NULL;
          }
        }
      }
      if (free_matrices == 1) {
        if (vec_free_D2[(D2->nrows - 1 - (map->pri[i] - map->npiv))/__GBLA_NROWS_MULTILINE] == 0) {
          vec_free_D2[(D2->nrows - 1 - (map->pri[i] - map->npiv))/__GBLA_NROWS_MULTILINE]  = 1;
        }
      }
    } else { /*  we are in A resp. B */
      block_row_idx = (B2->nrows - 1 - map->pri[i]) / B2->bheight;
      row_idx_block = ((B2->nrows - 1 - map->pri[i]) % B2->bheight) / __GBLA_NROWS_MULTILINE;

      line  = (B2->nrows - 1 - map->pri[i]) % __GBLA_NROWS_MULTILINE;

      /*  A should be always freed already! */
      if (A_freed == 0) {
      } else { /*  add 1 at pivot position for A part */
        M->rows[local_piv][M->rwidth[local_piv]] = (re_t)1;
        M->pos[local_piv][M->rwidth[local_piv]]  = map->pc[i];
        M->rwidth[local_piv]++;
      }
      /*  add B part to row in M */
      for (j=0; j<clB; ++j) {
        if (B2->blocks[block_row_idx][j] == NULL)
          continue;
        if (B2->blocks[block_row_idx][j][row_idx_block].val == NULL)
          continue;
        if (B2->blocks[block_row_idx][j][row_idx_block].sz == 0) {
          free(B2->blocks[block_row_idx][j][row_idx_block].val);
          B2->blocks[block_row_idx][j][row_idx_block].val  = NULL;
          free(B2->blocks[block_row_idx][j][row_idx_block].idx);
          B2->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
          continue;
        }

        block_start_idx = B2->bwidth * j;

        if (B2->blocks[block_row_idx][j][row_idx_block].dense == 0) {
          for (k=0; k<B2->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
            idx = B2->blocks[block_row_idx][j][row_idx_block].idx[k];
            val = B2->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
            if (val != 0) {
              M->rows[local_piv][M->rwidth[local_piv]] = val;
              M->pos[local_piv][M->rwidth[local_piv]]  = block_start_idx + idx;
              M->rwidth[local_piv]++;
            }
          }
        } else { /*  B2 block multiline is dense */
          for (k=0; k<B2->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
            val = B2->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
            if (val != 0) {
              M->rows[local_piv][M->rwidth[local_piv]] = val;
              M->pos[local_piv][M->rwidth[local_piv]]  = block_start_idx + k;
              M->rwidth[local_piv]++;
            }
          }
        }
        if (free_matrices == 1) {
          if (vec_free_B2[(B2->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 1) {
            if (B2->blocks[block_row_idx][j][row_idx_block].dense == 0) {
              free (B2->blocks[block_row_idx][j][row_idx_block].idx);
              B2->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
            }
            /* printf("frees block [%d][%d][%d]\n", block_row_idx,j,row_idx_block); */
            free (B2->blocks[block_row_idx][j][row_idx_block].val);
            B2->blocks[block_row_idx][j][row_idx_block].val  = NULL;
          }
        }
      }
      if (free_matrices == 1) {
        if (vec_free_B2[(B2->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 0) {
          vec_free_B2[(B2->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE]  = 1;
        }
      }
    }

    /*  adjust memory, rows in M are done */
    M->nnz  = M->nnz + M->rwidth[local_piv];
    M->rows[local_piv] = realloc(M->rows[local_piv],
        (M->rwidth[local_piv]) * sizeof(re_t));
    M->pos[local_piv] = realloc(M->pos[local_piv],
        (M->rwidth[local_piv]) * sizeof(ci_t));

    /*  next pivot row */
    local_piv++;
  }

  /*  free memory of B2 and D2 */
  if (free_matrices == 1) {
    free (vec_free_B2);
    free (vec_free_D2);
    for (i=0; i<rlB; ++i) {
      for (j=0; j<clB; ++j) {
        free(B2->blocks[i][j]);
        B2->blocks[i][j] = NULL;
      }
      free(B2->blocks[i]);
      B2->blocks[i]  = NULL;
    }
    free(B2->blocks);
    free (B2);
    B2 = NULL;

    for (i=0; i<rlD; ++i) {
      for (j=0; j<clD; ++j) {
        free(D2->blocks[i][j]);
        D2->blocks[i][j] = NULL;
      }
      free(D2->blocks[i]);
      D2->blocks[i]  = NULL;
    }
    free(D2->blocks);
    free (D2);
    D2 = NULL;
  }

  /*  free memory for map */
  free(map->pc);
  free(map->npc);
  free(map->pc_rev);
  free(map->npc_rev);
  free(map->pri);
  free(map->npri);
  free(map);
  map = NULL;
}


void reconstruct_matrix_no_multiline_keep_A(sm_t *M, sm_fl_t *A, sb_fl_t *B,
    dm_t *D, map_fl_t *map, const int nthrds)
{
  ri_t i, j, k, l;

  const ri_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const ci_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const ri_t rlD  = (uint32_t) ceil((float)D->rank / __GBLA_SIMD_BLOCK_SIZE);

  M->nnz    = 0;
  M->nrows  = map->npiv + D->rank;
  M->rows   = (re_t **)malloc(M->nrows * sizeof(re_t *));
  M->pos    = (ci_t **)malloc(M->nrows * sizeof(ci_t *));

  // reallocate memory for rows in M: they have a max length of
  // 1+D->ncols for rows from A+B, and
  // D->ncols for rows from D
#pragma omp parallel num_threads(nthrds)
{
  ri_t block_row_idx;
  ri_t line_idx;
  ci_t start_idx;
  ri_t min_range;
  ri_t min_range_blocks;
#pragma omp for private(i,j,k) schedule(dynamic)
  // each thread gets a whole block row
  for (l=0; l<rlB; ++l) {
    block_row_idx = l*__GBLA_SIMD_BLOCK_SIZE;
    // getting a row from A+B
    if (l<rlB-1) {
      min_range         = (l+1)*__GBLA_SIMD_BLOCK_SIZE;
      min_range_blocks  = __GBLA_SIMD_BLOCK_SIZE;
    } else {
      min_range         = map->npiv;
      min_range_blocks  = map->npiv - (rlB-1)*__GBLA_SIMD_BLOCK_SIZE;
    }

    // allocate memory and set first elements for the whole block
    for (i=block_row_idx; i<min_range; ++i) {
      M->rows[i]    = (re_t *)malloc((A->sz[i]+D->ncols) * sizeof(re_t));
      M->pos[i]     = (ci_t *)malloc((A->sz[i]+D->ncols) * sizeof(ci_t));
      M->rwidth[i]  = 0;

      // do A
      for (j=0; j<A->sz[i]; ++j) {
        M->rows[i][M->rwidth[i]]  = A->row[i][j];
        M->pos[i][M->rwidth[i]]   = map->pc_rev[A->pos[i][j]];
        M->rwidth[i]++;
      }
      free(A->row[i]);
      free(A->pos[i]);
    }

    // now write B part, the blocks of B are stored in column-first
    // representation!
    for (i=0; i<clB; ++i) {
      if (B->blocks[l][i].val ==  NULL)
        continue;
      start_idx = i * __GBLA_SIMD_BLOCK_SIZE;
      for ( j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
        if (B->blocks[l][i].val[j] == NULL)
          continue;
        for ( k=0; k<B->blocks[l][i].sz[j]; ++k) {
          M->rows[block_row_idx+B->blocks[l][i].pos[j][k]][M->rwidth[block_row_idx+B->blocks[l][i].pos[j][k]]]  =
            B->blocks[l][i].val[j][k];
          M->pos[block_row_idx+B->blocks[l][i].pos[j][k]][M->rwidth[block_row_idx+B->blocks[l][i].pos[j][k]]]   =
            map->npc_rev[map->npiv+j+start_idx];
          M->rwidth[block_row_idx+B->blocks[l][i].pos[j][k]]++;
        }
      }
    }
    // remove data for blocks of B
    for (j=0; j<clB; ++j) {
      free(B->blocks[l][j].val);
      B->blocks[l][j].val = NULL;
    }
    //printf("rwidth[%u] = %u\n",i,M->rwidth[i]);
    // relloc row size
    for (i=block_row_idx; i<min_range; ++i) {
      M->rows[i]  = realloc(M->rows[i], M->rwidth[i] * sizeof(re_t));
      M->pos[i]   = realloc(M->pos[i], M->rwidth[i] * sizeof(ci_t));
    }
    // if we had the last line in the block, we can remove it
#if __GBLA_BENCHMARKS
    for (i=block_row_idx; i<min_range; ++i) {
      free(M->rows[i]);
      free(M->pos[i]);
      M->rows[i]  = NULL;
      M->pos[i]   = NULL;
    }
#endif
  }
}

#pragma omp parallel num_threads(nthrds)
{
  ri_t block_row_idx;
  ri_t line_idx;
  ci_t start_idx;
  ri_t min_range;
  ri_t min_range_blocks;
#pragma omp for private(i,j,k) schedule(dynamic)
  for (l=0; l<rlD; ++l) {
    block_row_idx = l*__GBLA_SIMD_BLOCK_SIZE + map->npiv;
    // getting rows from D
    if (l<rlD-1) {
      min_range = (l+1)*__GBLA_SIMD_BLOCK_SIZE + map->npiv;
    } else {
      min_range = M->nrows;
    }

    // allocate memory and set first elements for the whole block
    for (i=block_row_idx; i<min_range; ++i) {
      line_idx      = i-map->npiv;
      M->rows[i]    = (re_t *)malloc((D->ncols-D->row[line_idx]->piv_lead+1) * sizeof(re_t));
      M->pos[i]     = (ci_t *)malloc((D->ncols-D->row[line_idx]->piv_lead+1) * sizeof(ci_t));
      M->rwidth[i]  = 0;
    }

    for (i=block_row_idx; i<min_range; ++i) {
      line_idx  = i-map->npiv;
      if (D->row[line_idx]->piv_val != NULL) {
        for (j=D->row[line_idx]->piv_lead; j<D->ncols; ++j) {
          if (D->row[line_idx]->piv_val[j] != 0) {
            M->rows[i][M->rwidth[i]]  = D->row[line_idx]->piv_val[j];
            M->pos[i][M->rwidth[i]]   = map->npc_rev[map->npiv+j];
            M->rwidth[i]++;
          }
        }
        // relloc row size
        M->rows[i]  = realloc(M->rows[i], M->rwidth[i] * sizeof(re_t));
        M->pos[i]   = realloc(M->pos[i], M->rwidth[i] * sizeof(ci_t));
        // free row in D
        free(D->row[line_idx]->piv_val);
        free(D->row[line_idx]);
        D->row[line_idx] = NULL;
      }
    }
#if __GBLA_BENCHMARKS
    for (i=block_row_idx; i<min_range; ++i) {
      free(M->rows[i]);
      free(M->pos[i]);
      M->rows[i]  = NULL;
      M->pos[i]   = NULL;
    }
#endif
  }
}

  for (i=0; i<M->nrows; ++i)
    M->nnz  +=  M->rwidth[i];

  // free A,  B and D
  free(A->row);
  free(A->pos);
  free(A->sz);
  free(A->buf);
  A = NULL;

  for (i=0; i<rlB; ++i) {
    for (j=0; j<clB; ++j) {
      free(B->blocks[i][j].sz);
      free(B->blocks[i][j].buf);
    }
  }
  free(B->blocks);
  free(B);
  B = NULL;

  free(D);
  D = NULL;

  // free map
  free(map->pc);
  free(map->npc);
  free(map->pc_rev);
  free(map->npc_rev);
  free(map->pri);
  free(map->npri);
  free(map);
  map = NULL;
}

void reconstruct_matrix_block_no_multiline(sm_t *M, sb_fl_t *A, dbm_fl_t *B, dm_t *D,
    map_fl_t *map, const int nthrds)
{
  /*  1 pivot from A and otherwise at most D->ncols = B->ncols elements */
  /*  ==> no reallocations! */
  /* ci_t buffer = (ci_t) ceil((float)M->ncols / 10); */

  ri_t i, j, k, l;

  const ri_t rlB  = (uint32_t) ceil((float)B->nrows / __GBLA_SIMD_BLOCK_SIZE);
  const ci_t clB  = (uint32_t) ceil((float)B->ncols / __GBLA_SIMD_BLOCK_SIZE);
  const ri_t rlD  = (uint32_t) ceil((float)D->rank / __GBLA_SIMD_BLOCK_SIZE);

  M->nnz    = 0;
  M->nrows  = map->npiv + D->rank;
  M->rows   = (re_t **)malloc(M->nrows * sizeof(re_t *));
  M->pos    = (ci_t **)malloc(M->nrows * sizeof(ci_t *));


  // reallocate memory for rows in M: they have a max length of
  // 1+D->ncols for rows from A+B, and
  // D->ncols for rows from D
#pragma omp parallel num_threads(nthrds)
{
  ri_t block_row_idx;
  ri_t line_idx;
  ci_t start_idx;
  ri_t min_range;
  ri_t min_range_blocks;
#pragma omp for private(i,j,k) schedule(dynamic)
  // each thread gets a whole block row
  for (l=0; l<rlB; ++l) {
    block_row_idx = l*__GBLA_SIMD_BLOCK_SIZE;
    // getting a row from A+B
    if (l<rlB-1) {
      min_range         = (l+1)*__GBLA_SIMD_BLOCK_SIZE;
      min_range_blocks  = __GBLA_SIMD_BLOCK_SIZE;
    } else {
      min_range         = map->npiv;
      min_range_blocks  = map->npiv - (rlB-1)*__GBLA_SIMD_BLOCK_SIZE;
    }

    // allocate memory and set first elements for the whole block
    for (i=block_row_idx; i<min_range; ++i) {
      //printf("%u || %u\n",i,D->ncols);
      M->rows[i]    = (re_t *)malloc((1+D->ncols) * sizeof(re_t));
      M->pos[i]     = (ci_t *)malloc((1+D->ncols) * sizeof(ci_t));
      M->rwidth[i]  = 0;

      // write 1 for A part
      M->rows[i][M->rwidth[i]]  = 1;
      M->pos[i][M->rwidth[i]]   = map->pc[i];
      M->rwidth[i]++;
    }

    // now write B part

    for (i=0; i<clB; ++i) {
      if (B->blocks[l][i].val ==  NULL)
        continue;
      start_idx = i * __GBLA_SIMD_BLOCK_SIZE;
      for ( j=0; j<min_range_blocks; ++j) {
        line_idx  = j * __GBLA_SIMD_BLOCK_SIZE;
        for ( k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          if (B->blocks[l][i].val[k+line_idx] != 0) {
            M->rows[block_row_idx+j][M->rwidth[block_row_idx+j]]  = B->blocks[l][i].val[k+line_idx];
            M->pos[block_row_idx+j][M->rwidth[block_row_idx+j]]   = map->npc_rev[map->npiv+k+start_idx];
            M->rwidth[block_row_idx+j]++;
          }
        }
        //printf("width[%u] = %u\n",block_row_idx+j,M->rwidth[block_row_idx+j]);
      }
    }
    // remove data for blocks of B
    for (j=0; j<clB; ++j) {
      free(B->blocks[l][j].val);
      B->blocks[l][j].val = NULL;
    }
    //printf("rwidth[%u] = %u\n",i,M->rwidth[i]);
    // relloc row size
    for (i=block_row_idx; i<min_range; ++i) {
      //printf("i %u (%u | %u)/ %u || %p %p ||| %u %u\n",i,M->rwidth[i],1+D->ncols,M->nrows, M->rows[i], M->pos[i], min_range, min_range_blocks);
      M->rows[i]  = realloc(M->rows[i], M->rwidth[i] * sizeof(re_t));
      M->pos[i]   = realloc(M->pos[i], M->rwidth[i] * sizeof(ci_t));
    }
    //printf("\n");
    // if we had the last line in the block, we can remove it
  }
}

#pragma omp parallel num_threads(nthrds)
{
  ri_t block_row_idx;
  ri_t line_idx;
  ci_t start_idx;
  ri_t min_range;
  ri_t min_range_blocks;
#pragma omp for private(i,j,k) schedule(dynamic)
  for (l=0; l<rlD; ++l) {
    block_row_idx = l*__GBLA_SIMD_BLOCK_SIZE + map->npiv;
    // getting rows from D
    if (l<rlD-1) {
      min_range = (l+1)*__GBLA_SIMD_BLOCK_SIZE + map->npiv;
    } else {
      min_range = M->nrows;
    }

    // allocate memory and set first elements for the whole block
    for (i=block_row_idx; i<min_range; ++i) {
      line_idx      = i-map->npiv;
      M->rows[i]    = (re_t *)malloc((D->ncols-D->row[line_idx]->piv_lead+1) * sizeof(re_t));
      M->pos[i]     = (ci_t *)malloc((D->ncols-D->row[line_idx]->piv_lead+1) * sizeof(ci_t));
      M->rwidth[i]  = 0;
    }

    for (i=block_row_idx; i<min_range; ++i) {
      line_idx  = i-map->npiv;
      if (D->row[line_idx]->piv_val != NULL) {
        for (j=D->row[line_idx]->piv_lead; j<D->ncols; ++j) {
          if (D->row[line_idx]->piv_val[j] != 0) {
            M->rows[i][M->rwidth[i]]  = D->row[line_idx]->piv_val[j];
            M->pos[i][M->rwidth[i]]   = map->npc_rev[map->npiv+j];
            M->rwidth[i]++;
          }
        }
        // relloc row size
        M->rows[i]  = realloc(M->rows[i], M->rwidth[i] * sizeof(re_t));
        M->pos[i]   = realloc(M->pos[i], M->rwidth[i] * sizeof(ci_t));
        // free row in D
        free(D->row[line_idx]->piv_val);
        free(D->row[line_idx]);
        D->row[line_idx] = NULL;
      }
    }
  }
}

  for (i=0; i<M->nrows; ++i)
    M->nnz  +=  M->rwidth[i];

  // free B and D
  free(B->blocks);
  free(B);
  B = NULL;
  free(D);
  D = NULL;

  // free map
  free(map->pc);
  free(map->npc);
  free(map->pc_rev);
  free(map->npc_rev);
  free(map->pri);
  free(map->npri);
  free(map);
  map = NULL;
}

void reconstruct_matrix_block(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sm_fl_ml_t *D,
    map_fl_t *map, const ci_t coldim, int free_matrices,
    int M_freed, int A_freed, int D_freed, int nthrds) {

  int ret;
  uint8_t *vec_free_AB, *vec_free_D;
  /*  1 pivot from A and otherwise at most D->ncols = B->ncols elements */
  /*  ==> no reallocations! */
  /* ci_t buffer = (ci_t) ceil((float)M->ncols / 10); */

  if (free_matrices == 1) {
    ret = posix_memalign((void *)&vec_free_AB, ALIGNT, (B->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    memset(vec_free_AB, 0, (B->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    ret = posix_memalign((void *)&vec_free_D, ALIGNT, (D->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    memset(vec_free_D, 0, (D->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
  }
  
  uint32_t i, j, k;
  uint32_t new_piv_to_row[coldim];
  uint32_t new_piv  = 0;
  for (i=0; i<coldim; ++i) {
    if (map->pri[i] == __GBLA_MINUS_ONE_32)
      continue;
    new_piv_to_row[i] = new_piv;
    new_piv++;
  }
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  M->nnz  = 0;
  ci_t block_row_idx, row_idx_block, block_start_idx, idx;
  re_t val;
  ri_t mline;
  uint32_t local_piv;
  uint16_t line;

#pragma omp parallel num_threads(1)
  {
#pragma omp for private(j,k) schedule(dynamic)
    for (i=0; i<coldim; ++i) {
      if (map->pri[i] == __GBLA_MINUS_ONE_32) {
        continue;
      }

      local_piv = new_piv_to_row[i];

      /*  clear old data */
      /* memset(*row_M, 0, tmp_width * sizeof(re_t)); */
      /* memset(*pos_M, 0, tmp_width * sizeof(ci_t)); */
      /*  clear data first */
      /* printf("realloc %d memory blocks\n",buffer); */
      M->rows[local_piv]    = realloc(M->rows[local_piv], (1 + D->ncols) * sizeof(re_t));
      M->pos[local_piv]     = realloc(M->pos[local_piv], (1 + D->ncols) * sizeof(ci_t));
      M->rwidth[local_piv]  = 0;

      if (map->pri[i] >= map->npiv) { /*  we are in D */
        mline = (map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE;
        line  = (map->pri[i]-map->npiv)%__GBLA_NROWS_MULTILINE;

				/* uint32_t idx_D; */
        if (D->ml[mline].sz > 0) {
          for (j=0; j<D->ml[mline].sz; ++j) {
            if (D->ml[mline].val[2*j+line] == 0)
              continue;

            M->rows[local_piv][M->rwidth[local_piv]]  = D->ml[mline].val[2*j+line];
            M->pos[local_piv][M->rwidth[local_piv]]   = map->npc_rev[j];
            M->rwidth[local_piv]++;
          }
          if (free_matrices == 1) {
            if (vec_free_D[(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE] == 0) {
              vec_free_D[(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE]  = 1;
            } else {
              /* printf("frees multiline %d\n",(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE); */
              free(D->ml[mline].val);
              D->ml[mline].val  = NULL;
            }
          }
        }
      } else { /*  we are in A resp. B */
        block_row_idx = (B->nrows - 1 - map->pri[i]) / B->bheight;
        row_idx_block = ((B->nrows - 1 - map->pri[i]) % B->bheight) / __GBLA_NROWS_MULTILINE;

        line  = (B->nrows - 1 - map->pri[i]) % __GBLA_NROWS_MULTILINE;

        /*  A should be always freed already! */
        if (A_freed == 0) {
        } else { /*  add 1 at pivot position for A part */
          M->rows[local_piv][M->rwidth[local_piv]]  = (re_t)1;
          M->pos[local_piv][M->rwidth[local_piv]]   = map->pc[i];
          M->rwidth[local_piv]++;
        }
        /*  add B part to row in M */
        for (j=0; j<clB; ++j) {
          if (B->blocks[block_row_idx][j] == NULL)
            continue;
          if (B->blocks[block_row_idx][j][row_idx_block].val == NULL)
            continue;
          if (B->blocks[block_row_idx][j][row_idx_block].sz == 0) {
            free(B->blocks[block_row_idx][j][row_idx_block].val);
            B->blocks[block_row_idx][j][row_idx_block].val  = NULL;
            free(B->blocks[block_row_idx][j][row_idx_block].idx);
            B->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
            continue;
          }

          block_start_idx = B->bwidth * j;

          if (B->blocks[block_row_idx][j][row_idx_block].dense == 0) {
            for (k=0; k<B->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
              idx = B->blocks[block_row_idx][j][row_idx_block].idx[k];
              val = B->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
              if (val != 0) {
                M->rows[local_piv][M->rwidth[local_piv]]  = val;
                M->pos[local_piv][M->rwidth[local_piv]]   = block_start_idx + idx;
                M->rwidth[local_piv]++;
              }
            }
          } else { /*  B block multiline is dense */
            for (k=0; k<B->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
              val = B->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
              if (val != 0) {
                M->rows[local_piv][M->rwidth[local_piv]]  = val;
                M->pos[local_piv][M->rwidth[local_piv]]   = block_start_idx + k;
                M->rwidth[local_piv]++;
              }
            }
          }
          if (free_matrices == 1) {
            if (vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 1) {
              if (B->blocks[block_row_idx][j][row_idx_block].dense == 0) {
                free (B->blocks[block_row_idx][j][row_idx_block].idx);
                B->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
              }
              free (B->blocks[block_row_idx][j][row_idx_block].val);
              B->blocks[block_row_idx][j][row_idx_block].val  = NULL;
            }
          }
        }
        if (free_matrices == 1) {
          if (vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 0) {
            vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE]  = 1;
          }
        }
      }
      /*  adjust memory, rows in M are done */
      M->nnz  = M->nnz + M->rwidth[local_piv];
      M->rows[local_piv] = realloc(M->rows[local_piv],
          M->rwidth[local_piv] * sizeof(re_t));
      M->pos[local_piv] = realloc(M->pos[local_piv],
          M->rwidth[local_piv] * sizeof(ci_t));
      /*  next pivot row */
#if __GBLA_BENCHMARKS
      free(M->rows[local_piv]);
      free(M->pos[local_piv]);
      M->rows[local_piv]  = NULL;
      M->pos[local_piv]   = NULL;
#endif
    }
  }
  if (free_matrices == 1) {
    free (vec_free_AB);
    free (vec_free_D);
    free (D->ml);
    free (D);
    D = NULL;
    for (i=0; i<rlB; ++i) {
      for (j=0; j<clB; ++j) {
        free(B->blocks[i][j]);
        B->blocks[i][j] = NULL;
      }
      free(B->blocks[i]);
      B->blocks[i]  = NULL;
    }
    free(B->blocks);
    free (B);
    B = NULL;
  }

  /*  free memory for map */
  free(map->pc);
  free(map->npc);
  free(map->pc_rev);
  free(map->npc_rev);
  free(map->pri);
  free(map->npri);
  free(map);
  map = NULL;
}

void reconstruct_matrix_ml(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, sm_fl_ml_t *D,
    map_fl_t *map, const ci_t coldim, int free_matrices,
    int M_freed, int A_freed, int D_freed, int nthrds) {

  int ret;
  uint8_t *vec_free_AB, *vec_free_D;
  /*  we need size of A->ml and otherwise at most D->ncols = B->ncols elements */
  /*  ==> no reallocations! */
  /* ci_t buffer = (ci_t) ceil((float)M->ncols / 10); */

  if (free_matrices == 1) {
    ret = posix_memalign((void *)&vec_free_AB, ALIGNT, (B->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    memset(vec_free_AB, 0, (B->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    ret = posix_memalign((void *)&vec_free_D, ALIGNT, (D->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
    memset(vec_free_D, 0, (D->nrows / __GBLA_NROWS_MULTILINE + 1)
        * sizeof(uint8_t));
  }

  uint32_t i, j, k;
  uint32_t new_piv_to_row[coldim];
  uint32_t new_piv  = 0;
  for (i=0; i<coldim; ++i) {
    if (map->pri[i] == __GBLA_MINUS_ONE_32)
      continue;
    new_piv_to_row[i] = new_piv;
    new_piv++;
  }
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  M->nnz  = 0;

#pragma omp parallel num_threads(1)
  {
    ml_t *row_A, *row_D;
		/* re_t **row_M; */
		/* ci_t **pos_M; */
		/* ci_t row_M_width; */
		/* ci_t tmp_width; */
    uint32_t local_piv;
    uint16_t line;
#pragma omp for private(i,j,k) schedule(dynamic)
    for (i=0; i<coldim; ++i) {
      if (map->pri[i] == __GBLA_MINUS_ONE_32)
        continue;

      local_piv = new_piv_to_row[i];

      if (map->pri[i] >= map->npiv) { /*  we are in D */
        /*  we need at most D->ncols = B->ncols elements */
        /*  ==> no reallocations! */
        M->rows[local_piv]    = realloc(M->rows[local_piv], D->ncols * sizeof(re_t));
        M->pos[local_piv]     = realloc(M->pos[local_piv], D->ncols * sizeof(ci_t));
        M->rwidth[local_piv]  = 0;

        row_D = &(D->ml[(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE]);
        line  = (map->pri[i]-map->npiv)%__GBLA_NROWS_MULTILINE;

				/* uint32_t idx_D; */
        if (row_D != NULL && row_D->sz > 0) {
          for (j=0; j<row_D->sz; ++j) {
            if (row_D->val[2*j+line] == 0)
              continue;

            M->rows[local_piv][M->rwidth[local_piv]]  = row_D->val[2*j+line];
            M->pos[local_piv][M->rwidth[local_piv]]   = map->npc_rev[j];
            /* printf("1 %d -- %u\n",map->npc_rev[j],row_D->val[2*j+line]); */
            M->rwidth[local_piv]++;
          }
          if (free_matrices == 1) {
            if (vec_free_D[(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE] == 0) {
              vec_free_D[(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE]  = 1;
            } else {
              /* printf("frees multiline %d\n",(map->pri[i]-map->npiv)/__GBLA_NROWS_MULTILINE); */
              free(row_D->val);
              row_D->val  = NULL;
              free(row_D->idx);
              row_D->idx  = NULL;
            }
          }
        }
      } else { /*  we are in A resp. B */
        /*  we need size of A->ml and otherwise at most D->ncols = B->ncols elements */
        /*  ==> no reallocations! */
        row_A = &(A->ml[(A->nrows - 1 - (map->pri[i]))/__GBLA_NROWS_MULTILINE]);
        M->rows[local_piv]    = realloc(M->rows[local_piv], (row_A->sz+B->ncols) * sizeof(re_t));
        M->pos[local_piv]     = realloc(M->pos[local_piv], (row_A->sz+B->ncols) * sizeof(ci_t));
        M->rwidth[local_piv]  = 0;
        ci_t block_row_idx, row_idx_block, block_start_idx, idx;
        re_t val;
        block_row_idx = (B->nrows - 1 - map->pri[i]) / B->bheight;
        row_idx_block = ((B->nrows - 1 - map->pri[i]) % B->bheight) / __GBLA_NROWS_MULTILINE;

        line  = (B->nrows - 1 - map->pri[i]) % __GBLA_NROWS_MULTILINE;

        if (row_A->dense == 0) {
          for (j=0; j<row_A->sz; ++j) {
            M->rows[local_piv][M->rwidth[local_piv]]  = row_A->val[2*j+line];
            M->pos[local_piv][M->rwidth[local_piv]]   = row_A->idx[j];
            M->rwidth[local_piv]++;
          }
        } else {
          for (j=0; j<row_A->sz; ++j) {
            M->rows[local_piv][M->rwidth[local_piv]]  = row_A->val[2*j+line];
            M->pos[local_piv][M->rwidth[local_piv]]   = j;
            M->rwidth[local_piv]++;
          }
        }
        /*  add B part to row in M */
        for (j=0; j<clB; ++j) {
          if (B->blocks[block_row_idx][j] == NULL)
            continue;
          if (B->blocks[block_row_idx][j][row_idx_block].val == NULL)
            continue;
          if (B->blocks[block_row_idx][j][row_idx_block].sz == 0) {
            free(B->blocks[block_row_idx][j][row_idx_block].val);
            B->blocks[block_row_idx][j][row_idx_block].val  = NULL;
            free(B->blocks[block_row_idx][j][row_idx_block].idx);
            B->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
            continue;
          }

          block_start_idx = B->bwidth * j;

          if (B->blocks[block_row_idx][j][row_idx_block].dense == 0) {
            for (k=0; k<B->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
              idx = B->blocks[block_row_idx][j][row_idx_block].idx[k];
              val = B->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
              if (val != 0) {
                M->rows[local_piv][M->rwidth[local_piv]]  = val;
                M->pos[local_piv][M->rwidth[local_piv]]   = block_start_idx + idx;
                M->rwidth[local_piv]++;
              }
            }
          } else { /*  B block multiline is dense */
            for (k=0; k<B->blocks[block_row_idx][j][row_idx_block].sz; ++k) {
              val = B->blocks[block_row_idx][j][row_idx_block].val[2*k+line];
              if (val != 0) {
                M->rows[local_piv][M->rwidth[local_piv]]  = val;
                M->pos[local_piv][M->rwidth[local_piv]]   = block_start_idx + k;
                M->rwidth[local_piv]++;
              }
            }
          }
          if (free_matrices == 1) {
            if (vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 1) {
              if (B->blocks[block_row_idx][j][row_idx_block].dense == 0) {
                free (B->blocks[block_row_idx][j][row_idx_block].idx);
                B->blocks[block_row_idx][j][row_idx_block].idx  = NULL;
              }
              free (B->blocks[block_row_idx][j][row_idx_block].val);
              B->blocks[block_row_idx][j][row_idx_block].val  = NULL;
            }
          }
        }
        if (free_matrices == 1) {
          if (vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 1) {
            free(row_A->val);
            row_A->val  = NULL;
            free(row_A->idx);
            row_A->idx  = NULL;
          }
          if (vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE] == 0) {
            vec_free_AB[(B->nrows - 1 - map->pri[i])/__GBLA_NROWS_MULTILINE]  = 1;
          }
        }
      }
      /*  adjust memory, rows in M are done */
      /* printf("final row width: %d for %d\n",row_M_width, local_piv); */
      M->rows[local_piv]  = (re_t*) realloc(M->rows[local_piv],
          M->rwidth[local_piv] * sizeof(re_t));
      M->pos[local_piv]   = (ci_t*) realloc(M->pos[local_piv],
          M->rwidth[local_piv] * sizeof(ci_t));
      M->nnz  = M->nnz + M->rwidth[local_piv];
      /*  next pivot row */
#if __GBLA_BENCHMARKS
      free(M->rows[local_piv]);
      free(M->pos[local_piv]);
      M->rows[local_piv]  = NULL;
      M->pos[local_piv]   = NULL;
#endif
    }
  }
  if (free_matrices == 1) {
    free (vec_free_AB);
    free (vec_free_D);
    free (D->ml);
    free (D);
    D = NULL;
    free (A->ml);
    free (A);
    A = NULL;
    for (i=0; i<rlB; ++i) {
      for (j=0; j<clB; ++j) {
        free(B->blocks[i][j]);
        B->blocks[i][j] = NULL;
      }
      free(B->blocks[i]);
      B->blocks[i]  = NULL;
    }
    free(B->blocks);
    free (B);
    B = NULL;
  }


  /*  free memory for map */
  free(map->pc);
  free(map->npc);
  free(map->pc_rev);
  free(map->npc_rev);
  free(map->pri);
  free(map->npri);
  free(map);
  map = NULL;
}

void splice_fl_matrix_reduced(sbm_fl_t *B, sm_fl_ml_t *D, sbm_fl_t *B1, sbm_fl_t *B2,
    sbm_fl_t *D1, sbm_fl_t *D2, map_fl_t *map, ri_t complete_nrows, ci_t complete_ncols,
    int block_dim, int rows_multiline, int nthrds, int destruct_input_matrix,
    int verbose) {

  int bdim  = block_dim / rows_multiline;
  /*  construct index map for Faugère-Lachartre decomposition of matrix M */
  /* construct_fl_map(M, map); */

  int  l/*, m*/;
	uint32_t i,j,k ;

  /*  initialize meta data for block submatrices */
  B1->nrows   = B->nrows;   /*  row dimension */
  B1->ncols   = map->npiv;  /*  col dimension */
  B1->bheight = bdim;       /*  block height */
  B1->bwidth  = bdim;       /*  block width */
  B1->ba      = dtrl;       /*  block alignment */
  B1->fe      = 0;          /*  fill empty blocks? */
  B1->hr      = 0;          /*  allow hybrid rows? */
  B1->nnz     = 0;          /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlB1 = (uint32_t) ceil((float)B1->nrows / B1->bheight);
  const uint32_t clB1 = (uint32_t) ceil((float)B1->ncols / B1->bwidth);

  B1->blocks = (mbl_t ***)malloc(rlB1 * sizeof(mbl_t **));
  for (i=0; i<rlB1; ++i) {
    B1->blocks[i] = (mbl_t **)malloc(clB1 * sizeof(mbl_t *));
    for (j=0; j<clB1; ++j) {
      B1->blocks[i][j] = (mbl_t *)malloc(
          (B1->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(B1->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        B1->blocks[i][j][k].val = NULL;
        B1->blocks[i][j][k].idx = NULL;
        B1->blocks[i][j][k].sz  = B1->blocks[i][j][k].dense = 0;
      }
    }
  }

  B2->nrows    = B->nrows;              /*  row dimension */
  B2->ncols    = B->ncols - map->npiv;  /*  col dimension */
  B2->bheight  = bdim;                  /*  block height */
  B2->bwidth   = bdim;                  /*  block width */
  B2->ba       = dtlr;                  /*  block alignment */
  B2->fe       = 1;                     /*  fill empty blocks? */
  B2->hr       = 1;                     /*  allow hybrid rows? */
  B2->nnz      = 0;                     /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlB2  = (uint32_t) ceil((float)B2->nrows / B2->bheight);
  const uint32_t clB2  = (uint32_t) ceil((float)B2->ncols / B2->bwidth);

  B2->blocks = (mbl_t ***)malloc(rlB2 * sizeof(mbl_t **));
  for (i=0; i<rlB2; ++i) {
    B2->blocks[i] = (mbl_t **)malloc(clB2 * sizeof(mbl_t *));
    for (j=0; j<clB2; ++j) {
      B2->blocks[i][j] = (mbl_t *)malloc(
          (B2->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(B2->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        B2->blocks[i][j][k].val = NULL;
        B2->blocks[i][j][k].idx = NULL;
        B2->blocks[i][j][k].sz  = B2->blocks[i][j][k].dense = 0;
      }
    }
  }

  /*  initialize meta data for block submatrices */
  D1->nrows   = map->npiv;            /*  row dimension */
  D1->ncols   = D->ncols - map->npiv; /*  col dimension */
  D1->bheight = bdim;                 /*  block height */
  D1->bwidth  = bdim;                 /*  block width */
  D1->ba      = dtrl;                 /*  block alignment */
  D1->fe      = 1;                    /*  fill empty blocks? */
  D1->hr      = 1;                    /*  allow hybrid rows? */
  D1->nnz     = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlD1 = (uint32_t) ceil((float)D1->nrows / D1->bheight);
  const uint32_t clD1 = (uint32_t) ceil((float)D1->ncols / D1->bwidth);

  D1->blocks = (mbl_t ***)malloc(rlD1 * sizeof(mbl_t **));
  for (i=0; i<rlD1; ++i) {
    D1->blocks[i] = (mbl_t **)malloc(clD1 * sizeof(mbl_t *));
    for (j=0; j<clD1; ++j) {
      D1->blocks[i][j] = (mbl_t *)malloc(
          (D1->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(D1->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        D1->blocks[i][j][k].val = NULL;
        D1->blocks[i][j][k].idx = NULL;
        D1->blocks[i][j][k].sz  = D1->blocks[i][j][k].dense = 0;
      }
    }
  }

  D2->nrows    = map->npiv;            /*  row dimension */
  D2->ncols    = complete_ncols - map->npiv; /*  col dimension */
  D2->bheight  = bdim;                 /*  block height */
  D2->bwidth   = bdim;                 /*  block width */
  D2->ba       = dtlr;                 /*  block alignment */
  D2->fe       = 1;                    /*  fill empty blocks? */
  D2->hr       = 1;                    /*  allow hybrid rows? */
  D2->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlD2  = (uint32_t) ceil((float)D2->nrows / D2->bheight);
  const uint32_t clD2  = (uint32_t) ceil((float)D2->ncols / D2->bwidth);

  D2->blocks = (mbl_t ***)malloc(rlD2 * sizeof(mbl_t **));
  for (i=0; i<rlD2; ++i) {
    D2->blocks[i] = (mbl_t **)malloc(clD2 * sizeof(mbl_t *));
    for (j=0; j<clD2; ++j) {
      D2->blocks[i][j] = (mbl_t *)malloc(
          (D2->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(D2->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        D2->blocks[i][j][k].val = NULL;
        D2->blocks[i][j][k].idx = NULL;
        D2->blocks[i][j][k].sz  = D2->blocks[i][j][k].dense = 0;
      }
    }
  }

#if __GBLA_DEBUG
  printf("B1[%dx%d]\n",B1->nrows,B1->ncols);
  printf("B2[%dx%d]\n",B2->nrows,B2->ncols);
  printf("D1[%dx%d]\n",D1->nrows,D1->ncols);
  printf("D2[%dx%d]\n",D2->nrows,D2->ncols);
#endif

  /*  copy B block-row-wise to B1 and B2 */

  /*  allocate memory for a sparse matrix of B->bheight rows and move B's entries */
  /*  to B1 and B2 correspondingly */
  sm_t *tmp_B = (sm_t *)malloc(sizeof(sm_t));
  tmp_B->rows   = (re_t **)malloc(B->bheight * sizeof(re_t *));
  tmp_B->pos    = (ci_t **)malloc(B->bheight * sizeof(ci_t *));
  tmp_B->rwidth = (ci_t *)malloc(B->bheight * sizeof(ci_t));
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i) nowait
    for (i=0; i<B->bheight; ++i) {
      tmp_B->rows[i]  = (re_t *)malloc(B->ncols * sizeof(re_t));
      tmp_B->pos[i]   = (ci_t *)malloc(B->ncols * sizeof(ci_t));
    }
  }
	/* bi_t rbi  = 0;          |+  row block index +| */
	/* bi_t cvb  = 0;          |+  current vector in block +| */
	/* ri_t rihb[B1->bheight]; |+  row indices horizontal blocks +| */

	assert((int)B->nrows>=0);
	assert((int)B->bheight>=0);
  for (l=(int)B->nrows - 1; l>=0; --l) {
    /*  clear tmp_B */
    memset(tmp_B->rwidth, 0, B->bheight * sizeof(ci_t));
    for (l=0; l<(int)B->bheight; ++l) {
      memset(tmp_B->rows[l], 0, tmp_B->ncols * sizeof(re_t));
      memset(tmp_B->pos[l], 0, B->ncols * sizeof(ci_t));
    }
    /*  block row index */
		/* ri_t bri  = (B->nrows - 1 - l) / B->bheight; */
    /*  row index in block */
		/* ri_t rib  = ((B->nrows - 1 - l) % B->bheight) / __GBLA_NROWS_MULTILINE; */
    /*  block start index */
		/* ri_t bsi; */
		/* re_t val; */
		/* ri_t idx; */
		/* bi_t line = (B->nrows - 1 - l) % __GBLA_NROWS_MULTILINE; */

  }
    /* for (j=0; j< */
}

void splice_fl_matrix_dense(sm_t *M, dbm_fl_t *A, dbm_fl_t *B, dbm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads) {

  // construct index map for Faugère-Lachartre decomposition of matrix M
  if (map_defined == 0)
    construct_fl_map(M, map);

  // initialize submatrices
  init_dbm(A, map->npiv, map->npiv);
  init_dbm(B, map->npiv, M->ncols - map->npiv);
  init_dbm(C, M->nrows - map->npiv, map->npiv);
  init_dbm(D, M->nrows - map->npiv, M->ncols - map->npiv);

#if GBLA_DEBUG
  check_dimensions(M, A, B, C, D);
#endif
  ri_t npiv = 0; // number pivots handled
  ri_t *piv_start_idx;

  // write data in A and B, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->pri, M->ncols, A->nrows);
  fill_dense_submatrices_diagonalize(M, A, B, map, map->pri, piv_start_idx, destruct_input_matrix,
      nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
  // write data in C and D, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->npri, M->nrows, C->nrows);
  fill_dense_submatrices(M, C, D, map, map->npri, piv_start_idx, destruct_input_matrix, nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
}

void splice_fl_matrix_hybrid(sm_t *M, hbm_fl_t *A, hbm_fl_t *B, hbm_fl_t *C,
    hbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads) {

  // construct index map for Faugère-Lachartre decomposition of matrix M
  if (map_defined == 0)
    construct_fl_map(M, map);

  // initialize submatrices
  init_hbm(A, map->npiv, map->npiv);
  init_hbm(B, map->npiv, M->ncols - map->npiv);
  init_hbm(C, M->nrows - map->npiv, map->npiv);
  init_hbm(D, M->nrows - map->npiv, M->ncols - map->npiv);

#if GBLA_DEBUG
  check_dimensions_hbm(M, A, B, C, D);
#endif
  ri_t npiv = 0; // number pivots handled
  ri_t *piv_start_idx;

  // write data in A and B, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->pri, M->ncols, A->nrows);
  fill_hybrid_submatrices(M, A, B, map, map->pri, piv_start_idx, destruct_input_matrix,
      nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
  // write data in C and D, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->npri, M->nrows, C->nrows);
  fill_hybrid_submatrices(M, C, D, map, map->npri, piv_start_idx, destruct_input_matrix, nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
}

void splice_fl_matrix_hybrid_dense(sm_t *M, hbm_fl_t *A, dbm_fl_t *B, hbm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads) {

  // construct index map for Faugère-Lachartre decomposition of matrix M
  if (map_defined == 0)
    construct_fl_map(M, map);

  // initialize submatrices
  init_hbm(A, map->npiv, map->npiv);
  init_dbm(B, map->npiv, M->ncols - map->npiv);
  init_hbm(C, M->nrows - map->npiv, map->npiv);
  init_dbm(D, M->nrows - map->npiv, M->ncols - map->npiv);

#if GBLA_DEBUG
  //check_dimensions_hbm(M, A, B, C, D);
#endif
  ri_t npiv = 0; // number pivots handled
  ri_t *piv_start_idx;

  // write data in A and B, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->pri, M->ncols, A->nrows);
  fill_hybrid_dense_submatrices(M, A, B, map, map->pri, piv_start_idx, destruct_input_matrix,
      nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
  // write data in C and D, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->npri, M->nrows, C->nrows);
  fill_hybrid_dense_submatrices(M, C, D, map, map->npri, piv_start_idx, destruct_input_matrix, nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
}

void splice_fl_matrix_sparse_dense_keep_A(sm_t *M, sm_fl_t *A, sb_fl_t *B, sm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads) {

  // construct index map for Faugère-Lachartre decomposition of matrix M
  if (map_defined == 0)
    construct_fl_map(M, map);

  // initialize submatrices
  init_sm(A, map->npiv, map->npiv);
  init_sb(B, map->npiv, M->ncols - map->npiv);
  init_sm(C, M->nrows - map->npiv, map->npiv);
  init_dbm(D, M->nrows - map->npiv, M->ncols - map->npiv);

#if GBLA_DEBUG
  //check_dimensions_hbm(M, A, B, C, D);
#endif
  ri_t npiv = 0; // number pivots handled
  ri_t *piv_start_idx;

  // write data in A and B, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->pri, M->ncols, A->nrows);
  fill_sparse_sparse_submatrices_keep_A(M, A, B, map, map->pri, piv_start_idx,
      destruct_input_matrix, nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
  // write data in C and D, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->npri, M->nrows, C->nrows);
  fill_sparse_dense_submatrices_keep_A(M, C, D, map, map->npri, piv_start_idx,
      destruct_input_matrix, nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
}

void splice_fl_matrix_sparse_dense_2(sm_t *M, sb_fl_t *A, dbm_fl_t *B, sb_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads) {

  // construct index map for Faugère-Lachartre decomposition of matrix M
  if (map_defined == 0)
    construct_fl_map(M, map);

  // initialize submatrices
  init_sb(A, map->npiv, map->npiv);
  init_dbm(B, map->npiv, M->ncols - map->npiv);
  init_sb(C, M->nrows - map->npiv, map->npiv);
  init_dbm(D, M->nrows - map->npiv, M->ncols - map->npiv);

#if GBLA_DEBUG
  //check_dimensions_hbm(M, A, B, C, D);
#endif
  ri_t npiv = 0; // number pivots handled
  ri_t *piv_start_idx;

  // write data in A and B, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->pri, M->ncols, A->nrows);
  fill_sparse_dense_submatrices(M, A, B, map, map->pri, piv_start_idx, destruct_input_matrix,
      nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
  // write data in C and D, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->npri, M->nrows, C->nrows);
  fill_sparse_dense_submatrices(M, C, D, map, map->npri, piv_start_idx, destruct_input_matrix,
      nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
}

void splice_fl_matrix_sparse_dense(sm_t *M, sb_fl_t *A, dbm_fl_t *B, dbm_fl_t *C,
    dbm_fl_t *D, map_fl_t *map, const int map_defined,
    const int destruct_input_matrix, const int verbose, const int nthreads) {

  // construct index map for Faugère-Lachartre decomposition of matrix M
  if (map_defined == 0)
    construct_fl_map(M, map);

  // initialize submatrices
  init_sb(A, map->npiv, map->npiv);
  init_dbm(B, map->npiv, M->ncols - map->npiv);
  init_dbm(C, M->nrows - map->npiv, map->npiv);
  init_dbm(D, M->nrows - map->npiv, M->ncols - map->npiv);

#if GBLA_DEBUG
  //check_dimensions_hbm(M, A, B, C, D);
#endif
  ri_t npiv = 0; // number pivots handled
  ri_t *piv_start_idx;

  // write data in A and B, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->pri, M->ncols, A->nrows);
  fill_sparse_dense_submatrices(M, A, B, map, map->pri, piv_start_idx, destruct_input_matrix,
      nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
  // write data in C and D, i.e. fill upper part of splicing
  init_pivot_block_start_indices(&piv_start_idx, &npiv, map->npri, M->nrows, C->nrows);
  fill_dense_submatrices(M, C, D, map, map->npri, piv_start_idx, destruct_input_matrix, nthreads);

  free(piv_start_idx);
  piv_start_idx = NULL;
}


void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
                    map_fl_t *map, ri_t complete_nrows, ci_t complete_ncols,
                    int block_dim, int rows_multiline,
                    int nthreads, int destruct_input_matrix, int verbose,
                    int map_defined) {

  int bdim  = block_dim / rows_multiline;
  /*  construct index map for Faugère-Lachartre decomposition of matrix M */
  if (map_defined == 0)
    construct_fl_map(M, map);

  int l/*, m*/;
	uint32_t i,j,k;

  /*  initialize meta data for block submatrices */
  A->nrows    = map->npiv;            /*  row dimension */
  A->ncols    = map->npiv;            /*  col dimension */
  A->bheight  = bdim;                 /*  block height */
  A->bwidth   = bdim;                 /*  block width */
  A->ba       = dtrl;                 /*  block alignment */
  A->fe       = 0;                    /*  fill empty blocks? */
  A->hr       = 0;                    /*  allow hybrid rows? */
  A->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / A->bheight);
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);

  A->blocks = (mbl_t ***)malloc(rlA * sizeof(mbl_t **));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (mbl_t **)malloc(clA * sizeof(mbl_t *));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j] = (mbl_t *)malloc(
          (A->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(A->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        A->blocks[i][j][k].val  = NULL;
        A->blocks[i][j][k].idx  = NULL;
        A->blocks[i][j][k].sz   = A->blocks[i][j][k].dense  = 0;
      }
    }
  }

  B->nrows    = map->npiv;            /*  row dimension */
  B->ncols    = complete_ncols - map->npiv; /*  col dimension */
  B->bheight  = bdim;                 /*  block height */
  B->bwidth   = bdim;                 /*  block width */
  B->ba       = dtlr;                 /*  block alignment */
  B->fe       = 1;                    /*  fill empty blocks? */
  B->hr       = 1;                    /*  allow hybrid rows? */
  B->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  B->blocks = (mbl_t ***)malloc(rlB * sizeof(mbl_t **));
  for (i=0; i<rlB; ++i) {
    B->blocks[i]  = (mbl_t **)malloc(clB * sizeof(mbl_t *));
    for (j=0; j<clB; ++j) {
      B->blocks[i][j] = (mbl_t *)malloc(
          (B->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(B->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        B->blocks[i][j][k].val  = NULL;
        B->blocks[i][j][k].idx  = NULL;
        B->blocks[i][j][k].sz   = B->blocks[i][j][k].dense  = 0;
      }
    }
  }

  C->nrows    = complete_nrows - map->npiv; /*  row dimension */
  C->ncols    = map->npiv;            /*  col dimension */
  C->bheight  = bdim;                 /*  block height */
  C->bwidth   = bdim;                 /*  block width */
  C->ba       = dtrl;                 /*  block alignment */
  C->fe       = 0;                    /*  fill empty blocks? */
  C->hr       = 0;                    /*  allow hybrid rows? */
  C->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlC  = (uint32_t) ceil((float)C->nrows / C->bheight);
  const uint32_t clC  = (uint32_t) ceil((float)C->ncols / C->bwidth);

  C->blocks = (mbl_t ***)malloc(rlC * sizeof(mbl_t **));
  for (i=0; i<rlC; ++i) {
    C->blocks[i]  = (mbl_t **)malloc(clC * sizeof(mbl_t *));
    for (j=0; j<clC; ++j) {
      C->blocks[i][j] = (mbl_t *)malloc(
          (C->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(C->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        C->blocks[i][j][k].val  = NULL;
        C->blocks[i][j][k].idx  = NULL;
        C->blocks[i][j][k].sz   = C->blocks[i][j][k].dense  = 0;
      }
    }
  }

  D->nrows    = complete_nrows - map->npiv; /*  row dimension */
  D->ncols    = complete_ncols - map->npiv; /*  col dimension */
  D->bheight  = bdim;                 /*  block height */
  D->bwidth   = bdim;                 /*  block width */
  D->ba       = dtlr;                 /*  block alignment */
  D->fe       = 1;                    /*  fill empty blocks? */
  D->hr       = 1;                    /*  allow hybrid rows? */
  D->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);

  D->blocks = (mbl_t ***)malloc(rlD * sizeof(mbl_t **));
  for (i=0; i<rlD; ++i) {
    D->blocks[i]  = (mbl_t **)malloc(clD * sizeof(mbl_t *));
    for (j=0; j<clD; ++j) {
      D->blocks[i][j] = (mbl_t *)malloc(
          (D->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(D->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        D->blocks[i][j][k].val  = NULL;
        D->blocks[i][j][k].idx  = NULL;
        D->blocks[i][j][k].sz   = D->blocks[i][j][k].dense  = 0;
      }
    }
  }

#if __GBLA_DEBUG
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

  uint32_t npiv = 0; /*  number pivots handled */
  /*  Take maximal number of rows in blocks to be able to use array */
  /*  piv_start_idx for construction of A & B as well as for the construction of */
  /*  C & D. */
  uint32_t max_nrows =  (A->nrows > C->nrows) ? A->nrows : C->nrows;
  uint32_t piv_start_idx[(max_nrows / B->bheight) + 2];
  uint32_t block_idx/*, block_idx_2*/;

  /*  find blocks for construction of A & B */
	assert((int)M->ncols>=0);
  for (l = (int)M->ncols-1; l > -1; --l) {
    if (map->pri[l] != __GBLA_MINUS_ONE_32) {
      npiv++;
    }
    if (npiv % B->bheight == 0) {
      piv_start_idx[npiv/B->bheight]  = l;
    }
  }
  /*  loop might overwrite piv_start_idx[0] with a wrong index; */
  /*  instead of checking "npiv > 0" in each if clause we just reset */
  /*  piv_start_idx[0] after the for loop */
  piv_start_idx[0]  = M->ncols;

#if __GBLA_NEW_SPLICER
  /*
   * For code in __GBLA_NEW_SPLICER
   * Trying to optimize splicing: We know that the number of elements in the
   * lefthand side are getting fewer and fewer since it is an upper triangular
   * matrix (at least for the upper left part). So we try to give each thread 2
   * tasks at a time, namely one from the top (more stuff to do) and one from
   * the bottom (less stuff to do). Thus we hope that the workload is better
   * balanced.
   */
  uint32_t rihb[B->bheight];  /*  rows indices horizontal block */
  uint16_t cvb  = 0;          /*  current vector in block */
  uint32_t nb   = npiv/B->bheight;
  if (npiv % B->bheight != 0)
    nb++;

  /*  set leftout entries to zero */
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  if (nb % 2) {
    nb--;
    for (i = ((int)piv_start_idx[nb]-1);
        i > (int)piv_start_idx[nb+1]-1; --i) {
      if (map->pri[i] != __GBLA_MINUS_ONE_32) {
        rihb[cvb] = map->pri[i];
        cvb++;
      }
      if (cvb == B->bheight || i == 0) {
        /*  multiline version for A */
        write_blocks_lr_matrix(M, A, B, map, rihb, cvb, nb, 0);

        /*  TODO: Destruct input matrix on the go */
        if (destruct_input_matrix) {
          for (j=0; j<cvb; ++j) {
            free(M->pos[rihb[j]]);
            M->pos[rihb[j]] = NULL;
            free(M->rows[rihb[j]]);
            M->rows[rihb[j]] = NULL;
          }
        /*  free memory */
        }
        cvb = 0;
      }
    }
  }
  omp_set_dynamic(0);
#pragma omp parallel private(cvb, rihb, block_idx, block_idx_2, i, j) num_threads(nthreads)
  {
    cvb  = 0;
#if __GBLA_DEBUG
    printf("nthreads %d\n",nthreads);
    printf("npiv %d -- bheight %d\n",map->npiv,B->bheight);
    printf("div %d\n", map->npiv/B->bheight);
#endif

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx < nb/2; ++block_idx) {
      /*  construct block submatrices A & B */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
      block_idx_2 = nb - block_idx - 1;
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (map->pri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        if (cvb == B->bheight || i == 0) {
          /*  multiline version for A */
          write_blocks_lr_matrix(M, A, B, map, rihb, cvb, block_idx, 0);

          /*  TODO: Destruct input matrix on the go */
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          /*  free memory */
          }
          cvb = 0;
        }
      }
      for (i = ((int)piv_start_idx[block_idx_2]-1);
          i > (int)piv_start_idx[block_idx_2+1]-1; --i) {
        if (map->pri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        if (cvb == B->bheight || i == 0) {
          /*  multiline version for A */
          write_blocks_lr_matrix(M, A, B, map, rihb, cvb, block_idx_2, 0);

          /*  TODO: Destruct input matrix on the go */
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          /*  free memory */
          }
          cvb = 0;
        }
      }
    }
  }

  /*  find blocks for construction of C & D */
  npiv  = 0;
  for (i = (int)M->ncols-1; i > -1; --i) {
    if (map->npri[i] != __GBLA_MINUS_ONE_32)
      npiv++;

    if (npiv % D->bheight == 0) {
      piv_start_idx[npiv/D->bheight]  = i;
    }
  }
  /*  loop might overwrite piv_start_idx[0] with a wrong index; */
  /*  instead of checking "npiv > 0" in each if clause we just reset */
  /*  piv_start_idx[0] after the for loop */
  piv_start_idx[0]  = M->ncols;
  /*  set leftout entries to zero */
  for (i=npiv/D->bheight+1; i < (max_nrows / D->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  nb  = npiv/D->bheight;
  if (npiv % D->bheight != 0)
    nb++;
  if (nb % 2) {
    nb--;
    for (i = ((int)piv_start_idx[nb]-1);
        i > (int)piv_start_idx[nb+1]-1; --i) {
      if (map->npri[i] != __GBLA_MINUS_ONE_32) {
        rihb[cvb] = map->npri[i];
        cvb++;
      }
      if (cvb == D->bheight || i == 0) {
        /*  multiline version for C */
        write_blocks_lr_matrix(M, C, D, map, rihb, cvb, nb, 1);

        if (destruct_input_matrix) {
          for (j=0; j<cvb; ++j) {
            free(M->pos[rihb[j]]);
            M->pos[rihb[j]] = NULL;
            free(M->rows[rihb[j]]);
            M->rows[rihb[j]] = NULL;
          }
        }
        cvb = 0;
      }
    }
  }
  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, block_idx_2, i, j) num_threads(nthreads)
  {
    cvb  = 0;
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx < nb/2; ++block_idx) {
      /*  construct block submatrices C & D */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
      block_idx_2 = nb - block_idx - 1;
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (map->npri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          /*  multiline version for C */
          write_blocks_lr_matrix(M, C, D, map, rihb, cvb, block_idx, 1);

          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          }
          cvb = 0;
        }
      }
      for (i = ((int)piv_start_idx[block_idx_2]-1);
          i > (int)piv_start_idx[block_idx_2+1]-1; --i) {
        if (map->npri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          /*  multiline version for C */
          write_blocks_lr_matrix(M, C, D, map, rihb, cvb, block_idx_2, 1);

          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          }
          cvb = 0;
        }
      }
    }
  }
#else

  /*  set leftout entries to zero */
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i, j, l) num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  /*  rows indices horizontal block */
    uint16_t cvb  = 0;          /*  current vector in block */

#if __GBLA_DEBUG
    printf("nthreads %d\n",nthreads);
    printf("npiv %d -- bheight %d\n",map->npiv,B->bheight);
    printf("div %d\n", map->npiv/B->bheight);
#endif
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= npiv/A->bheight; ++block_idx) {
#if __GBLA_DEBUG
      printf("bi %d\n", block_idx);
      printf("piv_idx[block] %d\n",piv_start_idx[block_idx]);
      printf("piv_idx[block+1] %d\n",piv_start_idx[block_idx+1]);
#endif
      /*  construct block submatrices A & B */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
			assert((int)piv_start_idx[block_idx]>=0);
      for (l = ((int)piv_start_idx[block_idx]-1);
          l > (int)piv_start_idx[block_idx+1]-1; --l) {
        if (map->pri[l] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->pri[l];
          cvb++;
        }
        if (cvb == B->bheight || l == 0) {
          write_blocks_lr_matrix(M, A, B, map, rihb, cvb, block_idx, 0);

          /*  TODO: Destruct input matrix on the go */
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          /*  free memory */
          }
          cvb = 0;
        }
      }
    }
  }

  /*  find blocks for construction of C & D */
  npiv  = 0;
	assert((int)M->nrows>=0);
  for (l = (int)M->nrows-1; l > -1; --l) {
    if (map->npri[l] != __GBLA_MINUS_ONE_32)
      npiv++;

    if (npiv % B->bheight == 0) {
      piv_start_idx[npiv/B->bheight]  = l;
    }
  }
  /*  loop might overwrite piv_start_idx[0] with a wrong index; */
  /*  instead of checking "npiv > 0" in each if clause we just reset */
  /*  piv_start_idx[0] after the for loop */
  piv_start_idx[0]  = M->nrows;

  /*  set leftout entries to zero */
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i, j, l) num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  /*  rows indices horizontal block */
    uint16_t cvb  = 0;          /*  current vector in block */
    int ctr = 0;
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= npiv/C->bheight; ++block_idx) {
      /*  construct block submatrices C & D */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
			assert((int)piv_start_idx[block_idx]>=0);
      for (l = ((int)piv_start_idx[block_idx]-1);
          l > (int)piv_start_idx[block_idx+1]-1; --l) {
        if (map->npri[l] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->npri[l];
          cvb++;
        }
        if (cvb == D->bheight || l == 0) {
          ctr++;
          write_blocks_lr_matrix(M, C, D, map, rihb, cvb, block_idx, 1);

          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          }
          cvb = 0;
        }
      }
    }
  }
#endif
}


void splice_fl_matrix_ml_A_C(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, sm_fl_ml_t *C,
                    sbm_fl_t *D, map_fl_t *map, int block_dim, int rows_multiline,
                    int nthreads, int destruct_input_matrix, int verbose) {

  int bdim  = block_dim / rows_multiline;
  /*  construct index map for Faugère-Lachartre decomposition of matrix M */
  construct_fl_map(M, map);

	int  l/*, m*/;
	uint32_t i,j,k;

  /*  initialize meta data for multiline sub matrices A and C */
  A->nrows    = map->npiv;            /*  row dimension */
  A->ncols    = map->npiv;            /*  col dimension */
  A->ba       = dtrl;                 /*  block alignment */
  A->fe       = 0;                    /*  fill empty blocks? */
  A->hr       = 0;                    /*  allow hybrid rows? */
  A->nnz      = 0;                    /*  number nonzero elements */

  ri_t rlA  = A->nrows / 2;
  if (A->nrows % 2)
    rlA++;

  A->ml = (ml_t *)malloc(rlA * sizeof(ml_t));
  for (i=0; i<rlA; ++i) {
    A->ml[i].val    = NULL;
    A->ml[i].idx    = NULL;
    A->ml[i].sz     = 0;
    A->ml[i].dense  = 0;
  }


  C->nrows    = M->nrows - map->npiv; /*  row dimension */
  C->ncols    = map->npiv;            /*  col dimension */
  C->ba       = dtrl;                 /*  block alignment */
  C->fe       = 0;                    /*  fill empty blocks? */
  C->hr       = 0;                    /*  allow hybrid rows? */
  C->nnz      = 0;                    /*  number nonzero elements */

  ri_t rlC  = C->nrows / 2;
  if (C->nrows % 2)
    rlC++;

  C->ml = (ml_t *)malloc(rlC * sizeof(ml_t));
  for (i=0; i<rlC; ++i) {
    C->ml[i].val    = NULL;
    C->ml[i].idx    = NULL;
    C->ml[i].sz     = 0;
    C->ml[i].dense  = 0;
  }


  /*  initialize meta data for block submatrices */
  B->nrows    = map->npiv;            /*  row dimension */
  B->ncols    = M->ncols - map->npiv; /*  col dimension */
  B->bheight  = bdim;                 /*  block height */
  B->bwidth   = bdim;                 /*  block width */
  B->ba       = dtlr;                 /*  block alignment */
  B->fe       = 1;                    /*  fill empty blocks? */
  B->hr       = 1;                    /*  allow hybrid rows? */
  B->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlB  = (uint32_t) ceil((float)B->nrows / B->bheight);
  const uint32_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  B->blocks = (mbl_t ***)malloc(rlB * sizeof(mbl_t **));
  for (i=0; i<rlB; ++i) {
    B->blocks[i]  = (mbl_t **)malloc(clB * sizeof(mbl_t *));
    for (j=0; j<clB; ++j) {
      B->blocks[i][j] = (mbl_t *)malloc(
          (B->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(B->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        B->blocks[i][j][k].val  = NULL;
        B->blocks[i][j][k].idx  = NULL;
        B->blocks[i][j][k].sz   = B->blocks[i][j][k].dense  = 0;
      }
    }
  }

  D->nrows    = M->nrows - map->npiv; /*  row dimension */
  D->ncols    = M->ncols - map->npiv; /*  col dimension */
  D->bheight  = bdim;                 /*  block height */
  D->bwidth   = bdim;                 /*  block width */
  D->ba       = dtlr;                 /*  block alignment */
  D->fe       = 1;                    /*  fill empty blocks? */
  D->hr       = 1;                    /*  allow hybrid rows? */
  D->nnz      = 0;                    /*  number nonzero elements */

  /*  allocate memory for blocks */

  /*  row and column loops */
  const uint32_t rlD  = (uint32_t) ceil((float)D->nrows / D->bheight);
  const uint32_t clD  = (uint32_t) ceil((float)D->ncols / D->bwidth);

  D->blocks = (mbl_t ***)malloc(rlD * sizeof(mbl_t **));
  for (i=0; i<rlD; ++i) {
    D->blocks[i]  = (mbl_t **)malloc(clD * sizeof(mbl_t *));
    for (j=0; j<clD; ++j) {
      D->blocks[i][j] = (mbl_t *)malloc(
          (D->bheight / __GBLA_NROWS_MULTILINE) * sizeof(mbl_t));
      for (k=0; k<(D->bheight / __GBLA_NROWS_MULTILINE); ++k) {
        D->blocks[i][j][k].val  = NULL;
        D->blocks[i][j][k].idx  = NULL;
        D->blocks[i][j][k].sz   = D->blocks[i][j][k].dense  = 0;
      }
    }
  }

#if __GBLA_DEBUG
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

  uint32_t npiv = 0; /*  number pivots handled */
  /*  Take maximal number of rows in blocks to be able to use array */
  /*  piv_start_idx for construction of A & B as well as for the construction of */
  /*  C & D. */
  uint32_t max_nrows =  (A->nrows > C->nrows) ? A->nrows : C->nrows;
  uint32_t piv_start_idx[(max_nrows / B->bheight) + 2];
  uint32_t block_idx/*, block_idx_2*/;

  /*  find blocks for construction of A & B */
	assert((int)M->ncols>=0);
  for (l = (int)M->ncols-1; l > -1; --l) {
    if (map->pri[l] != __GBLA_MINUS_ONE_32) {
      npiv++;
    }
    if (npiv % B->bheight == 0) {
      piv_start_idx[npiv/B->bheight]  = l;
    }
  }
  /*  loop might overwrite piv_start_idx[0] with a wrong index; */
  /*  instead of checking "npiv > 0" in each if clause we just reset */
  /*  piv_start_idx[0] after the for loop */
  piv_start_idx[0]  = M->ncols;

#if __GBLA_NEW_SPLICER
  /*
   * For code in __GBLA_NEW_SPLICER
   * Trying to optimize splicing: We know that the number of elements in the
   * lefthand side are getting fewer and fewer since it is an upper triangular
   * matrix (at least for the upper left part). So we try to give each thread 2
   * tasks at a time, namely one from the top (more stuff to do) and one from
   * the bottom (less stuff to do). Thus we hope that the workload is better
   * balanced.
   */
  uint32_t rihb[B->bheight];  /*  rows indices horizontal block */
  uint16_t cvb  = 0;          /*  current vector in block */
  uint32_t nb   = npiv/B->bheight;
  if (npiv % B->bheight != 0)
    nb++;

  /*  set leftout entries to zero */
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  if (nb % 2) {
    nb--;
    for (i = ((int)piv_start_idx[nb]-1);
        i > (int)piv_start_idx[nb+1]-1; --i) {
      if (map->pri[i] != __GBLA_MINUS_ONE_32) {
        rihb[cvb] = map->pri[i];
        cvb++;
      }
      if (cvb == B->bheight || i == 0) {
        /*  multiline version for A */
        write_lr_matrix_ml(M, A, B, map, rihb, cvb, nb, 0);

        /*  TODO: Destruct input matrix on the go */
        if (destruct_input_matrix) {
          for (j=0; j<cvb; ++j) {
            free(M->pos[rihb[j]]);
            M->pos[rihb[j]] = NULL;
            free(M->rows[rihb[j]]);
            M->rows[rihb[j]] = NULL;
          }
        /*  free memory */
        }
        cvb = 0;
      }
    }
  }
  omp_set_dynamic(0);
#pragma omp parallel private(cvb, rihb, block_idx, block_idx_2, i, j) num_threads(nthreads)
  {
    cvb  = 0;
#if __GBLA_DEBUG
    printf("nthreads %d\n",nthreads);
    printf("npiv %d -- bheight %d\n",map->npiv,B->bheight);
    printf("div %d\n", map->npiv/B->bheight);
#endif

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx < nb/2; ++block_idx) {
      /*  construct block submatrices A & B */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
      block_idx_2 = nb - block_idx - 1;
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (map->pri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        if (cvb == B->bheight || i == 0) {
          /*  multiline version for A */
          write_lr_matrix_ml(M, A, B, map, rihb, cvb, block_idx, 0);

          /*  TODO: Destruct input matrix on the go */
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          /*  free memory */
          }
          cvb = 0;
        }
      }
      for (i = ((int)piv_start_idx[block_idx_2]-1);
          i > (int)piv_start_idx[block_idx_2+1]-1; --i) {
        if (map->pri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->pri[i];
          cvb++;
        }
        if (cvb == B->bheight || i == 0) {
          /*  multiline version for A */
          write_lr_matrix_ml(M, A, B, map, rihb, cvb, block_idx_2, 0);

          /*  TODO: Destruct input matrix on the go */
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          /*  free memory */
          }
          cvb = 0;
        }
      }
    }
  }

  /*  find blocks for construction of C & D */
  npiv  = 0;
  for (i = (int)M->nrows-1; i > -1; --i) {
    if (map->npri[i] != __GBLA_MINUS_ONE_32)
      npiv++;

    if (npiv % D->bheight == 0) {
      piv_start_idx[npiv/D->bheight]  = i;
    }
  }
  /*  loop might overwrite piv_start_idx[0] with a wrong index; */
  /*  instead of checking "npiv > 0" in each if clause we just reset */
  /*  piv_start_idx[0] after the for loop */
  piv_start_idx[0]  = M->ncols;
  /*  set leftout entries to zero */
  for (i=npiv/D->bheight+1; i < (max_nrows / D->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  nb  = npiv/D->bheight;
  if (npiv % D->bheight != 0)
    nb++;
  if (nb % 2) {
    nb--;
    for (i = ((int)piv_start_idx[nb]-1);
        i > (int)piv_start_idx[nb+1]-1; --i) {
      if (map->npri[i] != __GBLA_MINUS_ONE_32) {
        rihb[cvb] = map->npri[i];
        cvb++;
      }
      if (cvb == D->bheight || i == 0) {
        /*  multiline version for C */
        write_lr_matrix_ml(M, C, D, map, rihb, cvb, nb, 1);

        if (destruct_input_matrix) {
          for (j=0; j<cvb; ++j) {
            free(M->pos[rihb[j]]);
            M->pos[rihb[j]] = NULL;
            free(M->rows[rihb[j]]);
            M->rows[rihb[j]] = NULL;
          }
        }
        cvb = 0;
      }
    }
  }
  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, block_idx_2, i, j) num_threads(nthreads)
  {
    cvb  = 0;
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx < nb/2; ++block_idx) {
      /*  construct block submatrices C & D */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
      block_idx_2 = nb - block_idx - 1;
      for (i = ((int)piv_start_idx[block_idx]-1);
          i > (int)piv_start_idx[block_idx+1]-1; --i) {
        if (map->npri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          /*  multiline version for C */
          write_lr_matrix_ml(M, C, D, map, rihb, cvb, block_idx, 1);

          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          }
          cvb = 0;
        }
      }
      for (i = ((int)piv_start_idx[block_idx_2]-1);
          i > (int)piv_start_idx[block_idx_2+1]-1; --i) {
        if (map->npri[i] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->npri[i];
          cvb++;
        }
        if (cvb == D->bheight || i == 0) {
          /*  multiline version for C */
          write_lr_matrix_ml(M, C, D, map, rihb, cvb, block_idx_2, 1);

          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          }
          cvb = 0;
        }
      }
    }
  }
#else
  /*  set leftout entries to zero */
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i, j, l) num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  /*  rows indices horizontal block */
    uint16_t cvb  = 0;          /*  current vector in block */

#if __GBLA_DEBUG
    printf("nthreads %d\n",nthreads);
    printf("npiv %d -- bheight %d\n",map->npiv,B->bheight);
    printf("div %d\n", map->npiv/B->bheight);
#endif
#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= npiv/B->bheight; ++block_idx) {
#if __GBLA_DEBUG
      printf("bi %d\n", block_idx);
      printf("piv_idx[block] %d\n",piv_start_idx[block_idx]);
      printf("piv_idx[block+1] %d\n",piv_start_idx[block_idx+1]);
#endif
      /*  construct block submatrices A & B */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
			assert((int)piv_start_idx[block_idx]>=0);
      for (l = ((int)piv_start_idx[block_idx]-1);
          l > (int)piv_start_idx[block_idx+1]-1; --l) {
        if (map->pri[l] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->pri[l];
          cvb++;
        }
        if (cvb == B->bheight || l == 0) {
          /*  multiline version for A */
          write_lr_matrix_ml(M, A, B, map, rihb, cvb, block_idx, 0);
          /*  TODO: Destruct input matrix on the go */
          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          /*  free memory */
          }
          cvb = 0;
        }
      }
    }
  }

  /*  find blocks for construction of C & D */
  npiv  = 0;
	assert((int)M->nrows>=0);
  for (l = (int)M->nrows-1; l > -1; --l) {
    if (map->npri[l] != __GBLA_MINUS_ONE_32)
      npiv++;

    if (npiv % B->bheight == 0) {
      piv_start_idx[npiv/B->bheight]  = l;
    }
  }
  /*  loop might overwrite piv_start_idx[0] with a wrong index; */
  /*  instead of checking "npiv > 0" in each if clause we just reset */
  /*  piv_start_idx[0] after the for loop */
  piv_start_idx[0]  = M->nrows;

  /*  set leftout entries to zero */
  for (i=npiv/B->bheight+1; i < (max_nrows / B->bheight) + 2; ++i)
    piv_start_idx[i] = 0;

  omp_set_dynamic(0);
#pragma omp parallel private(block_idx, i, j, l) num_threads(nthreads)
  {
    uint32_t rihb[B->bheight];  /*  rows indices horizontal block */
    uint16_t cvb  = 0;          /*  current vector in block */

#pragma omp for schedule(dynamic) nowait
    for (block_idx = 0; block_idx <= npiv/D->bheight; ++block_idx) {
      /*  construct block submatrices C & D */
      /*  Note: In the for loop we always construct block "block+1" and not block */
      /*  "block". */
      /*  TODO: Try to improve this rather strange looping. */
      for (l = ((int)piv_start_idx[block_idx]-1);
          l > (int)piv_start_idx[block_idx+1]-1; --l) {
        if (map->npri[l] != __GBLA_MINUS_ONE_32) {
          rihb[cvb] = map->npri[l];
          cvb++;
        }
        if (cvb == D->bheight || l == 0) {
          /*  multiline version for C */
          write_lr_matrix_ml(M, C, D, map, rihb, cvb, block_idx, 1);


          if (destruct_input_matrix) {
            for (j=0; j<cvb; ++j) {
              free(M->pos[rihb[j]]);
              M->pos[rihb[j]] = NULL;
              free(M->rows[rihb[j]]);
              M->rows[rihb[j]] = NULL;
            }
          }
          cvb = 0;
        }
      }
    }
  }
#endif
}



void write_blocks_lr_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, map_fl_t *map,
    ri_t *rihb, const ri_t cvb, const ri_t rbi, const bi_t density_level) {

  bi_t  bir;    /*  block index in row */
  bi_t  eil;    /*  element index in line */
  bi_t  lib;    /*  line index in block */
	/* bi_t  length; |+  local helper for block line length arithmetic +| */
  ci_t  it1, it2, i1, i2;

  /*  memory for block entries is already allocated in splice_fl_matrix() */

  /*  current loop variable i, block indices 1 (rihb[i]) and 2 (rihb[i+1]) */
  ri_t i, j, k, l, bi1, bi2;
  /* const uint32_t loop_size  = (uint32_t) ceil(cvb / __GBLA_NROWS_MULTILINE); */


  /*  column loops */
  const ci_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);
  const ci_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  /*  Usually blocks in B tend to be denser than blocks in A, thus we allocate */
  /*  already at the beginning more memory for those lines. */
  /*  We use a density threshold (which is globally defined in GBLA_config.h.in) to */
  /*  specify how much memory we allocate initially. */

  bi_t init_buffer_A, init_buffer_B;
  switch (density_level) {
    case 0: /*  splicing upper part (A & B) -- tends to be sparse */
      if (M->density > __GBLA_DENSITY_THRESHOLD) {
        init_buffer_A = (bi_t)(A->bwidth/16);
        init_buffer_B = (bi_t)(B->bwidth/4);
      } else {
        init_buffer_A = (bi_t)(A->bwidth/32);
        init_buffer_B = (bi_t)(B->bwidth/4);
      }
      break;
    case 1: /*  splicing lower part (C & D) -- tends to be denser */
      if (M->density > __GBLA_DENSITY_THRESHOLD) {
        init_buffer_A = (bi_t)(A->bwidth/4);
        init_buffer_B = (bi_t)(B->bwidth/2);
      } else {
        init_buffer_A = (bi_t)(A->bwidth/8);
        init_buffer_B = (bi_t)(B->bwidth/2);
      }
      break;
  }

  bi_t buffer_A[clA];
  bi_t buffer_B[clB];

  /*  NOTE: In LELA Martani uses a stack vector to store the values for sub */
  /*  block matrix A and first writes to B and this stack vector. Afterwards a */
  /*  looping is done over the stack vector and all values are written into A in */
  /*  one run. */
  /*  Here we do not use this feature, but write directly into A. */

  /*  ususally cvb is divisible by 2, but for the last row of blocks there might */
  /*  be only an odd number of lines in the blocks */
  uint16_t rounded_cvb  = cvb;
  if (cvb % 2)
    rounded_cvb = cvb-1;
  for (i=0; i<rounded_cvb; i=i+2) {
    bi1 = rihb[i];
    bi2 = rihb[i+1];
    i1  = 0;
    i2  = 0;
    lib = i/2;
    memset(buffer_A, 0, clA * sizeof(bi_t));
    memset(buffer_B, 0, clB * sizeof(bi_t));

    /*  loop over rows i and i+1 of M and splice correspondingly into A & B */
    while (i1 < M->rwidth[bi1] && i2 < M->rwidth[bi2]) {
      it1 = M->pos[bi1][i1];
      it2 = M->pos[bi2][i2];
      if (it1 < it2) {
        if (map->pc[it1] != __GBLA_MINUS_ONE_32) {
          bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
          eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
          /*  realloc memory if needed */
          if (A->blocks[rbi][bir][lib].sz == buffer_A[bir]) {
            realloc_block_rows(A, rbi, bir, lib, init_buffer_A, &buffer_A[bir]);
          }
          /*  set values */
          insert_block_row_data_ml_1_1(A, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG
          printf("1 %d -- %d -- %d\n", rbi,bir,lib);
          printf("1 eil %d\n",eil);
          printf("%d\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-2]);
          printf("%d\n\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-1]);
#endif
        } else {
          bir = map->npc[it1] / B->bwidth;
          eil = map->npc[it1] % B->bwidth;
          /*  realloc memory if needed */
          if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
            realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
          /*  set values */
          insert_block_row_data_ml_1_1(B, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG_LL
            printf("1 %d -- %d -- %d\n", rbi,bir,lib);
            printf("1 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
          printf("2 %d -- %d -- %d\n", rbi,bir,lib);
          printf("2 eil %d\n",eil);
#endif
        }
        i1++;
      } else {
        if (it1 > it2) {
          if (map->pc[it2] != __GBLA_MINUS_ONE_32) {
            bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
            eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
            /*  realloc memory if needed */
            if (A->blocks[rbi][bir][lib].sz == buffer_A[bir])
              realloc_block_rows(A, rbi, bir, lib, init_buffer_A, &buffer_A[bir]);
            /*  set values */
            insert_block_row_data_ml_1_2(A, M, rbi, bir, lib, eil, bi2, i2);
#if __GBLA_DEBUG
            printf("3 %d -- %d -- %d\n", rbi,bir,lib);
            printf("3 eil %d\n",eil);
          printf("%d\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-2]);
          printf("%d\n\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-1]);
#endif
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            /*  realloc memory if needed */
            if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
              realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
            /*  set values */
            insert_block_row_data_ml_1_2(B, M, rbi, bir, lib, eil, bi2, i2);
#if __GBLA_DEBUG_LL
            printf("2 %d -- %d -- %d\n", rbi,bir,lib);
            printf("2 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
            printf("4 %d -- %d -- %d\n", rbi,bir,lib);
            printf("4 eil %d\n",eil);
#endif
          }
          i2++;
        } else { /*  it1 == it2 */
          if (map->pc[it1] != __GBLA_MINUS_ONE_32) { /*  holds also for it2 */
            bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
            eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
            /*  realloc memory if needed */
            if (A->blocks[rbi][bir][lib].sz == buffer_A[bir])
              realloc_block_rows(A, rbi, bir, lib, init_buffer_A, &buffer_A[bir]);
            /*  set values */
            insert_block_row_data_ml_2(A, M, rbi, bir, lib, eil, bi1, i1, bi2, i2);
#if __GBLA_DEBUG
            printf("5 %d -- %d -- %d\n", rbi,bir,lib);
            printf("5 eil %d\n",eil);
          printf("%d\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-2]);
          printf("%d\n\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-1]);
#endif
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            /*  realloc memory if needed */
            if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
              realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
            /*  set values */
            insert_block_row_data_ml_2(B, M, rbi, bir, lib, eil, bi1, i1, bi2, i2);
#if __GBLA_DEBUG_LL
            printf("3 %d -- %d -- %d\n", rbi,bir,lib);
            printf("3 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
          }
          i1++;
          i2++;
        }
      }
    }
    /*  Now we splice rows i and i+1 separately (one of the rows we have already */
    /*  completely spliced since we are out of the above while loop. Thus maximal */
    /*  one of the following two while loops are executed (if rows i and i+1 have */
    /*  the same number of elements none of the following while loops is */
    /*  executed). */
    while (i1 < (M->rwidth[bi1])) {
      it1 = M->pos[bi1][i1];
      if (map->pc[it1] != __GBLA_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
        /*  realloc memory if needed */
        if (A->blocks[rbi][bir][lib].sz == buffer_A[bir])
          realloc_block_rows(A, rbi, bir, lib, init_buffer_A, &buffer_A[bir]);
        /*  set values */
        insert_block_row_data_ml_1_1(A, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG
        printf("7 %d -- %d -- %d\n", rbi,bir,lib);
        printf("7 eil %d\n",eil);
          printf("%d\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-2]);
          printf("%d\n\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-1]);
#endif
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        /*  set values */
        insert_block_row_data_ml_1_1(B, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG_LL
            printf("4 %d -- %d -- %d\n", rbi,bir,lib);
            printf("4 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
        printf("8 %d -- %d -- %d\n", rbi,bir,lib);
        printf("8 eil %d\n",eil);
#endif
      }
      i1++;
    }
    while (i2 < (M->rwidth[bi2])) {
      it2 = M->pos[bi2][i2];
      /* printf("it2 %d\n",it2); */
      if (map->pc[it2] != __GBLA_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it2]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it2]) % A->bwidth;
        /*  realloc memory if needed */
        if (A->blocks[rbi][bir][lib].sz == buffer_A[bir])
          realloc_block_rows(A, rbi, bir, lib, init_buffer_A, &buffer_A[bir]);
        /*  set values */
        insert_block_row_data_ml_1_2(A, M, rbi, bir, lib, eil, bi2, i2);
#if __GBLA_DEBUG
        printf("9 %d -- %d -- %d\n", rbi,bir,lib);
        printf("9 eil %d\n",eil);
          printf("%d\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-2]);
          printf("%d\n\n",A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz-1]);
#endif
      } else {
        bir = map->npc[it2] / B->bwidth;
        eil = map->npc[it2] % B->bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        /*  set values */
        insert_block_row_data_ml_1_2(B, M, rbi, bir, lib, eil, bi2, i2);
#if __GBLA_DEBUG_LL
            printf("5 %d -- %d -- %d\n", rbi,bir,lib);
            printf("5 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
        printf("10 %d -- %d -- %d\n", rbi,bir,lib);
        printf("10 eil %d\n",eil);
#endif
      }
      i2++;
    }
  }

#if __GBLA_DEBUG
  printf("i %d -- cvb %d\n",i,cvb);
#endif
  if (i<cvb) {
    /*  now we know the correct number of multilines per block */
    rounded_cvb = (rounded_cvb+2);
    bi1 = rihb[i];
    i1  = 0;
    lib = i/2;
    memset(buffer_A, 0, clA * sizeof(bi_t));
    memset(buffer_B, 0, clB * sizeof(bi_t));

    while (i1 < (M->rwidth[bi1])) {
      it1 = M->pos[bi1][i1];
      if (map->pc[it1] != __GBLA_MINUS_ONE_32) {
        bir = (A->ncols - 1 - map->pc[it1]) / A->bwidth;
        eil = (A->ncols - 1 - map->pc[it1]) % A->bwidth;
        /*  realloc memory if needed */
        if (A->blocks[rbi][bir][lib].sz == buffer_A[bir])
          realloc_block_rows(A, rbi, bir, lib, init_buffer_A, &buffer_A[bir]);
        /*  set values */
        insert_block_row_data_ml_1_1(A, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG
        printf("11 %d -- %d -- %d\n", rbi,bir,lib);
        printf("11 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        /*  set values */
        insert_block_row_data_ml_1_1(B, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG_LL
            printf("6 %d -- %d -- %d\n", rbi,bir,lib);
            printf("6 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
        printf("12 %d -- %d -- %d\n", rbi,bir,lib);
        printf("12 eil %d\n",eil);
#endif
      }
      i1++;
    }
  }
  /*  write data to A (from right to left) */
  swap_block_data(A, clA, rbi, cvb);

  /*  hybrid multirows for the righthand side block matrices? */
  if (B->hr) {
    uint32_t idx  = 0;
    /*  TODO: Implement hybrid stuff */
    for (i=0; i<clB; ++i) {
      for (j=0; j<rounded_cvb/__GBLA_NROWS_MULTILINE; ++j) {
        if ((float)B->blocks[rbi][i][j].sz / (float)B->bwidth
            < __GBLA_HYBRID_THRESHOLD) {
          B->blocks[rbi][i][j].idx =  realloc(
              B->blocks[rbi][i][j].idx,
              B->blocks[rbi][i][j].sz * sizeof(bi_t));
          B->blocks[rbi][i][j].val =  realloc(
              B->blocks[rbi][i][j].val,
              2 * B->blocks[rbi][i][j].sz * sizeof(re_t));
          continue;
        }
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
        B->blocks[rbi][i][j].idx    = NULL;
        free(B->blocks[rbi][i][j].val);
        B->blocks[rbi][i][j].val    = tmp_val_ptr;
        B->blocks[rbi][i][j].sz     = B->bwidth;
        B->blocks[rbi][i][j].dense  = 1;
      }
    }
  } else { /*  cut down memory usage */
    /*  Realloc memory usage: */
    /*  Note that A is reallocated during the swapping of the data, so we */
    /*  reallocate only B here. */
    int ctr;
    for (k=0; k<clB; ++k) {
      ctr = 0;
      for (l=0; l<rounded_cvb/2; ++l) {
        if (B->blocks[rbi][k][l].sz>0) {
          ctr = 1;
          B->blocks[rbi][k][l].idx = realloc(
              B->blocks[rbi][k][l].idx,
              B->blocks[rbi][k][l].sz * sizeof(bi_t));
          B->blocks[rbi][k][l].val = realloc(
              B->blocks[rbi][k][l].val,
              2 * B->blocks[rbi][k][l].sz  * sizeof(re_t));
        }
      }
      /*  if full block is empty, remove it! */
      if (ctr == 0) {
        free(B->blocks[rbi][k]);
        B->blocks[rbi][k] = NULL;
      }
    }
  }
}


void write_lr_matrix_ml(sm_t *M, sm_fl_ml_t *A, sbm_fl_t *B, map_fl_t *map,
    ri_t *rihb, const ri_t cvb, const ri_t rbi, const bi_t density_level) {

  bi_t  bir;    /*  block index in row */
  bi_t  eil;    /*  element index in line */
  bi_t  lib;    /*  line index in block */
	/* bi_t  length; |+  local helper for block line length arithmetic +| */
  ci_t  it1, it2, i1, i2;

  /*  memory for block entries is already allocated in splice_fl_matrix() */

  /*  current loop variable i, block indices 1 (rihb[i]) and 2 (rihb[i+1]) */
  uint32_t i, j, k, l, bi1, bi2;
  mli_t mli; /*  multiline index in A */

  /* const uint32_t loop_size  = (uint32_t) ceil(cvb / __GBLA_NROWS_MULTILINE); */


  /*  column loops */
  const ci_t clB  = (uint32_t) ceil((float)B->ncols / B->bwidth);

  /*  Usually blocks in B tend to be denser than blocks in A, but here A is just */
  /*  a multiline sub matrix without blocks, so we also allocate enough initial */
  /*  memory for the multilines in A. */
  /*  We use a density threshold (which is globally defined in GBLA_config.h.in) to */
  /*  specify how much memory we allocate initially. */
  bi_t init_buffer_A, init_buffer_B;
  switch (density_level) {
    case 0: /*  splicing upper part (A & B) -- tends to be sparse */
      if (M->density > __GBLA_DENSITY_THRESHOLD) {
        /* init_buffer_A = (bi_t)(M->density*(M->ncols/100/2)); */
        init_buffer_A = (bi_t)(2*B->bwidth);
        init_buffer_B = (bi_t)(B->bwidth/2);
      } else {
        /* init_buffer_A = (bi_t)(M->density*(M->ncols/100/8)); */
        init_buffer_A = (bi_t)(B->bwidth);
        init_buffer_B = (bi_t)(B->bwidth/2);
      }
      break;
    case 1: /*  splicing lower part (C & D) -- tends to be denser */
      if (M->density > __GBLA_DENSITY_THRESHOLD) {
        /* init_buffer_A = (bi_t)(M->density*(M->ncols/100/2)); */
        init_buffer_A = (bi_t)(4*B->bwidth);
        init_buffer_B = (bi_t)( B->bwidth/2);
      } else {
        /* init_buffer_A = (bi_t)(M->density*(M->ncols/100/4)); */
        init_buffer_A = (bi_t)(2*B->bwidth);
        init_buffer_B = (bi_t)(B->bwidth/2);
      }
      break;
  }

  mli_t buffer_A; /*  16bit is not enough, overflows appear for lager matrices! */
  bi_t buffer_B[clB];

  /*  NOTE: In LELA Martani uses a stack vector to store the values for sub */
  /*  block matrix A and first writes to B and this stack vector. Afterwards a */
  /*  looping is done over the stack vector and all values are written into A in */
  /*  one run. */
  /*  Here we do not use this feature, but write directly into A. */

  /*  ususally cvb is divisible by 2, but for the last row of blocks there might */
  /*  be only an odd number of lines in the blocks */
  ri_t rounded_cvb  = cvb;
  if (cvb % 2)
    rounded_cvb = cvb-1;
  for (i=0; i<rounded_cvb; i=i+2) {
    bi1 = rihb[i];
    bi2 = rihb[i+1];
    i1  = 0;
    i2  = 0;
    lib = i/2;

    /*  compute corresponding multiline index for A */
    mli = (rbi * B->bheight / __GBLA_NROWS_MULTILINE) + lib;

    buffer_A = 0;
    memset(buffer_B, 0, clB * sizeof(bi_t));

    /*  loop over rows i and i+1 of M and splice correspondingly into A & B */
    while (i1 < M->rwidth[bi1] && i2 < M->rwidth[bi2]) {
      it1 = M->pos[bi1][i1];
      it2 = M->pos[bi2][i2];
      if (it1 < it2) {
        if (map->pc[it1] != __GBLA_MINUS_ONE_32) {
          /*  realloc memory if needed */
          if (A->ml[mli].sz == buffer_A) {
            realloc_rows_ml(A, mli, init_buffer_A, &buffer_A);
          }
          /*  set values */
          insert_row_data_ml_1_1(A, M, mli, map->pc[it1], bi1, i1);
#if __GBLA_DEBUG
          printf("1 %d -- %d -- %d\n", rbi,bir,lib);
          printf("1 eil %d\n",eil);
#endif
        } else {
          bir = map->npc[it1] / B->bwidth;
          eil = map->npc[it1] % B->bwidth;
          /*  realloc memory if needed */
          if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
            realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
          /*  set values */
          insert_block_row_data_ml_1_1(B, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG_LL
            printf("1 %d -- %d -- %d\n", rbi,bir,lib);
            printf("1 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
          printf("2 %d -- %d -- %d\n", rbi,bir,lib);
          printf("2 eil %d\n",eil);
#endif
        }
        i1++;
      } else {
        if (it1 > it2) {
          if (map->pc[it2] != __GBLA_MINUS_ONE_32) {
          /*  realloc memory if needed */
          if (A->ml[mli].sz == buffer_A) {
            realloc_rows_ml(A, mli, init_buffer_A, &buffer_A);
          }
          /*  set values */
          insert_row_data_ml_1_2(A, M, mli, map->pc[it2], bi2, i2);
#if __GBLA_DEBUG
            printf("3 %d -- %d -- %d\n", rbi,bir,lib);
            printf("3 eil %d\n",eil);
#endif
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            /*  realloc memory if needed */
            if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
              realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
            /*  set values */
            insert_block_row_data_ml_1_2(B, M, rbi, bir, lib, eil, bi2, i2);
#if __GBLA_DEBUG_LL
            printf("2 %d -- %d -- %d\n", rbi,bir,lib);
            printf("2 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
            printf("4 %d -- %d -- %d\n", rbi,bir,lib);
            printf("4 eil %d\n",eil);
#endif
          }
          i2++;
        } else { /*  it1 == it2 */
          if (map->pc[it1] != __GBLA_MINUS_ONE_32) { /*  holds also for it2 */
            /*  realloc memory if needed */
            if (A->ml[mli].sz == buffer_A) {
              realloc_rows_ml(A, mli, init_buffer_A, &buffer_A);
            }
            /*  set values */
            insert_row_data_ml_2(A, M, mli, map->pc[it1], bi1, i1, bi2, i2);
#if __GBLA_DEBUG
            printf("5 %d -- %d -- %d\n", rbi,bir,lib);
            printf("5 eil %d\n",eil);
#endif
          } else {
            bir = map->npc[it2] / B->bwidth;
            eil = map->npc[it2] % B->bwidth;
            /*  realloc memory if needed */
            if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
              realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
            /*  set values */
            insert_block_row_data_ml_2(B, M, rbi, bir, lib, eil, bi1, i1, bi2, i2);
#if __GBLA_DEBUG_LL
            printf("3 %d -- %d -- %d\n", rbi,bir,lib);
            printf("3 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
          }
          i1++;
          i2++;
        }
      }
    }
    /*  Now we splice rows i and i+1 separately (one of the rows we have already */
    /*  completely spliced since we are out of the above while loop. Thus maximal */
    /*  one of the following two while loops are executed (if rows i and i+1 have */
    /*  the same number of elements none of the following while loops is */
    /*  executed). */
    while (i1 < (M->rwidth[bi1])) {
      it1 = M->pos[bi1][i1];
      if (map->pc[it1] != __GBLA_MINUS_ONE_32) {
        /*  realloc memory if needed */
        if (A->ml[mli].sz == buffer_A) {
          realloc_rows_ml(A, mli, init_buffer_A, &buffer_A);
        }
        /*  set values */
        insert_row_data_ml_1_1(A, M, mli, map->pc[it1], bi1, i1);
#if __GBLA_DEBUG
        printf("7 %d -- %d -- %d\n", rbi,bir,lib);
        printf("7 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        /*  set values */
        insert_block_row_data_ml_1_1(B, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG_LL
            printf("4 %d -- %d -- %d\n", rbi,bir,lib);
            printf("4 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
        printf("8 %d -- %d -- %d\n", rbi,bir,lib);
        printf("8 eil %d\n",eil);
#endif
      }
      i1++;
    }
    while (i2 < (M->rwidth[bi2])) {
      it2 = M->pos[bi2][i2];
      /* printf("it2 %d\n",it2); */
      if (map->pc[it2] != __GBLA_MINUS_ONE_32) {
        /*  realloc memory if needed */
        if (A->ml[mli].sz == buffer_A) {
          realloc_rows_ml(A, mli, init_buffer_A, &buffer_A);
        }
        /*  set values */
        insert_row_data_ml_1_2(A, M, mli, map->pc[it2], bi2, i2);
#if __GBLA_DEBUG
        printf("9 %d -- %d -- %d\n", rbi,bir,lib);
        printf("9 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it2] / B->bwidth;
        eil = map->npc[it2] % B->bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        /*  set values */
        insert_block_row_data_ml_1_2(B, M, rbi, bir, lib, eil, bi2, i2);
#if __GBLA_DEBUG_LL
            printf("5 %d -- %d -- %d\n", rbi,bir,lib);
            printf("5 eil %d\n",eil);
            printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
            printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
        printf("10 %d -- %d -- %d\n", rbi,bir,lib);
        printf("10 eil %d\n",eil);
#endif
      }
      i2++;
    }
    if (A->ml[mli].sz>0) {
      A->ml[mli].idx  = (mli_t*) realloc(A->ml[mli].idx, A->ml[mli].sz* sizeof(mli_t));
      A->ml[mli].val  = (re_t*)  realloc(A->ml[mli].val, 2 * A->ml[mli].sz* sizeof(re_t));
    }
  }

#if __GBLA_DEBUG
  printf("i %d -- cvb %d\n",i,cvb);
#endif
  if (i<cvb) {
    bi1 = rihb[i];
    i1  = 0;
    lib = i/2;
    buffer_A = 0;
    memset(buffer_B, 0, clB * sizeof(bi_t));

    /*  compute corresponding multiline index for A */
    mli = (rbi * B->bheight / __GBLA_NROWS_MULTILINE) + lib;

    while (i1 < (M->rwidth[bi1])) {
      it1 = M->pos[bi1][i1];
      if (map->pc[it1] != __GBLA_MINUS_ONE_32) {
        /*  realloc memory if needed */
        if (A->ml[mli].sz == buffer_A) {
          realloc_rows_ml(A, mli, init_buffer_A, &buffer_A);
        }
        /*  set values */
        insert_row_data_ml_1_1(A, M, mli, map->pc[it1], bi1, i1);
#if __GBLA_DEBUG
        printf("11 %d -- %d -- %d\n", rbi,bir,lib);
        printf("11 eil %d\n",eil);
#endif
      } else {
        bir = map->npc[it1] / B->bwidth;
        eil = map->npc[it1] % B->bwidth;
        /*  realloc memory if needed */
        if (B->blocks[rbi][bir][lib].sz == buffer_B[bir])
          realloc_block_rows(B, rbi, bir, lib, init_buffer_B, &buffer_B[bir]);
        /*  set values */
        insert_block_row_data_ml_1_1(B, M, rbi, bir, lib, eil, bi1, i1);
#if __GBLA_DEBUG_LL
        printf("6 %d -- %d -- %d\n", rbi,bir,lib);
        printf("6 eil %d\n",eil);
        printf("1 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-2]);
        printf("2 %d\n",B->blocks[rbi][bir][lib].val[2*B->blocks[rbi][bir][lib].sz-1]);
#endif
#if __GBLA_DEBUG
        printf("12 %d -- %d -- %d\n", rbi,bir,lib);
        printf("12 eil %d\n",eil);
#endif
      }
      i1++;
    }
    if (A->ml[mli].sz >0) {
      A->ml[mli].idx  = (mli_t*) realloc(A->ml[mli].idx, A->ml[mli].sz* sizeof(mli_t));
      A->ml[mli].val  = (re_t*)  realloc(A->ml[mli].val, 2 * A->ml[mli].sz* sizeof(re_t));
    }
  }

  /*  hybrid multirows for the righthand side block matrices? */
  if (B->hr) {
    uint32_t idx  = 0;
    /*  TODO: Implement hybrid stuff */
    for (i=0; i<clB; ++i) {
      for (j=0; j<rounded_cvb/__GBLA_NROWS_MULTILINE; ++j) {
        if ((float)B->blocks[rbi][i][j].sz / (float)B->bwidth
            < __GBLA_HYBRID_THRESHOLD) {
          B->blocks[rbi][i][j].idx =  realloc(
              B->blocks[rbi][i][j].idx,
              B->blocks[rbi][i][j].sz * sizeof(bi_t));
          B->blocks[rbi][i][j].val =  realloc(
              B->blocks[rbi][i][j].val,
              2 * B->blocks[rbi][i][j].sz * sizeof(re_t));
          continue;
        }
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
        B->blocks[rbi][i][j].idx    = NULL;
        free(B->blocks[rbi][i][j].val);
        B->blocks[rbi][i][j].val    = tmp_val_ptr;
        B->blocks[rbi][i][j].sz     = B->bwidth;
        B->blocks[rbi][i][j].dense  = 1;
      }
    }
  } else { /*  cut down memory usage */
    /*  Realloc memory usage: */
    /*  Note that A is reallocated during the swapping of the data, so we */
    /*  reallocate only B here. */
    int ctr;
    for (k=0; k<clB; ++k) {
      ctr = 0;
      for (l=0; l<rounded_cvb/2; ++l) {
        if (B->blocks[rbi][k][l].sz>0) {
          ctr = 1;
          B->blocks[rbi][k][l].idx = realloc(
              B->blocks[rbi][k][l].idx,
              B->blocks[rbi][k][l].sz * sizeof(bi_t));
          B->blocks[rbi][k][l].val = realloc(
              B->blocks[rbi][k][l].val,
              2 * B->blocks[rbi][k][l].sz  * sizeof(re_t));
        }
      }
      /*  if full block is empty, remove it! */
      if (ctr == 0) {
        free(B->blocks[rbi][k]);
        B->blocks[rbi][k] = NULL;
      }
    }
  }
}

/* vim:sts=2:sw=2:ts=2:
 */
