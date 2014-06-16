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
      if (map->pr_idxs_by_entry[entry] == __GB_MINUS_ONE_32) {
        map->pr_idxs_by_entry[entry]  = i;
        npiv++;
      } else { // check for a sparser pivot row (see ELAGB talk from Lachartre)
        if (M->rwidth[map->pr_idxs_by_entry[entry]] > M->rwidth[i]) {
          map->npr_idxs[map->pr_idxs_by_entry[entry]]  = map->pr_idxs_by_entry[entry];
          map->pr_idxs_by_entry[entry] = i;
        } else {
          map->npr_idxs[i]  = i;
        }
      }
    } else {
      map->npr_idxs[i]  = i;
    }
  }

  ci_t pc_idx = 0, npc_idx = 0, j;

  // construct pivot columns and non-pivot columns maps and the corresponding
  // reverse maps
  for (j=0; j<M->ncols; ++j) {
    if (map->pr_idxs_by_entry[j] !=  __GB_MINUS_ONE_32) {
      map->pc[j]           = pc_idx;
      map->pc_rev[pc_idx]  = j;
      pc_idx++;
    } else {
      map->npc[j]          = npc_idx;
      map->npc_rev[npc_idx] = j;
      npc_idx++;
    }
  }
  return map;
}

void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D, map_fl_t *map) {
  
  // construct index map for Faugère-Lachartre decomposition of matrix M
  map  = construct_fl_map(M);

  // 
}
