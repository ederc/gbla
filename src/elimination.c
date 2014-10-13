#include "elimination.h"

int elim_fl_A_block(sbm_fl_t *A, sbm_fl_t *B, int nthrds) {
  ci_t i, rc;
  const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth);
  const ri_t rlA  = (ci_t) ceil((float) A->nrows / A->bheight);
  #pragma omp parallel num_threads(nthrds)
  {
    #pragma omp for
      // each task takes one block column of B
      for (i=0; i<clB; ++i) {
        #pragma omp task
        {
          rc  = elim_fl_A_blocks_task(A, B, i, rlA);
        }
      }
    #pragma omp taskwait
  }
  return 0;
}

int elim_fl_A_ml_block(sm_fl_ml_t *A, sbm_fl_t *B, int nthrds) {
  return 0;
}

int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A) {
  bi_t i;
  ri_t j, k;
  //re_l_t *dense_block[B->bheight] __attribute__((aligned(0x1000)));
  re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *));
  ci_t size = B->bwidth * sizeof(re_l_t);
  for (i=0; i<B->bheight; ++i) {
    posix_memalign((void **)&dense_block[i], 16, size);
  }

  for (j=0; j<nbrows_A; ++j) {
    const ri_t first_block_idx  = 0;
    // TODO: Martani uses a minimum possibly for zero blocks at the end
    const ri_t last_block_idx   = nbrows_A- (nbrows_A - 1) +1;

    // set dense block entries to zero
    for (k=0; k<B->bheight; ++k)
      memset(dense_block[k], 0, size);

    // copy sparse block data to dense representation
    if (B->blocks[j][block_col_idx_B] != NULL)
      copy_sparse_to_dense_block(B->blocks[j][block_col_idx_B], dense_block, B->bheight, B->bwidth);

    for (k=0; k<last_block_idx; ++k)
      red_with_rectangular_block(A->blocks[j][k], B->blocks[j][k], dense_block, B->bheight);

    red_with_triangular_block(A->blocks[j][last_block_idx], B->blocks[j][k], dense_block, B->bheight);

  }
  for (i=0; i<B->bheight; ++i)
    free(dense_block[i]);

  return 0;
}

void red_with_triangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight) {
  bi_t i;
  if (block_A == NULL || block_B == NULL)
    return;

  for (i=0; i<bheight/2; ++i) {

  }
}
void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B, ri_t bheight) {
  bi_t i;
  if (block_A == NULL || block_B == NULL)
    return;

  for (i=0; i<bheight/2; ++i) {

  }
}
