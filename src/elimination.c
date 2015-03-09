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

#include "elimination.h"

#define DDEBUG 0
#define DDEBUG_C 0
#define DDEBUG_D 0
#define DDEBUG_DD 0
#define DDEBUG_DONE 0
#define DDEBUG_D_ONE 0
#define DENSE_MULTILINE_C_AFTER_ELIM  0

#define TASKS 0
#define GBLA_WITH_FFLAS 0

// global variables for echelonization of D
static omp_lock_t echelonize_lock;
static ri_t global_next_row_to_reduce;
static ri_t global_last_piv;
static wl_t waiting_global;

// global variables for reduction of C
static omp_lock_t reduce_C_lock;
static ri_t reduce_C_next_col_to_reduce;

// global variables for reduction of C
static omp_lock_t reduce_A_lock;
static ri_t reduce_A_next_col_to_reduce;

#if TASKS
int elim_fl_A_block(sbm_fl_t **A_in, sbm_fl_t *B, mod_t modulus, int nthrds) {
  sbm_fl_t *A = *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth);
  const ri_t rlA  = (ri_t) ceil((float) A->nrows / A->bheight);
  const ci_t clA  = (ci_t) ceil((float) A->ncols / A->bwidth);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<clB; ++i) {
#pragma omp task
      {
        rc  = elim_fl_A_blocks_task(A, B, i, rlA, modulus);
      }
    }
#pragma omp taskwait
  }
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i] != NULL) {
          for (k=0; k<A->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(A->blocks[j][i][k].idx);
            A->blocks[j][i][k].idx  = NULL;
            free(A->blocks[j][i][k].val);
            A->blocks[j][i][k].val  = NULL;
          }
          free(A->blocks[j][i]);
          A->blocks[j][i] = NULL;
        }
      }
      free(A->blocks[j]);
      A->blocks[j]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
  return 0;
}

int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A, mod_t modulus) {
bi_t i;
ri_t j, k;
re_l_t *dense_block[B->bheight] __attribute__((aligned(0x1000)));
//re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *));
uint64_t size = B->bwidth * sizeof(re_l_t);
for (i=0; i<B->bheight; ++i) {
  posix_memalign((void **)&dense_block[i], 16, size);
}
for (j=0; j<nbrows_A; ++j) {
  const ri_t first_block_idx  = 0;

  // set dense block entries to zero
  for (k=0; k<B->bheight; ++k)
    memset(dense_block[k], 0, size);

  // copy sparse block data to dense representation
  if (B->blocks[j][block_col_idx_B] != NULL)
    copy_sparse_to_dense_block(B->blocks[j][block_col_idx_B], dense_block,
        B->bheight, B->bwidth);

  for (k=0; k<j; ++k) {
    red_with_rectangular_block(A->blocks[j][k], B->blocks[k][block_col_idx_B],
        dense_block, B->bheight, 1, modulus);
    /*
        printf("RECTANGULAR DONE\n");
        for (int kk=0; kk<B->bheight; ++kk) {
        for (int ll=0; ll<B->bheight; ++ll) {
        printf("(%d,%d) %ld ",kk,ll,dense_block[kk][ll]);
        }
        printf("\n");
        }
        */
  }

  red_with_triangular_block(A->blocks[j][j], dense_block,
      B->bheight, 1, modulus);
  /*printf("TRIANGULAR DONE\n");
    for (int kk=0; kk<B->bheight; ++kk) {
    for (int ll=0; ll<B->bheight; ++ll) {
    printf("%ld ",dense_block[kk][ll]);
    }
    printf("\n");
    }
    */


  //printf("OUT BEFORE %p\n",B->blocks[j][block_col_idx_B]);
  copy_dense_block_to_sparse(
      dense_block, &B->blocks[j][block_col_idx_B],
      B->bheight, B->bwidth, modulus);
  //printf("OUT AFTERWARDS %p\n",B->blocks[j][block_col_idx_B]);
#if DDDEBUG
  printf("after copying\n");
  if (B->blocks[j][block_col_idx_B] != NULL) {
    for (int kk=0; kk<B->bheight/__GB_NROWS_MULTILINE; ++kk) {
      if (B->blocks[j][block_col_idx_B][kk].sz>0) {
        printf("%d\n",kk);
        for (int ll=0; ll<B->blocks[j][block_col_idx_B][kk].sz; ++ll) {
          printf("%d %d ",B->blocks[j][block_col_idx_B][kk].val[2*ll], B->blocks[j][block_col_idx_B][kk].val[2*ll+1]);
        }
        printf("\n");
      }
    }
  }
#endif
}
for (i=0; i<B->bheight; ++i) {
  free(dense_block[i]);
  dense_block[i]  = NULL;
}

return 0;
}
#else
int elim_fl_A_block(sbm_fl_t **A_in, sbm_fl_t *B, mod_t modulus, int nthrds) {
  sbm_fl_t *A = *A_in;
  ci_t i, rc;
  ri_t j, k;
  const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth);
  const ri_t rlA  = (ri_t) ceil((float) A->nrows / A->bheight);
  const ci_t clA  = (ci_t) ceil((float) A->ncols / A->bwidth);

  reduce_A_next_col_to_reduce = 0;
  omp_init_lock(&reduce_A_lock);

#pragma omp parallel shared(reduce_A_next_col_to_reduce) num_threads(nthrds)
    {
#pragma omp for nowait
      for (i=0; i<nthrds; ++i) {
        rc  = elim_fl_A_blocks_task(A, B, i, rlA, modulus);
      }
    }
  omp_destroy_lock(&reduce_A_lock);
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i] != NULL) {
          for (k=0; k<A->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(A->blocks[j][i][k].idx);
            A->blocks[j][i][k].idx  = NULL;
            free(A->blocks[j][i][k].val);
            A->blocks[j][i][k].val  = NULL;
          }
          free(A->blocks[j][i]);
          A->blocks[j][i] = NULL;
        }
      }
      free(A->blocks[j]);
      A->blocks[j]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
  return 0;
}

int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B, ci_t block_col_idx_B, ri_t nbrows_A, mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[B->bheight] __attribute__((aligned(0x1000)));
  //re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *));
  uint64_t size = B->bwidth * sizeof(re_l_t);
  for (i=0; i<B->bheight; ++i) {
    posix_memalign((void **)&dense_block[i], 16, size);
  }

  const ci_t clB  = (ci_t) ceil((float) B->ncols / B->bwidth);
  uint32_t lci; // local column index
  while (1) {
    omp_set_lock(&reduce_A_lock);
    if (reduce_A_next_col_to_reduce < clB) {
      lci = reduce_A_next_col_to_reduce;
      ++reduce_A_next_col_to_reduce;
    } else {
      omp_unset_lock(&reduce_A_lock);
      break;
    }
    omp_unset_lock(&reduce_A_lock);
    for (j=0; j<nbrows_A; ++j) {
      const ri_t first_block_idx  = 0;

      // set dense block entries to zero
      for (k=0; k<B->bheight; ++k) {
        memset(dense_block[k], 0, size);
      }
      // copy sparse block data to dense representation
      if (B->blocks[j][lci] != NULL)
        copy_sparse_to_dense_block(B->blocks[j][lci], dense_block,
            B->bheight, B->bwidth);

      for (k=0; k<j; ++k) {
        red_with_rectangular_block(A->blocks[j][k], B->blocks[k][lci],
            dense_block, B->bheight, 1, modulus);
      }

      red_with_triangular_block(A->blocks[j][j], dense_block,
          B->bheight, 1, modulus);


      copy_dense_block_to_sparse(
          dense_block, &B->blocks[j][lci],
          B->bheight, B->bwidth, modulus);
    }
  }
  for (i=0; i<B->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#endif
#if TASKS
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C_in, sbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  sbm_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = (ci_t) ceil((float) D->ncols / D->bwidth);
  const ri_t rlC  = (ri_t) ceil((float) C->nrows / C->bheight);
  const ri_t clC  = (ci_t) ceil((float) C->ncols / C->bwidth);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<clD; ++i) {
#pragma omp task
      {
        rc  = elim_fl_C_blocks_task(B, C, D, i, rlC, clC, inv_scalars, modulus);
      }
    }
#pragma omp taskwait
  }
  // free C
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlC; ++j) {
      for (i=0; i<clC; ++i) {
        if (C->blocks[j][i] != NULL) {
          for (k=0; k<C->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(C->blocks[j][i][k].idx);
            C->blocks[j][i][k].idx  = NULL;
            free(C->blocks[j][i][k].val);
            C->blocks[j][i][k].val  = NULL;
          }
          free(C->blocks[j][i]);
          C->blocks[j][i] = NULL;
        }
      }
      free(C->blocks[j]);
      C->blocks[j]  = NULL;
    }
  }
  free(C->blocks);
  C->blocks = NULL;
  free(C);
  C = NULL;
  *C_in = C;
  return 0;
}

int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[D->bheight] __attribute__((aligned(0x1000)));
  //re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *));
  uint64_t size = D->bwidth * sizeof(re_l_t);
  for (i=0; i<D->bheight; ++i) {
    posix_memalign((void **)&dense_block[i], 16, size);
  }

  // take maximum number and check if blocks are NULL in loop
  const ri_t last_block_idx = nbcols_C;

  for (j=0; j<nbrows_C; ++j) {
    const ri_t first_block_idx  = 0;

    // set dense block entries to zero
    for (k=0; k<D->bheight; ++k)
      memset(dense_block[k], 0, size);

    // copy sparse block data to dense representation
    if (D->blocks[j][block_col_idx_D] != NULL)
      copy_sparse_to_dense_block(D->blocks[j][block_col_idx_D], dense_block,
          D->bheight, D->bwidth);

    for (k=0; k<last_block_idx; ++k) {
#if DDEBUG_C
      printf("j %d -- k %d -- lci %d\n", j, k, block_col_idx_D);
#endif
      if (C->blocks[j][k] != NULL) {
        red_with_rectangular_block(C->blocks[j][k], B->blocks[k][block_col_idx_D],
            dense_block, B->bheight, inv_scalars, modulus);
      }
      /*
         printf("RECTANGULAR DONE\n");
         for (int kk=0; kk<B->bheight; ++kk) {
         for (int ll=0; ll<B->bheight; ++ll) {
         printf("(%d,%d) %ld ",kk,ll,dense_block[kk][ll]);
         }
         printf("\n");
         }
         */
    }
    //printf("OUT BEFORE %p\n",B->blocks[j][block_col_idx_B]);
    copy_dense_block_to_sparse(
        dense_block, &D->blocks[j][block_col_idx_D],
        D->bheight, D->bwidth, modulus);
    //printf("OUT AFTERWARDS %p\n",B->blocks[j][block_col_idx_B]);
#if DDEBUG_C
    printf("after copying\n");
    if (D->blocks[j][block_col_idx_D] != NULL) {
      for (int kk=0; kk<D->bheight/__GB_NROWS_MULTILINE; ++kk) {
        if (D->blocks[j][block_col_idx_D][kk].sz>0) {
          printf("%d\n",kk);
          for (int ll=0; ll<D->blocks[j][block_col_idx_D][kk].sz; ++ll) {
            printf("%d %d ",D->blocks[j][block_col_idx_D][kk].val[2*ll], D->blocks[j][block_col_idx_D][kk].val[2*ll+1]);
          }
          printf("\n");
        }
      }
    }
#endif
  }

  for (i=0; i<D->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#else
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C_in, sbm_fl_t *D,
    const int inv_scalars, const mod_t modulus, const int nthrds) {

  sbm_fl_t *C = *C_in;

  ci_t i, rc;
  ri_t j, k;
  const ci_t clD  = (ci_t) ceil((float) D->ncols / D->bwidth);
  const ri_t rlC  = (ri_t) ceil((float) C->nrows / C->bheight);
  const ri_t clC  = (ci_t) ceil((float) C->ncols / C->bwidth);

  reduce_C_next_col_to_reduce = 0;
  omp_init_lock(&reduce_C_lock);

#pragma omp parallel shared(reduce_C_next_col_to_reduce) num_threads(nthrds)
    {
#pragma omp for nowait
      for (i=0; i<nthrds; ++i) {
        rc  = elim_fl_C_blocks_task(B, C,  D, i, rlC, clC, inv_scalars, modulus);
      }
    }
  omp_destroy_lock(&reduce_C_lock);
  // free C
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j,k)
    for (j=0; j<rlC; ++j) {
      for (i=0; i<clC; ++i) {
        if (C->blocks[j][i] != NULL) {
          for (k=0; k<C->bheight/__GB_NROWS_MULTILINE; ++k) {
            free(C->blocks[j][i][k].idx);
            C->blocks[j][i][k].idx  = NULL;
            free(C->blocks[j][i][k].val);
            C->blocks[j][i][k].val  = NULL;
          }
          free(C->blocks[j][i]);
          C->blocks[j][i] = NULL;
        }
      }
      free(C->blocks[j]);
      C->blocks[j]  = NULL;
    }
  }
  free(C->blocks);
  C->blocks = NULL;
  free(C);
  C = NULL;
  *C_in = C;
  return 0;
}

int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
  const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
  const int inv_scalars, const mod_t modulus) {

  const ci_t clD  = (ci_t) ceil((float) D->ncols / D->bwidth);
  bi_t i;
  ri_t j, k;
  re_l_t *dense_block[D->bheight] __attribute__((aligned(0x1000)));
  //re_l_t **dense_block  = (re_l_t **)malloc(B->bheight * sizeof(re_l_t *));
  uint64_t size = D->bwidth * sizeof(re_l_t);
  for (i=0; i<D->bheight; ++i) {
    posix_memalign((void **)&dense_block[i], 16, size);
  }

  uint32_t lci; // local column index
  // take maximum number and check if blocks are NULL in loop
  const ri_t last_block_idx = nbcols_C;

  while (1) {
    omp_set_lock(&reduce_C_lock);
    if (reduce_C_next_col_to_reduce < clD) {
      lci = reduce_C_next_col_to_reduce;
      ++reduce_C_next_col_to_reduce;
    } else {
      omp_unset_lock(&reduce_C_lock);
      break;
    }
    omp_unset_lock(&reduce_C_lock);

    for (j=0; j<nbrows_C; ++j) {
      const ri_t first_block_idx  = 0;

      // set dense block entries to zero
      for (k=0; k<D->bheight; ++k)
        memset(dense_block[k], 0, size);

      // copy sparse block data to dense representation
      if (D->blocks[j][lci] != NULL) {
        copy_sparse_to_dense_block(D->blocks[j][lci], dense_block,
            D->bheight, D->bwidth);
      }
      for (k=0; k<last_block_idx; ++k) {
#if DDEBUG_C
        printf("j %d -- k %d -- lci %d\n", j, k, block_col_idx_D);
#endif
        if (C->blocks[j][k] != NULL) {
          red_with_rectangular_block(C->blocks[j][k], B->blocks[k][lci],
              dense_block, B->bheight, inv_scalars, modulus);
        }
      }
      copy_dense_block_to_sparse(
          dense_block, &D->blocks[j][lci],
          D->bheight, D->bwidth, modulus);
    }
  }

  for (i=0; i<D->bheight; ++i) {
    free(dense_block[i]);
    dense_block[i]  = NULL;
  }

  return 0;
}
#endif

void red_with_triangular_block(mbl_t *block_A, re_l_t **dense_block,
    const ri_t bheight, const int inv_scalars, const mod_t modulus) {
int i, j, k;

int last_idx;

for (i=0; i<bheight/2; ++i) {
  if (block_A[i].sz == 0)
    continue;

  last_idx  = -1;
  if (block_A[i].val[2*(block_A[i].sz-1)+1] == 0)
    last_idx  = block_A[i].sz-1;
  else
    last_idx  = block_A[i].sz-2;

  register uint32_t Av1_col1, Av2_col1;
  register uint32_t Av1_col2, Av2_col2;
  bi_t Ap1, Ap2;

#if DDDEBUG
  printf("lidx %d\n",last_idx);
#endif
  for (j=0; j<last_idx; ++j) {
    Ap1       = block_A[i].idx[j];
    Av1_col1  = block_A[i].val[2*j];
    Av2_col1  = block_A[i].val[2*j+1];

    if (inv_scalars == 1) {
      if (Av1_col1 != 0)
        Av1_col1  = (uint32_t)modulus - Av1_col1;
      if (Av2_col1 != 0)
        Av2_col1  = (uint32_t)modulus - Av2_col1;
    }

    if ((Ap1 % 2) == 0 && (j < last_idx-1)) {
      Ap2  = block_A[i].idx[j+1];
      if (Ap2 == Ap1+1) { // AXPY two rows
        Av1_col2  = block_A[i].val[2*(j+1)];
        Av2_col2  = block_A[i].val[2*(j+1)+1];

        if (inv_scalars == 1) {
          if (Av1_col2 != 0)
            Av1_col2  = (uint32_t)modulus - Av1_col2;
          if (Av2_col2 != 0)
            Av2_col2  = (uint32_t)modulus - Av2_col2;
        }
        ++j;

        dense_scal_mul_sub_2_rows_array_array(
            Av1_col1, Av2_col1, Av1_col2, Av2_col2, bheight,
            dense_block[Ap1], dense_block[Ap1+1],
            dense_block[2*i], dense_block[2*i+1]);
      } else { // AXPY one row
        dense_scal_mul_sub_1_row_array_array(
            Av1_col1, Av2_col1, bheight,
            dense_block[Ap1],
            dense_block[2*i], dense_block[2*i+1]);
      }
    } else { // AXPY one row
      dense_scal_mul_sub_1_row_array_array(
          Av1_col1, Av2_col1, bheight,
          dense_block[Ap1],
          dense_block[2*i], dense_block[2*i+1]);
    }
  }

  // do modular reduction on dense array row
  //printf("outside %d\n", 2*i);
  red_dense_array_modular(dense_block[2*i], bheight, modulus);

  // reduce lines within the same multiline
  if (block_A[i].sz > 1) {
    j         = block_A[i].sz-2;
    Av1_col1  = block_A[i].val[2*j+1];
    Ap1       = block_A[i].idx[j];

    if(Av1_col1 != 0) {
      if (inv_scalars == 1)
        if (Av1_col1 != 0)
          Av1_col1  = (uint32_t)modulus - Av1_col1;

      const bi_t offset1  = 2*i+1;
      const bi_t offset2  = 2*i;

      for (k=0; k<bheight; ++k)
        dense_block[offset1][k] +=  (uint32_t) Av1_col1 * dense_block[offset2][k];
    }
  }

  // do modular reduction on dense array row
  //printf("outside2 %d\n", 2*i+1);
  red_dense_array_modular(dense_block[2*i+1], bheight, modulus);
}
}

void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_block,
    const ri_t bheight, const int inv_scalars, const mod_t modulus) {
bi_t i, j;
if (block_A == NULL || block_B == NULL)
  return;

for (i=0; i<bheight/2; ++i) {
  const bi_t is_sparse  = block_A[i].dense == 0 ? 1 : 0;
  const bi_t N          = is_sparse == 1 ? block_A[i].sz : bheight;
  for (j=0; j<N; ++j) {
    //printf("%d && %d\n",i,j);
    /*
  for (int kk=0; kk<bheight; ++kk) {
    for (int ll=0; ll<bheight; ++ll) {
      printf("%ld ",dense_block[kk][ll]);
    }
    printf("\n");
  }
    */
    //printf("%d // %d // %d\n",i,j, N);
    const bi_t Ap1  = is_sparse == 1 ? block_A[i].idx[j] : j;
    register uint32_t Av1_col1  = block_A[i].val[2*j];
    register uint32_t Av2_col1  = block_A[i].val[2*j+1];

    if (inv_scalars == 1) {
      if (Av1_col1 != 0)
        Av1_col1  = (uint32_t)modulus - Av1_col1;
      if (Av2_col1 != 0)
        Av2_col1  = (uint32_t)modulus - Av2_col1;
    }
    if (((Ap1 % 2) == 0) && (j < (uint32_t)(N - 1))) {
      const bi_t Ap2  = (is_sparse == 1 ? block_A[i].idx[j+1] : j+1);
      if (Ap2 == Ap1+1) { // AXPY two rows
        register uint32_t Av1_col2  = block_A[i].val[2*(j+1)];
        register uint32_t Av2_col2  = block_A[i].val[2*(j+1)+1];

        if (inv_scalars == 1) {
          if (Av1_col2 != 0)
            Av1_col2  = (uint32_t)modulus - Av1_col2;
          if (Av2_col2 != 0)
            Av2_col2  = (uint32_t)modulus - Av2_col2;
        }
        ++j;
        if (block_B[Ap1 / __GB_NROWS_MULTILINE].dense == 0) {
          //printf("1S Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE);
          sparse_scal_mul_sub_2_rows_vect_array(
              Av1_col1, Av2_col1, Av1_col2, Av2_col2,
              block_B[Ap1 / __GB_NROWS_MULTILINE],
              dense_block[2*i], dense_block[2*i+1]);
          //printf("1 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]);
        } else {
          //printf("1D Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE);
          dense_scal_mul_sub_2_rows_vect_array(
              Av1_col1, Av2_col1, Av1_col2, Av2_col2,
              block_B[Ap1 / __GB_NROWS_MULTILINE], bheight,
              dense_block[2*i], dense_block[2*i+1]);
          //printf("2 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]);
        }
      } else { // AXPY one row
        if (block_B[Ap1 / __GB_NROWS_MULTILINE].dense == 0) {
          //printf("2S Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE);
          sparse_scal_mul_sub_1_row_vect_array(
              Av1_col1, Av2_col1,
              block_B[Ap1 / __GB_NROWS_MULTILINE],
              Ap1 % __GB_NROWS_MULTILINE,
              dense_block[2*i], dense_block[2*i+1]);
          //printf("3 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]);
        } else {
          //printf("2D Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE);
          dense_scal_mul_sub_1_row_vect_array(
              Av1_col1, Av2_col1,
              block_B[Ap1 / __GB_NROWS_MULTILINE],
              Ap1 % __GB_NROWS_MULTILINE, bheight,
              dense_block[2*i], dense_block[2*i+1]);
          //printf("4 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]);
        }
      }
    } else { // AXPY one row
      if (block_B[Ap1 / __GB_NROWS_MULTILINE].dense == 0) {
      //printf("3S Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE);
        sparse_scal_mul_sub_1_row_vect_array(
            Av1_col1, Av2_col1,
            block_B[Ap1 / __GB_NROWS_MULTILINE],
            Ap1 % __GB_NROWS_MULTILINE,
            dense_block[2*i], dense_block[2*i+1]);
          //printf("5 %ld|%ld  ",dense_block[2*i][255],dense_block[2*i+1][255]);
      } else {
      //printf("3D Ap1 %d -- Ap1/ML %d\n", Ap1, Ap1/__GB_NROWS_MULTILINE);
        dense_scal_mul_sub_1_row_vect_array(
            Av1_col1, Av2_col1,
            block_B[Ap1 / __GB_NROWS_MULTILINE],
            Ap1 % __GB_NROWS_MULTILINE, bheight,
            dense_block[2*i], dense_block[2*i+1]);
      }
    }
  }
}
}

void copy_dense_block_to_sparse(
  re_l_t **dense_block, mbl_t **sparse_block_in, int bheight, int bwidth, mod_t modulus) {

mbl_t *sparse_block = *sparse_block_in;
bi_t i,j, k, ctr, buffer, rows_empty;

// if block was empty in the beginning, reallocate memory
if (sparse_block == NULL) {
  sparse_block  = (mbl_t *)malloc(bheight / __GB_NROWS_MULTILINE * sizeof(mbl_t));
  for (i=0; i<(bheight / __GB_NROWS_MULTILINE); ++i) {
    sparse_block[i].val = NULL;
    sparse_block[i].idx = NULL;
    sparse_block[i].sz  = sparse_block[i].dense =  0;
  }
}
rows_empty  = 0;
for (i=0; i<bheight/2; ++i) {
  if (sparse_block[i].dense == 1) {
    sparse_block[i].idx = (bi_t *)malloc(bwidth * sizeof(bi_t));
  }
    ctr                 = 0;
    buffer              = sparse_block[i].sz;
    sparse_block[i].sz  = 0;
    sparse_block[i].dense = 0;
    for (j=0; j<bwidth; ++j) {
      if (dense_block[2*i][j] != 0) {
        dense_block[2*i][j] = (re_l_t)(dense_block[2*i][j] % modulus);
      }
      if (dense_block[2*i+1][j] != 0) {
        dense_block[2*i+1][j] = (re_l_t)(dense_block[2*i+1][j] % modulus);
      }
      if (dense_block[2*i][j] != 0 || dense_block[2*i+1][j] != 0) {
        if (ctr >= buffer) {
          // if this happens just allocate memory for full block multiline row
          sparse_block[i].idx = realloc(
              sparse_block[i].idx, bwidth * sizeof(bi_t));
          sparse_block[i].val = realloc(
              sparse_block[i].val, 2 * bwidth * sizeof(re_t));
          buffer  = bwidth;
        }
        sparse_block[i].idx[ctr]      = j;
        sparse_block[i].val[2*ctr]    = (re_t) dense_block[2*i][j];
        sparse_block[i].val[2*ctr+1]  = (re_t) dense_block[2*i+1][j];
        //printf("%d-vals %d .. %d   ",ctr, sparse_block[i].val[2*ctr], sparse_block[i].val[2*ctr+1]);
        ctr++;
      }
    }
    sparse_block[i].sz  = ctr;
    // try to get hybrid representation, i.e. try to make multiline row dense

    if ((float)sparse_block[i].sz / (float)bwidth < __GB_HYBRID_THRESHOLD) {
      // realloc memory, cut it down as much as possible
      if (sparse_block[i].sz>0) {
        sparse_block[i].idx = realloc(
            sparse_block[i].idx,
            sparse_block[i].sz * sizeof(bi_t));
        sparse_block[i].val = realloc(
            sparse_block[i].val,
            2 * sparse_block[i].sz  * sizeof(re_t));
      } else {
        free(sparse_block[i].idx);
        sparse_block[i].idx = NULL;
        free(sparse_block[i].val);
        sparse_block[i].val = NULL;
        rows_empty++;
      }
    } else { // dense multiline row
      re_t *tmp  = (re_t *)malloc(2 * bwidth * sizeof(re_t));
      ctr  = 0;
      //for (k=0; k<bwidth; ++k) {
      k = 0;
      while (ctr<sparse_block[i].sz) {
        if (sparse_block[i].idx[ctr] == k) {
          tmp[2*k]   = sparse_block[i].val[2*ctr];
          tmp[2*k+1] = sparse_block[i].val[2*ctr+1];
          ctr++;
          k++;
        } else {
          tmp[2*k]   = 0;
          tmp[2*k+1] = 0;
          k++;
        }
      }
      for (; k<bwidth; ++k) {
        tmp[2*k]   = 0;
        tmp[2*k+1] = 0;
      }
      free(sparse_block[i].idx);
      sparse_block[i].idx   = NULL;
      free(sparse_block[i].val);
      sparse_block[i].val   = tmp;
      sparse_block[i].sz    = bwidth;
      sparse_block[i].dense = 1;
    }
  //}
}
// if block is completely empty remove memory
if (rows_empty == bheight/2) {
  free(sparse_block);
  sparse_block  = NULL;
}
*sparse_block_in  = sparse_block;
}


ri_t elim_fl_D_fflas_ffpack(sbm_fl_t *D_old, mod_t modulus, int nthrds) {
  
  // generate DNS matrix D out of D_old
  SAFE_MALLOC_DECL(D, 1, DNS);
  initDenseUnit(D);
  copyMetaData(D, D_old, modulus);

  // dense representation of D, alloc memory at once and set entries to zero
  SAFE_CALLOC(D->ptr, (index_t)D->row * (index_t)D->ld, elemt_t);

  // copy multiline block matrix to DNS format
  copy_block_ml_matrix_to_dns_matrix(&D_old, &D);

  // row reduce D with FFLAS-FFPACK
  ri_t rank;
#if GBLA_WITH_FFLAS
  rank  = Mjoin(RowReduce,elemt_t)(D->mod,D->ptr,D->row,D->col,D->ld, nthrds); 
#endif
  /*
  size_t i, j, k;
  ri_t ctr = 0;
  for (i=0; i<D->row; ++i) {
    ctr = 0;
    printf("\n%d || %d\n",i,ctr/256);
    for (j=0; j<D->col; ++j) {
      printf("%.1f ",D->ptr[i*D->ld+j]);
      ctr++;
      if (ctr % 256 == 0)
        printf("\n%d || %d\n",i,ctr/256);
    }
  }
  */

  return rank;
}

ri_t elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, mod_t modulus, int nthrds) {

  ri_t i;
  ci_t rc;

  // row indices for subdividing echelonization parts in D_red
  global_next_row_to_reduce  = nthrds * 2;
  global_last_piv            = global_next_row_to_reduce - 1;

  // meta data for the computation of the rank of D_red at the end
  int head_line_1       = -1;
  int head_line_2       = -1;
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;
  const ci_t coldim     = D->ncols;
  ri_t wl_dim; // waiting list dimension
  re_t h_a1;

  wl_dim  = (D->nrows/__GB_NROWS_MULTILINE);
  if (D->nrows%__GB_NROWS_MULTILINE)
    wl_dim++;

  // global waiting list
  waiting_global.list = (wle_t *)malloc(wl_dim * sizeof(wle_t));
  waiting_global.sidx = 0;
  waiting_global.slp  = 0;
  waiting_global.sz   = 0;

  // copy D to D_red and delete D
  D_red = copy_block_matrix_to_multiline_matrix(&D, D_red, 1, nthrds);
#if DDEBUG_DD
  printf("BEFORE\n");
  const uint32_t rlD  = (uint32_t) ceil((float)D_red->nrows / __GB_NROWS_MULTILINE);
  int ii,jj,kk,ll;
  for (ii=0; ii<rlD; ++ii) {
    printf("%d .. \n",ii);
    //printf("size %d\n", D_red->ml[ii].sz);
    if (D_red->ml[ii].sz>0) {
      for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
        if (D_red->ml[ii].idx != NULL)
          printf("%d -- ", D_red->ml[ii].idx[ll]);
        else
          printf("%d -- ", ll);
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
      printf("\n");
    } else {
      printf("ml %d is zero! %p\n", ii, D_red->ml[ii]);
    }
  }
#endif
  echelonize_rows_sequential(D_red, 0, global_last_piv, modulus);
#if DDEBUG_DD
  printf("AFTER\n");
  for (int ii=0; ii<rlD; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", D_red->ml[ii].sz);
    if (D_red->ml[ii].sz>0) {
      for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
        /*
           if (D_red->ml[ii].idx != NULL)
           printf("%d -- ", D_red->ml[ii].idx[ll]);
           else
           printf("%d -- ", ll);
           */
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
    } else {
    }
    printf("\n");
  }
#endif

  const ri_t ml_nrows_D_red = (D_red->nrows % __GB_NROWS_MULTILINE == 0) ?
    D_red->nrows / __GB_NROWS_MULTILINE :
    D_red->nrows / __GB_NROWS_MULTILINE + 1;

  // define lock
  //omp_lock_t echelonize_lock;
  omp_init_lock(&echelonize_lock);

  // if there are rows left do elimination in parallel
  if (ml_nrows_D_red >= global_next_row_to_reduce) {
    // TODO: parallel elimination with OpenMP
#pragma omp parallel shared(D_red, waiting_global, global_next_row_to_reduce, global_last_piv) num_threads(nthrds)
    {
#pragma omp for nowait
      for (i=0; i<nthrds; ++i) {
        rc  = echelonize_rows_task(D_red, ml_nrows_D_red,
            //global_next_row_to_reduce, global_last_piv,
            //&waiting_global,
            modulus
            //, echelonize_lock
            );
      }
    }
  }
  omp_destroy_lock(&echelonize_lock);

  ri_t rank = 0;

  for (i=0; i<ml_nrows_D_red; ++i) {
    if (D_red->ml[i].sz == 0) {
      continue;
    }


    head_line_1 = get_head_multiline_hybrid(&(D_red->ml[i]), 0,
        &h_a1, &head_line_1_idx, coldim);
    head_line_2 = get_head_multiline_hybrid(&(D_red->ml[i]), 1,
        &h_a1, &head_line_2_idx, coldim);

    if (head_line_1 != -1)
      rank++;
    if (head_line_2 != -1)
      rank++;
  }
#if DDEBUG_DONE
  const uint32_t rlD  = (uint32_t) ceil((float)D_red->nrows / __GB_NROWS_MULTILINE);
  int ii,jj,kk,ll;
  printf("AFTER ALL\n");
  for (ii=0; ii<rlD; ++ii) {
    printf("%d .. \n",ii);
    printf("size %d\n", D_red->ml[ii].sz * 2);
    if (D_red->ml[ii].sz>0) {
      for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
        /*
           if (D_red->ml[ii].idx != NULL)
           printf("%d -- ", D_red->ml[ii].idx[ll]);
           else
           printf("%d -- ", ll);
           */
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
      printf("\n");
    }
  }
#endif

  return rank;
}

ri_t echelonize_rows_sequential(sm_fl_ml_t *A, const ri_t from, const ri_t to,
    const mod_t modulus) {

  if (A->nrows == 0)
    return 0;

  ri_t npiv_real  = 0;
  ri_t N          = A->nrows / __GB_NROWS_MULTILINE +
                    A->nrows % __GB_NROWS_MULTILINE;
  const ci_t coldim = A->ncols;

  ml_t *ml_row;
  re_l_t *dense_array_1, *dense_array_2;
  posix_memalign((void **)&dense_array_1, 16, coldim * sizeof(re_l_t));
  posix_memalign((void **)&dense_array_2, 16, coldim * sizeof(re_l_t));

  ri_t i;
  ci_t j, k;
  int head_line_1       = -1;
  int head_line_2       = -1;
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;

  re_t intermediate_val;
  uint32_t tmp_val;

  normalize_multiline(&A->ml[from], coldim, modulus);

  ri_t min_loop = to > N-1 ? N-1 : to;
  for (i=from; i<=min_loop; ++i) {
    if (A->ml[i].val != NULL) {
      memset(dense_array_1, 0, coldim * sizeof(re_l_t));
      memset(dense_array_2, 0, coldim * sizeof(re_l_t));
      copy_multiline_to_dense_array(A->ml[i], dense_array_1, dense_array_2, coldim);
#if DDEBUG_D
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("3-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("4-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
#endif
      re_t h_a1 = 1, h_a2 = 1;

      re_t v1_col1 = 0, v2_col1 = 0, v1_col2 = 0, v2_col2 = 0;
      ri_t tmp  = 0;

      for (j=from; j<i; ++j) {
        ml_row  = &(A->ml[j]);
        if (ml_row->sz == 0)
          continue;

        head_line_1 = get_head_multiline_hybrid(ml_row, 0,
            &h_a1, &head_line_1_idx, coldim);
        head_line_2 = get_head_multiline_hybrid(ml_row, 1,
            &h_a2, &head_line_2_idx, coldim);

        if (head_line_1 != -1) {
          v1_col1 = dense_array_1[head_line_1] % modulus;
          v2_col1 = dense_array_2[head_line_1] % modulus;

          if (v1_col1 != 0)
            v1_col1 = modulus - v1_col1;
          if (v2_col1 != 0)
            v2_col1 = modulus - v2_col1;
        } else {
          v1_col1 = 0;
          v2_col1 = 0;
        }
#if DDEBUG_D_ONE
        printf("v11 %d\n",v1_col1);
        printf("v21 %d\n",v2_col1);
#endif
        if (head_line_2 != -1) {
          v1_col2 = dense_array_1[head_line_2] % modulus;
          v2_col2 = dense_array_2[head_line_2] % modulus;

          intermediate_val  = ml_row->val[2*head_line_2_idx];
          tmp_val = v1_col2 + (uint32_t)v1_col1 * intermediate_val;
          v1_col2 = tmp_val % modulus;
          tmp_val = v2_col2 + (uint32_t)v2_col1 * intermediate_val;
          v2_col2 = tmp_val % modulus;

          if (v1_col2 != 0)
            v1_col2 = modulus - v1_col2;
          if (v2_col2 != 0)
            v2_col2 = modulus - v2_col2;
        } else {
          v1_col2 = 0;
          v2_col2 = 0;
        }
#if DDEBUG_D_ONE
        printf("v12 %d\n",v1_col2);
        printf("v22 %d\n",v2_col2);
#endif
        if (ml_row->sz < coldim) {
          sparse_scal_mul_sub_2_rows_vect_array_multiline(
              v1_col1, v2_col1,
              v1_col2, v2_col2,
              *ml_row,
              dense_array_1,
              dense_array_2);
        } else {
          dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
              v1_col1, v2_col1,
              v1_col2, v2_col2,
              *ml_row,
              dense_array_1,
              dense_array_2,
              head_line_1,
              head_line_2);
        }
      }
#if DDEBUG_D
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("5-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("6-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
#endif
      // normalize dense arrays
      head_line_1 = normalize_dense_array(dense_array_1, coldim, modulus);
      head_line_2 = normalize_dense_array(dense_array_2, coldim, modulus);
#if DDEBUG_D
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("7-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("8-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
#endif
      // reduce by same multiline
      if (head_line_1 >= head_line_2 && head_line_1 != -1 && head_line_2 != -1) {
        dense_array_2[head_line_1] %= modulus;
        if (dense_array_2[head_line_1] != 0) {
          register uint32_t h = modulus - dense_array_2[head_line_1];
          register uint32_t v__;
          for (k=head_line_1; k<coldim; ++k) {
            v__               =   dense_array_1[k] & 0x000000000000ffff;
            dense_array_2[k]  +=  h * v__;
          }
        }
      }

      head_line_2 = get_head_dense_array(dense_array_2, &h_a2, coldim, modulus);

      // make A->ml[i] dense, i.e. free memory for idx and enlarge memory for val
      // initialize all values in val with 0, set size of A->ml[i] to coldim
      if (A->ml[i].dense == 0) {
        free (A->ml[i].idx);
        A->ml[i].idx  = NULL;
        A->ml[i].val  = realloc(A->ml[i].val, 2 * coldim * sizeof(re_t));
        A->ml[i].sz   = coldim;
      }
      memset(A->ml[i].val, 0, 2 * coldim * sizeof(re_t));

#if DDEBUG_Dd
      printf("DENSE1\n");
      for (int uu=0; uu<coldim; ++uu)
        printf("%u :: ",dense_array_1[uu]);
      printf("DENSE2\n");
      for (int uu=0; uu<coldim; ++uu)
        printf("%u :: ",dense_array_2[uu]);
#endif

      // save the line with the smallest column entry first
      if (head_line_1 == -1) {
        copy_dense_array_to_zero_dense_multiline(
            dense_array_2, head_line_2, &A->ml[i], coldim, modulus);
      } else {
        if (head_line_2 == -1) {
          copy_dense_array_to_zero_dense_multiline(
              dense_array_1, head_line_1, &A->ml[i], coldim, modulus);
        } else { // both are not empty
          if (head_line_1 > head_line_2) {
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_2, dense_array_1, head_line_2, &A->ml[i], coldim, modulus);
          } else {
            //copy_dense_arrays_to_dense_multiline(
            //    dense_array_1, dense_array_2, ml, coldim, modulus);
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_1, dense_array_2, head_line_1, &A->ml[i], coldim, modulus);
          }
        }
      }
#if DDEBUG_D
      printf("BEFORE NORMALIZE\n");
      const uint32_t rlD  = (uint32_t) ceil((float)A->nrows / __GB_NROWS_MULTILINE);
      int ii,jj,kk,ll;
      for (ii=0; ii<rlD; ++ii) {
        printf("%d .. \n",ii);
        printf("size %d\n", A->ml[ii].sz * 2);
        if (A->ml[ii].sz>0) {
          for (ll=0; ll<A->ml[ii].sz; ++ll) {
            /*
               if (D_red->ml[ii].idx != NULL)
               printf("%d -- ", D_red->ml[ii].idx[ll]);
               else
               printf("%d -- ", ll);
               */
            printf("%d %d ", A->ml[ii].val[2*ll], A->ml[ii].val[2*ll+1]);
          }
          printf("\n");
        }
      }
#endif
      // normalize multiline
      normalize_multiline(&A->ml[i], coldim, modulus);
      if (head_line_1 != -1)
        npiv_real++;
      if (head_line_2 != -1)
        npiv_real++;
    }
  }
  free(dense_array_1);
  dense_array_1 = NULL;
  free(dense_array_2);
  dense_array_2 = NULL;

  return npiv_real;
}

int echelonize_rows_task(sm_fl_ml_t *A, const ri_t N,
    //ri_t global_next_row_to_reduce, ri_t global_last_piv,
    //wl_t *waiting_global,
    const mod_t modulus
    //, omp_lock_t echelonize_lock
    ) {

  const ci_t coldim = A->ncols;
  ri_t curr_row_to_reduce;
  ri_t local_last_piv;

  int ready_for_waiting_list = 0;
  int curr_row_fully_reduced = 0;

  ri_t from_row;
  int nreduced_consecutively  = 0;
  // local waiting list entries
  ri_t wl_idx  = 0;
  ri_t wl_lp   = 0;

  re_l_t *dense_array_1, *dense_array_2;
  posix_memalign((void **)&dense_array_1, 16, coldim * sizeof(re_l_t));
  posix_memalign((void **)&dense_array_2, 16, coldim * sizeof(re_l_t));
  int tid = omp_get_thread_num();

  // do the computations
  while (1) {
    local_last_piv = global_last_piv;
#if DEBUG_ECHELONIZE
    printf("seeting llp\n");
    printf("llp %d\n",local_last_piv);
    printf("grr %d\n",global_next_row_to_reduce);
#endif
    if (global_last_piv >= N) {
      omp_unset_lock(&echelonize_lock);
      break;
    }
    omp_set_lock(&echelonize_lock);
    if ((ready_for_waiting_list == 0) && (global_next_row_to_reduce < N)) {
      curr_row_to_reduce = global_next_row_to_reduce;
#if DEBUG_ECHELONIZE
          printf("%d -- inc grr\n",tid);
#endif
      global_next_row_to_reduce++;
      from_row = 0;
    } else {
      ready_for_waiting_list = 1;
    }
    if (ready_for_waiting_list == 1) {
      if (get_smallest_waiting_row(&waiting_global, &wl_idx, &wl_lp) == 0) { // no waiting rows
        if (global_next_row_to_reduce >= N) { // no more rows to reduce
          omp_unset_lock(&echelonize_lock);
          break;
        } else {
          if (local_last_piv >= N) { // we are also done
            omp_unset_lock(&echelonize_lock);
            break;
          } else {
            ready_for_waiting_list  = 0;
            omp_unset_lock(&echelonize_lock);
            continue;
          }
        }
      }
      from_row            = wl_lp + 1;
      curr_row_to_reduce  = wl_idx;
#if DEBUG_ECHELONIZE
      printf("(%d) from row %d -- crr %d -- gllp %d\n", tid,from_row, curr_row_to_reduce, global_last_piv);
#endif
    }

    omp_unset_lock(&echelonize_lock);
    if (A->ml[curr_row_to_reduce].val != NULL) {
    // set zero
    memset(dense_array_1, 0, coldim * sizeof(re_l_t));
    memset(dense_array_2, 0, coldim * sizeof(re_l_t));

    copy_multiline_to_dense_array(A->ml[curr_row_to_reduce],
        dense_array_1, dense_array_2, coldim);
    // echelonize one row
#if DEBUG_ECHELONIZE
    printf("thread %d reduces %d with rows %d -- %d\n",tid, curr_row_to_reduce, from_row, local_last_piv);
#endif
    echelonize_one_row(A, dense_array_1, dense_array_2, from_row,
        local_last_piv, modulus); // TODO
#if DEBUG_ECHELONIZE
    printf("thread %d done with rows %d -- %d\n",tid, from_row, local_last_piv);
#endif
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("13-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("14-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    }
    if (curr_row_to_reduce == local_last_piv + 1) {
      curr_row_fully_reduced  = 1;
      ready_for_waiting_list  = 0;
    } else {
      curr_row_fully_reduced  = 0;
      ready_for_waiting_list  = 0;
      ++nreduced_consecutively;
      if (nreduced_consecutively > 4) {
        nreduced_consecutively  = 0;
        ready_for_waiting_list  = 1;
      }
    }

    // save back to multiline
    /*
    printf("1 %p -- 2 %p\n",dense_array_1,dense_array_2);
    printf("BEFORE A.val %p\n",A->ml[curr_row_to_reduce].val);
    */
    if (A->ml[curr_row_to_reduce].val != NULL) {
      save_back_and_reduce(&(A->ml[curr_row_to_reduce]), dense_array_1,
          dense_array_2, coldim, modulus, curr_row_fully_reduced, curr_row_to_reduce);
    }
    /*
    printf("AFTER A.val %p\n",A->ml[curr_row_to_reduce].val);
    printf("1 %p -- 2 %p\n",dense_array_1,dense_array_2);
    printf("3-mlidx %p -- mlval %p\n",A->ml[curr_row_to_reduce].idx,A->ml[curr_row_to_reduce].val);
      printf("MULTILINE AFTER SAVE\n");
      for (int kk=0; kk<A->ml[curr_row_to_reduce].sz; ++kk) {
        printf("%d : %u -- %u ",kk,A->ml[curr_row_to_reduce].val[2*kk],A->ml[curr_row_to_reduce].val[2*kk+1]);
      }
      printf("\n");
      */

    if (curr_row_fully_reduced == 1) {
      omp_set_lock(&echelonize_lock);
      ++global_last_piv;
#if DEBUG_ECHELONIZE
      printf("%d -- inc glp: %d\n", tid, global_last_piv);
#endif
      omp_unset_lock(&echelonize_lock);
    } else {
      omp_set_lock(&echelonize_lock);
#if DEBUG_ECHELONIZE
      printf("%d -- pushes %d / %d\n",tid, curr_row_to_reduce, local_last_piv);
#endif
      push_row_to_waiting_list(&waiting_global, curr_row_to_reduce, local_last_piv);
      omp_unset_lock(&echelonize_lock);
    }
  }
  // free memory
  free(dense_array_1);
  dense_array_1 = NULL;
  free(dense_array_2);
  dense_array_2 = NULL;

  return 0;
}

void echelonize_one_row(sm_fl_ml_t *A,
    re_l_t *dense_array_1, re_l_t *dense_array_2,
    const ri_t first_piv, const ri_t last_piv,
    const mod_t modulus) {

  ci_t j;
  int head_line_1       = -1;
  int head_line_2       = -1;
  int nzc               = 0; // tracks if there is a nonzero coeff at all
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;

  re_t intermediate_val;
  uint32_t tmp_val;
  const ci_t coldim = A->ncols;

  ml_t *ml_row;

  re_t h_a1 = 1, h_a2 = 1;

  re_t v1_col1 = 0, v2_col1 = 0, v1_col2 = 0, v2_col2 = 0;
  ri_t tmp  = 0;

#if DEBUG_ECHELONIZE
  printf("fp %d -- lp %d\n",first_piv, last_piv);
#endif
  for (j=first_piv; j<last_piv+1; ++j) {
    nzc = 0;
    ml_row  = &(A->ml[j]);
#if DEBUG_ECHELONIZE
    printf("j %d\n",j);
#endif
    if (ml_row->val == NULL || ml_row->sz == 0)
      continue;
#if DDEBUG_D
    if (ml_row != NULL) {
      printf("MULTILINE\n");
      for (int kk=0; kk<ml_row->sz; ++kk) {
        printf("%d : %u -- %u ",kk,ml_row->val[2*kk],ml_row->val[2*kk+1]);
      }
      printf("\n");
    }
#endif


    head_line_1 = get_head_multiline_hybrid(ml_row, 0,
        &h_a1, &head_line_1_idx, coldim);
    head_line_2 = get_head_multiline_hybrid(ml_row, 1,
        &h_a2, &head_line_2_idx, coldim);

#if DEBUG_ECHELONIZE
    printf("hl1 %d || h2 %d\n",head_line_1,head_line_2);
#endif

    if (head_line_1 != -1) {
#if DDEBUG_D_ONE
      printf("d11 %lu (%p)\n",dense_array_1[head_line_1],dense_array_1[head_line_1]);
      printf("d21 %lu\n",dense_array_2[head_line_1]);
      printf("d11 %lu\n",dense_array_1[head_line_1] % modulus);
      printf("d21 %lu\n",dense_array_2[head_line_1] % modulus);
#endif
      v1_col1 = dense_array_1[head_line_1] % modulus;
      v2_col1 = dense_array_2[head_line_1] % modulus;
#if DDEBUG_D_ONE
      printf("v11 %d\n",v1_col1);
      printf("v21 %d\n",v2_col1);
#endif

      if (v1_col1 != 0) {
        v1_col1 = modulus - v1_col1;
        nzc = 1;
      }
      if (v2_col1 != 0) {
        v2_col1 = modulus - v2_col1;
        nzc = 1;
      }
    } else {
      v1_col1 = 0;
      v2_col1 = 0;
    }
#if DDEBUG_D_ONE
    printf("v11 %d\n",v1_col1);
    printf("v21 %d\n",v2_col1);
#endif
    //printf("hl2 %d => %d ?\n",head_line_2,head_line_2 != -1);
    if (head_line_2 != -1) {
      v1_col2 = dense_array_1[head_line_2] % modulus;
      v2_col2 = dense_array_2[head_line_2] % modulus;
#if DDEBUG_D_ONE
      printf("v12 %d\n",v1_col2);
      printf("v22 %d\n",v2_col2);
#endif

      intermediate_val  = ml_row->val[2*head_line_2_idx];
      tmp_val = v1_col2 + (uint32_t)v1_col1 * intermediate_val;
      v1_col2 = tmp_val % modulus;
      tmp_val = v2_col2 + (uint32_t)v2_col1 * intermediate_val;
      v2_col2 = tmp_val % modulus;
#if DDEBUG_D_ONE
      printf("v12 %d\n",v1_col2);
      printf("v22 %d\n",v2_col2);
#endif

      if (v1_col2 != 0) {
        v1_col2 = modulus - v1_col2;
        nzc = 1;
      }
      if (v2_col2 != 0) {
        v2_col2 = modulus - v2_col2;
        nzc = 1;
      }
    } else {
      v1_col2 = 0;
      v2_col2 = 0;
    }
#if DDEBUG_D_ONE
    printf("v12 %d\n",v1_col2);
    printf("v22 %d\n",v2_col2);
#endif
    if (nzc == 1) {
      if (ml_row->sz < coldim) {
        sparse_scal_mul_sub_2_rows_vect_array_multiline(
            v1_col1, v2_col1,
            v1_col2, v2_col2,
            *ml_row,
            dense_array_1,
            dense_array_2);
      } else {
        dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
            v1_col1, v2_col1,
            v1_col2, v2_col2,
            *ml_row,
            dense_array_1,
            dense_array_2,
            head_line_1,
            head_line_2);
      }
    }
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("15-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("16-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
  }
}

void save_back_and_reduce(ml_t *ml, re_l_t *dense_array_1,
    re_l_t *dense_array_2, const ci_t coldim, const mod_t modulus,
    const int reduce, const ri_t curr_row_to_reduce) {


  if (reduce == 1) {
    int head_line_1 = -1;
    int head_line_2 = -1;
    re_t h_a1 = 1, h_a2 = 1;
    ci_t k;
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("17-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("18-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    // normalize dense arrays
    head_line_1 = normalize_dense_array(dense_array_1, coldim, modulus);
    head_line_2 = normalize_dense_array(dense_array_2, coldim, modulus);
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("19-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("20-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    // reduce by same multiline
    if (head_line_1 >= head_line_2 && head_line_1 != -1 && head_line_2 != -1) {
      dense_array_2[head_line_1] %= modulus;
      if (dense_array_2[head_line_1] != 0) {
        register uint32_t h = modulus - dense_array_2[head_line_1];
        register uint32_t v__;
        for (k=head_line_1; k<coldim; ++k) {
          v__               =   dense_array_1[k] & 0x000000000000ffff;
          dense_array_2[k]  +=  h * v__;
        }
      }
    }
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("21-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("22-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
    head_line_2 = get_head_dense_array(dense_array_2, &h_a2, coldim, modulus);

#if DEBUG_ECHELONIZE
    printf("headl1 %d | headl2 %d\n",head_line_1,head_line_2);
    printf("mlidx %p -- mlval %p : sz %d\n",ml->idx,ml->val, ml->sz);
#endif
    // possibly we have a zero row
    if (head_line_1 == -1 && head_line_2 == -1) {
      free(ml->val);
      ml->val  = NULL;
      ml->sz   = 0;
    } else {
      memset(ml->val, 0, 2 * coldim * sizeof(re_t));

      if (head_line_1 == -1) {
        copy_dense_array_to_zero_dense_multiline(
            dense_array_2, head_line_2, ml, coldim, modulus);
      } else {
        if (head_line_2 == -1) {
          copy_dense_array_to_zero_dense_multiline(
              dense_array_1, head_line_1, ml, coldim, modulus);
        } else { // both are not empty
          if (head_line_1 > head_line_2) {
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_2, dense_array_1, head_line_2, ml, coldim, modulus);
          } else {
            copy_dense_arrays_to_zero_dense_multiline(
                dense_array_1, dense_array_2, head_line_1, ml, coldim, modulus);
          }
        }
      }
      // normalize multiline
      normalize_multiline(ml, coldim, modulus);
    }
  } else { // do not reduce
    // make A->ml[i] dense, i.e. free memory for idx and enlarge memory for val
    // initialize all values in val with 0, set size of A->ml[i] to coldim
    //free(ml->idx);
    //ml->idx  = NULL;
    //ml->sz   = coldim;
    memset(ml->val, 0, 2 * coldim * sizeof(re_t));
    //copy_dense_arrays_to_dense_multiline(
       // dense_array_1, dense_array_2, ml, coldim, modulus);
    copy_dense_arrays_to_zero_dense_multiline(
        dense_array_1, dense_array_2, 0, ml, coldim, modulus);
  }
}


/******************************************************************************
 * MULTILINE VARIANT OF ELIMINATION PROCESS
 *****************************************************************************/

int elim_fl_C_ml(sm_fl_ml_t *C, sm_fl_ml_t *A, mod_t modulus, int nthrds) {
  ci_t i;
  ri_t j, k, rc;

  const ri_t rlC  = (ri_t) ceil((float) C->nrows / __GB_NROWS_MULTILINE);
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for
    // each task takes one block column of B
    for (i=0; i<rlC; ++i) {
#pragma omp task
      {
        rc  = elim_fl_C_ml_task(C, A, i, modulus);
      }
    }
#pragma omp taskwait
  }

  return 0;
}

int elim_fl_C_ml_task(sm_fl_ml_t *C, sm_fl_ml_t *A, ri_t row_idx, mod_t modulus) {

  ci_t i;
  const ci_t coldim = C->ncols;

  ci_t start_idx;
  // Note that all our multilines are stored in a sparse fashion at the moment.
  // For compatibility and later changes we keep this check in the code.
  if (C->ml[row_idx].dense == 0)
    start_idx = C->ml[row_idx].idx[0];
  else
    start_idx = 0;

  uint32_t Cv1_col1, Cv2_col1;
  uint32_t Cv1_col2, Cv2_col2;
  uint32_t Cp1, Cp2;
  uint32_t tmp  = 0;
  ri_t row_in_A_mod, row_in_A;

  re_l_t *dense_array_C_1, *dense_array_C_2;
  posix_memalign((void **)&dense_array_C_1, 16, coldim * sizeof(re_l_t));
  posix_memalign((void **)&dense_array_C_2, 16, coldim * sizeof(re_l_t));

  memset(dense_array_C_1, 0, coldim * sizeof(re_l_t));
  memset(dense_array_C_2, 0, coldim * sizeof(re_l_t));

  copy_multiline_to_dense_array(C->ml[row_idx], dense_array_C_1,
      dense_array_C_2, coldim);

  for (i=start_idx; i<coldim; ++i) {
    Cp1 = i;
    Cv1_col1  = dense_array_C_1[i] % modulus;
    Cv2_col1  = dense_array_C_2[i] % modulus;

    if (Cv1_col1 == 0 && Cv2_col1 == 0)
      continue;

    if (Cv1_col1 != 0)
      Cv1_col1  = (uint32_t)modulus - Cv1_col1;
    if (Cv2_col1 != 0)
      Cv2_col1  = (uint32_t)modulus - Cv2_col1;

    row_in_A      = (A->nrows - 1 - Cp1) / __GB_NROWS_MULTILINE;
    row_in_A_mod  = (A->nrows - 1 - Cp1) % __GB_NROWS_MULTILINE;

    if (row_in_A_mod == 1) {
      Cp2 = i+1;
      Cv1_col2  = dense_array_C_1[i+1] % modulus;
      Cv2_col2  = dense_array_C_2[i+1] % modulus;

      register re_t v__ = A->ml[row_in_A].val[2*1+1];

      tmp = Cv1_col2 + (uint32_t)Cv1_col1 * v__;
      Cv1_col2 = tmp % modulus;
      tmp = Cv2_col2 + (uint32_t)Cv2_col1 * v__;
      Cv2_col2 = tmp % modulus;

      if (Cv1_col2 != 0)
        Cv1_col2  = (uint32_t)modulus - Cv1_col2;
      if (Cv2_col2 != 0)
        Cv2_col2  = (uint32_t)modulus - Cv2_col2;

      sparse_scal_mul_sub_2_rows_vect_array_multiline(
          Cv1_col2, Cv2_col2,
          Cv1_col1, Cv2_col1,
          A->ml[row_in_A],
          dense_array_C_1,
          dense_array_C_2);

      dense_array_C_1[Cp1]  = Cv1_col1;
      dense_array_C_1[Cp2]  = Cv1_col2;

      dense_array_C_2[Cp1]  = Cv2_col1;
      dense_array_C_2[Cp2]  = Cv2_col2;

      // since we have done two-in-one
      ++i;
    } else {
      sparse_scal_mul_sub_1_row_vect_array_multiline(
          Cv1_col1, Cv2_col1,
          A->ml[row_in_A],
          row_in_A_mod,
          dense_array_C_1,
          dense_array_C_2);

      dense_array_C_1[Cp1]  = Cv1_col1;
      dense_array_C_2[Cp1]  = Cv2_col1;
    }
  }

  ml_t *ml = &C->ml[row_idx];
#if DENSE_MULTILINE_C_AFTER_ELIM
  // make multiline in C dense
  free(ml->idx);
  ml->idx   = NULL;
  ml->dense = 1;
  ml->sz    = coldim;
  ml->val   = realloc(ml->val, 2 * coldim * sizeof(re_t));
  memset(ml->val, 0, 2 * coldim * sizeof(re_t));

  copy_dense_arrays_to_zero_dense_multiline(
      dense_array_C_1, dense_array_C_2, 0, ml, coldim, modulus);
#else
  copy_dense_arrays_to_multiline(
      dense_array_C_1, dense_array_C_2, 0, ml, coldim, modulus);
#endif

  free(dense_array_C_1);
  free(dense_array_C_2);

  return 0;
}
