#include "elimination.h"

#define DDEBUG 0
#define DDEBUG_C 0
#define DDEBUG_D 0
#define DDEBUG_DD 0
#define DDEBUG_DONE 1
#define DDEBUG_D_ONE 0

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
          if (A->blocks[j][i][k].dense == 0) {
            if (A->blocks[j][i][k].idx != NULL) {
              free(A->blocks[j][i][k].idx);
              A->blocks[j][i][k].idx  = NULL;
            }
          }
          if (A->blocks[j][i][k].val != NULL) {
            free(A->blocks[j][i][k].val);
            A->blocks[j][i][k].val  = NULL;
          }
        }
        if (A->blocks[j][i] != NULL) {
          free(A->blocks[j][i]);
          A->blocks[j][i] = NULL;
        }
      }
    }
    if (A->blocks[j] != NULL) {
      free(A->blocks[j]);
      A->blocks[j]  = NULL;
    }
  }
}
free(A->blocks);
A->blocks = NULL;
free(A);
A = NULL;
*A_in  = A;
return 0;
}

int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C_in, sbm_fl_t *D, mod_t modulus, int nthrds) {

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
        rc  = elim_fl_C_blocks_task(B, C, D, i, rlC, clC, modulus);
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
          if (C->blocks[j][i][k].dense == 0) {
            if (C->blocks[j][i][k].idx != NULL) {
              free(C->blocks[j][i][k].idx);
              C->blocks[j][i][k].idx  = NULL;
            }
          }
          if (C->blocks[j][i][k].val != NULL) {
            free(C->blocks[j][i][k].val);
            C->blocks[j][i][k].val  = NULL;
          }
        }
        if (C->blocks[j][i] != NULL) {
          free(C->blocks[j][i]);
          C->blocks[j][i] = NULL;
        }
      }
    }
    if (C->blocks[j] != NULL) {
      free(C->blocks[j]);
      C->blocks[j]  = NULL;
    }
  }
}
free(C->blocks);
C->blocks = NULL;
free(C);
C = NULL;
*C_in = C;
return 0;
}

int elim_fl_A_ml_block(sm_fl_ml_t *A, sbm_fl_t *B, int nthrds) {
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
  // TODO: Martani uses a minimum possibly for zero blocks at the end
  const ri_t last_block_idx   = j;

  // set dense block entries to zero
  for (k=0; k<B->bheight; ++k)
    memset(dense_block[k], 0, size);

  // copy sparse block data to dense representation
  if (B->blocks[j][block_col_idx_B] != NULL)
    copy_sparse_to_dense_block(B->blocks[j][block_col_idx_B], dense_block,
        B->bheight, B->bwidth);

  for (k=0; k<last_block_idx; ++k) {
#if DDDEBUG
    printf("j %d -- k %d -- lci %d\n", j, k, block_col_idx_B);
#endif
    red_with_rectangular_block(A->blocks[j][k], B->blocks[k][block_col_idx_B],
        dense_block, B->bheight, modulus, 1);
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


#if DDDEBUG
  printf("j %d -- lbi %d\n", j, last_block_idx);
#endif
  red_with_triangular_block(A->blocks[j][last_block_idx], dense_block,
      B->bheight, modulus, 1);
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

int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
  ci_t block_col_idx_D, ri_t nbrows_C, ci_t nbcols_C, mod_t modulus) {
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
          dense_block, B->bheight, modulus, 1);
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

void red_with_triangular_block(mbl_t *block_A, re_l_t **dense_block, ri_t bheight, mod_t modulus, int inv_scalars) {
bi_t i, j, k;

bi_t last_idx;

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

    if (inv_scalars) {
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

        if (inv_scalars) {
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
      if (inv_scalars)
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

void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_block, ri_t bheight, mod_t modulus, int inv_scalars) {
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

    if (inv_scalars) {
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

        if (inv_scalars) {
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

ri_t elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, mod_t modulus, int nthrds) {

ri_t i;
ci_t rc;

// row indices for subdividing echelonization parts in D_red
ri_t global_next_row_to_reduce  = nthrds * 2;
ri_t global_last_piv            = global_next_row_to_reduce - 1;

// meta data for the computation of the rank of D_red at the end
ri_t rank             = 0;
int head_line_1       = -1;
int head_line_2       = -1;
ci_t head_line_1_idx  = 0;
ci_t head_line_2_idx  = 0;
const ci_t coldim     = D->ncols;
re_t h_a1;

wl_t waiting;
waiting.list  = (wle_t *)malloc(coldim * sizeof(wle_t));
waiting.sidx  = 0;
waiting.slp   = 0;
waiting.sz    = 0;

// copy D to D_red and delete D
printf("D %p\n",D);
D_red = copy_block_matrix_to_multiline_matrix(&D, D_red, 1, nthrds);
printf("D %p\n",D);
printf("BEFORE\n");
const uint32_t rlD  = (uint32_t) ceil((float)D_red->nrows / __GB_NROWS_MULTILINE);
int ii,jj,kk,ll;
for (ii=0; ii<rlD; ++ii) {
  printf("%d .. \n",ii);
  //printf("size %d\n", D_red->ml[ii].sz);
  if (D_red->ml[ii].sz>0) {
    for (ll=0; ll<D_red->ml[ii].sz; ++ll) {
      /*
      if (D_red->ml[ii].idx != NULL)
        printf("%d -- ", D_red->ml[ii].idx[ll]);
      else
        printf("%d -- ", ll);
        */
      /*
      if (ii==270) {
        printf("%u -- %u == 0?\n",
            D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
      }
      */
        printf("%d %d ", D_red->ml[ii].val[2*ll], D_red->ml[ii].val[2*ll+1]);
    }
    printf("\n");
  } else {
    //printf("ml %d is zero! %p\n", ii, D_red->ml[ii]);
  }
}
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

// if there are rows left do elimination in parallel
printf("MLNDRED %d / %d\n",ml_nrows_D_red, global_next_row_to_reduce);
if (ml_nrows_D_red >= global_next_row_to_reduce) {
  // TODO: parallel elimination with OpenMP
  #pragma omp parallel shared(waiting, global_next_row_to_reduce, global_last_piv) num_threads(nthrds)
  {
    #pragma omp for
      for (i=0; i<nthrds; ++i) {
        rc  = echelonize_rows_task(D_red, ml_nrows_D_red,
            global_next_row_to_reduce, global_last_piv,
            &waiting, modulus);
      }
    }
  }

  for (i=0; i<ml_nrows_D_red; ++i) {
    if (D_red->ml[i].sz == 0)
      continue;

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
  //const uint32_t rlD  = (uint32_t) ceil((float)D_red->nrows / __GB_NROWS_MULTILINE);
  //int ii,jj,kk,ll;
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
  printf("COLUMN DIM %u\n",coldim);
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
#if DDEBUG_D
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
        if (v2_col1 != 0)
          v2_col2 = modulus - v2_col2;
      } else {
        v1_col2 = 0;
        v2_col2 = 0;
      }
#if DDEBUG_D
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
            head_line_1);
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
    /*
    free (A->ml[i].idx);
    A->ml[i].idx  = NULL;
    A->ml[i].val  = realloc(A->ml[i].val, 2 * coldim * sizeof(re_t));
    A->ml[i].sz   = coldim;
    */
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
        printf("hl2 == -1\n");
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
  free(dense_array_1);
  dense_array_1 = NULL;
  free(dense_array_2);
  dense_array_2 = NULL;

  return npiv_real;
}

int echelonize_rows_task(sm_fl_ml_t *A, const ri_t N,
    ri_t global_next_row_to_reduce, ri_t global_last_piv,
    wl_t *waiting, const mod_t modulus) {

  const ci_t coldim = A->ncols;
  ri_t curr_row_to_reduce;
  ri_t local_last_piv;

  int ready_for_waiting_list = 0;
  int curr_row_fully_reduced = 0;

  ri_t from_row;
  int nreduced_consecutively  = 0;

  printf("COLUMN DIM %d\n",coldim);
  re_l_t *dense_array_1, *dense_array_2;
  posix_memalign((void **)&dense_array_1, 16, coldim * sizeof(re_l_t));
  posix_memalign((void **)&dense_array_2, 16, coldim * sizeof(re_l_t));
  printf("addresses after allocation %p -- %p\n",dense_array_1,dense_array_2);

  // do the computations
  while (1) {
    local_last_piv = global_last_piv;
    printf("llp %d\n",local_last_piv);
    printf("grr %d\n",global_next_row_to_reduce);
    if (global_last_piv >= N)
      break;

    if (!ready_for_waiting_list && (global_next_row_to_reduce < N)) {
      curr_row_to_reduce = global_next_row_to_reduce;
      ++global_next_row_to_reduce;
      from_row = 0;
    } else {
      ready_for_waiting_list = 1;
    }

    if (ready_for_waiting_list) {
      if (get_smallest_waiting_row(waiting) == 0) { // no waiting rows
        if (global_next_row_to_reduce >= N) { // no more rows to reduce
          break;
        } else {
          if (local_last_piv >= N) { // we are also done
            break;
          } else {
            ready_for_waiting_list  = 0;
            continue;
          }
        }
      }
      from_row            = waiting->slp + 1;
      curr_row_to_reduce  = waiting->sidx;
    }
    // set zero
    memset(dense_array_1, 0, coldim * sizeof(re_l_t));
    memset(dense_array_2, 0, coldim * sizeof(re_l_t));

    // copy multiline to dense representation
    printf("crtr %d 1st element %u (%p)\n",curr_row_to_reduce, A->ml[curr_row_to_reduce].val[0], A->ml[curr_row_to_reduce].val);
    copy_multiline_to_dense_array(A->ml[curr_row_to_reduce],
        dense_array_1, dense_array_2, coldim);
    printf("crtr %d 1st element %u mlval(%p)\n",curr_row_to_reduce, A->ml[curr_row_to_reduce].val[0], A->ml[curr_row_to_reduce].val);
    printf("crtr %d 1st element %u (%p)\n",curr_row_to_reduce, dense_array_1[0], dense_array_1);
    if (local_last_piv == 414) {
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("11-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
      }
      for (int kk = 0; kk< coldim/2; ++kk) {
        printf("12-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
      }
    }
    // echelonize one row
    echelonize_one_row(A, dense_array_1, dense_array_2, from_row,
        local_last_piv, modulus); // TODO
#if DDEBUG_D
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("13-%d ,, %lu || %lu\n",kk,dense_array_1[2*kk], dense_array_1[2*kk+1]);
    }
    for (int kk = 0; kk< coldim/2; ++kk) {
      printf("14-%d ,, %lu || %lu\n",kk,dense_array_2[2*kk], dense_array_2[2*kk+1]);
    }
#endif
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
    printf("0-mlidx %p -- mlval %p\n",A->ml[curr_row_to_reduce].idx,A->ml[curr_row_to_reduce].val);
    printf("1 %p -- 2 %p\n",dense_array_1,dense_array_2);
    printf("BEFORE A.val %p\n",A->ml[curr_row_to_reduce].val);
    save_back_and_reduce(&(A->ml[curr_row_to_reduce]), dense_array_1,
        dense_array_2, coldim, modulus, curr_row_fully_reduced);
    printf("AFTER A.val %p\n",A->ml[curr_row_to_reduce].val);
    printf("1 %p -- 2 %p\n",dense_array_1,dense_array_2);
    printf("3-mlidx %p -- mlval %p\n",A->ml[curr_row_to_reduce].idx,A->ml[curr_row_to_reduce].val);
      printf("MULTILINE AFTER SAVE\n");
      for (int kk=0; kk<A->ml[curr_row_to_reduce].sz; ++kk) {
        printf("%d : %u -- %u ",kk,A->ml[curr_row_to_reduce].val[2*kk],A->ml[curr_row_to_reduce].val[2*kk+1]);
      }
      printf("\n");

    if (curr_row_fully_reduced)
      ++global_last_piv;
    else
      push_row_to_waiting_list(waiting, curr_row_to_reduce, local_last_piv);

  }
  // free memory
  printf("addresses before freeing %p -- %p\n",dense_array_1,dense_array_2);
  free(dense_array_1);
  printf("1\n");
  dense_array_1 = NULL;
  printf("2\n");
  free(dense_array_2);
  printf("3\n");
  dense_array_2 = NULL;
  printf("4\n");

  return 0;
}

void echelonize_one_row(sm_fl_ml_t *A,
    re_l_t *dense_array_1, re_l_t *dense_array_2,
    const ri_t first_piv, const ri_t last_piv,
    const mod_t modulus) {

  ci_t j;
  int head_line_1       = -1;
  int head_line_2       = -1;
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;

  re_t intermediate_val;
  uint32_t tmp_val;
  const ci_t coldim = A->ncols;

  ml_t *ml_row;

  re_t h_a1 = 1, h_a2 = 1;

  re_t v1_col1 = 0, v2_col1 = 0, v1_col2 = 0, v2_col2 = 0;
  ri_t tmp  = 0;

  printf("fp %d -- lp %d\n",first_piv, last_piv);
  for (j=first_piv; j<last_piv+1; ++j) {
    ml_row  = &(A->ml[j]);
    printf("j %d\n",j);
    if (ml_row->sz == 0)
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

    printf("hl1 %d || h2 %d\n",head_line_1,head_line_2);

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
          head_line_1);
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
    const int reduce) {

printf("COLDIMSBAR %d\n",coldim);

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

    printf("headl1 %d | headl2 %d\n",head_line_1,head_line_2);
    printf("mlidx %p -- mlval %p : sz %d\n",ml->idx,ml->val, ml->sz);
    // possibly we have a zero row
    if (head_line_1 == -1 && head_line_2 == -1) {
      printf("mlv %p\n",ml->val);
      free(ml->val);
      printf("mlv %p\n",ml->val);
      ml->val  = NULL;
      printf("mlv %p\n",ml->val);
      ml->sz   = 0;
    } else {
      memset(ml->val, 0, 2 * coldim * sizeof(re_t));

      if (head_line_1 == -1) {
        copy_dense_array_to_zero_dense_multiline(
            dense_array_2, head_line_2, ml, coldim, modulus);
      } else {
        if (head_line_2 == -1) {
          printf("hl2 == -1\n");
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
      /*
      // save the line with the smallest column entry first
      if ((head_line_1 == -1) || ((head_line_2 != -1) && (head_line_1 > head_line_2))) {
    //copy_dense_arrays_to_dense_multiline(
    //    dense_array_2, dense_array_1, ml, coldim, modulus);
        copy_dense_arrays_to_zero_dense_multiline(
            dense_array_2, dense_array_1, ml, coldim, modulus);
      } else {
    //copy_dense_arrays_to_dense_multiline(
    //    dense_array_1, dense_array_2, ml, coldim, modulus);
        copy_dense_arrays_to_zero_dense_multiline(
            dense_array_1, dense_array_2, ml, coldim, modulus);
      }
      */
      /*
      printf("MULTILINE in SAVE BEFORE NORMALIZE\n");
      for (int kk=0; kk<ml->sz; ++kk) {
        printf("%d : %u -- %u ",kk,ml->val[2*kk],ml->val[2*kk+1]);
      }
      printf("\n");
      */
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
    printf("2-mlidx %p -- mlval %p\n",ml->idx,ml->val);
    /*
      printf("MULTILINE in SAVE\n");
      for (int kk=0; kk<ml->sz; ++kk) {
        printf("%d : %u -- %u ",kk,ml->val[2*kk],ml->val[2*kk+1]);
      }
      printf("\n");
      */
}
