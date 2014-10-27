#include "elimination.h"

#define DDEBUG 0
#define DDEBUG_C 0

int elim_fl_A_block(sbm_fl_t *A, sbm_fl_t *B, mod_t modulus, int nthrds) {
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
    #pragma omp for private(i,k)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i] != NULL) {
          for (k=0; k<A->bheight/__GB_NROWS_MULTILINE; ++k) {
            if (A->blocks[j][i][k].dense == 0)
              free(A->blocks[j][i][k].idx);
            free(A->blocks[j][i][k].val);
          }
          free(A->blocks[j][i]);
        }
      }
      free(A->blocks[j]);
    }
  }
  free(A);
  A = NULL;
  return 0;
}

int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D, mod_t modulus, int nthrds) {
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
    #pragma omp for private(i,k)
    for (j=0; j<rlC; ++j) {
      for (i=0; i<clC; ++i) {
        if (C->blocks[j][i] != NULL) {
          for (k=0; k<C->bheight/__GB_NROWS_MULTILINE; ++k) {
            if (C->blocks[j][i][k].dense == 0)
              free(C->blocks[j][i][k].idx);
            free(C->blocks[j][i][k].val);
          }
          free(C->blocks[j][i]);
        }
      }
      free(C->blocks[j]);
    }
  }
  free(C);
  C = NULL;
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
    B->blocks[j][block_col_idx_B] = copy_dense_block_to_sparse(
        dense_block, B->blocks[j][block_col_idx_B],
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
    D->blocks[j][block_col_idx_D] = copy_dense_block_to_sparse(
        dense_block, D->blocks[j][block_col_idx_D],
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
    bi_t Ap1,  Ap2;
    
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

mbl_t *copy_dense_block_to_sparse(
    re_l_t **dense_block, mbl_t *sparse_block, int bheight, int bwidth, mod_t modulus) {

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
          if (sparse_block[i].idx != NULL) {
            free(sparse_block[i].idx);
            sparse_block[i].idx = NULL;
          }
          if (sparse_block[i].val != NULL) {
            free(sparse_block[i].val);
            sparse_block[i].val = NULL;
          }
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
  return sparse_block;
}

int elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, mod_t modulus, int nthrds) {
  // copy D to D_red and delete D
  copyBlockMatrixToMultilineMatrix(D, D_red, 1, nthrds);

  ri_t global_next_row_to_reduce  = nthrds * 2;
  ri_t global_last_pivot          = global_next_row_to_reduc - 1;
}

ri_t echelonizeRowsSequential(sm_fl_ml_t *A, ri_t from, ri_t to, mod_t modulus) {
  if (A->nrows == 0)
    return 0;

  ri_t npiv_real  = 0;
  ri_t N          = A->nrows / __GB_NROWS_MULTILINE + 
                    A->nrows % __GB_NROWS_MULTILINE;

  ml_t *ml_row;
  re_l_t *dense_array_1, dense_array_2;
  posix_memalign((void **)&dense_array_1, 16, A->ncols * sizeof(re_l_t));
  posix_memalign((void **)&dense_array_2, 16, A->ncols * sizeof(re_l_t));

  ri_t i;
  long head_line_1      = -1;
  long head_line_2      = -1;
  ci_t head_line_1_idx  = 0;
  ci_t head_line_2_idx  = 0;

  normalizeMultiline(A->ml[from]);

  for (i=from; i<=min(to, N-1); ++i) {
    memset(dense_array_1, 0, A->ncols * sizeof(re_l_t));
    memset(dense_array_2, 0, A->ncols * sizeof(re_l_t));
  }
}
