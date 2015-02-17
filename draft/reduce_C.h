/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Brice Boyer <brice.boyer@lip6.fr>
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

#ifndef __GB_reduce_C_H
#define __GB_reduce_C_H

#include "sparseops.h"
#include "sparseops_block.h"

void sparse_copy(
		dimen_t row
		, index_t * start
		, dimen_t * colid
		, elemt_t * data
		, elemt_t * temp_C
		, dimen_t ldc
		, int conv_c
		)
{
	if (conv_c == 1) {
		sparse_dcopy_block(row,start,colid,data,temp_C,ldc);
	}
	else {
		assert(conv_c == 0);
		sparse_dcopy(row,start,colid,data,temp_C,ldc);
	}
}



void reduce_chunk_1(
		dimen_t blk_i
		, GBMatrix_t * A
		, int conv_a
		, GBMatrix_t * B
		, int conv_b
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p)
{
	dimen_t j, ii ;
	for ( j = 0 ; j < A->col ; ++j) { /* A->col == C->col */
		/* XXX invert loops ? */
		for ( ii = 0 ; ii < blk_i ; ii += 1 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			if (tc != 0.) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = A->sub + kk ;
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				if (conv_a == 1)
					spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
				else
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k = B->sub + kk;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				if (conv_b == 1)
					spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);
				else
					spaxpy(tc,B_k->data+B_k->start[jj],
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);
			}
		}
	}
}



void reduce_chunk_2(
		dimen_t blk_i
		, GBMatrix_t * A
		, int conv_a
		, GBMatrix_t * B
		, int conv_b
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		)
{
	if (blk_i == 1) {
		reduce_chunk_1(blk_i,A,conv_a,B,conv_b,temp_C,ldc,temp_D,ldd,p);
		return ;
	}

	dimen_t j, ii ;
	for ( j = 0 ; j < A->col ; ++j) {
		/* XXX invert loops ? */
		for ( ii = 0 ; ii < blk_i/2*2 ; ii += 2 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			elemt_t td = Mjoin(fmod,elemt_t)(-temp_C[(ii+1)*ldc+j],p) ;
			/* temp_C -= temp_C[j] * A[j] */
			if (tc != 0.) {
				if (td != 0.) {
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					CSR * A_k = A->sub + kk ;
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
					if (conv_a == 1)
						spaxpy2_block(tc,td,A_k->data+A_k->start[jj]*UNRL,
								sz,
								A_k->colid+A_k->start[jj]
								,temp_C+ii*ldc,ldc);
					else
						spaxpy2(tc,td,A_k->data+A_k->start[jj],
								sz,
								A_k->colid+A_k->start[jj]
								,temp_C+ii*ldc,ldc);

					/* temp_D -= temp_C[j] * B[j] */

					CSR * B_k = B->sub + kk ;
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
					if (conv_b == 1)
						spaxpy2_block(tc,td,B_k->data+B_k->start[jj]*UNRL,
								sz,
								B_k->colid+B_k->start[jj],
								temp_D+ii*ldd,ldd);
					else
						spaxpy2(tc,td,B_k->data+B_k->start[jj],
								sz,
								B_k->colid+B_k->start[jj],
								temp_D+ii*ldd,ldd);
				}
				else {
					/* temp_C -= temp_C[j] * A[j] */
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					CSR * A_k = A->sub + kk ;
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
					if (conv_a == 1)
						spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
								sz,
								A_k->colid+A_k->start[jj]
								,temp_C+ii*ldc);
					else
						spaxpy(tc,A_k->data+A_k->start[jj],
								sz,
								A_k->colid+A_k->start[jj]
								,temp_C+ii*ldc);

					/* temp_D -= temp_C[j] * B[j] */

					CSR * B_k =  B->sub + kk;
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
					if (conv_b == 1)
						spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
								sz,
								B_k->colid+B_k->start[jj],
								temp_D+ii*ldd);
					else
						spaxpy(tc,B_k->data+B_k->start[jj],
								sz,
								B_k->colid+B_k->start[jj],
								temp_D+ii*ldd);

				}
			}
			else if (td != 0) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = A->sub + kk;
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				if (conv_a == 1)
					spaxpy_block(td,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+(ii+1)*ldc);
				else
					spaxpy(td,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+(ii+1)*ldc);

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k = B->sub + kk;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				if (conv_b == 1)
					spaxpy_block(td,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+(ii+1)*ldd);
				else
					spaxpy(td,B_k->data+B_k->start[jj],
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+(ii+1)*ldd);


			}
		}
		for (  ; ii < blk_i ; ii += 1 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			/* temp_C -= temp_C[j] * A[j] */
			if (tc != 0.) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = A->sub + kk ;
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				if (conv_a == 1)
					spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
				else
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k = B->sub + kk ;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				if (conv_b == 1)
					spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);
				else
					spaxpy(tc,B_k->data+B_k->start[jj],
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);

			}
		}
	}
}

void reduce_chunk(
		dimen_t blk_i
		, GBMatrix_t * A
		, int conv_a
		, GBMatrix_t * B
		, int conv_b
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		, int algo
		)
{
	switch (algo) {
		case 1 :
			reduce_chunk_1(blk_i,A,conv_a,B,conv_b,temp_C,ldc,temp_D,ldd,p);
			break;
		case 2 :
			reduce_chunk_2(blk_i,A,conv_a,B,conv_b,temp_C,ldc,temp_D,ldd,p);
			break;
		default:
			fprintf(stderr,"bad reduction algo\n");
			exit(-5);
	}
}




#ifndef _OPENMP

void reduce_C(
		GBMatrix_t      * A
		, int conv_a
		, GBMatrix_t    * B
		, int conv_b
		, GBMatrix_t    * C
		, int conv_c
		, DNS * D
		, int algo_red
		, int nth  )
{
	struct timeval tic,tac ;
	gettimeofday(&tic,NULL);
	GBMatrix_t * AH, * BH, * CH  ;

	if (conv_b == 1) {
		SAFE_MALLOC(BH,1,GBMatrix_t);
		initSparse(BH);
		convert_CSR_2_CSR_block(BH,B);
	}
	else
		BH = B ;

	if (conv_a == 1) {
		SAFE_MALLOC(AH,1,GBMatrix_t);
		initSparse(AH);
		convert_CSR_2_CSR_block(AH,A);
	}
	else
		AH = A ;

	if (conv_c == 1) {
		SAFE_MALLOC(CH,1,GBMatrix_t);
		initSparse(CH);
		convert_CSR_2_CSR_block(CH,C);
	}
	else
		CH = C ;

	gettimeofday(&tac,NULL);
	fprintf(stderr,"   -- convert time    : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

	dimen_t ldd = D->ld ;
	dimen_t ldc = ALIGN(CH->col) ;


	elemt_t p = AH->mod ;
	elemt_t * Dd = D->ptr;


	dimen_t blk = MAT_SUB_BLK ;

	assert((index_t)ldd*(index_t)blk < UINT32_MAX);
	assert((index_t)ldc*(index_t)blk < UINT32_MAX);


	assert(CH->sub_col == 1);

	if (CH->sub_col == 1) {

		/* XXX no need for temp_D */
		SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk,elemt_t);
		SAFE_MALLOC_DECL(temp_C,(index_t)ldc*(index_t)blk,elemt_t);

		dimen_t k ;
		assert(D->col == B->col);


		for (k = 0 ; k < CH->sub_row ; ++k) {
			CSR * C_k = CH->sub + k  ;
			dimen_t i ;
			for ( i = 0 ; i < C_k->row ; i += blk ) {
				dimen_t blk_i = min(blk,C_k->row - i);

				index_t i_offset = k*MAT_ROW_BLK + i ;
				cblas_dscal(ldc*blk_i,0.,temp_C,1);
				/* memset(temp_C,0,(index_t)ldc*(index_t)blk_i*sizeof(elemt_t)); */
				sparse_copy( blk_i, C_k->start+i, C_k->colid, C_k->data , temp_C, ldc, conv_c);
				/* elemt_t * temp_D = Dd+i_offset*(index_t)ldd ; */
				cblas_dcopy(blk_i*ldd,Dd+i_offset*(index_t)ldd,1,temp_D,1);

				reduce_chunk(blk_i,AH,conv_a,BH,conv_b,temp_C,ldc,temp_D,ldd,p,algo_red);

				cblas_dcopy(ldd*blk_i,temp_D,1,Dd+i_offset*(index_t)ldd,1);

			}
		}

		free(temp_D);
		free(temp_C);
		Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;
	}
#if 0
	else { /* col division */
		for (l = 0 ; l < C->sub_col ; ++l) {
			for (k = 0 ; k < C->sub_row ; ++k) {
				CSR * C_0 = CH->sub + k*C->sub_col+l  ;
				dimen_t rows = C_k->row ;
				dimen_t cols = C_k->col ;

				SAFE_MALLOC_DECL(temp_D,rows*cols,elemt_t);
				SAFE_CALLOC_DECL(temp_C,rows*cols,elemt_t);
				SAFE_CALLOC_DECL(Cinv,rows*cols,elemt_t);
				index_t i_offset = k*MAT_ROW_BLK + i ;
				index_t j_offset = l*MAT_COL_BLK  ;

				sparse_copy( rows, C_k->start+i, C_k->colid, C_k->data , temp_C, ldc, conv_c);
				reduce_leader(rows,AH,temp_C,invC,ldc,p,algo_red);

				for (r = l+1 ; r < C->sub_col ; ++r) {
					reduce_other(rows,AH,temp_C,invC,ldc,p,algo_red);
				}
			}
		}

		free(temp_D);
		free(temp_C);
		Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;

	}
#endif

	if (conv_b) {
		freeMat(BH);
		free(BH);
	}

	if (conv_a) {
		freeMat(AH);
		free(AH);
	}

	if (conv_c) {
		freeMat(CH);
		free(CH);
	}
}


#else /* _OPENMP */

void reduce_C(
		GBMatrix_t      * A
		, GBMatrix_t    * B
		, GBMatrix_t    * C
		, DNS * D
		, int nth  )
{
	struct timeval tic,tac ;
	gettimeofday(&tic,NULL);

#ifdef CONV_B
	SAFE_MALLOC_DECL(BH,1,GBMatrix_t);
	initSparse(BH);
	convert_CSR_2_CSR_block(BH,B);
#else
	GBMatrix_t * BH = B ;
#endif

#ifdef CONV_A
	SAFE_MALLOC_DECL(AH,1,GBMatrix_t);
	initSparse(AH);
	convert_CSR_2_CSR_block(AH,A);
#else
	GBMatrix_t * AH = A ;
#endif

#ifdef CONV_C
	SAFE_MALLOC_DECL(CH,1,GBMatrix_t);
	initSparse(CH);
	convert_CSR_2_CSR_block(CH,C);
#else
	GBMatrix_t * CH = C ;
#endif

	gettimeofday(&tac,NULL);
	fprintf(stderr,"   -- convert time    : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

	dimen_t ldd = D->ld ;
	dimen_t ldc = ALIGN(C->col) ;


	elemt_t p = AH->mod ;
	elemt_t * Dd = D->ptr;
	dimen_t i ;


	dimen_t blk = MAT_SUB_BLK ;

	assert((index_t)ldd*(index_t)blk < UINT32_MAX);
	assert((index_t)ldc*(index_t)blk < UINT32_MAX);



	dimen_t k ;
	assert(D->col == B->col);


	dimen_t l;
	for (l = 0 ; l < C->sub_col ; ++l) {
		CSR * C_k = CH->sub + k*C->sub_col+l  ;

#pragma omp parallel
#pragma omp single nowait
		{
			dimen_t i ;
			for ( i = 0 ; i < C_k->row ; i += blk ) {
#pragma omp task
				{
					dimen_t blk_i = min(blk,C_k->row - i);
					SAFE_MALLOC_DECL(temp_D,((index_t)ldd*(index_t)blk_i),elemt_t);
					SAFE_CALLOC_DECL(temp_C,((index_t)ldc*(index_t)blk_i),elemt_t);

					index_t i_offset = k*MAT_ROW_BLK + i ;
#ifdef CONV_C
					sparse_copy_block( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);
#else
					sparse_copy( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);
#endif

					cblas_dcopy(blk_i*ldd,Dd+i_offset*(index_t)ldd,1,temp_D,1);


#ifdef USE_SAXPY
#if  defined(USE_SAXPY2)
#error "make a choice"
#endif
					reduce_chunk_1(blk_i,AH,BH,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPY */

#ifdef USE_SAXPY2
#if defined(USE_SAXPY)
#error "make a choice"
#endif
					reduce_chunk_2(blk_i,AH,BH,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPY2 */


					cblas_dcopy(ldd*blk_i,temp_D,1,Dd+i_offset*(index_t)ldd,1);

					free(temp_D);
					free(temp_C);
				} /* task */
			} /* for */
		} /* single */
	} /* for */

	Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;

#ifdef CONV_B
	freeMat(BH);
	free(BH);
#endif

#ifdef CONV_A
	freeMat(AH);
	free(AH);
#endif
#ifdef CONV_C
	freeMat(CH);
	free(CH);
#endif
}

#endif /* _OPENMP */

#endif /* __GB_reduce_C_H */

/* vim: set ft=c: */
