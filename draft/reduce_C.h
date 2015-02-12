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

void reduce_chunk_1(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
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
#ifdef CONV_A
				spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);
#else
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);
#endif

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k = B->sub + kk;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
#ifdef CONV_B
				spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);
#else
				spaxpy(tc,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);
#endif
			}
		}
	}
}



void reduce_chunk_2(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		)
{
	if (blk_i == 1) {
		reduce_chunk_1(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
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
#ifdef CONV_A
					spaxpy2_block(tc,td,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc,ldc);
#else
					spaxpy2(tc,td,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc,ldc);

#endif
					/* temp_D -= temp_C[j] * B[j] */

					CSR * B_k = B->sub + kk ;
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
#ifdef CONV_B
					spaxpy2_block(tc,td,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd,ldd);
#else
					spaxpy2(tc,td,B_k->data+B_k->start[jj],
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd,ldd);
#endif
				}
				else {
					/* temp_C -= temp_C[j] * A[j] */
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					CSR * A_k = A->sub + kk ;
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
					spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
#else
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
#endif

					/* temp_D -= temp_C[j] * B[j] */

					CSR * B_k =  B->sub + kk;
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
#ifdef CONV_B
					spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);
#else
					spaxpy(tc,B_k->data+B_k->start[jj],
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);
#endif

				}
			}
			else if (td != 0) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = A->sub + kk;
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
				spaxpy_block(td,A_k->data+A_k->start[jj]*UNRL,
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+(ii+1)*ldc);
#else
				spaxpy(td,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+(ii+1)*ldc);
#endif

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k = B->sub + kk;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
#ifdef CONV_B
				spaxpy_block(td,B_k->data+B_k->start[jj]*UNRL,
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+(ii+1)*ldd);
#else
				spaxpy(td,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+(ii+1)*ldd);

#endif

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
#ifdef CONV_A
				spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);
#else
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);
#endif

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k = B->sub + kk ;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
#ifdef CONV_B
				spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);
#else
				spaxpy(tc,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);
#endif

			}
		}
	}
}



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


#ifndef _OPENMP
	SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk,elemt_t);
	SAFE_MALLOC_DECL(temp_C,(index_t)ldc*(index_t)blk,elemt_t);
#endif

	dimen_t k ;
	assert(D->col == B->col);


	for (k = 0 ; k < C->sub_row ; ++k) {
		dimen_t l;
		for (l = 0 ; l < C->sub_col ; ++l) {
			CSR * C_k = CH->sub + k*C->sub_col+l  ;

#ifdef _OPENMP
			/* #pragma omp parallel for */
#pragma omp parallel
#pragma omp single nowait
			{
				dimen_t i ;
				/* #pragma omp for */
#endif
				for ( i = 0 ; i < C_k->row ; i += blk ) {
#ifdef _OPENMP
#pragma omp task
					{
#endif
						dimen_t blk_i = min(blk,C_k->row - i);
#ifdef _OPENMP
						SAFE_MALLOC_DECL(temp_D,((index_t)ldd*(index_t)blk_i),elemt_t);
						SAFE_CALLOC_DECL(temp_C,((index_t)ldc*(index_t)blk_i),elemt_t);
#endif

						index_t i_offset = k*MAT_ROW_BLK + i ;
#ifndef _OPENMP
						cblas_dscal(ldc*blk_i,0.,temp_C,1); 						/* memset(temp_C,0,(index_t)ldc*(index_t)blk_i*sizeof(elemt_t)); */
#endif
#ifdef CONV_C
						sparse_dcopy_block( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);
#else
						sparse_dcopy( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);
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

#ifdef _OPENMP
						free(temp_D);
						free(temp_C);
					} /* task */
#endif /* _OPENMP */
				} /* for */
#ifdef _OPENMP
			} /* single */
#endif
		} /* for */

#ifndef _OPENMP
		free(temp_D);
		free(temp_C);
#endif /* _OPENMP */
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


#endif /* __GB_reduce_C_H */

/* vim: set ft=c: */


