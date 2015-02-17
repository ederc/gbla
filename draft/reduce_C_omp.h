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

#ifndef __GB_reduce_C_omp_H
#define __GB_reduce_C_omp_H

#include "reduce_C_tools.h"

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


		dimen_t k ;
		assert(D->col == B->col);


		for (k = 0 ; k < CH->sub_row ; ++k) {
			CSR * C_k = CH->sub + k  ;
#pragma omp parallel
#pragma omp single nowait
			{
				dimen_t i ;
				for ( i = 0 ; i < C_k->row ; i += blk ) {
#pragma omp task
					{
						dimen_t blk_i = min(blk,C_k->row - i);
						index_t i_offset = k*MAT_ROW_BLK + i ;
						/* SAFE_MALLOC_DECL(temp_D,((index_t)ldd*(index_t)blk_i),elemt_t); */
						SAFE_CALLOC_DECL(temp_C,((index_t)ldc*(index_t)blk_i),elemt_t);

						sparse_copy( blk_i, C_k->start+i, C_k->colid, C_k->data , temp_C, ldc, conv_c);
						elemt_t * temp_D = Dd+i_offset*(index_t)ldd ;
						/* cblas_dcopy(blk_i*ldd,Dd+i_offset*(index_t)ldd,1,temp_D,1); */

						reduce_chunk(blk_i,AH,conv_a,BH,conv_b,temp_C,ldc,temp_D,ldd,p,algo_red);

						/* cblas_dcopy(ldd*blk_i,temp_D,1,Dd+i_offset*(index_t)ldd,1); */

						/* free(temp_D); */
						free(temp_C);
					} /* task */
				} /* for */
			} /* single */
		} /* for */
		Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;
	}
#if 0
	else { /* col division */

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

#endif /* __GB_reduce_C_omp_H */

/* vim: set ft=c: */
