/* gbla: Gröbner Basis Linear Algebra
 * This file is part of gbla.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA, or see  <http://www.gnu.org/licenses/>. */

#ifndef __GBLA_reduce_C_H
#define __GBLA_reduce_C_H

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

		/* SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk,elemt_t); */
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
				elemt_t * temp_D = Dd+i_offset*(index_t)ldd ;

				reduce_chunk(blk_i,AH,conv_a,BH,conv_b,temp_C,ldc,temp_D,ldd,p,algo_red);

			}
		}

		/* free(temp_D); */
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

#endif /* __GBLA_reduce_C_H */

/* vim: set ft=c: */
