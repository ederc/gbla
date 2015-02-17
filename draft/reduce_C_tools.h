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

#ifndef __GB_reduce_C_tools_H
#define __GB_reduce_C_tools_H

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

#endif /* __GB_reduce_C_tools_H */
