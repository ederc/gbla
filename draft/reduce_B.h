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

#ifndef __GB_reduce_B_H
#define __GB_reduce_B_H

void     Freduce_double(double p, double * A, index_t n);
int cblas_daxpy(const int N, const double alpha, const double * X, const int incX, double * Y, const int incY);
int cblas_dscal(const int N, const double alpha, const double * X, const int incX);
int cblas_dcopy(const int N, const double * X, const int incX, double * Y, const int incY);



void reduce_B(
		GBMatrix_t      * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D )
{
	/*
	 * // y = U x
	 * for i from n-1 to 0 do
	 *   x[i] = y[i] ;
	 *   for j from i+1 to n-1 do
	 *     x[i,k] -= U[i,j] * x [j,k] ;
	 *   x[i,k] /=  U[i,i] ;
	 *
	 *
	 */
	dimen_t k  ;
	dimen_t ldb = B->ld ;
	dimen_t N = B->col ;
	elemt_t p = A->mod ;
	dimen_t i ;
	SAFE_CALLOC_DECL(row_beg,B->row,dimen_t);
	for ( k = 0 ; k < B->row ; ++k) {
		dimen_t ii = 0 ;
		while ( *(B->ptr+k*ldb+ii) == 0) {
			row_beg[k] += 1 ;
			++ii ;
		}
	}


	for ( k = A->sub_nb ; k --  ; ) {
		CSR * Ad = & (A->sub[k]);
		dimen_t M = Ad->row ;
		/* B = A^(-1) B */

		for ( i = M ; i--    ; ) {
			index_t i_offset = k * MAT_ROW_BLK + i;
			assert( (elemt_t)-1<1); /* unsigned will fail */
			index_t jz  ;
			for ( jz = Ad->start[i]+1 ; jz < Ad->start[i+1] ; ++jz)
			{
				dimen_t kz = Ad->colid[jz];
				dimen_t rs = min(row_beg[i_offset],row_beg[kz]);
				cblas_daxpy(N-rs,-Ad->data[jz],B->ptr+kz*(index_t)ldb+rs,1,B->ptr+i_offset*(index_t)ldb+rs,1);
				row_beg[i_offset] = rs ;
			}
			Mjoin(Freduce,elemt_t)(p,B->ptr+i_offset*(index_t)ldb,N);
			assert(Ad->data[Ad->start[i]] == 1);
		}
	}

	for ( k = 0 ; k < C->sub_nb  ; ++k ) {
		CSR * Cd = C->sub + k ;
		dimen_t ldd = D->ld ;
		/* D = D - C . B */
		assert( (elemt_t)-1<1); /* unsigned will fail */
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for ( i = 0 ; i < Cd->row ;  ++i) {
			index_t i_offset = k * MAT_ROW_BLK + i;
			index_t jz  ;
			for ( jz = Cd->start[i]; jz < Cd->start[i+1] ; ++jz ) {
				dimen_t kz = Cd->colid[jz];
				dimen_t rs = row_beg[kz] ;
				cblas_daxpy(N-rs,-Cd->data[jz],B->ptr+kz*(index_t)ldb+rs,1,D->ptr+i_offset*(index_t)ldd+rs,1);
			}
			Mjoin(Freduce,elemt_t)(p,D->ptr+i_offset*(index_t)ldd, N) ;
		}
	}

	free(row_beg);
}

#endif /*  __GB_reduce_B_H*/

/* vim: set ft=c: */
