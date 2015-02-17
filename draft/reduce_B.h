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



#ifndef _OPENMP
void reduce_B(
		GBMatrix_t      * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D
		, int nth
		)
{
	assert( (elemt_t)-1<1); /* unsigned will fail */
	dimen_t l,k  ;
	dimen_t ldb = B->ld ;
	dimen_t N = B->col ;
	elemt_t p = A->mod ;
	dimen_t i ;
	SAFE_CALLOC_DECL(row_beg,B->row,dimen_t);
	for ( k = 0 ; k < B->row ; ++k) {
		dimen_t ii = 0 ;
		elemt_t * B_off =B->ptr+(index_t)k*(index_t)ldb ;
		while ( *(B_off+ii) == 0) {
			row_beg[k] += 1 ;
			++ii ;
		}
	}


	{ /* B = A^(-1) B */
		for ( k = A->sub_row ; k --  ; ) {
			for ( l = 0 ; l < A->sub_col ; ++l) { /* XXX triangular */
				CSR * Ad = A->sub + k * A->sub_col + l ;
				dimen_t M = Ad->row ;
				for ( i = M ; i--    ; ) {
					dimen_t i_offset = k * MAT_ROW_BLK + i;
					elemt_t * B_off = B->ptr+(index_t)i_offset*(index_t)ldb;
					dimen_t rs ;

					index_t jz = Ad->start[i]+1 ;
					for (  ; jz < Ad->start[i+1] ; ++jz)
					{
						dimen_t kz = Ad->colid[jz];
						rs = min(row_beg[i_offset],row_beg[kz]);
						elemt_t * B_offset = B->ptr+(index_t)kz*(index_t)ldb ;
						cblas_daxpy(N-rs,-Ad->data[jz],B_offset+rs,1,B_off+rs,1);
						row_beg[i_offset] = rs ;
					}
					rs = row_beg[i_offset] ;
					Mjoin(Freduce,elemt_t)(p,B_off+rs,N-rs);
					assert(Ad->data[Ad->start[i]] == 1);
				}
			}
		}
	}

	{ /* D = D - C . B */
		dimen_t ldd = D->ld ;
		for ( k = 0 ; k < C->sub_row  ; ++k ) {
			for ( l = 0 ; l < C->sub_col  ; ++l ) {
				CSR * Cd = C->sub + k*C->sub_col + l ;
				index_t j_offset = l * MAT_COL_BLK ;
				for ( i = 0 ; i < Cd->row ;  ++i) {
					index_t i_offset = k * MAT_ROW_BLK + i;
					index_t jz  ;
					elemt_t * D_off = D->ptr+(index_t)i_offset*(index_t)ldd;
					elemt_t * B_off = B->ptr+(index_t)j_offset*(index_t)ldb;
					assert((Cd->start[i+1]-Cd->start[i])*(p-1)*(p-1) < REDLIM);
					for ( jz = Cd->start[i]; jz < Cd->start[i+1] ; ++jz ) {
						dimen_t kz = Cd->colid[jz];
						dimen_t rs = row_beg[kz] ;
						elemt_t * B_offset = B_off+(index_t)kz*(index_t)ldb ;
						cblas_daxpy(N-rs,-Cd->data[jz],B_offset+rs,1,D_off+rs,1);
					}
				}
				Mjoin(Freduce,elemt_t)(p,D->ptr+k*MAT_ROW_BLK, N*Cd->row) ;
			}
		}
	}

	free(row_beg);
}

#else /* _OPENMP */

void reduce_B(
		GBMatrix_t      * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D
		, int nth
		)
{
	exit(-1);
	assert( (elemt_t)-1<1); /* unsigned will fail */
	dimen_t ldb = B->ld ;
	dimen_t ldd = D->ld ;
	dimen_t N = B->col ;
	elemt_t p = A->mod ;

	{ /* B = A^(-1) B */
		dimen_t l  ;
#pragma omp parallel for
		for ( l = 0 ; l < A->sub_col ; ++l) {
			dimen_t k  ;
			for ( k = A->sub_row ; k --  ; ) {
				CSR * Ad = A->sub + k*A->sub_col + l;
				dimen_t M = Ad->row ;
				dimen_t j_offset = l * MAT_COL_BLK ;
				elemt_t * B_off = B->ptr + j_offset ;
				dimen_t col_width = min((dimen_t)MAT_COL_BLK,B->col-j_offset);
				for ( i = M ; i--    ; ) {
					dimen_t i_offset = k * MAT_ROW_BLK + i;
					index_t jz  ;
					elemt_t * B_offset = B_off+(index_t)i_offset*(index_t)ldb;
					for ( jz = Ad->start[i]+1 ; jz < Ad->start[i+1] ; ++jz)
					{
						dimen_t kz = Ad->colid[jz];
						cblas_daxpy(col_width,-Ad->data[jz],B_off+(index_t)kz*(index_t)ldb,1,B_offset,1);
					}
					Mjoin(Freduce,elemt_t)(p,B_offset,col_width);
				}
			}
		}
	}

	{ /* D = D - C . B */
		dimen_t k ;
#pragma omp parallel for
		for ( k = 0 ; k < C->sub_row  ; ++k ) {
			dimen_t l ;
			for ( l = 0 ; l < C->sub_col  ; ++l ) {
				CSR * Cd = C->sub + k * C->sub_col + l ;
				index_t j_offset = l * MAT_COL_BLK ;
				dimen_t i;
				for ( i = 0 ; i < Cd->row ;  ++i) {
					index_t i_offset = k * MAT_ROW_BLK + i;
					index_t jz  ;
					elemt_t * D_off = D->ptr+(index_t)i_offset*(index_t)ldd;
					elemt_t * B_off = B->ptr+(index_t)j_offset*(index_t)ldb;
					for ( jz = Cd->start[i]; jz < Cd->start[i+1] ; ++jz ) {
						dimen_t kz = Cd->colid[jz];
						cblas_daxpy(N,-Cd->data[jz],B_off+(index_t)kz*(index_t)ldb,1,D_off,1);
					}
					Mjoin(Freduce,elemt_t)(p,D_off, N) ;
				}
			}
		}
	}

}

#endif /* _OPENMP */


#endif /*  __GB_reduce_B_H*/

/* vim: set ft=c: */
