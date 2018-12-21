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

#ifndef __GBLA_reduce_B_omp_H
#define __GBLA_reduce_B_omp_H

void     Freduce_double(double p, double * A, index_t n);
int cblas_daxpy(const int N, const double alpha, const double * X, const int incX, double * Y, const int incY);
int cblas_dscal(const int N, const double alpha, const double * X, const int incX);
int cblas_dcopy(const int N, const double * X, const int incX, double * Y, const int incY);

void reduce_B(
		GBMatrix_t      * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D
		, int nth
		)
{
	assert( (elemt_t)-1<1); /* unsigned will fail */
	dimen_t ldb = B->ld ;
	dimen_t N = B->col ;
	elemt_t p = A->mod ;

	{ /* B = A^(-1) B */
		dimen_t k  ;
		for ( k = A->sub_row ; k --  ; ) {
			CSR * Ad = A->sub + k ;
			dimen_t M = Ad->row ;
			dimen_t band ;
#pragma omp parallel for
			for (band = 0 ; band < DIVIDE_INTO(B->col,MAT_COL_BLK) ; ++band) {
				dimen_t col_start = band*MAT_COL_BLK;
				dimen_t col_width = min ((dimen_t)MAT_COL_BLK,B->col-col_start);
				dimen_t i;

				for ( i = M ; i--    ; ) {
					dimen_t i_offset = k * MAT_ROW_BLK + i;
					index_t jz  ;
					elemt_t * B_off = B->ptr+(index_t)i_offset*(index_t)ldb+col_start;
					for ( jz = Ad->start[i]+1 ; jz < Ad->start[i+1] ; ++jz)
					{
						dimen_t kz = Ad->colid[jz];
						cblas_daxpy(col_width,-Ad->data[jz],B->ptr+(index_t)kz*(index_t)ldb+col_start,1,B_off,1);
					}
					Mjoin(Freduce,elemt_t)(p,B_off,col_width);
				}
			}
		}
	}

	{ /* D = D - C . B */
		dimen_t k ;
		for ( k = 0 ; k < C->sub_row  ; ++k ) {
			CSR * Cd = C->sub + k ;
			dimen_t ldd = D->ld ;
			dimen_t band ;
#pragma omp parallel for
			for (band = 0 ; band < DIVIDE_INTO(Cd->row,MAT_SUB_BLK) ; ++band) {
				dimen_t i;
				for ( i = band*MAT_SUB_BLK ; i < min((band+1)*MAT_SUB_BLK,Cd->row) ;  ++i) {
					index_t i_offset = k * MAT_ROW_BLK + i;
					index_t jz  ;
					elemt_t * D_off = D->ptr+(index_t)i_offset*(index_t)ldd;
					for ( jz = Cd->start[i]; jz < Cd->start[i+1] ; ++jz ) {
						dimen_t kz = Cd->colid[jz];
						cblas_daxpy(N,-Cd->data[jz],B->ptr+(index_t)kz*(index_t)ldb,1,D_off,1);
					}
					Mjoin(Freduce,elemt_t)(p,D_off, N) ;
				}
			}
		}
	}

}


#endif /* __GBLA_reduce_B_omp_H */

/* vim: set ft=c: */
