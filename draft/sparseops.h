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

#ifndef __GBLA_sparseops_H
#define __GBLA_sparseops_H

void     Freduce_double(double p, double * A, index_t n);
int cblas_daxpy(const int N, const double alpha, const double * X, const int incX, double * Y, const int incY);
int cblas_dscal(const int N, const double alpha, const double * X, const int incX);
int cblas_dcopy(const int N, const double * X, const int incX, double * Y, const int incY);

void spaxpy(
		elemt_t         tc,
		const elemt_t * A,
		const dimen_t   nb,
		const dimen_t * colid,
		elemt_t       * B)
{
	dimen_t jz = 0 ;

#ifndef DEROULE
	for ( jz = 0 ; jz < nb ; ++jz) {
		B[colid[jz]] += tc * A[jz] ;
	}
#else /* UNROLL */
	for ( jz = 0 ; jz < (nb/UNRL)*UNRL ; jz += UNRL) {
			B[colid[jz]]   += tc * A[jz] ;
			B[colid[jz+1]] += tc * A[jz+1] ;
#if (UNRL>2)
			B[colid[jz+2]] += tc * A[jz+2] ;
			B[colid[jz+3]] += tc * A[jz+3] ;
#endif
#if (UNRL>4)
			B[colid[jz+4]] += tc * A[jz+4] ;
			B[colid[jz+5]] += tc * A[jz+5] ;
#endif
#if (UNRL>6)
			B[colid[jz+6]] += tc * A[jz+6] ;
			B[colid[jz+7]] += tc * A[jz+7] ;
#endif
	}
	assert((int64_t)nb - (int64_t)jz < (int64_t)UNRL);
	for (  ; jz < nb ; ++jz) {
		B[colid[jz]] += tc * A[jz] ;
	}
#endif
}

void spaxpy2(
		elemt_t         tc,
		elemt_t         td,
		const elemt_t * A,
		const dimen_t   nb,
		const dimen_t * colid,
		elemt_t       * B,
		const dimen_t   ld
	    )
{
	dimen_t jz = 0;

#ifdef DEROULE
	for ( jz = 0 ; jz < (nb/UNRL)*UNRL ; jz +=UNRL) {
			B[colid[jz  ]]    += tc * A[jz  ] ;
			B[colid[jz+1]]    += tc * A[jz+1] ;
#if (UNRL>2)
			B[colid[jz+2]]    += tc * A[jz+2] ;
			B[colid[jz+3]]    += tc * A[jz+3] ;
#endif
#if (UNRL>4)
			B[colid[jz+4]]    += tc * A[jz+4] ;
			B[colid[jz+5]]    += tc * A[jz+5] ;
#endif
#if (UNRL>6)
			B[colid[jz+6]]    += tc * A[jz+6] ;
			B[colid[jz+7]]    += tc * A[jz+7] ;
#endif

			B[colid[jz  ]+ld] += td * A[jz  ] ;
			B[colid[jz+1]+ld] += td * A[jz+1] ;
#if (UNRL>2)
			B[colid[jz+2]+ld] += td * A[jz+2] ;
			B[colid[jz+3]+ld] += td * A[jz+3] ;
#endif
#if (UNRL>4)
			B[colid[jz+4]+ld] += td * A[jz+4] ;
			B[colid[jz+5]+ld] += td * A[jz+5] ;
#endif
#if (UNRL>6)
			B[colid[jz+6]+ld] += td * A[jz+6] ;
			B[colid[jz+7]+ld] += td * A[jz+7] ;
#endif
	}
	for (  ; jz < nb ; ++jz) {
		B[colid[jz]]    += tc * A[jz] ;
		B[colid[jz]+ld] += td * A[jz] ;
	}
#else /* UNROLL */
	for ( jz = 0 ; jz < nb ; ++jz) {
		B[colid[jz]]    += tc * A[jz] ;
		B[colid[jz]+ld] += td * A[jz] ;
	}
#endif
}

void sparse_dcopy(
		 dimen_t row
		, index_t * start
		, dimen_t * colid
		, elemt_t * data
		, elemt_t * temp_C
		, dimen_t ldc
		)
{
	dimen_t ii = 0 ;
	for ( ii = 0 ; ii < row ;  ++ii) {
		elemt_t * C_off = temp_C+(index_t)ii*(index_t)ldc;
		index_t jz = start[ii] ;
#ifdef DEROULE
		index_t last = start[ii] + (start[ii+1]-start[ii])/UNRL*UNRL ;
		assert(start[ii+1]-last < UNRL);
		for ( ; jz <  last ; jz += UNRL) {
			C_off[colid[jz  ]] = data[jz  ] ;
			C_off[colid[jz+1]] = data[jz+1] ;
#if (UNRL>2)
			C_off[colid[jz+2]] = data[jz+2] ;
			C_off[colid[jz+3]] = data[jz+3] ;
#endif
#if (UNRL>4)
			C_off[colid[jz+4]] = data[jz+4] ;
			C_off[colid[jz+5]] = data[jz+5] ;
#endif
#if (UNRL>6)
			C_off[colid[jz+6]] = data[jz+6] ;
			C_off[colid[jz+7]] = data[jz+7] ;
#endif
		}
#endif /*  DEROULE */

		for (  ; jz < start[ii+1] ; ++jz) {
			C_off[colid[jz]] = data[jz] ;
		}
	}
}


#endif /*  __GBLA_sparseops_H */

/* vim: set ft=c: */
