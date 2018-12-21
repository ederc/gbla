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

#ifndef __GBLA_sparseops_block_H
#define __GBLA_sparseops_block_H

void     Freduce_double(double p, double * A, index_t n);
int cblas_daxpy(const int N, const double alpha, const double * X, const int incX, double * Y, const int incY);
int cblas_dscal(const int N, const double alpha, const double * X, const int incX);
int cblas_dcopy(const int N, const double * X, const int incX, double * Y, const int incY);

void spaxpy_block(
		elemt_t         tc,
		const elemt_t * A,
		const dimen_t   nb,
		const dimen_t * colid,
		elemt_t       * B)
{
	dimen_t jz = 0 ;

	assert(!((uintptr_t)A % 32u));
	assert(!((uintptr_t)B % 32u));

#ifdef SIMD
	SET1_SIMD(dc,tc);
#endif
	for ( jz = 0 ; jz < nb ; ++jz) {
		dimen_t col = colid[jz];
		dimen_t jzu = jz*UNRL ;
#ifdef SIMD
		AXPY_SIMD(B+col,dc,A+jzu);
#else
		B[col  ] += tc * A[jzu  ] ;
		B[col+1] += tc * A[jzu+1] ;
#if (UNRL>2)
		B[col+2] += tc * A[jzu+2] ;
		B[col+3] += tc * A[jzu+3] ;
#endif
#if (UNRL>4)
		B[col+4] += tc * A[jzu+4] ;
		B[col+5] += tc * A[jzu+5] ;
#endif
#if (UNRL>6)
		B[col+6] += tc * A[jzu+6] ;
		B[col+7] += tc * A[jzu+7] ;
#endif
#endif /* SIMD */

	}
}

void spaxpy2_block(
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

	assert(!((uintptr_t)A % 32u));
	assert(!((uintptr_t)B % 32u));


#ifdef SIMD
	SET1_SIMD(dc,tc);
	SET1_SIMD(dd,td);
#endif
	for ( jz = 0 ; jz < nb ; ++jz) {
		dimen_t col = colid[jz];
		dimen_t jzu = jz*UNRL;
#ifdef SIMD
		AXPY_SIMD(B+col,dc,A+jzu);
#else
		B[col  ]  += tc * A[jzu  ] ;
		B[col+1]  += tc * A[jzu+1] ;
#if (UNRL>2)
		B[col+2]  += tc * A[jzu+2] ;
		B[col+3]  += tc * A[jzu+3] ;
#endif
#if (UNRL>4)
		B[col+4]  += tc * A[jzu+4] ;
		B[col+5]  += tc * A[jzu+5] ;
#endif
#if (UNRL>6)
		B[col+6]  += tc * A[jzu+6] ;
		B[col+7]  += tc * A[jzu+7] ;
#endif
#endif /* SIMD */

#ifdef SIMD
		AXPY_SIMD(B+col+ld,dd,A+jzu);
#else
		B[col+ld  ] += td * A[jzu  ] ;
		B[col+ld+1] += td * A[jzu+1] ;
#if (UNRL>2)
		B[col+ld+2] += td * A[jzu+2] ;
		B[col+ld+3] += td * A[jzu+3] ;
#endif
#if (UNRL>4)
		B[col+ld+4] += td * A[jzu+4] ;
		B[col+ld+5] += td * A[jzu+5] ;
#endif
#if (UNRL>6)
		B[col+ld+6] += td * A[jzu+6] ;
		B[col+ld+7] += td * A[jzu+7] ;
#endif
#endif /* SIMD */
	}
}

void sparse_dcopy_block(
		 dimen_t row
		, index_t * start
		, dimen_t * colid
		, elemt_t * data
		, elemt_t * temp_C
		, dimen_t ldc
		)
{
	assert(!((uintptr_t)data   % 32u));
	assert(!((uintptr_t)temp_C % 32u));

	dimen_t ii = 0 ;
	for ( ii = 0 ; ii < row ;  ++ii) {
		elemt_t * C_off = temp_C+(index_t)ii*(index_t)ldc;
		index_t jz = start[ii] ;
		for ( ; jz < start[ii+1] ; ++jz) {
			dimen_t col = colid[jz] ;
#ifdef SIMD
			COPY_SIMD(C_off+col,data+UNRL*jz);
#else
			C_off[col  ] = data[UNRL*jz  ] ;
			C_off[col+1] = data[UNRL*jz+1] ;
#if (UNRL>2)
			C_off[col+2] = data[UNRL*jz+2] ;
			C_off[col+3] = data[UNRL*jz+3] ;
#endif
#if (UNRL>4)
			C_off[col+4] = data[UNRL*jz+4] ;
			C_off[col+5] = data[UNRL*jz+5] ;
#endif
#if (UNRL>6)
			C_off[col+6] = data[UNRL*jz+6] ;
			C_off[col+7] = data[UNRL*jz+7] ;
#endif
#endif /* SIMD */

		}
	}
}


#endif /*  __GBLA_sparseops_block_H */

/* vim: set ft=c: */
