#ifndef __GB_io_H
#define __GB_io_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "matrix.h"
#include "field_ops.h"

/**
 * the matrix we have as input is "almost upper triangular", ie below every first
 * element in a row, there may be non zero elements on the next rows below.
 * There is no empty line.
 *
 */







/* buffer:
*/
uint32_t getSparsestRows(uint32_t * colid
		, uint64_t * start
		, uint32_t row
		, uint32_t col
		, uint32_t * pivots_data /* pivots_data[j] is the sparsest row with j as pivot */
		)
{
	/* SAFE_MALLOC_DECL(sparse,col,uint32_t); */
	/* SAFE_MALLOC_DECL(pivot_row,col,uint32_t); */
	uint32_t i = 0;
	uint32_t k_dim = 0;
	/* for ( i = 0 ; i < row ; ++i ) { */
		/* sparse[i] = (uint32_t)(-1); */
	/* } */
	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (uint32_t)(-1);
	}
	for ( i = 0 ; i < row ; ++i ) {
		uint32_t pivot_j = colid[start[i]] ;     /* first column row i */
		uint32_t creux   = start[i+1]-start[i] ; /* length of row i */
		/* sparse[i] = creux ; */
		/* if (pivot_j <= i) break ; */
		/* assert(k_dim < row); */
		assert(pivot_j < col);
		uint32_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (uint32_t)(-1)) {
			pivots_data[pivot_j] = i ;
			++k_dim;
		}
		else  {
			/* XXX no old_creux does not need to be recomputed */
			uint32_t old_creux = start[old_i+1]-start[old_i];
			if (old_creux > creux) { /* this row is sparser */
				pivots_data[pivot_j] = i ;
			}
		}
	}

	/* for ( i = 0 ; i < col; ++i ) {
		if (pivots_data[i] == (uint32_t)(-1)) {
			fprintf(stderr,"column %u: not pivot\n",i);
		}
	}                                                              */

	fprintf(stderr,"k dim will be : %u\n",k_dim);
	return k_dim ;
}

uint32_t createPivots(uint32_t * pivots /* list of pivot columns */
		, uint32_t * pivots_size
		, uint32_t * nonpiv /* list of not pivot columns */
		, uint32_t * nonpiv_size
		, uint32_t   k_dim
		, uint32_t * data
		, uint32_t   n)
{
	uint32_t j ;
	/* uint32_t cur_pi = 0 ; */
	/* uint32_t cur_pe = 0 ; */
	*nonpiv_size =  0;
	*pivots_size = 0 ;
	for ( j = 0 ; j < n ; ++j) {
		if (data[j] == (uint32_t) -1) {
			nonpiv[(*nonpiv_size)++] = j ;
			assert(*nonpiv_size <= n-k_dim);
		}
		else {
			pivots[(*pivots_size)++] = data[j] ;
			if (k_dim == *pivots_size)
				break;
		}
	}
	 /* XXX k_dim == pivots_size */

	uint32_t k_split =  *pivots_size + *nonpiv_size ;

	return   k_split ;
}

void splitHorizontal(
		GBMatrix_t * A
		, GBMatrix_t * C
		, uint32_t * colid_zo
		, uint64_t * start_zo
		, uint32_t row
		, uint32_t nnz
		, uint32_t * pivots
		, uint32_t pivots_size
		, uint32_t * map_zo_pol
		)
{
	uint32_t last_pivot=0;
	uint32_t i = 0 ;
	assert(start_zo[row] == nnz);
	for ( ; i < row ; ++i)  {
		assert(start_zo[i] < nnz);
		if ( (last_pivot < pivots_size) && (pivots[last_pivot] == i) ) {
			appendRow(A,colid_zo+start_zo[i],start_zo[i+1]-start_zo[i],map_zo_pol[i]);
			++last_pivot;
		}
		else {
			appendRow(C,colid_zo+start_zo[i],start_zo[i+1]-start_zo[i],map_zo_pol[i]);
		}
	}
	return ;
}


void expandColid( const uint32_t * compress, uint32_t size_compressed
		, uint32_t * expand
		, uint32_t size_expand)
{
	uint32_t mask = (1<<31);
	uint32_t i = 0,j=0 ;
	uint32_t col ;
	for ( ; i < size_compressed ;) {
		col = compress[i++] ;
		if (col & mask) {
			expand[j++] = col ^ mask ;
		}
		else {
			uint32_t k = 0 ;
			for (; k < compress[i] ;++k) {
				expand[j++] = col + k;
			}
			++i;
		}
	}
	assert(j == size_expand);
	assert(i == size_compressed);
}


void splitVerticalUnit(
		const GBMatrix_t * A_init,
		const uint32_t * nonpiv,
		const uint32_t nonpiv_size,
		const CSR_pol * polys,
		GBMatrix_t * A,
		DenseMatrix_t * B
		)
{

	uint32_t max_col = A->col + nonpiv_size ;

	/* Init dense */

	SAFE_CALLOC(B->data,B->row*B->col,elem_t);


	/* Init sparse */

	const CSR_zo * A_k = getLastMatrixConst(A_init);
	appendMatrix(A);
	CSR_zo * Ad = getLastMatrix(A);
	Ad->row = A->row ;
	assert(A->row == A_k->row);
	Ad->col = A->col ;
	Ad->nnz = 0 ;
	/* Ad->mod = A_k->mod ; */

	SAFE_REALLOC(Ad->start_zo,Ad->row+1,uint64_t);
	Ad->start_zo[0] = 0 ;
	SAFE_MALLOC(Ad->colid_zo,A_k->nnz,uint32_t);
	SAFE_MALLOC(Ad->data    ,A_k->nnz,elem_t);
	assert(Ad->map_zo_pol==NULL);


	uint32_t ldb = B->col ;
	uint32_t here  = 0; /* shift */
	uint32_t there = 0;
	here = 0;
	/* uint32_t a_off = 0; */
	uint32_t i,j ;
	for (i = 0 ; i <= Ad->row ; ++i) {
		Ad->start_zo[i] = 0 ;
	}
	for (i = 0 ; i < A_k->row ; ++i) {
		uint32_t start_p = polys->start_pol[ A_k->map_zo_pol[i] ] ;
		elem_t * d = polys->vals_pol+start_p ;
		for (j = A_k->start_zo[i] ; j < A_k->start_zo[i+1] ; ++j) {
			uint32_t k = A_k->colid_zo[j]  ;
			/* int b_touch = 0 ; */
			if ( k  < max_col ) {
				while (here < nonpiv_size && k > nonpiv[here]) {
					here ++ ;
				}
				if (here < nonpiv_size && k == nonpiv[here]) {
					B->data[i*ldb+here] = *d; /* XXX */
				}
				else {
					Ad->start_zo[i+1] += 1 ;
					Ad->colid_zo[there] = k-here ;
					Ad->data[there++] = *d ;
				}
			}
			else {
				assert(nonpiv_size+k-max_col < B->col);
				B->data[i*ldb+(nonpiv_size+k-max_col)] = *d ; /* XXX */
			}
			++d;
		}
		assert(here <= nonpiv_size);
		here = 0 ;
	}
	for ( i = 1 ; i < Ad->row ; ++i) {
		Ad->start_zo[i+1] += Ad->start_zo[i] ;
	}
	A->nnz = Ad->nnz = Ad->start_zo[A->row] ;
	assert(A->nnz);
	assert(A->nnz == there);
	/* realloc to smaller */
	SAFE_REALLOC(Ad->colid_zo,Ad->nnz,uint32_t);
	SAFE_REALLOC(Ad->data    ,Ad->nnz,elem_t);

}

void splitVertical(
		const GBMatrix_t * A_init
		, const GBMatrix_t * C_init
		, const uint32_t * nonpiv
		, const uint32_t   nonpiv_size
		, const CSR_pol * polys
		, GBMatrix_t * A
		, DenseMatrix_t * B
		/* , GBMatrix_t * Bt */
		, GBMatrix_t * C
		, DenseMatrix_t * D
		)
{
	A->row = B->row = A_init->row ;
	C->row = D->row = C_init->row ;
	A->mod = B->mod = C->mod = D->mod = A_init->mod;
	A->col = C->col = A_init->row ;
	B->col = D->col = A_init->col - A_init->row ;

	/* uint32_t nonpiv_size = A_init->col - A_init->row ; |+ could be smaller +| */


	assert(A->row);
	assert(C->row);

	splitVerticalUnit(A_init,nonpiv,nonpiv_size,polys,A,B);
	splitVerticalUnit(C_init,nonpiv,nonpiv_size,polys,C,D);

}


/* matrix reader and row splitter
 * A_init [out] the top part with upper triangular left
 * B_init [out] the bottom part
 * polys [out] the polynomials used in the matrices (shared)
 * fh [in] the file in appropriate format
 */
uint32_t * readFileSplit(
		GBMatrix_t * A,
		DenseMatrix_t * B,
		GBMatrix_t * C,
		DenseMatrix_t * D
		/* GBMatrix_t * A_init */
		/* , GBMatrix_t * C_init */
		/* , CSR_pol * polys */
		, FILE * fh
	     )
{
	SAFE_MALLOC_DECL(A_init,1,GBMatrix_t);
	SAFE_MALLOC_DECL(C_init,1,GBMatrix_t);
	init(A_init);
	init(C_init);
	SAFE_MALLOC_DECL(polys,1,CSR_pol);

	/* uint32_t i = 0 ; */

	/* sanity */
	assert(fh);

	/* format */
	SAFE_READ_DECL_V(b,uint32_t,fh);
	/* XXX set elem_s here and C++-ise*/

	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,uint32_t,fh);
	SAFE_READ_DECL_V(n,uint32_t,fh);
	SAFE_READ_DECL_V(mod,elem_s,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,uint64_t,fh);

	fprintf(stderr," Mat is %u x %u (sparsity : %f) mod %lu\n",m,n,(double)(nnz)/(double)(m*n),(int64_t)mod);

	A_init->col = C_init->col = n ;
	A_init->mod = C_init->mod = mod ;

	/* READ in ZO start */
	SAFE_READ_DECL_P(start_zo,m+1,uint64_t,fh);

	/* largest row */
	/* uint32_t big_row = 0 ;
	   for (i = 0 ; i < m ; ++i)
	   big_row = max(big_row,start_zo[i+1]-start_zo[i]);    */

	/* pol/zo correspondance */
	SAFE_READ_DECL_P(map_zo_pol,m,uint32_t,fh);


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,uint64_t,fh);
	/* fprintf(stderr,"colid_size %u\n",colid_size); */
	SAFE_READ_DECL_P(buffer,colid_size,uint32_t,fh); /* buffer has the matrix */
	SAFE_MALLOC_DECL(colid_zo,nnz,uint32_t); /* colid expands buffer */

	expandColid(buffer,colid_size,colid_zo,nnz);
	free(buffer);

	/* uint32_t last_start = 0 ; */
	/* i = 0 ; */
	SAFE_CALLOC_DECL(pivots_data,n,uint32_t);

	uint32_t k_dim = getSparsestRows(colid_zo,start_zo,m,n, pivots_data);
	/* i =  0 ; for( ; i < m ; ++i) fprintf(stderr, "%u ", pivots[i]); fprintf(stderr,"\n"); i = 0; */
	assert(k_dim);
	SAFE_MALLOC_DECL(pivots,k_dim,uint32_t);
	SAFE_MALLOC_DECL(nonpiv,n-k_dim,uint32_t);

	uint32_t pivots_size, nonpiv_size ;

	uint32_t k_split = createPivots(pivots,&pivots_size,nonpiv,&nonpiv_size,k_dim,pivots_data,n);

	assert(k_dim >= pivots_size);
	assert(n-k_dim >= nonpiv_size);

	SAFE_REALLOC(pivots,pivots_size,uint32_t);
	SAFE_REALLOC(nonpiv,nonpiv_size,uint32_t);

	splitHorizontal(A_init,C_init,colid_zo,start_zo,m,nnz,pivots,pivots_size,map_zo_pol);

	free(pivots);

	A_init->row = k_dim;
	C_init->row = m - k_dim;

	assert(nnz == A_init->nnz + C_init->nnz);


	/* size of nnz for pols: */
	SAFE_READ_V(polys->nb,uint32_t,fh);

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_P(polys->start_pol,polys->nb+1,uint32_t,fh);

	/* XXX what if elem_s == elem_t ??? */
	SAFE_READ_DECL_P(polys_vals_pol,polys->start_pol[polys->nb],elem_s,fh);
	SAFE_MEMCPY_CVT(polys->vals_pol,elem_t,polys_vals_pol,polys->start_pol[polys->nb]);

	splitVertical(A_init,C_init,nonpiv,nonpiv_size,polys,A,B,C,D);

	return nonpiv ;
}


uint32_t RowReduce_int32_t ( int32_t p, int32_t * A, uint32_t m, uint32_t n, uint32_t lda) ;
uint32_t RowReduce_double  ( double p, double * A, uint32_t m, uint32_t n, uint32_t lda) ;
void     Freduce_double(double p, double * A, uint32_t n);
void     Finit_double  (double p, const double * A, double * B, uint32_t n);
int cblas_daxpy(const int N, const double alpha, const double * X, const int incX, double * Y, const int incY);

void reduce( GBMatrix_t * A
		, DenseMatrix_t * B
		/* , GBMatrix_t * Bt */
		, GBMatrix_t * C
		, DenseMatrix_t * D )
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
	int32_t blk = A->matrix_nb-1 ;
	uint32_t ldb = B->col ;
	uint32_t N = B->col ;
	elem_t p = A->mod ;
	for ( ; blk >=0  ; --blk) {
		CSR_zo * Ad = & (A->matrix_zo[blk]);
		uint32_t M = Ad->row ;
		int32_t i = M -1;
		/* B = A^(-1) B */
		uint64_t i_offset = blk * MAT_ROW_BLOCK;
		for ( ; i >= 0 ; --i) {
			/* fprintf(stderr,"row %u \n",i); */
			uint64_t jz   ;
			uint64_t k  ;
			/* TODO invert jz/k ? */
			assert( (elem_t)-1<1); /* unsigned will fail */
#if 1
			for ( jz = Ad->start_zo[i]+1 ; jz < Ad->start_zo[i+1] ; ++jz) {
				uint64_t j = Ad->colid_zo[jz];
				cblas_daxpy(N,-Ad->data[jz],B->data+j*ldb,1,B->data+(i_offset+i)*ldb,1);
				/* for ( k = 0 ; k < N ; ++k) {
					B->data[(i_offset+i)*ldb+k] -=  Ad->data[jz] * B->data[j*ldb+k] ;;
					|+ Mjoin(fmodin,elem_t)(B->data+(i_offset+i)*ldb+k, p) ; +|
				}                                                                             */
			}
#else
			for ( k = 0 ; k < N ; ++k) {
				for ( jz = Ad->start_zo[i]+1 ; jz < Ad->start_zo[i+1] ; ++jz) {
					uint64_t j = Ad->colid_zo[jz];
					B->data[(i_offset+i)*ldb+k] -=  Ad->data[jz] * B->data[j*ldb+k] ;;
					/* Mjoin(fmodin,elem_t)(B->data+(i_offset+i)*ldb+k, p) ; */
				}
			}

#endif
			Mjoin(Freduce,elem_t)(p,B->data+(i_offset+i)*ldb,N);
			/* for ( k = 0 ; k < N ; ++k) {
				Mjoin(fmodin,elem_t)(B->data+(i_offset+i)*ldb+k, p) ;
			}                                                                */
#if 1
			assert(Ad->data[Ad->start_zo[i]] == 1);
#else
			elem_t d = Ad->data[Ad->start_zo[i]] ;
			elem_t di = Mjoin(invert,double)(d,p);
			/* fprintf(stderr,"diag inv elt "); */
			/* Mjoin(print,elem_t)(d); */
			/* fprintf(stderr," -> "); */
			/* Mjoin(print,elem_t)(di); */
			/* fprintf(stderr,"\n"); */
			k = 0;
			for ( ; k < N ; ++k) {
				B->data[(i_offset+i)*ldb+k] *= di ;
				/* B->data[(i_offset+i)*ldb+k] %= p ; */
				Mjoin(fmodin,elem_t)(B->data+(i_offset+i)*ldb+k, p) ;
			}
#endif
		}
	}

	blk = C->matrix_nb-1 ;
	for ( ; blk >=0  ; --blk) {
		CSR_zo * Cd = &(C->matrix_zo[blk]);
		uint32_t i_offset = blk * MAT_ROW_BLOCK;
		int32_t i  ;
		uint32_t ldd = D->col ;
		/* D = D - C . B */
		assert( (elem_t)-1<1); /* unsigned will fail */
		uint64_t jz,k  ;
		for ( i = 0 ; i < (int32_t)Cd->row ;  ++i) {
			uint32_t j = 0 ;
#if 1
			for ( jz = Cd->start_zo[i]; jz < Cd->start_zo[i+1] ; ++jz ) {
				k = Cd->colid_zo[jz];
				/* cblas_daxpy */
				cblas_daxpy(B->col,-Cd->data[jz],B->data+k*ldb,1,D->data+(i_offset+i)*ldd,1);
				/* for ( j = 0 ; j < B->col ; ++j) {
					D->data[(i_offset+i)*ldd+j] -= ( Cd->data[jz] * B->data[k*ldb+j] ) ;
				}                                                                               */
			}
#else
			for ( j = 0 ; j < B->col ; ++j) {
				for ( jz = Cd->start_zo[i]; jz < Cd->start_zo[i+1] ; ++jz ) {
				k = Cd->colid_zo[jz];
					D->data[(i_offset+i)*ldd+j] -= ( Cd->data[jz] * B->data[k*ldb+j] ) ;
				}
			}
#endif
			/* for ( j = 0 ; j < B->col ; ++j) {
				Mjoin(fmodin,elem_t)(D->data+(i_offset+i)*ldd+j, p) ;
			}                                                                */
			Mjoin(Freduce,elem_t)(p,D->data+(i_offset+i)*ldd, B->col) ;
		}
	}
}

void reduce_fast( GBMatrix_t * A
		, DenseMatrix_t * B
		/* , GBMatrix_t * Bt */
		, GBMatrix_t * C
		, DenseMatrix_t * D )
{
	int32_t blk = A->matrix_nb-1 ;
	assert(blk == 0);
	uint32_t ldb = B->col ;
	uint32_t ldd = D->col ;
	uint32_t N = B->col ;
	elem_t p = A->mod ;
	CSR_zo * Ad = &(A->matrix_zo[blk]);
	CSR_zo * Cd = &(C->matrix_zo[blk]);
	elem_t * Bd = B->data ;
	elem_t * Dd = D->data ;

	SAFE_MALLOC_DECL(temp_D,D->col,elem_t);
	uint32_t i ;
	assert(Cd->row == D->row);
	assert(Cd->col == A->col);
	assert(D->col == B->col);
	SAFE_MALLOC_DECL(temp_C,Cd->col,elem_t);
	for (i = 0 ; i < Cd->row ; ++i) {
		uint64_t jz  ;
		for ( jz = 0 ; jz < Cd->col ; ++jz)
			temp_C[jz] = 0. ;
		for ( jz = Cd->start_zo[i] ; jz < Cd->start_zo[i+1] ; ++jz) {
			assert(Cd->colid_zo[jz]<Cd->col);
			temp_C[Cd->colid_zo[jz]] = Cd->data[jz] ;
		}
		cblas_dcopy(D->col,Dd+i*ldd,1,temp_D,1);

		uint32_t j ;
		for ( j = 0 ; j < Cd->col ; ++j) {
			elem_t tc = Mjoin(fmod,elem_t)(-temp_C[j],p) ;
			if (tc != 0) {
				/* temp_C -= temp_C[j] * A[j] */
				for ( jz = Ad->start_zo[j] ; jz < Ad->start_zo[j+1] ; ++jz) {
					assert(Ad->colid_zo[jz]<Ad->col);
					temp_C[Ad->colid_zo[jz]] += tc * Ad->data[jz] ;
				}

				/* temp_D -= temp_C[j] * B[j] */
				cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1);
			}
		}
		Mjoin(Finit,elem_t)(p, temp_D, Dd+i*ldd, D->col) ;
	}
}

void echelonD( GBMatrix_t * A
	     , DenseMatrix_t * D)
{

	uint32_t r = Mjoin(RowReduce,elem_t)(D->mod,D->data,D->row,D->col,D->col);

	fprintf(stderr,"result : %u\n",r+A->row);
}

#endif /* __GB_io_H */

/* vim: set ft=c: */
