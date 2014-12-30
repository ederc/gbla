#ifndef __GB_io_H
#define __GB_io_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "matrix.h"
#include "field_ops.h"
#include "selecter.h"

/**
 * The matrix we have as input is "almost upper triangular", ie below every first
 * element in a row, there may be non zero elements on the next rows below.
 * A non pivot column can be permuted.
 * There is no empty line.
 *
 */



/* buffer:
*/

taille_t getSparsestRows_fast(
		taille_t   * colid
		, index_t  * start
		, taille_t   row
		, taille_t   col
		, taille_t * pivots_data /* pivots_data[j] is the sparsest row with j as pivot */
		)
{
	taille_t i = 0;
	taille_t k_dim = 0;

	SAFE_MALLOC_DECL(creux_v,row,taille_t);
	for ( i = 0 ; i < row ; ++i)
		creux_v[i] = start[i+1]-start[i];

	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (taille_t)(-1);
	}

	for ( i = 0 ; i < row ; ++i ) {
		taille_t pivot_j = colid[start[i]] ;     /* first column row i */
		taille_t creux   = creux_v[i] ; /* length of row i */
		assert(pivot_j < col);
		taille_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (taille_t)(-1)) {
			pivots_data[pivot_j] = i ;
			++k_dim;
		}
		else  {
			taille_t old_creux = creux_v[old_i];
			if (old_creux > creux) { /* this row is sparser */
				pivots_data[pivot_j] = i ;
			}
		}
	}

	free(creux_v);
	fprintf(stderr,"  -- number of pivots : %u\n",k_dim);
	return k_dim ;
}

taille_t getSparsestRows(
		taille_t   * colid
		, index_t  * start
		, taille_t   row
		, taille_t   col
		, taille_t * pivots_data /* pivots_data[j] is the sparsest row with j as pivot */
		)
{
	taille_t i = 0;
	taille_t k_dim = 0;

	/*  get k */

	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (taille_t)(-1);
	}

	for ( i = 0 ; i < row ; ++i ) {
		taille_t pivot_j = colid[start[i]] ;     /* first column row i */
		assert(pivot_j < col);
		taille_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (taille_t)(-1)) {
			pivots_data[pivot_j] = i ;
			++k_dim;
		}
	}
	printf("good k : %u\n",k_dim);

	/*  get number of columms to pivot */

	taille_t goodk = 0 ;
	taille_t badk = 0;
	for ( i = 0 ; i < col ; ++i) {
		if (pivots_data[i] != (taille_t)(-1)) {
			++goodk ;
			if (goodk == k_dim) {
				break;
			}
		}
		else {
			++badk ;
		}
	}
	printf("bad k : %u\n",badk);

	/*  get best sparsity  for A */

	SAFE_CALLOC_DECL(creux_v,row,taille_t);

	for ( i = 0 ; i < row ; ++i) {
		index_t jz ;
		for (jz = start[i] ; jz < start[i+1] ; ++jz) {
			if (colid[jz] >= (k_dim+badk))
					break;
			if (pivots_data[colid[jz]] != (taille_t)(-1)) {
				creux_v[i] += 1 ;
			}
		}
	}

	/*  reinit (not necessary ? */


	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (taille_t)(-1);
	}

	for ( i = 0 ; i < row ; ++i ) {
		taille_t pivot_j = colid[start[i]] ;     /* first column row i */
		taille_t creux   = creux_v[i] ; /* length of row i */
		assert(pivot_j < col);
		taille_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (taille_t)(-1)) {
			pivots_data[pivot_j] = i ;
			/* ++k_dim; */
		}
		else  {
			taille_t old_creux = creux_v[old_i];
			if (old_creux > creux) { /* this row is sparser */
				pivots_data[pivot_j] = i ;
			}
		}
	}


	free(creux_v);
	fprintf(stderr,"  -- number of pivots : %u\n",k_dim);
	return k_dim ;
}

void createPivots(
		taille_t   * pivots /* list of pivot columns */
		, taille_t * pivots_size
		, taille_t * nonpiv /* list of not pivot columns */
		, taille_t * nonpiv_size
		, taille_t   k_dim
		, taille_t * pivots_data
		, taille_t   n)
{
	taille_t j ;
	*nonpiv_size =  0;
	*pivots_size = 0 ;
	for ( j = 0 ; j < n ; ++j) {
		if (pivots_data[j] == (taille_t) -1) {
			nonpiv[(*nonpiv_size)++] = j ;
			assert(*nonpiv_size <= n-k_dim);
		}
		else {
			pivots[(*pivots_size)++] = pivots_data[j] ;
			if (k_dim == *pivots_size)
				break;
		}
	}
	assert(k_dim == *pivots_size);

	return;
}

void splitHorizontal(
		GBMatrix_t   * A
		, GBMatrix_t * C
		, taille_t   * colid
		, index_t    * start
		, taille_t    row
#ifndef NDEBUG
		, index_t     nnz
#endif
		, taille_t   * pivots
		, taille_t     pivots_size
		, taille_t   * map_zo_pol
		)
{
	taille_t last_pivot=0;
	taille_t i = 0 ;
	assert(start[row] == nnz);
	for ( ; i < row ; ++i)  {
		assert(start[i] < nnz);
		if ( (last_pivot < pivots_size) && (pivots[last_pivot] == i) ) {
			appendRow(A,colid+start[i],start[i+1]-start[i],map_zo_pol[i]);
			++last_pivot;
		}
		else {
			appendRow(C,colid+start[i],start[i+1]-start[i],map_zo_pol[i]);
		}
	}

#ifndef NDEBUG
	assert(A->nnz + C->nnz == nnz);
#endif
	return ;
}


void expandColid(
	       	const uint32_t   * compress
		, index_t          size_compressed
		, taille_t       * expand
#ifndef NDEBUG
		, index_t          size_expand
		, taille_t         n
#endif
		)
{
	uint32_t mask = (1<<31);
	index_t i = 0 ;
	index_t j = 0 ;
	taille_t col ;
	for ( ; i < size_compressed ;) {
		col = compress[i++] ;
		if (col & mask) {
			expand[j++] = col ^ mask ;
		}
		else {
			taille_t k = 0 ;
			for (; k < compress[i] ;++k) {
				expand[j++] = col + k;
			}
			++i;
		}
	}
	assert(j == size_expand);
	assert(i == size_compressed);
#ifndef NDEBUG
	for ( i = 0 ; i < size_expand ; ++i) {
		assert(expand[i] < n);
	}
#endif
}


void splitVerticalUnit(
		const GBMatrix_t * A_init,
		const taille_t   * nonpiv,
		const taille_t     nonpiv_size,
		const CSR_pol    * polys,
		GBMatrix_t       * A,
		DNS    * B
		)
{

	taille_t max_col = A->col + nonpiv_size ;

	/* Init dense */

	SAFE_CALLOC(B->ptr,(index_t)B->row*(index_t)B->col,elem_t);


	/* Init sparse */

	taille_t j ;
	for ( j = 0 ; j < A_init->sub_nb ; ++j) {
		const CSR * A_k = &(A_init->sub[j]) ;
#ifndef NDEBUG
		checkMatUnit(A_k);
#endif
		appendMatrix(A);
		CSR * Ad = getLastMatrix(A);
		Ad->row = A_k->row ;
		Ad->col = A->col ;
		Ad->nnz = 0 ;
		/* A->row += Ad->row ; */
		assert(A->col == Ad->col);

		SAFE_REALLOC(Ad->start,Ad->row+1,index_t);
		Ad->start[0] = 0 ;
		SAFE_MALLOC(Ad->colid,A_k->nnz,taille_t);
		SAFE_MALLOC(Ad->data ,A_k->nnz,elem_t);
		assert(Ad->map_zo_pol==NULL);


		taille_t ldb = B->col ;
		taille_t here  = 0; /* shift */
		index_t there = 0;
		here = 0;
		taille_t i ;
		index_t jz ;
		for (i = 0 ; i <= Ad->row ; ++i) {
			Ad->start[i] = 0 ;
		}

		for (i = 0 ; i < A_k->row ; ++i) {
			index_t i_offset =  j * MAT_ROW_BLOCK + i ;
			taille_t start_p = polys->start_pol[ A_k->map_zo_pol[i] ] ;
			elem_t * d = polys->data_pol+start_p ;
			for (jz = A_k->start[i] ; jz < A_k->start[i+1] ; ++jz) {
				taille_t k = A_k->colid[jz]  ;
				assert(k < A_k->col);
				if ( k  < max_col ) {
					while (here < nonpiv_size && k > nonpiv[here]) {
						here ++ ;
					}
					if (here < nonpiv_size && k == nonpiv[here]) {
						B->ptr[i_offset*(index_t)ldb+(index_t)here] = *d; /* XXX */
						B->nnz += 1 ;
					}
					else {
						Ad->start[i+1] += 1 ;
						Ad->colid[there] = k-here ;
						Ad->data[there++] = *d ;
					}
				}
				else {
					assert(nonpiv_size+k-max_col < B->col);
					B->ptr[i_offset*(index_t)ldb+(index_t)(nonpiv_size+k-max_col)] = *d ; /* XXX */
					B->nnz += 1 ;
				}
				++d;
			}
			assert(here <= nonpiv_size);
			here = 0 ;
		}
		for ( i = 1 ; i < Ad->row ; ++i) {
			Ad->start[i+1] += Ad->start[i] ;
		}
		Ad->nnz = Ad->start[Ad->row] ;
		A->nnz += Ad->nnz ;
		assert(A->nnz);
		assert(Ad->nnz == there);
		/* realloc to smaller */
		SAFE_REALLOC(Ad->colid,Ad->nnz,taille_t);
		SAFE_REALLOC(Ad->data ,Ad->nnz,elem_t);

	}

	assert(A_init->nnz == A->nnz + B->nnz);
}

void splitVertical(
		const GBMatrix_t   * A_init
		, const GBMatrix_t * C_init
		, const taille_t   * nonpiv
		, const taille_t     nonpiv_size
		, const CSR_pol    * polys
		, GBMatrix_t    * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D
		)
{
	A->row = B->row = A_init->row ;
	C->row = D->row = C_init->row ;
	A->mod = B->mod = C->mod = D->mod = A_init->mod;
	A->col = C->col = A_init->row ;
	B->col = D->col = A_init->col - A_init->row ;

	assert(A->row);
	assert(C->row);

	splitVerticalUnit(A_init,nonpiv,nonpiv_size,polys,A,B);
	splitVerticalUnit(C_init,nonpiv,nonpiv_size,polys,C,D);

#ifndef NDEBUG
	checkMat(A);
	checkMat(C);
#endif
}


/* matrix reader and row splitter
 * A_init [out] the top part with upper triangular left
 * B_init [out] the bottom part
 * polys [out] the polynomials used in the matrices (shared)
 * fh [in] the file in appropriate format
 */
taille_t * readFileSplit(
		GBMatrix_t    * A,
		DNS * B,
		GBMatrix_t    * C,
		DNS * D
		, FILE        * fh
		)
{
	SAFE_MALLOC_DECL(A_init,1,GBMatrix_t);
	SAFE_MALLOC_DECL(C_init,1,GBMatrix_t);
	initSparse(A_init);
	initSparse(C_init);
	SAFE_MALLOC_DECL(polys,1,CSR_pol);


	/* sanity */
	assert(fh);

	/* format */
	SAFE_READ_DECL_V(b,uint32_t,fh);
	/* XXX set elem_s here and C++-ise*/
	assert(b== Mjoin(select,elem_s)());

	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,uint32_t,fh);
	/* assert(m < MAT_ROW_BLOCK); */
	SAFE_READ_DECL_V(n,uint32_t,fh);
	SAFE_READ_DECL_V(mod,elem_s,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,uint64_t,fh);

	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %f) mod %lu\n",m,n,nnz,(double)(nnz)/((double)m*(double)n),(int64_t)mod);

	A_init->col = C_init->col = n ;
	A_init->mod = C_init->mod = mod ;

	/* READ in ZO start */
	SAFE_READ_DECL_P(start,m+1,uint64_t,fh);

	/* pol/zo correspondance */
	SAFE_READ_DECL_P(map_zo_pol,m,uint32_t,fh);


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,uint64_t,fh);
	SAFE_READ_DECL_P(buffer,colid_size,uint32_t,fh); /* buffer has the matrix */
	SAFE_MALLOC_DECL(colid,nnz,taille_t); /* colid expands buffer */

	struct timeval tic,tac;
	gettimeofday(&tic,NULL);


	expandColid(buffer,colid_size,colid
#ifndef NDEBUG
			,nnz,n
#endif
		   );
	free(buffer);

	SAFE_CALLOC_DECL(pivots_data,n,taille_t);

	/* taille_t k_dim = getSparsestRows(colid,start,m,n, pivots_data); */
	taille_t k_dim = getSparsestRows_fast(colid,start,m,n, pivots_data);
	assert(k_dim);
	SAFE_MALLOC_DECL(pivots,k_dim,taille_t);
	SAFE_MALLOC_DECL(nonpiv,n-k_dim,taille_t);

	taille_t pivots_size, nonpiv_size ;

	/* taille_t k_split = */ createPivots(pivots,&pivots_size,nonpiv,&nonpiv_size,k_dim,pivots_data,n);

	free(pivots_data);

	assert(k_dim >= pivots_size);
	assert(n-k_dim >= nonpiv_size);

	SAFE_REALLOC(pivots,pivots_size,taille_t);
	SAFE_REALLOC(nonpiv,nonpiv_size,taille_t);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> introspect       : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

	gettimeofday(&tic,NULL);

	splitHorizontal(A_init,C_init,colid,start,m
#ifndef NDEBUG
			,nnz
#endif
			,pivots,pivots_size,map_zo_pol);
#ifndef NDEBUG
	checkMat(A_init);
	checkMat(C_init);
#endif

	free(pivots);
	free(start);
	free(map_zo_pol);
	free(colid);

	A_init->row = k_dim;
	C_init->row = m - k_dim;

	assert(nnz == A_init->nnz + C_init->nnz);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> split horizontal : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));



	/* size of nnz for pols: */
	SAFE_READ_V(polys->nb,uint32_t,fh);

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_P(polys->start_pol,polys->nb+1,uint32_t,fh);

	/* XXX what if elem_s == elem_t ??? */
	SAFE_READ_DECL_P(polys_data_pol,polys->start_pol[polys->nb],elem_s,fh);

	gettimeofday(&tic,NULL);
	SAFE_MEMCPY_CVT(polys->data_pol,elem_t,polys_data_pol,polys->start_pol[polys->nb]);
	free(polys_data_pol);

	splitVertical(A_init,C_init,nonpiv,nonpiv_size,polys,A,B,C,D);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> split vertical   : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));


	freeMat(A_init);
	free(A_init);
	freeMat(C_init);
	free(C_init);
	freePol(polys);
	free(polys);
	return nonpiv ;
}


taille_t RowReduce_int32_t ( int32_t p, int32_t * A, taille_t m, taille_t n, taille_t lda) ;
taille_t RowReduce_double  ( double p, double * A, taille_t m, taille_t n, taille_t lda) ;
void     Freduce_double(double p, double * A, index_t n);
void     Finit_double  (double p, const double * A, double * B, index_t n);
int cblas_daxpy(const int N, const double alpha, const double * X, const int incX, double * Y, const int incY);
int cblas_dscal(const int N, const double alpha, const double * X, const int incX);
int cblas_dcopy(const int N, const double * X, const int incX, double * Y, const int incY);

void reduce(
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
	taille_t k  ;
	taille_t ldb = B->col ;
	taille_t N = B->col ;
	elem_t p = A->mod ;
	taille_t i ;
	index_t jz  ;
	index_t kz ;
	for ( k = A->sub_nb ; k --  ; ) {
		CSR * Ad = & (A->sub[k]);
		taille_t M = Ad->row ;
		/* B = A^(-1) B */

		for ( i = M ; i--    ; ) {
			index_t i_offset = k * MAT_ROW_BLOCK + i;
			assert( (elem_t)-1<1); /* unsigned will fail */
/* #ifdef _OPENMP */
/* #pragma omp parallel for private (kz) schedule (static,5) */
/* #endif */
			for ( jz = Ad->start[i]+1 ; jz < Ad->start[i+1] ; ++jz) {
				kz = Ad->colid[jz];
				cblas_daxpy(N,-Ad->data[jz],B->ptr+kz*(index_t)ldb,1,B->ptr+i_offset*(index_t)ldb,1);
			}
			Mjoin(Freduce,elem_t)(p,B->ptr+i_offset*(index_t)ldb,N);
			assert(Ad->data[Ad->start[i]] == 1);
		}
	}


	for ( k = 0 ; k < C->sub_nb  ; ++k ) {
		CSR * Cd = &(C->sub[k]);
		taille_t ldd = D->col ;
		/* D = D - C . B */
		assert( (elem_t)-1<1); /* unsigned will fail */
#ifdef _OPENMP
#pragma omp parallel for private (kz,jz) schedule (static,5)
#endif

		for ( i = 0 ; i < Cd->row ;  ++i) {
			index_t i_offset = k * MAT_ROW_BLOCK + i;
			for ( jz = Cd->start[i]; jz < Cd->start[i+1] ; ++jz ) {
				kz = Cd->colid[jz];
				cblas_daxpy(N,-Cd->data[jz],B->ptr+kz*(index_t)ldb,1,D->ptr+i_offset*(index_t)ldd,1);
			}
			Mjoin(Freduce,elem_t)(p,D->ptr+i_offset*(index_t)ldd, N) ;
		}
	}
}


void spaxpy(
		elem_t           tc,
		const elem_t   * A,
		const taille_t   nb,
		const taille_t * colid,
		elem_t         * B)
{
	taille_t jz ;
	for (jz = 0 ; jz < nb ; ++jz) {
		B[colid[jz]] += tc * A[jz] ;
	}
}

void spaxpy2(
		elem_t         tc,
		elem_t         td,
		const elem_t   * A,
		const taille_t   nb,
		const taille_t * colid,
		elem_t         * B,
		const taille_t   ld
	    )
{
	taille_t jz ;
	for (jz = 0 ; jz < nb ; ++jz) {
		B[colid[jz]]    += tc * A[jz] ;
		B[colid[jz]+ld] += td * A[jz] ;
	}
}

void spaxpyn(
		elem_t         *  diag,
		const elem_t   * A,
		const taille_t   nb,
		const taille_t * colid,
		elem_t         * B,
		const taille_t * jump,
		const taille_t   js,
		const taille_t   ld
		)
{
	taille_t jz,ii ;
	for (jz = 0 ; jz < nb ; ++jz) {
		for ( ii = 0 ; ii < js ; ++ii)
			B[colid[jz]+ld*jump[ii]] += diag[ii] * A[jz] ;
	}
}



void reduce_fast(
		GBMatrix_t      * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D )
{
	taille_t ldb = B->col ;
	taille_t ldd = D->col ;
	elem_t p = A->mod ;
	elem_t * Bd = B->ptr;
	elem_t * Dd = D->ptr;
	taille_t i ;

	/* XXX convert B to sparse */
	SAFE_MALLOC_DECL(Bsp,1,GBMatrix_t);
	initSparse(Bsp);

	Bsp->mod = B->mod;
	Bsp->col = B->col;

	SAFE_MALLOC_DECL(col_buf,B->col,taille_t);
	SAFE_MALLOC_DECL(dat_buf,B->col,elem_t);

	for ( i = 0 ; i < B ->row ; ++i) {
		taille_t j ;
		taille_t length = 0 ;
		for (j = 0 ; j < B->col ; ++j) {
			if (Bd[i*ldb+j] != 0) {
				col_buf[length]   = j ;
				dat_buf[length++] = Bd[(index_t)i*(index_t)ldb+(index_t)j];
			}
		}
		appendRowData(Bsp,col_buf,length,dat_buf);
	}

	assert(Bsp->nnz == B->nnz);
	assert(Bsp->row == B->row);
	assert(Bsp->col == B->col);

	free(col_buf);
	free(dat_buf);

#ifndef NDEBUG
	checkMat(Bsp);
#endif

	taille_t blk = 16 ;

	index_t jz  ;
#ifndef _OPENMP
	SAFE_MALLOC_DECL(temp_D,(index_t)D->col*(index_t)blk,elem_t);
	SAFE_MALLOC_DECL(temp_C,(index_t)C->col*(index_t)blk,elem_t);
#endif

	/* SAFE_MALLOC_DECL(jump,blk,taille_t); */
	/* SAFE_MALLOC_DECL(diag,blk,elem_t); */
	taille_t k ;
	assert(D->col == B->col);
	for (k = 0 ; k < C->sub_nb ; ++k) {
		CSR * C_k = &(C->sub[k]) ;
#ifdef _OPENMP
#pragma omp parallel for private (jz) schedule (static,5)
#endif

		for ( i = 0 ; i < C_k->row ; i += blk ) {
#ifdef _OPENMP
			SAFE_MALLOC_DECL(temp_D,(index_t)D->col*(index_t)blk,elem_t);
			SAFE_MALLOC_DECL(temp_C,(index_t)C->col*(index_t)blk,elem_t);
#endif

			taille_t blk_i = min(blk,C_k->row - i);
			index_t i_offset = k*MAT_ROW_BLOCK + i ;
			cblas_dscal(C->col*blk_i,0.,temp_C,1); /* XXX blk * col < 2^32 */
			taille_t ii = 0 ;
			for ( ii = 0 ; ii < blk_i ; ++ii) {
				for ( jz = C_k->start[i+ii] ; jz < C_k->start[i+ii+1] ; ++jz) {
					assert(C_k->colid[jz] < C_k->col);
					temp_C[ii*(C->col)+C_k->colid[jz]] = C_k->data[jz] ;
				}
			}
			cblas_dcopy(D->col*blk_i,Dd+i_offset*(index_t)ldd,1,temp_D,1);

			taille_t j ;
#if 0
			for ( j = 0 ; j < C->col ; ++j) {
				/* XXX invert loops ? */
				for ( ii = 0 ; ii < blk_i ; ii += 2 ) {
					elem_t tc = Mjoin(fmod,elem_t)(-temp_C[ii*C->col+j],p) ;
					if (tc != 0.) {
						/* temp_C -= temp_C[j] * A[j] */
						taille_t jj = j%MAT_ROW_BLOCK ;
						taille_t kk = j/MAT_ROW_BLOCK ;
						CSR * A_k = &(A->sub[kk]);
						taille_t sz = (taille_t)(A_k->start[jj+1]-A_k->start[jj]);
						assert(kk*MAT_ROW_BLOCK+jj == j);
						spaxpy(tc,A_k->data+A_k->start[jj],
								sz,
								A_k->colid+A_k->start[jj]
								,temp_C+ii*(C->col));

						/* temp_D -= temp_C[j] * B[j] */

						/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
						CSR * B_k = &(Bsp->sub[kk]);
						sz = (taille_t)(B_k->start[jj+1]-B_k->start[jj]) ;
						spaxpy(tc,B_k->data+B_k->start[jj],
								sz,
								B_k->colid+B_k->start[jj],
								temp_D+ii*(D->col));
					}
				}
			}

#endif

#if 1
			for ( j = 0 ; j < C->col ; ++j) {
				/* XXX invert loops ? */
				for ( ii = 0 ; ii < blk_i ; ii += 2 ) {
					elem_t tc = Mjoin(fmod,elem_t)(-temp_C[ii*C->col+j],p) ;
					elem_t td = Mjoin(fmod,elem_t)(-temp_C[(ii+1)*C->col+j],p) ;
					/* temp_C -= temp_C[j] * A[j] */
					if (tc != 0.) {
						if (td != 0.) {
							taille_t jj = j%MAT_ROW_BLOCK ;
							taille_t kk = j/MAT_ROW_BLOCK ;
							CSR * A_k = &(A->sub[kk]);
							taille_t sz = (taille_t)(A_k->start[jj+1]-A_k->start[jj]);
							assert(kk*MAT_ROW_BLOCK+jj == j);
							spaxpy2(tc,td,A_k->data+A_k->start[jj],
									sz,
									A_k->colid+A_k->start[jj]
									,temp_C+ii*(C->col),C->col);

							/* temp_D -= temp_C[j] * B[j] */

							/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
							CSR * B_k = &(Bsp->sub[kk]);
							sz = (taille_t)(B_k->start[jj+1]-B_k->start[jj]) ;
							spaxpy2(tc,td,B_k->data+B_k->start[jj],
									sz,
									B_k->colid+B_k->start[jj],
									temp_D+ii*(D->col),D->col);
						}
						else {
							/* temp_C -= temp_C[j] * A[j] */
							taille_t jj = j%MAT_ROW_BLOCK ;
							taille_t kk = j/MAT_ROW_BLOCK ;
							CSR * A_k = &(A->sub[kk]);
							taille_t sz = (taille_t)(A_k->start[jj+1]-A_k->start[jj]);
							assert(kk*MAT_ROW_BLOCK+jj == j);
							spaxpy(tc,A_k->data+A_k->start[jj],
									sz,
									A_k->colid+A_k->start[jj]
									,temp_C+ii*(C->col));

							/* temp_D -= temp_C[j] * B[j] */

							/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
							CSR * B_k = &(Bsp->sub[kk]);
							sz = (taille_t)(B_k->start[jj+1]-B_k->start[jj]) ;
							spaxpy(tc,B_k->data+B_k->start[jj],
									sz,
									B_k->colid+B_k->start[jj],
									temp_D+ii*(D->col));

						}
					}
					else if (td != 0) {
						/* temp_C -= temp_C[j] * A[j] */
						taille_t jj = j%MAT_ROW_BLOCK ;
						taille_t kk = j/MAT_ROW_BLOCK ;
						CSR * A_k = &(A->sub[kk]);
						taille_t sz = (taille_t)(A_k->start[jj+1]-A_k->start[jj]);
						assert(kk*MAT_ROW_BLOCK+jj == j);
						spaxpy(td,A_k->data+A_k->start[jj],
								sz,
								A_k->colid+A_k->start[jj]
								,temp_C+(ii+1)*(C->col));

						/* temp_D -= temp_C[j] * B[j] */

						/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
						CSR * B_k = &(Bsp->sub[kk]);
						sz = (taille_t)(B_k->start[jj+1]-B_k->start[jj]) ;
						spaxpy(td,B_k->data+B_k->start[jj],
								sz,
								B_k->colid+B_k->start[jj],
								temp_D+(ii+1)*(D->col));

					}
				}
			}
#endif

#if 0

			for ( j = 0 ; j < C->col ; ++j) {
				taille_t nbnz = 0;
				for ( ii = 0 ; ii < blk_i ; ++ii) {
					elem_t tc = Mjoin(fmod,elem_t)(-temp_C[ii*C->col+j],p) ;
					if (tc != 0.) {
						diag[nbnz]   = tc ;
						jump[nbnz++] = ii ;
					}
				}
				taille_t jj = j%MAT_ROW_BLOCK ;
				taille_t kk = j/MAT_ROW_BLOCK ;
				CSR * A_k = &(A->sub[kk]);
				taille_t sz = (taille_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLOCK+jj == j);
				spaxpyn(diag,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C,jump,nbnz,C->col);

				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(Bsp->sub[kk]);
				sz = (taille_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				spaxpyn(diag,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D,jump,nbnz,D->col);
			}


#endif
			/* Mjoin(Finit,elem_t)(p, temp_D, Dd+i*ldd, D->col) ; */
			cblas_dcopy(D->col*blk_i,temp_D,1,Dd+i_offset*(index_t)ldd,1);

#ifdef _OPENMP
		free(temp_D);
		free(temp_C);
#endif
		}

	}
	freeMat(Bsp);
	free(Bsp);
#ifndef _OPENMP
	free(temp_D);
	free(temp_C);
#endif
	Mjoin(Freduce,elem_t)(p, Dd, (index_t)D->col*(index_t)D->row) ;
}

taille_t echelonD(
		GBMatrix_t      * A
		, DNS * D)
{

	taille_t r = Mjoin(RowReduce,elem_t)(D->mod,D->ptr,D->row,D->col,D->col);

	fprintf(stderr,"  -- residual rank    : %u\n",r);

	return r + A->row;
}

#endif /* __GB_io_H */

/* vim: set ft=c: */
