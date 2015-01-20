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

dimen_t getSparsestRows_fast(
		dimen_t   * colid
		, index_t  * start
		, dimen_t   row
		, dimen_t   col
		, dimen_t * pivots_data /* pivots_data[j] is the sparsest row with j as pivot */
		)
{
	dimen_t i = 0;
	dimen_t k_dim = 0;

	SAFE_MALLOC_DECL(creux_v,row,dimen_t);
	for ( i = 0 ; i < row ; ++i)
		creux_v[i] = start[i+1]-start[i];

	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (dimen_t)(-1);
	}

	for ( i = 0 ; i < row ; ++i ) {
		dimen_t pivot_j = colid[start[i]] ;     /* first column row i */
		dimen_t creux   = creux_v[i] ; /* length of row i */
		assert(pivot_j < col);
		dimen_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (dimen_t)(-1)) {
			pivots_data[pivot_j] = i ;
			++k_dim;
		}
		else  {
			dimen_t old_creux = creux_v[old_i];
			if (old_creux > creux) { /* this row is sparser */
				pivots_data[pivot_j] = i ;
			}
		}
	}

	free(creux_v);
	fprintf(stderr,"  -- number of pivots : %u\n",k_dim);
	return k_dim ;
}

dimen_t getSparsestRows(
		dimen_t   * colid
		, index_t  * start
		, dimen_t   row
		, dimen_t   col
		, dimen_t * pivots_data /* pivots_data[j] is the sparsest row with j as pivot */
		)
{
	dimen_t i = 0;
	dimen_t k_dim = 0;

	/*  get k */

	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (dimen_t)(-1);
	}

	for ( i = 0 ; i < row ; ++i ) {
		dimen_t pivot_j = colid[start[i]] ;     /* first column row i */
		assert(pivot_j < col);
		dimen_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (dimen_t)(-1)) {
			pivots_data[pivot_j] = i ;
			++k_dim;
		}
	}
	printf("good k : %u\n",k_dim);

	/*  get number of columms to pivot */

	dimen_t goodk = 0 ;
	dimen_t badk = 0;
	for ( i = 0 ; i < col ; ++i) {
		if (pivots_data[i] != (dimen_t)(-1)) {
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

	SAFE_CALLOC_DECL(creux_v,row,dimen_t);

	for ( i = 0 ; i < row ; ++i) {
		index_t jz ;
		for (jz = start[i] ; jz < start[i+1] ; ++jz) {
			if (colid[jz] >= (k_dim+badk))
					break;
			if (pivots_data[colid[jz]] != (dimen_t)(-1)) {
				creux_v[i] += 1 ;
			}
		}
	}

	/*  reinit (not necessary ? */


	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = (dimen_t)(-1);
	}

	for ( i = 0 ; i < row ; ++i ) {
		dimen_t pivot_j = colid[start[i]] ;     /* first column row i */
		dimen_t creux   = creux_v[i] ; /* length of row i */
		assert(pivot_j < col);
		dimen_t old_i = pivots_data[pivot_j] ; /* last row for pivot column */
		if (old_i == (dimen_t)(-1)) {
			pivots_data[pivot_j] = i ;
			/* ++k_dim; */
		}
		else  {
			dimen_t old_creux = creux_v[old_i];
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
		dimen_t   * pivots /* list of pivot columns */
		, dimen_t * pivots_size
		, dimen_t * nonpiv /* list of not pivot columns */
		, dimen_t * nonpiv_size
		, dimen_t   k_dim
		, dimen_t * pivots_data
		, dimen_t   n)
{
	dimen_t j ;
	*nonpiv_size =  0;
	*pivots_size = 0 ;
	for ( j = 0 ; j < n ; ++j) {
		if (pivots_data[j] == (dimen_t) -1) {
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
		, dimen_t   * colid
		, index_t    * start
		, dimen_t    row
#ifndef NDEBUG
		, index_t     nnz
#endif
		, dimen_t   * pivots
		, dimen_t     pivots_size
		, dimen_t   * map_zo_pol
		)
{
	dimen_t last_pivot=0;
	dimen_t i ;
	assert(start[row] == nnz);
#if 0
	dimen_t A_row = 0 ;
	index_t A_nnz = 0 ;
	for ( i = 0 ; i < row ; ++i)  {
		assert(start[i] < nnz);
		if ( (last_pivot < pivots_size) && (pivots[last_pivot] == i) ) {
			A_nnz += start[i+1] - start[i] ;
			A_row += 1 ;
			++last_pivot;
		}
	}
	A->nb_sub = DIVIDE_INTO(A_row,MAT_ROW_BLK);
	A->row = A_row ;
	A->nnz = A_nnz ;
	B->nb_sub = DIVIDE_INTO(A_row,MAT_ROW_BLK);
	B->row = row - A_row ;
	B->nnz = start[row]  - A_nnz ;
#endif

	last_pivot = 0 ;
	for ( i = 0 ; i < row ; ++i)  {
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
	       	const dimen_t * compress
		, index_t          size_compressed
		, dimen_t       * expand
#ifndef NDEBUG
		, index_t          size_expand
		, dimen_t         n
#endif
		)
{
	uint32_t mask = (1<<31);
	index_t i = 0 ;
	index_t j = 0 ;
	dimen_t col ;
	for ( ; i < size_compressed ;) {
		col = compress[i++] ;
		if (col & mask) {
			expand[j++] = col ^ mask ;
		}
		else {
			dimen_t k = 0 ;
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
		const dimen_t   * nonpiv,
		const dimen_t     nonpiv_size,
		const CSR_pol    * polys,
		GBMatrix_t       * A,
		DNS    * B
		)
{

	dimen_t max_col = A->col + nonpiv_size ;

	/* Init dense */

	dimen_t ldb = B->ld;
	SAFE_CALLOC(B->ptr,(index_t)B->row*(index_t)ldb,elemt_t);


	/* Init sparse */

	dimen_t j ;
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
		SAFE_MALLOC(Ad->colid,A_k->nnz,dimen_t);
		SAFE_MALLOC(Ad->data ,A_k->nnz,elemt_t);
		assert(Ad->map_zo_pol==NULL);


		dimen_t here  = 0; /* shift */
		index_t there = 0;
		here = 0;
		dimen_t i ;
		index_t jz ;
		for (i = 0 ; i <= Ad->row ; ++i) {
			Ad->start[i] = 0 ;
		}

		for (i = 0 ; i < A_k->row ; ++i) {
			index_t i_offset =  j * MAT_ROW_BLK + i ;
			dimen_t start_p = polys->start_pol[ A_k->map_zo_pol[i] ] ;
			elemt_t * d = polys->data_pol+start_p ;
			for (jz = A_k->start[i] ; jz < A_k->start[i+1] ; ++jz) {
				dimen_t k = A_k->colid[jz]  ;
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
		SAFE_REALLOC(Ad->colid,Ad->nnz,dimen_t);
		SAFE_REALLOC(Ad->data ,Ad->nnz,elemt_t);

	}

	assert(A_init->nnz == A->nnz + B->nnz);
}

void splitVerticalUnitSparse(
		const GBMatrix_t * A_init,
		const dimen_t   * nonpiv,
		const dimen_t     nonpiv_size,
		const CSR_pol    * polys,
		GBMatrix_t       * A,
		GBMatrix_t       * B
		)
{

	dimen_t max_col = A->col + nonpiv_size ;


	/* Init sparse */

	dimen_t j ;
	for ( j = 0 ; j < A_init->sub_nb ; ++j) {
		const CSR * A_k = &(A_init->sub[j]) ;
#ifndef NDEBUG
		checkMatUnit(A_k);
#endif
		appendMatrix(A);
		appendMatrix(B);
		CSR * Ad = getLastMatrix(A);
		CSR * Bd = getLastMatrix(B);
		Ad->row = A_k->row ;
		Bd->row = A_k->row ;
		Ad->col = A->col ;
		Bd->col = B->col ;
		Ad->nnz = 0 ;
		Bd->nnz = 0 ;
		/* A->row += Ad->row ; */
		assert(A->col == Ad->col);
		assert(B->col == Bd->col);

		SAFE_REALLOC(Ad->start,Ad->row+1,index_t);
		SAFE_REALLOC(Bd->start,Bd->row+1,index_t);
		assert(Ad->map_zo_pol==NULL);


		/* index_t there = 0; */
		dimen_t i ;
		index_t jz ;
		for (i = 0 ; i <= Ad->row ; ++i) {
			Ad->start[i] = 0 ;
		}
		for (i = 0 ; i <= Bd->row ; ++i) {
			Bd->start[i] = 0 ;
		}

		SAFE_CALLOC_DECL(b_row_off,Bd->row,dimen_t);
		SAFE_CALLOC_DECL(b_there,(Bd->row+1),index_t);
		SAFE_CALLOC_DECL(b_la,(Bd->row+1),index_t);


		for (i = 0 ; i < A_k->row ; ++i) {
			assert(b_there[i+1] == 0);
			assert(b_la[i+1] == 0);
			dimen_t here  = 0; /* shift */
			b_there[i+1] += b_there[i] ;
			for (jz = A_k->start[i] ; jz < A_k->start[i+1] ; ++jz) {
				dimen_t k = A_k->colid[jz]  ;
				assert(k < A_k->col);
				if ( k  < max_col ) {
					while (here < nonpiv_size && k > nonpiv[here]) {
						here ++ ;
					}
					if (here < nonpiv_size && k == nonpiv[here]) {
						b_row_off[i] += 1 ;
						Bd->nnz      += 1 ;
					}
					else {
						b_there[i+1]  += 1 ;
					}
				}
				else {
					assert(nonpiv_size+k-max_col < B->col);
					Bd->nnz   += 1 ;
					b_la[i+1] += 1 ;
				}
			}
			assert(here <= nonpiv_size);
		}

		Ad->nnz = A_k->nnz - Bd->nnz ;

		SAFE_MALLOC(Ad->colid,Ad->nnz,dimen_t);
		SAFE_MALLOC(Bd->colid,Bd->nnz,dimen_t);
		SAFE_MALLOC(Ad->data ,Ad->nnz,elemt_t);
		SAFE_MALLOC(Bd->data ,Bd->nnz,elemt_t);

		for (i = 0 ; i < A_k->row ; ++i) {
			b_la[i+1]  =  b_la[i+1] + b_row_off[i] ;
		}
		for (i = 0 ; i < A_k->row ; ++i) {
			b_la[i+1] +=  b_la[i]  ;
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
		/* dimen_t la = 0 ; */
		/* there = 0; */
		for (i = 0 ; i < A_k->row ; ++i) {
			dimen_t here  = 0; /* shift */
			index_t there = b_there[i] ;
			index_t ici =  b_la[i] ;
			index_t la  =  ici + b_row_off[i] ;
			dimen_t start_p = polys->start_pol[ A_k->map_zo_pol[i] ] ;
			elemt_t * d = polys->data_pol+start_p ;
			for (jz = A_k->start[i] ; jz < A_k->start[i+1] ; ++jz) {
				dimen_t k = A_k->colid[jz]  ;
				assert(k < A_k->col);
				if ( k  < max_col ) {
					while (here < nonpiv_size && k > nonpiv[here]) {
						here ++ ;
					}
					if (here < nonpiv_size && k == nonpiv[here]) {
						assert(ici < Bd->nnz);
						assert(here < Bd->col);
						Bd->start[i+1] += 1 ;
						Bd->colid[ici] = here ;
						Bd->data[ici++] = *d ;
						/* B->ptr[i_offset*(index_t)ldb+(index_t)here] = *d; |+ XXX +| */
					}
					else {
						Ad->start[i+1] += 1 ;
						Ad->colid[there] = k-here ;
						Ad->data[there++] = *d ;
					}
				}
				else {
					assert(nonpiv_size+k-max_col < B->col);
					assert(la < Bd->nnz);
					Bd->start[i+1] += 1 ;
					Bd->colid[la] = nonpiv_size+k-max_col ;
					Bd->data[la++] = *d ;
					/* B->ptr[i_offset*(index_t)ldb+(index_t)(nonpiv_size+k-max_col)] = *d ; |+ XXX +| */
				}
				++d;
			}
			assert(here <= nonpiv_size);
		}
		assert(Ad->nnz == b_there[Ad->row]);

		free(b_row_off);
		free(b_there);
		free(b_la);

		for ( i = 1 ; i < Ad->row ; ++i) {
			Ad->start[i+1] += Ad->start[i] ;
		}
		for ( i = 1 ; i < Bd->row ; ++i) {
			Bd->start[i+1] += Bd->start[i] ;
		}

		assert(Ad->nnz ==  Ad->start[Ad->row] );
		assert(Bd->nnz ==  Bd->start[Bd->row] );

		A->nnz += Ad->nnz ;
		B->nnz += Bd->nnz ;
		assert(A->nnz);

	}

	assert(A_init->nnz == A->nnz + B->nnz);
}

void splitVertical(
		const GBMatrix_t   * A_init
		, const GBMatrix_t * C_init
		, const dimen_t   * nonpiv
		, const dimen_t     nonpiv_size
		, const CSR_pol    * polys
		, GBMatrix_t    * A
		, GBMatrix_t    * B
		, GBMatrix_t    * C
		, DNS * D
		)
{
	A->row = B->row = A_init->row ;
	C->row = D->row = C_init->row ;
	A->mod = B->mod = C->mod = D->mod = A_init->mod;
	A->col = C->col = A_init->row ;
	B->col = D->col = A_init->col - A_init->row ;
	D->ld  = ALIGN(D->col);


	assert(A->row);
	assert(C->row);

	splitVerticalUnitSparse(A_init,nonpiv,nonpiv_size,polys,A,B);
	splitVerticalUnit(C_init,nonpiv,nonpiv_size,polys,C,D);

#ifndef NDEBUG
	checkMat(A);
	checkMat(B);
	checkMat(C);
#endif
}


/* matrix reader and row splitter
 * A_init [out] the top part with upper triangular left
 * B_init [out] the bottom part
 * polys [out] the polynomials used in the matrices (shared)
 * fh [in] the file in appropriate format
 */
dimen_t * readFileSplit(
		GBMatrix_t    * A,
		GBMatrix_t    * B,
		GBMatrix_t    * C,
		DNS * D
		, FILE        * fh
		)
{

	struct timeval tic,tac;
	gettimeofday(&tic,NULL);
	SAFE_MALLOC_DECL(A_init,1,GBMatrix_t);
	SAFE_MALLOC_DECL(C_init,1,GBMatrix_t);
	initSparse(A_init);
	initSparse(C_init);
	SAFE_MALLOC_DECL(polys,1,CSR_pol);


	/* sanity */
	assert(fh);

	/* format */
	SAFE_READ_DECL_V(b,uint32_t,fh);
	/* XXX set elemt_s here and C++-ise*/
	assert((b ^ VERMASK) == Mjoin(select,elemt_s)());

	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,uint32_t,fh);
	/* assert(m < MAT_ROW_BLK); */
	SAFE_READ_DECL_V(n,uint32_t,fh);
	SAFE_READ_DECL_V(mod,elemt_s,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,uint64_t,fh);

	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %.3f%%) mod %lu\n",m,n,nnz,(double)(nnz)/((double)m*(double)n)*100.,(int64_t)mod);

	A_init->col = C_init->col = n ;
	A_init->mod = C_init->mod = mod ;

	/* READ in ZO start */
	SAFE_READ_DECL_P(rows,m,uint32_t,fh);

	/* pol/zo correspondance */
	SAFE_READ_DECL_P(map_zo_pol,m,uint32_t,fh);


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,uint64_t,fh);
	SAFE_READ_DECL_P(buffer,colid_size,uint32_t,fh); /* buffer has the matrix */

	/* size of nnz for pols: */
	SAFE_READ_V(polys->nb,uint32_t,fh);

	SAFE_READ_DECL_V(pol_nnz,uint64_t,fh);

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_DECL_P(pol_rows,polys->nb,uint32_t,fh);

	/* XXX what if elemt_s == elemt_t ??? */
	SAFE_READ_DECL_P(polys_data_pol,pol_nnz,elemt_s,fh);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> end reading      : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));


	gettimeofday(&tic,NULL);


	SAFE_MALLOC_DECL(start,m+1,uint64_t);
	start[0] = 0 ;
	dimen_t i ;
	for ( i = 0 ; i < m ; ++i)
		start[i+1] = start[i] + rows[i];
	free(rows);

	SAFE_MALLOC(polys->start_pol,polys->nb+1,uint32_t);
	polys->start_pol[0] = 0 ;
	for ( i = 0 ; i < polys->nb ; ++i)
		polys->start_pol[i+1] = polys->start_pol[i]+pol_rows[i];
	free(pol_rows);

	assert(polys->start_pol[polys->nb] == pol_nnz);
	SAFE_MEMCPY_CVT(polys->data_pol,elemt_t,polys_data_pol,polys->start_pol[polys->nb]);
	free(polys_data_pol);


	SAFE_MALLOC_DECL(colid,nnz,dimen_t); /* colid expands buffer */


	expandColid(buffer,colid_size,colid
#ifndef NDEBUG
			,nnz,n
#endif
		   );
	free(buffer);

#ifdef ONLY_FFLAS
	D->row = m ;
	D->col = n ;
	D->ld  = ALIGN(n) ;
	D->nnz = nnz ;
	D->mod = mod ;
	SAFE_CALLOC(D->ptr,(index_t)D->row*(index_t)D->ld,elemt_t);


	for (i = 0 ; i < m  ; ++i) {
		index_t jz ;
		elemt_t * e = polys->data_pol + polys->start_pol[map_zo_pol[i]] ;
		for ( jz = start[i] ; jz < start[i+1] ; ++jz) {
			dimen_t j = colid[jz] ;
			D->ptr[i*D->ld+j] = *e  ;
			++e;
		}
	}

	gettimeofday(&tac,NULL);
	fprintf(stderr,"  >> fflas            : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));


	return NULL ;

#endif


	SAFE_CALLOC_DECL(pivots_data,n,dimen_t);

	/* dimen_t k_dim = getSparsestRows(colid,start,m,n, pivots_data); */
	dimen_t k_dim = getSparsestRows_fast(colid,start,m,n, pivots_data);
	assert(k_dim);
	SAFE_MALLOC_DECL(pivots,k_dim,dimen_t);
	SAFE_MALLOC_DECL(nonpiv,n-k_dim,dimen_t);

	dimen_t pivots_size, nonpiv_size ;

	/* dimen_t k_split = */ createPivots(pivots,&pivots_size,nonpiv,&nonpiv_size,k_dim,pivots_data,n);

	free(pivots_data);

	assert(k_dim >= pivots_size);
	assert(n-k_dim >= nonpiv_size);

	SAFE_REALLOC(pivots,pivots_size,dimen_t);
	SAFE_REALLOC(nonpiv,nonpiv_size,dimen_t);

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




	gettimeofday(&tic,NULL);

	splitVertical(A_init,C_init,nonpiv,nonpiv_size,polys,A,B,C,D);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> split vertical   : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));


#ifdef STATS
	{

		index_t acc ;

		acc = occupancySparse(A);
		fprintf(stderr,"A occupancy by %u : %lu - %lu : %.3f (gain %.3f)\n",UNRL,acc,A->nnz, (double)(acc)/(double)A->nnz,(double)(acc-A->nnz)/(double)A->nnz);

		acc = occupancySparse(B);
		fprintf(stderr,"B occupancy by %u : %lu - %lu : %.3f (gain %.3f)\n",UNRL,acc,B->nnz, (double)(acc)/(double)B->nnz,(double)(acc-B->nnz)/(double)B->nnz);

		acc = occupancySparse(C);
		fprintf(stderr,"C occupancy by %u : %lu - %lu : %.3f (gain %.3f)\n",UNRL,acc,C->nnz, (double)(acc)/(double)C->nnz,(double)(acc-C->nnz)/(double)C->nnz);


		acc = occupancyDense(D);
		fprintf(stderr,"D occupancy by %u : %lu - %lu %.3f (gain %.3f)\n",UNRL,acc,D->nnz, (double)(acc)/(double)D->nnz, (double)(acc-D->nnz)/(double)D->nnz);
	}
#endif

	freeMat(A_init);
	free(A_init);
	freeMat(C_init);
	free(C_init);
	freePol(polys);
	free(polys);

	fprintf(stderr,"   -- sparsity of A   : %.3f%% (%u x %u - %lu)\n",(double)A->nnz/(double)A->row/(double)A->col*100.,A->row,A->col,A->nnz);
	fprintf(stderr,"   -- sparsity of B   : %.3f%% (%u x %u - %lu)\n",(double)B->nnz/(double)B->row/(double)B->col*100.,B->row,B->col,B->nnz);
	fprintf(stderr,"   -- sparsity of C   : %.3f%% (%u x %u - %lu)\n",(double)C->nnz/(double)C->row/(double)C->col*100.,C->row,C->col,C->nnz);
	fprintf(stderr,"   -- sparsity of D   : %.3f%% (%u x %u - %lu)\n",(double)D->nnz/(double)D->row/(double)D->col*100.,D->row,D->col,D->nnz);



	return nonpiv ;
}


dimen_t RowReduce_int32_t ( int32_t p, int32_t * A, dimen_t m, dimen_t n, dimen_t lda) ;
dimen_t RowReduce_double  ( double p, double * A, dimen_t m, dimen_t n, dimen_t lda, uint32_t nt) ;
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
		CSR * Cd = &(C->sub[k]);
		dimen_t ldd = D->ld ;
		/* D = D - C . B */
		assert( (elemt_t)-1<1); /* unsigned will fail */
#ifdef _OPENMP
#pragma omp parallel for /* schedule (static,8) */
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
	assert((int64_t)nb - (int64_t)jz < UNRL);
	for (  ; jz < nb ; ++jz) {
		B[colid[jz]] += tc * A[jz] ;
	}
#endif
}

#ifdef BLOCK_CSR
void spaxpy_block(
		elemt_t         tc,
		const elemt_t * A,
		const dimen_t   nb,
		const dimen_t * colid,
		elemt_t       * B)
{
	assert(!((uintptr_t)A % 32u));
	assert(!((uintptr_t)B % 32u));

#ifdef SIMD
	SET1_SIMD(dc,tc);
#endif
	dimen_t jz  ;
	for ( jz = 0 ; jz < nb ; ++jz) {
		dimen_t col = colid[jz];
#ifdef SIMD
		AXPY_SIMD(B+col,dc,A+jz*UNRL);
#else
		B[col  ] += tc * A[jz*UNRL  ] ;
		B[col+1] += tc * A[jz*UNRL+1] ;
#if (UNRL>2)
		B[col+2] += tc * A[jz*UNRL+2] ;
		B[col+3] += tc * A[jz*UNRL+3] ;
#endif
#if (UNRL>4)
		B[col+4] += tc * A[jz*UNRL+4] ;
		B[col+5] += tc * A[jz*UNRL+5] ;
#endif
#if (UNRL>6)
		B[col+6] += tc * A[jz*UNRL+6] ;
		B[col+7] += tc * A[jz*UNRL+7] ;
#endif
#endif /* SIMD */

	}
}
#endif

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
	dimen_t jz = 0 ;
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

#ifdef BLOCK_CSR
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
	assert(!((uintptr_t)A % 32u));
	assert(!((uintptr_t)B % 32u));


	dimen_t jz  ;
#ifdef SIMD
	SET1_SIMD(dc,tc);
	SET1_SIMD(dd,td);
#endif
	for ( jz = 0 ; jz < nb ; ++jz) {
		dimen_t col = colid[jz];
#ifdef SIMD
		AXPY_SIMD(B+col,dc,A+jz*UNRL);
#else
		B[col  ]  += tc * A[jz*UNRL  ] ;
		B[col+1]  += tc * A[jz*UNRL+1] ;
#if (UNRL>2)
		B[col+2]  += tc * A[jz*UNRL+2] ;
		B[col+3]  += tc * A[jz*UNRL+3] ;
#endif
#if (UNRL>4)
		B[col+4]  += tc * A[jz*UNRL+4] ;
		B[col+5]  += tc * A[jz*UNRL+5] ;
#endif
#if (UNRL>6)
		B[col+6]  += tc * A[jz*UNRL+6] ;
		B[col+7]  += tc * A[jz*UNRL+7] ;
#endif
#endif /* SIMD */

#ifdef SIMD
		AXPY_SIMD(B+col+ld,dd,A+jz*UNRL);
#else
		B[col+ld  ] += td * A[jz*UNRL  ] ;
		B[col+ld+1] += td * A[jz*UNRL+1] ;
#if (UNRL>2)
		B[col+ld+2] += td * A[jz*UNRL+2] ;
		B[col+ld+3] += td * A[jz*UNRL+3] ;
#endif
#if (UNRL>4)
		B[col+ld+4] += td * A[jz*UNRL+4] ;
		B[col+ld+5] += td * A[jz*UNRL+5] ;
#endif
#if (UNRL>6)
		B[col+ld+6] += td * A[jz*UNRL+6] ;
		B[col+ld+7] += td * A[jz*UNRL+7] ;
#endif
#endif /* SIMD */

	}
}

#endif

void spaxpyn(
		elemt_t         *  diag,
		const elemt_t   * A,
		const dimen_t   nb,
		const dimen_t * colid,
		elemt_t         * B,
		const dimen_t * jump,
		const dimen_t   js,
		const dimen_t   ld
		)
{
	dimen_t jz,ii ;
	for (jz = 0 ; jz < nb ; ++jz) {
		for ( ii = 0 ; ii < js ; ++ii)
			B[colid[jz]+ld*jump[ii]] += diag[ii] * A[jz] ;
	}
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
	index_t jz;
	for ( ii = 0 ; ii < row ;  ++ii) {
		for ( jz = start[ii] ; jz < start[ii+1] ; ++jz) {
			temp_C[ii*ldc+colid[jz]] = data[jz] ;
		}
	}
}

#ifdef BLOCK_CSR
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
	index_t jz;
	for ( ii = 0 ; ii < row ;  ++ii) {
		elemt_t * C_off = temp_C+(index_t)ii*(index_t)ldc;
		for ( jz = start[ii] ; jz < start[ii+1] ; ++jz) {
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
#endif


void reduce_chunk_1(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
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
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);

				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(B->sub[kk]);
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				spaxpy(tc,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);
			}
		}
	}
}

#ifdef BLOCK_CSR
void reduce_chunk_1_block(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p)
{
	assert(A->col);
	dimen_t j, ii ;
	for ( j = 0 ; j < A->col ; ++j) { /* A->col == C->col */
		/* XXX invert loops ? */
		for ( ii = 0 ; ii < blk_i ; ii += 1 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			if (tc != 0.) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
				spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);
#else
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);
#endif


				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(B->sub[kk]);
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				/* assert(sz); */
				spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);
			}
		}
	}
}
#endif

void reduce_chunk_dense_1(
		dimen_t blk_i
		, GBMatrix_t * A
		, DNS * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		, dimen_t * row_beg
	       	)
{
	dimen_t j, ii ;
	for ( j = 0 ; j < A->col ; ++j) {
		for ( ii = 0 ; ii < blk_i ; ii += 1 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			if (tc != (elemt_t)0.) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				assert(kk*MAT_ROW_BLK+jj == j);
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);

				/* temp_D -= temp_C[j] * B[j] */
				dimen_t rs = row_beg[j] ;
				cblas_daxpy(B->col-rs, tc,B->ptr+(index_t)j*(index_t)B->ld+rs,1,temp_D+ii*ldd+rs,1);
			}
		}
	}
}


void reduce_chunk_2(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		)
{
	if (blk_i == 1) {
		reduce_chunk_1(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
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
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
					spaxpy2(tc,td,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc,ldc);

					/* temp_D -= temp_C[j] * B[j] */

					/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
					CSR * B_k = &(B->sub[kk]);
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
					spaxpy2(tc,td,B_k->data+B_k->start[jj],
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd,ldd);
				}
				else {
					/* temp_C -= temp_C[j] * A[j] */
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);

					/* temp_D -= temp_C[j] * B[j] */

					/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
					CSR * B_k = &(B->sub[kk]);
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
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
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpy(td,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+(ii+1)*ldc);

				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(B->sub[kk]);
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
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
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);

				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(B->sub[kk]);
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				spaxpy(tc,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D+ii*ldd);

			}
		}
	}
}

#ifdef BLOCK_CSR
void reduce_chunk_2_block(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		)
{
	if (blk_i == 1) {
		reduce_chunk_1_block(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
		return ;
	}
	dimen_t j, ii ;
	assert(A->col);
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
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
					spaxpy2_block(tc,td,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc,ldc);
#else
					spaxpy2(tc,td,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc,ldc);

#endif

					/* temp_D -= temp_C[j] * B[j] */

					/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
					CSR * B_k = &(B->sub[kk]);
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
					spaxpy2_block(tc,td,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd,ldd);
				}
				else {
					/* temp_C -= temp_C[j] * A[j] */
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
					spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
#else
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
#endif

					/* temp_D -= temp_C[j] * B[j] */

					/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
					CSR * B_k = &(B->sub[kk]);
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
					spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);

				}
			}
			else if (td != 0) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
				spaxpy_block(td,A_k->data+A_k->start[jj]*UNRL,
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+(ii+1)*ldc);
#else
				spaxpy(td,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+(ii+1)*ldc);
#endif

				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(B->sub[kk]);
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				spaxpy_block(td,B_k->data+B_k->start[jj]*UNRL,
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
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
#ifdef CONV_A
					spaxpy_block(tc,A_k->data+A_k->start[jj]*UNRL,
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
#else
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);
#endif

					/* temp_D -= temp_C[j] * B[j] */

					/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
					CSR * B_k = &(B->sub[kk]);
					sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
					spaxpy_block(tc,B_k->data+B_k->start[jj]*UNRL,
							sz,
							B_k->colid+B_k->start[jj],
							temp_D+ii*ldd);

				}
		}
	}
}

#endif

void reduce_chunk_dense_2(
		dimen_t blk_i
		, GBMatrix_t * A
		, DNS * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p
		, dimen_t * row_beg
		)
{
	if (blk_i == 1) {
		reduce_chunk_dense_1(blk_i,A,B,temp_C,ldc,temp_D,ldd,p,row_beg);
		return ;
	}

	elemt_t * Bd = B->ptr ;
	dimen_t ldb = B->ld ;

	dimen_t j, ii ;
	for ( j = 0 ; j < A->col ; ++j) {
		for ( ii = 0 ; ii < blk_i/2*2 ; ii += 2 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			elemt_t td = Mjoin(fmod,elemt_t)(-temp_C[(ii+1)*ldc+j],p) ;
			if (tc != (elemt_t)0.) {
				if (td != (elemt_t)0.) {
					/* temp_C -= temp_C[j] * A[j] */
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					assert(kk*MAT_ROW_BLK+jj == j);
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
					spaxpy2(tc,td,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc,ldc);

					/* temp_D -= temp_C[j] * B[j] */
					dimen_t rs = row_beg[j] ;
					cblas_daxpy(B->col-rs, tc,Bd+(index_t)j*(index_t)ldb+rs,1,temp_D+ii*ldd+rs,1);
					cblas_daxpy(B->col-rs, td,Bd+(index_t)j*(index_t)ldb+rs,1,temp_D+(ii+1)*ldd+rs,1);
				}
				else {
					/* temp_C -= temp_C[j] * A[j] */
					dimen_t jj = j%MAT_ROW_BLK ;
					dimen_t kk = j/MAT_ROW_BLK ;
					assert(kk*MAT_ROW_BLK+jj == j);
					CSR * A_k = &(A->sub[kk]);
					dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
					assert(kk*MAT_ROW_BLK+jj == j);
					spaxpy(tc,A_k->data+A_k->start[jj],
							sz,
							A_k->colid+A_k->start[jj]
							,temp_C+ii*ldc);

					/* temp_D -= temp_C[j] * B[j] */
					dimen_t rs = row_beg[j] ;
					cblas_daxpy(B->col-rs, tc,Bd+(index_t)j*(index_t)ldb+rs,1,temp_D+ii*ldd+rs,1);
				}
			}
			else if (td != (elemt_t)0.) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				assert(kk*MAT_ROW_BLK+jj == j);
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpy(td,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+(ii+1)*ldc);

				/* temp_D -= temp_C[j] * B[j] */
				dimen_t rs = row_beg[j] ;
				cblas_daxpy(B->col-rs, td,Bd+(index_t)j*(index_t)ldb+rs,1,temp_D+(ii+1)*ldd+rs,1);
			}
		}
		for (  ; ii < blk_i ; ii += 1 ) {
			elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
			if (tc != (elemt_t)0.) {
				/* temp_C -= temp_C[j] * A[j] */
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				assert(kk*MAT_ROW_BLK+jj == j);
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpy(tc,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C+ii*ldc);

				/* temp_D -= temp_C[j] * B[j] */
				dimen_t rs = row_beg[j] ;
				cblas_daxpy(B->col-rs, tc,Bd+(index_t)j*(index_t)ldb+rs,1,temp_D+ii*ldd+rs,1);
			}
		}
	}
}

void reduce_chunk_n(
		dimen_t blk_i
		, GBMatrix_t * A
		, GBMatrix_t * B
		, elemt_t * temp_C
		, dimen_t ldc
		, elemt_t * temp_D
		, dimen_t ldd
		, elemt_t p)
{
	dimen_t blk = MAT_SUB_BLK ;
	dimen_t j, ii ;
	SAFE_MALLOC_DECL(jump,blk,dimen_t);
	SAFE_MALLOC_DECL(diag,blk,elemt_t);
			for ( j = 0 ; j < A->col ; ++j) {
				dimen_t nbnz = 0;
				for ( ii = 0 ; ii < blk_i ; ++ii) {
					elemt_t tc = Mjoin(fmod,elemt_t)(-temp_C[ii*ldc+j],p) ;
					if (tc != 0.) {
						diag[nbnz]   = tc ;
						jump[nbnz++] = ii ;
					}
				}
				dimen_t jj = j%MAT_ROW_BLK ;
				dimen_t kk = j/MAT_ROW_BLK ;
				CSR * A_k = &(A->sub[kk]);
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpyn(diag,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C,jump,nbnz,ldc);

				/* temp_D -= temp_C[j] * B[j] */

				/* cblas_daxpy(D->col, tc,Bd+j*ldb,1,temp_D,1); */
				CSR * B_k = &(B->sub[kk]);
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				spaxpyn(diag,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D,jump,nbnz,ldd);
			}
}


void reduce_fast(
		GBMatrix_t      * A
		, GBMatrix_t    * B
		, GBMatrix_t    * C
		, DNS * D )
{
	dimen_t ldd = D->ld ;
	dimen_t ldc = C->col ;


	elemt_t p = A->mod ;
	/* elemt_t * Bd = B->ptr; */
	elemt_t * Dd = D->ptr;
	dimen_t i ;


	dimen_t blk = MAT_SUB_BLK ;

#ifndef _OPENMP
	SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk,elemt_t);
	SAFE_MALLOC_DECL(temp_C,(index_t)ldc*(index_t)blk,elemt_t);
#endif

	dimen_t k ;
	assert(D->col == B->col);


	for (k = 0 ; k < C->sub_nb ; ++k) {
		CSR * C_k = &(C->sub[k]) ;
#ifdef _OPENMP
#pragma omp parallel for
#endif

		for ( i = 0 ; i < C_k->row ; i += blk ) {
			dimen_t blk_i = min(blk,C_k->row - i);
#ifdef _OPENMP
			SAFE_MALLOC_DECL(temp_D,((index_t)ldd*(index_t)blk_i),elemt_t);
			SAFE_CALLOC_DECL(temp_C,((index_t)ldc*(index_t)blk_i),elemt_t);
#endif

			index_t i_offset = k*MAT_ROW_BLK + i ;
#ifndef _OPENMP
			/* cblas_dscal(ldc*blk_i,0.,temp_C,1); |+ XXX blk * col < 2^32 +| */
			memset(temp_C,0,(index_t)ldc*(index_t)blk_i*sizeof(elemt_t));
#endif
			sparse_dcopy( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);

			cblas_dcopy(ldd*blk_i,Dd+i_offset*(index_t)ldd,1,temp_D,1);


#ifdef USE_SAXPY
#if defined(USE_SAXPYn) || defined(USE_SAXPY2)
#error "make a choice"
#endif
			reduce_chunk_1(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPY */

#ifdef USE_SAXPY2
#if defined(USE_SAXPY) || defined(USE_SAXPYn)
#error "make a choice"
#endif
			reduce_chunk_2(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPY2 */

#ifdef USE_SAXPYn
#if defined(USE_SAXPY) || defined(USE_SAXPY2)
#error "make a choice"
#endif
			reduce_chunk_n(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPYn */

			/* Mjoin(Finit,elemt_t)(p, temp_D, Dd+i*ldd, D->col) ; */
			cblas_dcopy(ldd*blk_i,temp_D,1,Dd+i_offset*(index_t)ldd,1);

#ifdef _OPENMP
		free(temp_D);
		free(temp_C);
#endif /* _OPENMP */
		}

	}

#ifndef _OPENMP
	free(temp_D);
	free(temp_C);
#endif /* _OPENMP */
	Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;
}

#ifdef BLOCK_CSR
void reduce_fast_block(
		GBMatrix_t      * A
		, GBMatrix_t    * B
		, GBMatrix_t    * C
		, DNS * D )
{

	struct timeval tic,tac ;
	gettimeofday(&tic,NULL);

	SAFE_MALLOC_DECL(BH,1,GBMatrix_t);
	initSparse(BH);
	convert_CSR_2_CSR_block(BH,B);

#ifdef CONV_A
	SAFE_MALLOC_DECL(AH,1,GBMatrix_t);
	initSparse(AH);
	convert_CSR_2_CSR_block(AH,A);
#else
	GBMatrix_t * AH = A ;
#endif

#ifdef CONV_C
	SAFE_MALLOC_DECL(CH,1,GBMatrix_t);
	initSparse(CH);
	convert_CSR_2_CSR_block(CH,C);
#else
	GBMatrix_t * CH = C ;
#endif

	gettimeofday(&tac,NULL);
	fprintf(stderr,"   -- convert time    : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

	dimen_t ldd = D->ld ;
	dimen_t ldc = ALIGN(C->col) ;


	elemt_t p = AH->mod ;
	elemt_t * Dd = D->ptr;
	dimen_t i ;


	dimen_t blk = MAT_SUB_BLK ;

#ifndef _OPENMP
	SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk,elemt_t);
	SAFE_MALLOC_DECL(temp_C,(index_t)ldc*(index_t)blk,elemt_t);
#endif

	dimen_t k ;
	assert(D->col == B->col);


	for (k = 0 ; k < C->sub_nb ; ++k) {
		CSR * C_k = &(CH->sub[k]) ;
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for ( i = 0 ; i < C_k->row ; i += blk ) {
			dimen_t blk_i = min(blk,C_k->row - i);
#ifdef _OPENMP
			SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk_i,elemt_t);
			SAFE_CALLOC_DECL(temp_C,((index_t)ldc*(index_t)blk_i),elemt_t);
#endif

			index_t i_offset = k*MAT_ROW_BLK + i ;
#ifndef _OPENMP
			cblas_dscal(ldc*blk_i,0.,temp_C,1); /* XXX blk * col < 2^32 */
			/* memset(temp_C,0,(index_t)ldc*(index_t)blk_i*sizeof(elemt_t)); */
#endif
#ifdef CONV_C
			sparse_dcopy_block( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);
#else
			sparse_dcopy( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);
#endif

			cblas_dcopy(blk_i*ldd,Dd+i_offset*(index_t)ldd,1,temp_D,1);


#ifdef USE_SAXPY
#if defined(USE_SAXPYn) || defined(USE_SAXPY2)
#error "make a choice"
#endif
			reduce_chunk_1_block(blk_i,AH,BH,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPY */

#ifdef USE_SAXPY2
#if defined(USE_SAXPY) || defined(USE_SAXPYn)
#error "make a choice"
#endif
			reduce_chunk_2_block(blk_i,AH,BH,temp_C,ldc,temp_D,ldd,p);
#endif /* USE_SAXPY2 */

			/* #ifdef USE_SAXPYn
#if defined(USE_SAXPY) || defined(USE_SAXPY2)
#error "make a choice"
#endif
reduce_chunk_n_block(blk_i,A,B,temp_C,ldc,temp_D,ldd,p);
#endif |+ USE_SAXPYn +|                                                             */

			cblas_dcopy(blk_i*ldd,temp_D,1,Dd+i_offset*(index_t)ldd,1);

#ifdef _OPENMP
			free(temp_D);
			free(temp_C);
#endif /* _OPENMP */
		}

	}

#ifndef _OPENMP
	free(temp_D);
	free(temp_C);
#endif /* _OPENMP */
	Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;

	freeMat(BH);
	free(BH);
#ifdef CONV_A
	freeMat(AH);
	free(AH);
#endif
#ifdef CONV_C
	freeMat(CH);
	free(CH);
#endif
}
#endif


void reduce_fast_dense(
		GBMatrix_t      * A
		, DNS * B
		, GBMatrix_t    * C
		, DNS * D )
{
	dimen_t ldd = D->ld ;
	dimen_t ldb = B->ld ;
	dimen_t ldc = ALIGN(C->col) ;

	SAFE_CALLOC_DECL(row_beg,B->row,dimen_t);
	dimen_t i ;
	for ( i = 0 ; i < B->row ; ++i) {
		dimen_t ii = 0 ;
		while ( *(B->ptr+i*ldb+ii) == 0) {
			row_beg[i] += 1 ;
			++ii ;
		}
	}


	elemt_t p = A->mod ;
	elemt_t * Dd = D->ptr;


	dimen_t blk = MAT_SUB_BLK ;

#ifndef _OPENMP
	SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk,elemt_t);
	SAFE_MALLOC_DECL(temp_C,(index_t)ldc*(index_t)blk,elemt_t);
#endif

	dimen_t k ;
	assert(D->col == B->col);


	for (k = 0 ; k < C->sub_nb ; ++k) {
		CSR * C_k = &(C->sub[k]) ;
#ifdef _OPENMP
#pragma omp parallel for
#endif

		for ( i = 0 ; i < C_k->row ; i += blk ) {
			dimen_t blk_i = min(blk,C_k->row - i);
#ifdef _OPENMP
			SAFE_MALLOC_DECL(temp_D,(index_t)ldd*(index_t)blk_i,elemt_t);
			SAFE_CALLOC_DECL(temp_C,((index_t)ldc*(index_t)blk_i),elemt_t);
#endif

			index_t i_offset = k*MAT_ROW_BLK + i ;
#ifndef _OPENMP
			/* cblas_dscal(ldc*blk_i,0.,temp_C,1); |+ XXX blk * col < 2^32 +| */
			memset(temp_C,0,(index_t)ldc*(index_t)blk_i*sizeof(elemt_t));
#endif
			sparse_dcopy( blk_i, C_k->start+i, C_k->colid, C_k->data    , temp_C, ldc);

			cblas_dcopy(ldd*blk_i,Dd+i_offset*(index_t)ldd,1,temp_D,1);


#ifdef USE_SAXPY
#if defined(USE_SAXPY2) || defined(USE_SAXPYn)
#error "make a choice"
#endif
			reduce_chunk_dense_1(blk_i,A,B,temp_C,ldc,temp_D,ldd,p,row_beg);
#endif /* USE_SAXPY */

#ifdef USE_SAXPY2
#if defined(USE_SAXPY) || defined(USE_SAXPYn)
#error "make a choice"
#endif
			reduce_chunk_dense_2(blk_i,A,B,temp_C,ldc,temp_D,ldd,p,row_beg);
#endif /* USE_SAXPY2 */

			/* Mjoin(Finit,elemt_t)(p, temp_D, Dd+i*ldd, D->col) ; */
			cblas_dcopy(D->col*blk_i,temp_D,1,Dd+i_offset*(index_t)ldd,1);

#ifdef _OPENMP
			free(temp_D);
			free(temp_C);
#endif /* _OPENMP */
		}

	}
	free(row_beg);

#ifndef _OPENMP
	free(temp_D);
	free(temp_C);
#endif /* _OPENMP */
	Mjoin(Freduce,elemt_t)(p, Dd, (index_t)ldd*(index_t)D->row) ;
}




dimen_t echelonD(
		GBMatrix_t      * A
		, DNS * D
		, uint32_t nt)
{

	dimen_t r = Mjoin(RowReduce,elemt_t)(D->mod,D->ptr,D->row,D->col,D->ld, nt);

	fprintf(stderr,"  -- residual rank    : %u\n",r);

	return r + A->row;
}

#endif /* __GB_io_H */

/* vim: set ft=c: */
