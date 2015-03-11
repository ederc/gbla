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
#include "tools.h"

/**
 * The matrix we have as input is "almost upper triangular", ie below every first
 * element in a row, there may be non zero elements on the next rows below.
 * A non pivot column can be permuted.
 * There is no empty line.
 *
 */



/* buffer:
*/

dimen_t nextcol( dimen_t * colid,
		index_t * start,
		dimen_t i)
{
	assert(start[i+1]>start[i]);
	if (start[i+1]-start[i] == 1) { /* no next element in row */
		return 0 ;
	}
	return colid[start[i]+1];
}

dimen_t getSparsestRows_fast(
		dimen_t   * colid
		, index_t  * start
		, dimen_t   row
		, dimen_t   col
		, dimen_t * pivots_data /* pivots_data[j] is the sparsest row with j as pivot */
		)
{


	dimen_t i = 0;
#ifdef STATS
	dimen_t cnt = 0 ;
	for ( i = 0 ; i< row ; ++i) {
		dimen_t j =  nextcol(colid,start,i) ;
		if ( j == 0 ||  (j-i) >= (4))
			cnt ++ ;
	}
	printf("%u %u\n",cnt,row);
#endif


	dimen_t k_dim = 0;

	SAFE_MALLOC_DECL(creux_v,row,dimen_t);
	for ( i = 0 ; i < row ; ++i)
		creux_v[i] = (dimen_t)(start[i+1]-start[i]);

	for ( i = 0 ; i < col ; ++i ) {
		pivots_data[i] = MARKER32 ;
	}

	dimen_t new_i ;
	for ( new_i = 0 ; new_i < row ; ++new_i ) {
		dimen_t pivot_j = colid[start[new_i]] ;     /* first column row new_i */
		dimen_t creux   = creux_v[new_i] ;          /* length of row new_i */
		assert(pivot_j < col);
		dimen_t old_i = pivots_data[pivot_j] ;      /* last row for pivot column */
		if (old_i == MARKER32) {
			pivots_data[pivot_j] = new_i ;
			++k_dim;
		}
		else  {
			dimen_t old_creux = creux_v[old_i];
			if (old_creux == creux) {           /* favour zeros after initial 1 */
				dimen_t old_j_next = nextcol(colid,start,old_i);
				if (old_j_next > 0) {
					dimen_t new_j_next = nextcol(colid,start,new_i);
					if ( (new_j_next == 0) || (new_j_next > old_j_next))
						pivots_data[pivot_j] = new_i ;
				}
			}
			if (old_creux > creux) {            /* this row is sparser */
				pivots_data[pivot_j] = new_i ;
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


	assert(A->row == pivots_size );
	assert(C->row == row - pivots_size );
	assert(C->row); /* return here and terminate */

	A->sub_row = DIVIDE_INTO(A->row,MAT_ROW_BLK);
	A->sub_col = 1 ;
	assert(A->sub_row);
	C->sub_row = DIVIDE_INTO(C->row,MAT_ROW_BLK);
	C->sub_col = 1 ;
	assert(C->sub_row);


	SAFE_CALLOC_DECL(Asub_nnz,A->sub_row,index_t);
	SAFE_CALLOC_DECL(Csub_nnz,C->sub_row,index_t);
	SAFE_CALLOC_DECL(ami,row,uint8_t);
	SAFE_CALLOC_DECL(qui,row,dimen_t);

	dimen_t i_c = 0 , i_a = 0 ;
	for ( i = 0 ; i < row ; ++i)  {
		assert(start[i] < nnz);
		if ( (last_pivot < pivots_size) && (pivots[last_pivot] == i) ) {
			dimen_t i_loc = i_a / MAT_ROW_BLK ;
			Asub_nnz[i_loc] += (start[i+1] - start[i]) ;
			++last_pivot;
			ami[i] = 1 ;
			qui[i] = i_a  ;
			i_a ++ ;
		}
		else {
			dimen_t i_loc = i_c / MAT_ROW_BLK ;
			Csub_nnz[i_loc] += (start[i+1] - start[i]) ;
			qui[i] = i_c   ;
			i_c ++ ;
		}
	}

	SAFE_MALLOC(A->sub,A->sub_row,CSR);
	SAFE_MALLOC(C->sub,C->sub_row,CSR);

	index_t A_nnz = 0 ;

	for ( i = 0 ; i < A->sub_row ; ++i) {
		CSR * A_sub = A->sub+i ;
		A_sub->row = min((dimen_t)MAT_ROW_BLK,A->row-i*MAT_ROW_BLK);
		A_sub->col = A->col ;
		A_sub->nnz = Asub_nnz[i] ;
		A_nnz += Asub_nnz[i] ;
		SAFE_CALLOC(A_sub->start,(A_sub->row+1),index_t);
		SAFE_MALLOC(A_sub->colid,A_sub->nnz,dimen_t);
		SAFE_MALLOC(A_sub->map_zo_pol,A_sub->row,dimen_t);
		A_sub->data = NULL ;
	}

	for ( i = 0 ; i < C->sub_row; ++i) {
		CSR * C_sub = C->sub+i ;
		C_sub->row = min((dimen_t)MAT_ROW_BLK,C->row-i*MAT_ROW_BLK);
		C_sub->col = C->col ;
		/* C_sub->mod = C->mod ; */
		C_sub->nnz = Csub_nnz[i] ;
		SAFE_CALLOC(C_sub->start,(C_sub->row+1),index_t);
		SAFE_MALLOC(C_sub->colid,C_sub->nnz,dimen_t);
		SAFE_MALLOC(C_sub->map_zo_pol,C_sub->row,dimen_t);
		C_sub->data = NULL ;
	}

	free (Asub_nnz);
	free (Csub_nnz);

	A->nnz = A_nnz ;
	C->nnz = start[row]  - A_nnz ;

#ifdef _OPENMP
/* #pragma omp parallel for */
#endif
	for ( i = 0 ; i < row ; ++i)  {
		assert(start[i] < nnz);
		dimen_t s_loc = qui[i] / MAT_ROW_BLK ;
		dimen_t i_loc = qui[i] % MAT_ROW_BLK ;
		if ( ami[i] == 1 ) {
			assert(s_loc < A->sub_row);
			CSR * A_sub = A->sub+s_loc ;
			assert(qui[i_loc] < A_sub->row);
			setRow(A_sub,i_loc,colid+start[i],start[i+1]-start[i],map_zo_pol[i]);
		}
		else {
			CSR * C_sub = C->sub+s_loc ;
			setRow(C_sub,i_loc,colid+start[i],start[i+1]-start[i],map_zo_pol[i]);
		}


	}

	free(ami);
	free(qui);



#ifndef NDEBUG
	assert(A->nnz + C->nnz == nnz);
#endif
	return ;
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

	A->sub_row = A_init->sub_row ;
	A->sub_col = A_init->sub_col ;
	SAFE_MALLOC(A->sub,A->sub_row,CSR);

	/* Init sparse */

	dimen_t j ;
	for ( j = 0 ; j < A_init->sub_row; ++j) {
		const CSR * A_k = A_init->sub + j ;
#ifndef NDEBUG
		checkMatUnit(A_k);
#endif
		CSR * Ad = A->sub + j;
		Ad->row = A_k->row ;
		Ad->col = A->col ;
		Ad->nnz = 0 ;
		/* A->row += Ad->row ; */
		assert(A->col == Ad->col);

		SAFE_CALLOC(Ad->start,(Ad->row+1),index_t);
		SAFE_MALLOC(Ad->colid,A_k->nnz,dimen_t);
		SAFE_MALLOC(Ad->data ,A_k->nnz,elemt_t);
		Ad->map_zo_pol=NULL;


		SAFE_CALLOC_DECL(b_there,(A_k->row+1),index_t);
		dimen_t i ;

		for (i = 0 ; i < A_k->row ; ++i) {
			dimen_t here  = 0; /* shift */
			index_t jz ;
			b_there[i+1] += b_there[i] ;
			for (jz = A_k->start[i] ; jz < A_k->start[i+1] ; ++jz) {
				dimen_t k = A_k->colid[jz]  ;
				assert(k < A_k->col);
				if ( k  < max_col ) {
					while (here < nonpiv_size && k > nonpiv[here]) {
						here ++ ;
					}
					if (here < nonpiv_size && k == nonpiv[here]) {
						B->nnz += 1 ;
					}
					else {
						b_there[i+1] += 1 ;
					}
				}
				else {
					assert(nonpiv_size+k-max_col < B->col);
					B->nnz += 1 ;
				}
			}
			assert(here <= nonpiv_size);
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0 ; i < A_k->row ; ++i) {
			dimen_t here  = 0; /* shift */
			index_t there = b_there[i];
			index_t i_offset =  j * MAT_ROW_BLK + i ;
			index_t start_p = polys->start_pol[ A_k->map_zo_pol[i] ] ;
			elemt_t * d = polys->data_pol+start_p ;
			index_t jz ;
			for (jz = A_k->start[i] ; jz < A_k->start[i+1] ; ++jz) {
				dimen_t k = A_k->colid[jz]  ;
				assert(k < A_k->col);
				if ( k  < max_col ) {
					while (here < nonpiv_size && k > nonpiv[here]) {
						here ++ ;
					}
					if (here < nonpiv_size && k == nonpiv[here]) {
						B->ptr[i_offset*(index_t)ldb+(index_t)here] = *d; /* XXX */
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
				}
				++d;
			}
			assert(here <= nonpiv_size);
		}

		free(b_there);

		for ( i = 1 ; i < Ad->row ; ++i) {
			Ad->start[i+1] += Ad->start[i] ;
		}
		Ad->nnz = Ad->start[Ad->row] ;
		A->nnz += Ad->nnz ;
		assert(A->nnz);
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

	A->sub_row = B->sub_row = A_init->sub_row ;
	A->sub_col = B->sub_col = A_init->sub_col ;
	assert(A->sub_col);

	SAFE_MALLOC(A->sub,A->sub_row,CSR);
	SAFE_MALLOC(B->sub,B->sub_row,CSR);
	dimen_t j ;
	for ( j = 0 ; j < A_init->sub_row ; ++j) {
		const CSR * A_k = A_init->sub + j ;
#ifndef NDEBUG
		checkMatUnit(A_k);
#endif
		CSR * Ad = A->sub+j;
		CSR * Bd = B->sub+j;
		Ad->row = A_k->row ;
		Ad->col = A->col ;
		Ad->nnz = 0 ;
		Bd->row = A_k->row ;
		Bd->col = B->col ;
		Bd->nnz = 0 ;
		assert(A->col == Ad->col);
		assert(B->col == Bd->col);

		SAFE_CALLOC(Ad->start,(Ad->row+1),index_t);
		SAFE_CALLOC(Bd->start,(Bd->row+1),index_t);
		Ad->map_zo_pol=NULL;
		Bd->map_zo_pol=NULL;


		/* index_t there = 0; */
		dimen_t i ;
		index_t jz ;

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

		/* dimen_t la = 0 ; */
		/* there = 0; */
#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0 ; i < A_k->row ; ++i) {
			dimen_t here  = 0; /* shift */
			index_t there = b_there[i] ;
			index_t ici =  b_la[i] ;
			index_t la  =  ici + b_row_off[i] ;
			index_t start_p = polys->start_pol[ A_k->map_zo_pol[i] ] ;
			elemt_t * d = polys->data_pol+start_p ;
			index_t jz ;
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
		, struct timeval * toc
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
	SAFE_READ_DECL_V(b,dimen_t,fh);
	/* XXX set elemt_s here and C++-ise*/
	assert((b ^ VERMASK) == Mjoin(select,elemt_s)());

	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,dimen_t,fh);
	/* assert(m < MAT_ROW_BLK); */
	SAFE_READ_DECL_V(n,dimen_t,fh);
	SAFE_READ_DECL_V(mod,elemt_m,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,index_t,fh);

	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %5.2f%%) mod %lu\n",m,n,nnz,(double)(nnz)/((double)m*(double)n)*100.,(int64_t)mod);

	A_init->col = C_init->col = n ;
	A_init->mod = C_init->mod = mod ;

	/* READ in ZO start */
	SAFE_READ_DECL_P(rows,m,dimen_t,fh);

	/* pol/zo correspondance */
	SAFE_READ_DECL_P(map_zo_pol,m,dimen_t,fh);


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,index_t,fh);
	SAFE_READ_DECL_P(buffer,colid_size,dimen_t,fh); /* buffer has the matrix */

	/* size of nnz for pols: */
	SAFE_READ_V(polys->nb,dimen_t,fh);

	SAFE_READ_DECL_V(pol_nnz,index_t,fh);

	if (pol_nnz > UINT32_MAX) {
		fprintf(stderr," ** warning, there are really many polynomials in this matrix **\n");
	}

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_DECL_P(pol_rows,polys->nb,dimen_t,fh);

	/* XXX what if elemt_s == elemt_t ??? */
	SAFE_READ_DECL_P(polys_data_pol,pol_nnz,elemt_s,fh);

	gettimeofday(&tac,NULL);

	toc->tv_sec  = tac.tv_sec - tic.tv_sec ;
	toc->tv_usec = tac.tv_usec- tic.tv_usec ;


	fprintf(stderr,"  >> end reading      : %.3f s\n", ((double)(toc->tv_sec)
				+(double)(toc->tv_usec)/1e6));


	gettimeofday(&tic,NULL);


	SAFE_MALLOC_DECL(start,(m+1),index_t);
	start[0] = 0 ;
	dimen_t i ;
	for ( i = 0 ; i < m ; ++i)
		start[i+1] = start[i] + rows[i];
	free(rows);

	SAFE_MALLOC(polys->start_pol,(polys->nb+1),index_t);
	polys->start_pol[0] = 0 ;
	for ( i = 0 ; i < polys->nb ; ++i)
		polys->start_pol[i+1] = polys->start_pol[i]+pol_rows[i];
	free(pol_rows);

	assert(polys->start_pol[polys->nb] == pol_nnz);
	SAFE_MEMCPY_CVT(polys->data_pol,elemt_t,polys_data_pol,pol_nnz);
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

	A_init->row = k_dim;
	C_init->row = m - k_dim;


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

	assert(nnz == A_init->nnz + C_init->nnz);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> split horizontal : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));




	gettimeofday(&tic,NULL);

	SAFE_MALLOC_DECL(AH,1,GBMatrix_t);
	SAFE_MALLOC_DECL(BH,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CH,1,GBMatrix_t);

	initSparse(AH);
	initSparse(BH);
	initSparse(CH);

	splitVertical(A_init,C_init,nonpiv,nonpiv_size,polys,AH,BH,CH,D);

	freeMat(A_init);
	free(A_init);
	freeMat(C_init);
	free(C_init);
	freePol(polys);
	free(polys);

	gettimeofday(&tac,NULL);

	fprintf(stderr,"  >> split vertical   : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

	gettimeofday(&tic,NULL);

	createSubmatrices(A,AH);

	freeMat(AH);
	free(AH);

	createSubmatrices(B,BH);
	freeMat(BH);
	free(BH);

	createSubmatrices(C,CH);
	freeMat(CH);
	free(CH);


	gettimeofday(&tac,NULL);
	fprintf(stderr,"  >> submatrices      : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));
#ifdef STATS
	dimen_t cnt = 0 ;
	assert(A->sub_col == 1);
	for ( i = 0 ; i< A->sub_row ; ++i) {
		CSR * Ai = A->sub + i ;
		dimen_t i_off = i*MAT_ROW_BLK ;
		dimen_t ii = 0 ;
		for ( ii = 0 ; ii < Ai->row ; ++ii) {
			dimen_t j =  nextcol(Ai->colid,Ai->start,ii) ;
			if ( j == 0 ||  (j-ii-i_off*MAT_ROW_BLK) >= (4))
			cnt ++ ;
		}
	}
	printf("%u %u\n",cnt,A->row);
#endif



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



	fprintf(stderr,"   -- sparsity of A   : %5.2f%% (%8u x %8u - %12lu)\n",(double)A->nnz/(double)A->row/(double)A->col*100.,A->row,A->col,A->nnz);
	fprintf(stderr,"   -- sparsity of B   : %5.2f%% (%8u x %8u - %12lu)\n",(double)B->nnz/(double)B->row/(double)B->col*100.,B->row,B->col,B->nnz);
	fprintf(stderr,"   -- sparsity of C   : %5.2f%% (%8u x %8u - %12lu)\n",(double)C->nnz/(double)C->row/(double)C->col*100.,C->row,C->col,C->nnz);
	fprintf(stderr,"   -- sparsity of D   : %5.2f%% (%8u x %8u - %12lu)\n",(double)D->nnz/(double)D->row/(double)D->col*100.,D->row,D->col,D->nnz);



	return nonpiv ;
}



#endif /* __GB_io_H */

/* vim: set ft=c: */
