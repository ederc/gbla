#ifndef __GB_io_H
#define __GB_io_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include "matrix.h"

#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })

#define ERROR_READING \
	if (verbose > 0) { \
		printf("Error while reading file '%s'\n",fn); \
		fclose(fh); \
		return -1 ; \
	}

#define SWAP(a,b)  \
	t = a ; \
a = b ; \
b = t

#if 0
void insert_sort(uint32_t * liste, uint32_t  size)
{
	uint32_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && liste[d] < liste[d-1]) {
			SWAP(liste[d],liste[d-1]);
			d--;
		}
	}
}
#endif


void insert_sort_duo(uint32_t * liste, uint32_t  size, uint32_t * copain)
{
	uint32_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && (liste)[d] < (liste)[d-1]) {
			SWAP((liste)[d],(liste)[d-1]);
			SWAP((copain)[d],(copain)[d-1]);
			d--;
		}
	}
}

void insert_sort_duo_data(uint32_t * liste, uint32_t  size, TYPE * copain)
{
	uint32_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && liste[d] < liste[d-1]) {
			fprintf(stderr,"permuting %u and %u, that is %u and %u, and side %u and %u\n",
					d,d-1,liste[d],liste[d-1],copain[d],copain[d-1]);
			SWAP(liste[d],liste[d-1]);
			SWAP(copain[d],copain[d-1]);
			d--;
		}
	}
}

#undef SWAP


/* genericity */
#define copy(to,from,size) \
{ \
	uint32_t iii = 0 ; \
	for ( ; iii < size ; ++iii) {  \
		(to)[iii] = (from)[iii] ; \
	} \
}

/* buffer:
*/
uint32_t getSparsestRows(uint32_t * colid
		, uint64_t * start
		, uint32_t row
		, uint32_t * pivots /* pivots[j] is the sparsest row with j as pivot */
		)
{
	SAFE_MALLOC_DECL(sparse,row,uint32_t);
	uint32_t j = 0;
	uint32_t pivot_nb = 0;
	for ( ; j < row ; ++j) {
		sparse[j] = (uint32_t)(-1);
	}
	j = 0 ;
	for ( ; j < row ; ++j) {
		uint32_t pivot = colid[start[j]] ; /* first column in each row */
		uint32_t creux = start[j+1]-start[j] ; /* length of the row */
		if (sparse[pivot] > creux) { /* this row is sparsest */
			sparse[pivot] = creux ;
			pivots[pivot]  = j ;
			++ pivot_nb ;
		}
	}

	return pivot_nb ;
}

void splitRowsUpBottom(
		GBMatrix_t * A
		, GBMatrix_t * C
		, uint32_t  * colid
		, uint64_t * start
		, uint32_t row
		, uint32_t * pivots
		, uint32_t pivots_size
		, uint32_t * pols
		)
{
	uint32_t last_pivot=0;
	uint32_t i = 0 ;
	for ( ; i < row ; ++i)  {
		if ( (last_pivot < pivots_size) && (pivots[last_pivot] == i) ) {
			appendRow(A,colid+start[i],start[i+1]-start[i],pols[i]);
			++last_pivot;
		}
		else {
			appendRow(C,colid+start[i],start[i+1]-start[i],pols[i]);
		}
	}
	return ;
}


void expandColid( const uint32_t * compress, uint32_t size_compressed
		, uint32_t * expand  , uint32_t size_expand)
{
	uint32_t mask = (1<<31);
	uint32_t i = 0,j=0 ;
	uint32_t col ;
	/* for (i = 0 ; i < size_compressed ; ++i) { */
		/* fprintf(stderr,"%u ",compress[i]); */
	/* } */
		/* fprintf(stderr,"\n"); i = 0 ; */
	for ( ; i < size_compressed ;) {
		col = compress[i++] ;
		/* fprintf(stderr,"in col : %u\n",col); */
		if (col & mask) {
			expand[j++] = col ^ mask ;
		/* fprintf(stderr,"sing col : %u\n",col^mask); */
		}
		else {
			uint32_t k = 0 ;
			for (; k < compress[i] ;++k) {
				expand[j++] = col + k;
				/* fprintf(stderr,"cons col : %u\n",col+k); */
			}
			++i;
		}
	}
	assert(j == size_expand);
	assert(i == size_compressed);
}

/* matrix reader and row splitter
 * A_init [out] the top part with upper triangular left
 * B_init [out] the bottom part
 * polys [out] the polynomials used in the matrices (shared)
 * fh [in] the file in appropriate format
 */
int read_file(GBMatrix_t * A_init
		, GBMatrix_t * C_init
		, CSR_pol * polys
		, FILE * fh
	     )
{
	uint32_t i = 0 ;

	/* sanity */
	assert(fh);

	/* format */
	SAFE_READ_DECL_V(b,uint32_t,fh);
	/* XXX set TYPE here and C++-ise*/

	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,uint32_t,fh);
	SAFE_READ_DECL_V(n,uint32_t,fh);
	SAFE_READ_DECL_V(mod,TYPE,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,uint64_t,fh);

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

	/* SAFE_MALLOC_DECL(map_pol_zo_A,m,uint32_t); */
	/* uint32_t map_pol_zo_A_size = 0 ; */

	/* SAFE_MALLOC_DECL(map_pol_zo_B,m,uint32_t); */
	/* uint32_t map_pol_zo_B_size = 0 ; */


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,uint64_t,fh);
	/* fprintf(stderr,"colid_size %u\n",colid_size); */
	SAFE_READ_DECL_P(buffer,colid_size,uint32_t,fh); /* buffer has the matrix */
	SAFE_MALLOC_DECL(colid_zo,nnz,uint32_t); /* colid expands buffer */

	expandColid(buffer,colid_size,colid_zo,nnz);
	free(buffer);

	/* uint32_t last_start = 0 ; */
	i = 0 ;
	SAFE_CALLOC_DECL(pivots,m,uint32_t);

	uint32_t pivots_size = getSparsestRows(colid_zo,start_zo,m, pivots);
	/* i =  0 ; for( ; i < m ; ++i) fprintf(stderr, "%u ", pivots[i]); fprintf(stderr,"\n"); i = 0; */
	splitRowsUpBottom(A_init,C_init,colid_zo,start_zo,m,pivots,pivots_size,map_zo_pol);
	free(pivots);



	/* size of nnz for pols: */
	SAFE_READ_V(polys->nb,uint32_t,fh);

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_P(polys->start_pol,polys->nb+1,uint32_t,fh);

	SAFE_READ_P(polys->vals_pol,polys->start_pol[polys->nb],TYPE,fh);

	printMat(A_init);
	fprintf(stderr,"--------------\n");
	printMat(C_init);
	fprintf(stderr,"--------------\n");
	printPoly(polys);

	return 0 ;
}

uint32_t get_permute_columns(
		GBMatrix_t * A
		, uint32_t *  perm
		)
{
	uint32_t trans = 0 ;
	SAFE_CALLOC_DECL(col_size,A->col,uint32_t); /* sparsity of the columns */
	SAFE_CALLOC_DECL(last_elt,A->col-A->row,uint32_t); /* row id for the last element of each column out of the first square matrix*/
	uint32_t k = 0 ;
	for ( ; k < A->matrix_nb ; ++k ) {
		CSR_zo * A_k = &(A->matrix_zo[k]) ;
		uint32_t i = 0 ;
		for (; i < A_k->row ; ++i) {
			uint32_t j = A_k->start_zo[i] ;
			for ( ; j < A_k->start_zo[i+1] ; ++j) {
				col_size[A_k->colid_zo[j]] += 1 ;
				if (A_k->colid_zo[j] >= A->row )
					last_elt[A_k->colid_zo[j]-A->row] = (A->matrix_nb-1)*MAT_ROW_BLOCK+i ; /* zero based, but 0 is an index if column not empty :-) */
			}
		}
	}
	/* uint32_t i =  0 ; for( ; i < A->col ; ++i) fprintf(stderr, "%u ", col_size[i]); fprintf(stderr,"\n"); */
	/* i =  0 ; for( ; i < A->col-A->row ; ++i) fprintf(stderr, "%u ", last_elt[i]); fprintf(stderr,"\n"); */
	uint32_t j = A->row ;
	for ( ; j < A->col ; ++j) {
		if (col_size[j] == 0) continue ; /* empty column */
		uint32_t k = last_elt[j-A->row] ;
		if ( col_size[k] > col_size[j] ) { /* sparser column found */
			perm[k] = j ;
			++trans;
		}
	}
	/* fprintf(stderr,"trans : %u\n",trans); */
	return trans ;
}

void permuteCSR( CSR_zo * A_k , GBMatrix_t * A, uint32_t * start_b, uint32_t k
		, uint32_t * perm_i
		, uint32_t * perm_j
		, uint32_t perm_size
		, CSR_pol * polys
		)
{
	/* apply permutations (the rows keep the same length : in place) */
	uint32_t i = 0 ;
	for ( ; i < A_k->row ; ++i) {
		uint32_t here = 0 ;
		uint32_t j = A_k->start_zo[i] ;
		for (; j < A_k->start_zo[i+1] && here < perm_size ; ++j) {
			if (A_k->colid_zo[j] == perm_i[here])
				A_k->colid_zo[j] = perm_j[here++] ;
			else if (A_k->colid_zo[j] > perm_i[here]) {
				while ( here < perm_size && A_k->colid_zo[j] > perm_i[here])
					here++;
				if (here == perm_size)
					break ;
				if (A_k->colid_zo[j] == perm_i[here])
					A_k->colid_zo[j] = perm_j[here++] ;
			}
		}
	}
	assert(A_k->start_zo[A_k->row] == A_k->nnz);
	SAFE_MALLOC(A_k->data,A_k->nnz,TYPE);
	/* sort rows for we have messed them up */
	i = 0;
	for ( ; i < A_k->row ; ++i) {
		uint32_t j0 = A_k->start_zo[i] ;
		uint32_t j1 = A_k->start_zo[i+1] ;
		uint32_t start_idx = polys->start_pol[ A_k->map_zo_pol[i] ] ;
		fprintf(stderr,"row %u, poly %u, start %u\n",i,A_k->map_zo_pol[i],start_idx);
		copy(A_k->data+j0,polys->vals_pol+start_idx,j1-j0);
		/* uint32_t jj = 0 ; for ( ; jj < j1-j0 ; ++jj) fprintf(stderr,"%u ",A_k->data[j0+jj]) ; fprintf(stderr,"\n"); */
		/* insert_sort(&A_k->colid_zo[j0], j1-j0) ; */
		insert_sort_duo_data(A_k->colid_zo+j0, j1-j0, A_k->data+j0) ;
		/* jj = 0 ; for ( ; jj < j1-j0 ; ++jj) fprintf(stderr,"%u ",A_k->data[j0+jj]) ; fprintf(stderr,"\n"); */
	}

	fprintf(stderr,"permuted : ");
	i =  0 ; for( ; i < A_k->nnz ; ++i) fprintf(stderr, "%u ", A_k->colid_zo[i]); fprintf(stderr,"\n"); i = 0;

	/* where B starts in A_init */
	i = 0;
	for ( ; i < A_k->row ; ++i) {
		uint32_t j0 = A_k->start_zo[i] ;
		uint32_t j1 = A_k->start_zo[i+1] ;
		while(j0 < j1  && A_k->colid_zo[j0] < k)
			++j0;
		start_b[i] = j0 ;
	}
	fprintf(stderr,"col start: ");
	i =  0 ; for( ; i < A_k->row ; ++i) fprintf(stderr, "%u ", start_b[i]); fprintf(stderr,"\n"); i = 0;

	/* copy A */

	appendMatrix(A);
	CSR_zo * Ad = getLastMatrix(A);
	SAFE_MALLOC(Ad->start_zo,A_k->row+1,uint64_t);
	Ad->row = A_k->row ;
	Ad->col = k ;
	SAFE_MALLOC(Ad->map_zo_pol ,A_k->row,uint32_t);
	copy(Ad->map_zo_pol,A_k->map_zo_pol,A_k->row);

	i = 0;
	for ( ; i < A_k->row ; ++i) {
		Ad->nnz += start_b[i]-A_k->start_zo[i];
	}
	SAFE_MALLOC(Ad->colid_zo,Ad->nnz,uint32_t);
	SAFE_MALLOC(Ad->data,Ad->nnz,TYPE);

	i = 0;
	for ( ; i < A_k->row ; ++i) {
		uint32_t b = start_b[i]-A_k->start_zo[i];
		fprintf(stderr,"b : %u\n",b);
		Ad->start_zo[i+1]= Ad->start_zo[i]+b;
		fprintf(stderr,"start %u+1=%lu\n",i,Ad->start_zo[i+1]);
		uint32_t j = A_k->start_zo[i] ;
		copy(Ad->colid_zo+Ad->start_zo[i], A_k->colid_zo+A_k->start_zo[i], b);
		copy(Ad->data+Ad->start_zo[i], A_k->data+A_k->start_zo[i], b);
	}
	assert(Ad->start_zo[Ad->row] == Ad->nnz);


	fprintf(stderr,"Ad:");
	printMatUnit(Ad);
}

void do_permute_columns_up(
		GBMatrix_t * A_init,
		GBMatrix_t * A,
		GBMatrix_t * Bt
		, uint32_t * perm_i
		, uint32_t * perm_j
		, uint32_t perm_size
		, CSR_pol * polys
		)
{
	uint32_t k = A_init->row;

	A->row = Bt->col = k ;
	A->mod = Bt->mod = A_init->mod;
	A->col = k ;
	Bt->row= A_init->col - k ;

	uint32_t blk = 0;
	for ( ; blk < A_init->matrix_nb ; ++blk) {
		uint32_t i = 0 ;

		CSR_zo * A_k = &(A_init->matrix_zo[blk]); /* A_k band of A_init */

		SAFE_MALLOC_DECL(start_b,A_k->row,uint32_t);

		/*  permute columns and create A  */

		permuteCSR(A_k,A,start_b,k,perm_i,perm_j,perm_size,polys);


		/* Transpose B */

		appendMatrix(Bt);
		CSR_zo * Bd = getLastMatrix(Bt);

		Bd->col = A_k->row ;
		Bd->row = A_init->col - A_init->row ;
		SAFE_CALLOC(Bd->start_zo,Bd->row+1,uint64_t);
		uint64_t other_nnz = getLastMatrix(A)->nnz ;
		Bd->nnz = A_k->nnz - other_nnz ;
		SAFE_MALLOC(Bd->colid_zo,Bd->nnz,uint32_t);
		SAFE_MALLOC(Bd->map_zo_pol ,A_k->row,uint32_t);
		copy(Bd->map_zo_pol,A_k->map_zo_pol,A_k->row); /* useless */

		i = 0 ;
		for ( ; i < A_k->row ; ++i) {
			uint32_t j = start_b[i] ;
			for (; j < A_k->start_zo[i+1] ; ++j) {
				Bd->start_zo[A_k->colid_zo[j]-k+1] += 1 ;

			}
		}
		fprintf(stderr,"B row start: ");
		/* i =  0 ; for( ; i < Bd->row+1 ; ++i) fprintf(stderr, "%u ", Bd->start_zo[i]); fprintf(stderr,"\n"); i = 0; */


		i = 0 ;
		for ( ; i < Bd->row+1 ; ++i) {
			Bd->start_zo[i+1] += Bd->start_zo[i] ;
		}

		fprintf(stderr,"B row start: ");
		/* i =  0 ; for( ; i < Bd->row+1 ; ++i) fprintf(stderr, "%u ", Bd->start_zo[i]); fprintf(stderr,"\n"); i = 0; */


		SAFE_CALLOC_DECL(done_col,Bd->row,uint32_t);
		SAFE_MALLOC(Bd->data, Bd->nnz, TYPE);

		uint32_t j = 0 ;
		for (; j < A_k->row ; ++j) {
			i = start_b[j];
			while (i < A_k->start_zo[j+1]){
				uint32_t cur_place = Bd->start_zo[A_k->colid_zo[i]-k];
				cur_place  += done_col[A_k->colid_zo[i]-k] ;
				Bd->data [ cur_place ] = A_k->data[i] ;
				Bd->colid_zo[ cur_place ] =  j ;
				done_col[A_k->colid_zo[i]-k] += 1 ;
				++i;
			}
		}

		fprintf(stderr,"Bd:");
		printMatUnit(Bd);


		A->nnz += other_nnz ;
		Bt->nnz += Bd->nnz ;
	}

	fprintf(stderr,"A:");
	printMat(A);

	fprintf(stderr,"Bt:");
	printMat(Bt);
}


#if 1
void do_permute_columns_lo(
		uint32_t k,
		GBMatrix_t * C_init,
		GBMatrix_t * C,
		DenseMatrix_t * D
		, uint32_t * perm_i
		, uint32_t * perm_j
		, uint32_t perm_size
		, CSR_pol * polys
		)
{
	C->row = D->row = C_init->row ;
	C->mod = D->mod = C_init->mod;
	C->col = k ;
	D->col = C_init->col - k ;

	SAFE_CALLOC(D->data,D->row*D->col,TYPE);

	uint32_t blk = 0;
	for ( ; blk < C_init->matrix_nb ; ++blk) {
		uint32_t i = 0 ;

		CSR_zo * C_k = &(C_init->matrix_zo[blk]); /* A_k band of A_init */
		SAFE_MALLOC_DECL(start_b,C_k->row,uint32_t);

		/* C */
		permuteCSR(C_k,C,start_b,k,perm_i,perm_j,perm_size,polys);
		uint64_t other_nnz = getLastMatrix(C)->nnz ;

		C->nnz += other_nnz  ;

		/*  D */
		uint32_t row_offset = blk*MAT_ROW_BLOCK ;
		TYPE * data = D->data + (row_offset * D->col) ;

		i = 0 ;
		for ( ; i < C_k->row ; ++i) {
			uint32_t j = start_b[i] ;
			for (; j < C_k->start_zo[i+1] ; ++j) {
				uint32_t jj =  C_k->colid_zo[j]-k;
				data[D->col*i+jj] = C_k->data[j];
			}
		}
	}

	fprintf(stderr,"C:");
	printMat(C);

	fprintf(stderr,"D:");
	printMatDense(D);


}
#endif





/* transforms perm where perm[i] = j into two lists
 * such that ther permutations are (perm_i[k],perm_j[k]) for perm and perm^(-1).
 * perm_i is increasing.
 */
uint32_t split_permutations(
		const uint32_t * perm
		, uint32_t perm_size
		, uint32_t * perm_i
		, uint32_t * perm_j)
{
	uint32_t here = 0;
	uint32_t i = 0;
	for ( ; i < perm_size ; ++i)
		if (perm[i] != i) {
			perm_i[here] = i ;
			perm_j[here] = perm[i] ;
			++here;
		}
	i = 0 ;
	for ( ; i < perm_size ; ++i)
		if (perm[i] != i) {
			perm_i[here] = perm[i] ;
			perm_j[here] = i ;
			++here;
		}

	/* the last ones are not necessarily in order */
	insert_sort_duo(perm_i+here/2,here/2,perm_j+here/2);
	return here ;
}

void split_columns(
		GBMatrix_t * A_init
		, GBMatrix_t * C_init
		, CSR_pol * polys
		, GBMatrix_t * A
		, GBMatrix_t * Bt
		, GBMatrix_t * C
		, DenseMatrix_t * D
		)
{
	/* search for columns to permute to make A sparser */
	SAFE_MALLOC_DECL(perm,A->row,uint32_t); /* what columns the first row ones should be exchanged with */
	uint32_t i = 0;
	for ( ; i < A_init->row ; ++i)
		perm[i] = i ;
	uint32_t trans = get_permute_columns(A_init,perm);

	/* fprintf(stderr,"--------------\n"); */
	/* i =  0 ; for( ; i < A_init->row ; ++i) fprintf(stderr, "%u ", perm[i]); fprintf(stderr,"\n"); i = 0; */
	/* fprintf(stderr,"--------------\n"); */


	SAFE_MALLOC_DECL(perm_i,2*trans,uint32_t);
	SAFE_MALLOC_DECL(perm_j,2*trans,uint32_t);
	split_permutations(perm,A_init->row,perm_i,perm_j);
	/* i =  0 ; for( ; i < 2*trans ; ++i) fprintf(stderr, "%u ", perm_i[i]); fprintf(stderr,"\n"); i = 0; */
	/* i =  0 ; for( ; i < 2*trans ; ++i) fprintf(stderr, "%u ", perm_j[i]); fprintf(stderr,"\n"); i = 0; */

	/* copy last columns of A,C to B,D in proper format */
	do_permute_columns_up(A_init,A,Bt,perm_i, perm_j,2*trans,polys);
	do_permute_columns_lo(A_init->row,C_init,C,D,perm_i,perm_j,2*trans,polys);

	return ;
}


#endif /* __GB_io_H */
/* vim: set ft=c: */
