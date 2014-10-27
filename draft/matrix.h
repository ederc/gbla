#ifndef __GB_matrix_H
#define __GB_matrix_H

#include "types.h"

/*
   struct matrixFile_t {

   uint32_t row ;
   uint32_t col ;
   uint64_t nnz ;
   TYPE     mod ;
   uint32_t * start_zo ;
   uint32_t * map_zo_pol ;
   uint32_t * colid_zo ;
   uint32_t * start_pol;
   TYPE     * vals_pol;
   } ;
   */

typedef struct CSR_pol {
	uint32_t nb ;
	uint32_t * start_pol ;
	TYPE     * vals_pol ;
} CSR_pol;

typedef struct CSR_zo {
	uint32_t row;
	uint32_t col;
	uint64_t nnz ; /* this is start_zo[row], I only it were c++ */
	uint64_t * start_zo ;
	uint32_t * colid_zo ;
	uint32_t * map_zo_pol ;

} CSR_zo;

uint32_t * getRow(CSR_zo * mat, uint32_t i)
{
	return mat->colid_zo + mat->start_zo[i];
}

uint64_t size(CSR_zo * mat)
{
	return mat->start_zo[mat->row] ;
}

void initUnit(CSR_zo * mat)
{
	if (!mat) return ;
	mat->row=0;
	mat->col=0;
	SAFE_MALLOC(mat->start_zo,1,uint64_t);
	mat->start_zo[0]=0;
	SAFE_MALLOC(mat->colid_zo,0,uint32_t);
	SAFE_MALLOC(mat->map_zo_pol,0,uint32_t);
}

void appendRowUnit(CSR_zo * mat
		, uint32_t * colid
		, uint64_t size
		, uint32_t pol
		)
{
	uint64_t old = mat->start_zo[mat->row] ;
	mat->nnz = old + size;
	SAFE_REALLOC(mat->colid_zo,mat->nnz,uint32_t);
	uint64_t i = 0 ;
	for ( ; i < size ; ++i) {
		mat->colid_zo[old+i] = colid[i] ;
	}
	mat->row ++ ;
	SAFE_REALLOC(mat->start_zo,mat->row+1,uint64_t);
	mat->start_zo[mat->row] = mat->nnz  ;
	SAFE_REALLOC(mat->map_zo_pol,mat->row,uint32_t);
	mat->map_zo_pol[mat->row-1] = pol;
}

typedef struct GBMatrix_t {

	uint32_t row ;
	uint32_t col ;
	uint64_t nnz ;
	uint32_t mod ;
	/* uint32_t block_size ; */
	uint32_t matrix_nb ;  /* nb of 0/1 matrices */
	CSR_zo  * matrix_zo ; /* 0/1 matrices reprensenting positions */
} GBMatrix_t;


CSR_zo * getLastMatrix(GBMatrix_t * A)
{
	return &(A->matrix_zo[A->matrix_nb-1]);
}

typedef struct DenseMatrix_t {
	uint32_t row ;
	uint32_t col ;
	TYPE *  data ;
} DenseMatrix_t ;

void appendRow(GBMatrix_t * A
		, uint32_t * colid
		, uint64_t size
		, uint32_t pol
		)
{
	A->row += 1;
	A->nnz += size ;

	if (A->matrix_nb == 0) {
		A->matrix_nb++;
		SAFE_MALLOC(A->matrix_zo,A->matrix_nb,CSR_zo);
		initUnit(&(A->matrix_zo[A->matrix_nb-1]));
		(A->matrix_zo[A->matrix_nb-1]).col = A->col ;
	}
	else if ((A->matrix_zo[A->matrix_nb-1]).row == MAT_ROW_BLOCK) {
		A->matrix_nb++;
		SAFE_REALLOC(A->matrix_zo,A->matrix_nb,CSR_zo);
		initUnit(&(A->matrix_zo[A->matrix_nb-1]));
		(A->matrix_zo[A->matrix_nb-1]).col = A->col ;
	}
	appendRowUnit(&(A->matrix_zo[A->matrix_nb-1]),colid,size,pol);
}
void init(GBMatrix_t * A)
{
	if (!A) return;
	A->row = 0 ;
	A->col = 0 ;
	A->nnz = 0 ;
	A->matrix_nb = 0 ;
	SAFE_MALLOC(A->matrix_zo,A->matrix_nb,CSR_zo);
	initUnit(A->matrix_zo);
}

void printMatUnit(CSR_zo * A)
{
	fprintf(stderr,"block %u x %u - %lu",A->row, A->col,A->nnz);
	fprintf(stderr,"\nstart:\n");
	uint32_t i = 0 ;
	for ( ; i < A->row+1 ; ++i) {
		fprintf(stderr,"%lu ", A->start_zo[i]);
	}
	fprintf(stderr,"\ndata:\n");
	i = 0 ;
	for ( ; i < A->row ; ++i) {
		uint32_t j = A->start_zo[i] ;
		for ( ; j < A->start_zo[i+1] ; ++j) {
			fprintf(stderr,"%u ", A->colid_zo[j]);
		}
		fprintf(stderr,"|");
	}
	fprintf(stderr,"\npolys:\n");
	i = 0 ;
	for ( ; i < A->row ; ++i) {
		fprintf(stderr,"%u ", A->map_zo_pol[i]);
	}
	fprintf(stderr,"\n");
}

void printMat(GBMatrix_t * A)
{
	fprintf(stderr,"matrix %u x %u - %lu\n",A->row, A->col,A->nnz);
	fprintf(stderr,"mod %u\n",A->mod);
	uint32_t k = 0 ;
	for (  ; k < A->matrix_nb ; ++k ) {
		printMatUnit(&(A->matrix_zo[k]));
	}
}

void printPoly(CSR_pol * P)
{
	fprintf(stderr,"polys (%u)\n",P->nb);
	fprintf(stderr,"start\n");
	uint32_t i = 0 ;
	for ( ; i < P->nb+1 ; ++i) {
		fprintf(stderr,"%u ", P->start_pol[i]);
	}
	fprintf(stderr,"\ndata:\n");
	i = 0 ;
	for ( ; i < P->nb ; ++i) {
		uint32_t j = P->start_pol[i] ;
		for ( ; j < P->start_pol[i+1] ; ++j) {
			fprintf(stderr,"%u ", P->vals_pol[j]);
		}
		fprintf(stderr,"\n");
	}
}

#endif /* __GB_matrix_H */
/* vim: set ft=c: */
