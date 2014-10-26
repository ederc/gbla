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

void init(CSR_zo * mat)
{
	mat->row=0;
	mat->col=0;
	SAFE_MALLOC(mat->start_zo,1,uint64_t);
	mat->start_zo[0]=0;
	SAFE_MALLOC(mat->colid_zo,0,uint32_t);
	SAFE_MALLOC(mat->map_zo_pol,0,uint32_t);
}

void appendRow(CSR_zo * mat
		, uint32_t * colid
		, uint64_t size
		, uint32_t pol
		)
{
	uint64_t old = mat->start_zo[mat->row] ;
	uint64_t nnz = old + size;
	SAFE_REALLOC(mat->colid_zo,nnz,uint32_t);
	uint64_t i = 0 ;
	for ( ; i < size ; ++i) {
		mat->colid_zo[old+i] = colid[i] ;
	}
	mat->row ++ ;
	SAFE_REALLOC(mat->start_zo,mat->row+1,uint64_t);
	mat->start_zo[mat->row] = nnz  ;
	SAFE_REALLOC(mat->map_zo_pol,mat->row,uint32_t);
	mat->map_zo_pol[mat->row-1] = pol;
}

typedef struct GBMatrix_t {

	uint32_t row ;
	uint32_t col ;
	uint32_t nnz ;
	uint32_t mod ;
	/* uint32_t block_size ; */
	uint32_t matrix_nb ;  /* nb of 0/1 matrices */
	CSR_zo  * matrix_zo ; /* 0/1 matrices reprensenting positions */
} GBMatrix_t;


CSR_zo * getLastMatrix(GBMatrix_t * A)
{
	return &(A->matrix_zo[A->matrix_nb]);
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
	if (A->matrix_zo->row == MAT_ROW_BLOCK) {
		A->matrix_nb++;
		SAFE_REALLOC(A->matrix_zo,A->matrix_nb,CSR_zo);
		init(&(A->matrix_zo[A->matrix_nb]));
	}
	appendRow(&(A->matrix_zo[A->matrix_nb]),colid,size,pol);
}

#endif /* __GB_matrix_H */
