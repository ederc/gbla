#ifndef __GB_matrix_H
#define __GB_matrix_H

#include "types.h"
#include "tools.h"
#include "printer.h"



typedef struct CSR_pol {
	taille_t nb ;
	taille_t * start_pol ;
	elem_t   * data_pol ;
} CSR_pol;

typedef struct CSR_zo {
	taille_t row;
	taille_t col;
	index_t nnz ; /* this is start_zo[row] */
	index_t  * start_zo ;
	taille_t * colid_zo ;
	taille_t * map_zo_pol ;
	elem_t   * data ;

} CSR_zo;

taille_t * getRow(CSR_zo * mat, taille_t i)
{
	return mat->colid_zo + mat->start_zo[i];
}

index_t size(CSR_zo * mat)
{
	return mat->start_zo[mat->row] ;
}

void initUnit(CSR_zo * mat)
{
	/* if (!mat) return ; */
	assert(mat);
	mat->row=0;
	mat->col=0;
	mat->nnz=0;
	/* mat->mod=0; */
	SAFE_MALLOC(mat->start_zo,1,index_t);
	mat->start_zo[0]=0;
	/* SAFE_MALLOC(mat->colid_zo,0,taille_t); */
	/* SAFE_MALLOC(mat->map_zo_pol,0,taille_t); */
	/* SAFE_MALLOC(mat->data,0,elem_t); */
	mat->colid_zo = NULL ;
	mat->map_zo_pol = NULL ;
	mat->data = NULL ;
}

void appendRowUnit(CSR_zo * mat
		, taille_t * colid
		, index_t size
		, taille_t pol
		)
{
	index_t old = mat->start_zo[mat->row] ;
	mat->nnz = old + size;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->colid_zo,mat->nnz,taille_t);
	index_t i = 0 ;
	for ( ; i < size ; ++i) {
		mat->colid_zo[old+i] = colid[i] ;
	}
	mat->row ++ ;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->start_zo,mat->row+1,index_t);
	mat->start_zo[mat->row] = mat->nnz  ;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->map_zo_pol,mat->row,taille_t);
	mat->map_zo_pol[mat->row-1] = pol;
}

typedef struct GBMatrix_t {

	taille_t row ;
	taille_t col ;
	index_t nnz ;
	elem_t   mod ;
	/* taille_t block_size ; */
	taille_t matrix_nb ;  /* nb of 0/1 matrices */
	CSR_zo  * matrix_zo ; /* 0/1 matrices reprensenting positions */
} GBMatrix_t;


const CSR_zo * getLastMatrixConst(const GBMatrix_t * A)
{
	return &(A->matrix_zo[A->matrix_nb-1]);
}

CSR_zo * getLastMatrix(GBMatrix_t * A)
{
	return &(A->matrix_zo[A->matrix_nb-1]);
}


typedef struct DenseMatrix_t {
	taille_t row ;
	taille_t col ;
	elem_t   mod ;
	elem_t *  data ;
} DenseMatrix_t ;

void appendMatrix(GBMatrix_t * A)
{
	A->matrix_nb++;
	if (A->matrix_zo == NULL) {
		assert(A->matrix_nb == 1);
		SAFE_MALLOC(A->matrix_zo,A->matrix_nb,CSR_zo);
	}
	else {
		assert(A->matrix_nb > 1);
		SAFE_REALLOC(A->matrix_zo,A->matrix_nb,CSR_zo);
	}
	initUnit(&(A->matrix_zo[A->matrix_nb-1]));
	(A->matrix_zo[A->matrix_nb-1]).col = A->col ;
}

void appendRow(GBMatrix_t * A
		, taille_t * colid
		, index_t size
		, taille_t pol
		)
{
	A->row += 1;
	A->nnz += size ;

	if ( ( A->matrix_nb == 0) || ((A->matrix_zo[A->matrix_nb-1]).row == MAT_ROW_BLOCK) ){
		appendMatrix(A);
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
	/* SAFE_MALLOC(A->matrix_zo,A->matrix_nb,CSR_zo); */
	/* initUnit(A->matrix_zo); */
	A->matrix_zo = NULL ;
}

void printMatUnit(CSR_zo * A)
{
	fprintf(stderr,"block %u x %u - %lu",A->row, A->col,A->nnz);
	fprintf(stderr,"\nstart:\n<");
	taille_t i = 0 ;
	for ( ; i < A->row+1 ; ++i) {
		fprintf(stderr,"%lu ", A->start_zo[i]);
	}
	fprintf(stderr,">\ncolumns\n<");
	i = 0 ;
	for ( ; i < A->row ; ++i) {
		index_t j = A->start_zo[i] ;
		for ( ; j < A->start_zo[i+1] ; ++j) {
			fprintf(stderr,"%u ", A->colid_zo[j]);
		}
		fprintf(stderr,"|");
	}
	fprintf(stderr,">\npolys:\n<");
	i = 0 ;
	for ( ; i < A->row ; ++i) {
		fprintf(stderr,"%u ", A->map_zo_pol[i]);
	}

	if (A->data != NULL) {
		fprintf(stderr,">\nDATA:\n<");
		i = 0;
		for ( ; i < A->row ; ++i) {
			index_t j = A->start_zo[i] ;
			for ( ; j < A->start_zo[i+1] ; ++j) {
				Mjoin(print,elem_t)(A->data[j]);
				fprintf(stderr," ");
			}
			fprintf(stderr,"|");
		}
	}
	fprintf(stderr,">\n");
}

void printMat(GBMatrix_t * A)
{
	fprintf(stderr,"matrix %u x %u - %lu\n",A->row, A->col,A->nnz);
	fprintf(stderr,"mod ");
	Mjoin(print,elem_t)(A->mod);
	fprintf(stderr,"\n");
	taille_t k = 0 ;
	for (  ; k < A->matrix_nb ; ++k ) {
		printMatUnit(&(A->matrix_zo[k]));
	}
}

void printMatDense(DenseMatrix_t * A)
{
	fprintf(stderr,"matrix %u x %u - %lu\n",A->row, A->col, (index_t)A->row*(index_t)A->col);
	fprintf(stderr,"mod ");
	Mjoin(print,elem_t)(A->mod);
	fprintf(stderr,"\n");

	taille_t i = 0 ;
	for (  ; i < A->row ; ++i ) {
		taille_t j = 0 ;
		for (  ; j < A->col ; ++j ) {
			Mjoin(print,elem_t)(A->data[A->col*i+j]);
			fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
	}
}

void printPoly(CSR_pol * P)
{
	fprintf(stderr,"polys (%u)\n",P->nb);
	fprintf(stderr,"start\n");
	taille_t i = 0 ;
	for ( ; i < P->nb+1 ; ++i) {
		fprintf(stderr,"%u ", P->start_pol[i]);
	}
	fprintf(stderr,"\ndata:\n");
	i = 0 ;
	for ( ; i < P->nb ; ++i) {
		taille_t j = P->start_pol[i] ;
		for ( ; j < P->start_pol[i+1] ; ++j) {
			Mjoin(print,elem_t)(P->data_pol[j]);
			fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
	}
}


void checkMatUnit(const CSR_zo *Ak)
{

	if (Ak == NULL) exit(-1);
#ifndef NDEBUG
	taille_t i = 0 ;
	index_t jz = 0 ;
	assert(Ak->start_zo[0] == 0);
	for ( i = 0 ; i < Ak->row ; ++i) {
		assert(Ak->start_zo[i+1] <= Ak->nnz);
		for (jz = Ak->start_zo[i] ; jz < Ak->start_zo[i+1] ; ++jz) {
			taille_t k = Ak->colid_zo[jz];
			assert(k < Ak->col);
		}
	}
	assert(Ak->start_zo[Ak->row] == Ak->nnz);

#endif
	/* fprintf(stderr,"ok\n"); */
}

void checkMat(const GBMatrix_t *A)
{
	if (A == NULL) exit(-1);
#ifndef NDEBUG
	taille_t row = 0 ;
	index_t nnz = 0 ;
	taille_t i = 0 ;
	for ( ; i < A->matrix_nb ; ++i) {
		row += A->matrix_zo[i].row;
		nnz += A->matrix_zo[i].nnz;
	}
	assert (nnz == A->nnz);
	assert (row == A->row);

	assert( A->matrix_nb == 1);

	const CSR_zo * Ak = &(A->matrix_zo[0]);
	checkMatUnit(Ak);
#endif

}

void freeMatDense(DenseMatrix_t * A)
{
	free(A->data);
}

void freeMatUnit(CSR_zo * A)
{
	free(A->start_zo);
	free(A->colid_zo);
	free(A->data);
	free(A->map_zo_pol);
	/* free(A); */
}

void freePol( CSR_pol * A)
{
	free (A->start_pol);
	free (A->data_pol);
}


void freeMat(GBMatrix_t * A)
{
	taille_t i;
	for (i=0 ; i < A->matrix_nb ; ++i) {
		freeMatUnit(&(A->matrix_zo[i])) ;
	}
	free(A->matrix_zo);
}
#endif /* __GB_matrix_H */
/* vim: set ft=c: */
