#ifndef __GB_matrix_H
#define __GB_matrix_H

#include "config.h"
#include "types.h"
#include "macros.h"
#include "printer.h"


typedef struct CSR_pol {
	dimen_t nb ;
	dimen_t * start_pol ;
	elemt_t   * data_pol ;
} CSR_pol;

typedef struct CSR {
	dimen_t row;
	dimen_t col;
	index_t  nnz ; /* this is start[row] */
	index_t  * start ;
	dimen_t * colid ;
	dimen_t * map_zo_pol ;
	elemt_t   * data ;

} CSR;

typedef struct GBMatrix_t {

	dimen_t   row ;
	dimen_t   col ;
	index_t    nnz ;
	elemt_t     mod ;
	dimen_t   sub_nb ;  /* nb of 0/1 matrices */
	CSR      * sub ; /* 0/1 matrices reprensenting positions */
} GBMatrix_t;

typedef struct DNS {
	dimen_t  row ;
	dimen_t  col ;
	dimen_t  ld  ;
	elemt_t    mod ;
	index_t   nnz ;
	elemt_t  * ptr ;
} DNS ;

typedef struct DenseMatrix_t {
	dimen_t   row ;
	dimen_t   col ;
	elemt_t     mod ;
	index_t    nnz ;
	dimen_t   blk_nb ;
	DNS      * blk ;
} DenseMatrix_t ;

/* ACCESS */

index_t size(CSR * mat)
{
	return mat->start[mat->row] ;
}

dimen_t * getRow(CSR * mat, dimen_t i)
{
	return mat->colid + mat->start[i];
}

const CSR * getLastMatrixConst(const GBMatrix_t * A)
{
	return &(A->sub[A->sub_nb-1]);
}

CSR * getLastMatrix(GBMatrix_t * A)
{
	return &(A->sub[A->sub_nb-1]);
}

/* INIT */

void initSparseUnit(CSR * mat)
{
	/* if (!mat) return ; */
	assert(mat);
	mat->row=0;
	mat->col=0;
	mat->nnz=0;
	/* mat->mod=0; */
	SAFE_MALLOC(mat->start,1,index_t);
	mat->start[0]=0;
	/* SAFE_MALLOC(mat->colid,0,dimen_t); */
	/* SAFE_MALLOC(mat->map_zo_pol,0,dimen_t); */
	/* SAFE_MALLOC(mat->data,0,elemt_t); */
	mat->colid = NULL ;
	mat->map_zo_pol = NULL ;
	mat->data = NULL ;
}

void initSparse(GBMatrix_t * A)
{
	if (!A) return;
	A->row = 0 ;
	A->col = 0 ;
	A->nnz = 0 ;
	A->sub_nb = 0 ;
	/* SAFE_MALLOC(A->sub,A->sub_nb,CSR); */
	/* initSparseUnit(A->sub); */
	A->sub = NULL ;
}

void initDenseUnit (DNS * A)
{
	A->row = 0 ;
	A->col = 0 ;
	A->ld  = 0 ;
	A->mod = 0 ;
	A->nnz = 0 ;
	A->ptr = NULL ;
}

/* GROW */

void appendRowUnit(CSR * mat
		, dimen_t * colid
		, index_t size
		, dimen_t pol
		)
{
	index_t old = mat->start[mat->row] ;
	mat->nnz = old + size;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->colid,mat->nnz,dimen_t);
	index_t i = 0 ;
	for ( ; i < size ; ++i) {
		mat->colid[old+i] = colid[i] ;
	}
	mat->row ++ ;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->start,mat->row+1,index_t);
	mat->start[mat->row] = mat->nnz  ;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->map_zo_pol,mat->row,dimen_t);
	mat->map_zo_pol[mat->row-1] = pol;
}


void appendMatrix(GBMatrix_t * A)
{
	A->sub_nb++;
	if (A->sub == NULL) {
		assert(A->sub_nb == 1);
		SAFE_MALLOC(A->sub,A->sub_nb,CSR);
	}
	else {
		assert(A->sub_nb > 1);
		SAFE_REALLOC(A->sub,A->sub_nb,CSR);
	}
	initSparseUnit(&(A->sub[A->sub_nb-1]));
	(A->sub[A->sub_nb-1]).col = A->col ;
}

void appendRow(GBMatrix_t * A
		, dimen_t * colid
		, index_t size
		, dimen_t pol
		)
{
	A->row += 1;
	A->nnz += size ;

	if ( ( A->sub_nb == 0) || ((A->sub[A->sub_nb-1]).row == MAT_ROW_BLK) ){
		appendMatrix(A);
	}

	appendRowUnit(&(A->sub[A->sub_nb-1]),colid,size,pol);
}

/* PRINT */

void printMatUnit(CSR * A)
{
	fprintf(stderr,"block %u x %u - %lu",A->row, A->col,A->nnz);
	fprintf(stderr,"\nstart:\n<");
	dimen_t i = 0 ;
	for ( ; i < A->row+1 ; ++i) {
		fprintf(stderr,"%lu ", A->start[i]);
	}
	fprintf(stderr,">\ncolumns\n<");
	i = 0 ;
	for ( ; i < A->row ; ++i) {
		index_t j = A->start[i] ;
		for ( ; j < A->start[i+1] ; ++j) {
			fprintf(stderr,"%u ", A->colid[j]);
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
			index_t j = A->start[i] ;
			for ( ; j < A->start[i+1] ; ++j) {
				Mjoin(print,elemt_t)(A->data[j]);
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
	Mjoin(print,elemt_t)(A->mod);
	fprintf(stderr,"\n");
	dimen_t k = 0 ;
	for (  ; k < A->sub_nb ; ++k ) {
		printMatUnit(&(A->sub[k]));
	}
}

void printMatDense(DNS * A)
{
	fprintf(stderr,"matrix %u x %u - %lu\n",A->row, A->col, A->nnz);
	fprintf(stderr,"mod ");
	Mjoin(print,elemt_t)(A->mod);
	fprintf(stderr,"\n");

	dimen_t i = 0 ;
	for (  ; i < A->row ; ++i ) {
		dimen_t j = 0 ;
		for (  ; j < A->col ; ++j ) {
			Mjoin(print,elemt_t)(A->ptr[A->ld*i+j]);
			fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
	}
}

void printPoly(CSR_pol * P)
{
	fprintf(stderr,"polys (%u)\n",P->nb);
	fprintf(stderr,"start\n");
	dimen_t i = 0 ;
	for ( ; i < P->nb+1 ; ++i) {
		fprintf(stderr,"%u ", P->start_pol[i]);
	}
	fprintf(stderr,"\ndata:\n");
	i = 0 ;
	for ( ; i < P->nb ; ++i) {
		dimen_t j = P->start_pol[i] ;
		for ( ; j < P->start_pol[i+1] ; ++j) {
			Mjoin(print,elemt_t)(P->data_pol[j]);
			fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
	}
}

/* STATS */

index_t occupancySparse(GBMatrix_t * A)
{
	dimen_t k ;
	index_t  acc = 0;
	for ( k = 0 ; k < A->sub_nb  ; ++k ) {
		CSR * Ad = &(A->sub[k]);
		index_t j ;
		dimen_t i ;
		for (i = 0 ; i < Ad->row ; ++i) {
			SAFE_CALLOC_DECL(occup,(Ad->col/UNRL+1),dimen_t);
			for ( j = Ad->start[i] ; j < Ad->start[i+1] ; ++j) {
				occup[Ad->colid[j]/UNRL] += 1 ;
			}
			for ( j = 0 ; j < Ad->col/UNRL+1 ; ++j) {
				if (occup[j] > 0)
					acc += 1 ;
			}
			free(occup);

		}
	}
	acc *= UNRL ;

	return acc ;
}

index_t occupancyDense ( DNS * D)
{
	dimen_t i,j ;
	index_t  acc = 0;
	for (i = 0 ; i < D->row ; ++i) {
		SAFE_CALLOC_DECL(occup,(D->col/UNRL+1),dimen_t);
		for (j = 0 ; j < D->col ; ++j) {
			index_t k = (index_t)i*(index_t)D->col+j;
			if (D->ptr[k] != 0) {
				occup[j/UNRL] += 1 ;
			}
		}
		dimen_t k ;
		for ( k = 0 ; k < D->col/UNRL+1 ; ++k) {
			if (occup[k] > 0)
				acc += 1 ;
		}
		free(occup);

	}
	acc *= UNRL ;
	return acc ;
}


/* CHECK */

void checkMatUnit(const CSR *Ak)
{

	if (Ak == NULL) exit(-1);
#ifndef NDEBUG
	dimen_t i = 0 ;
	index_t jz = 0 ;
	assert(Ak->start[0] == 0);
	for ( i = 0 ; i < Ak->row ; ++i) {
		assert(Ak->start[i+1] <= Ak->nnz);
		for (jz = Ak->start[i] ; jz < Ak->start[i+1] ; ++jz) {
			dimen_t k = Ak->colid[jz];
			assert(k < Ak->col);
		}
	}
	assert(Ak->start[Ak->row] == Ak->nnz);

#endif
	/* fprintf(stderr,"ok\n"); */
}

void checkMat(const GBMatrix_t *A)
{
	if (A == NULL) exit(-1);
#ifndef NDEBUG
	dimen_t row = 0 ;
	index_t nnz = 0 ;
	dimen_t i = 0 ;
	for ( ; i < A->sub_nb ; ++i) {
		row += A->sub[i].row;
		nnz += A->sub[i].nnz;
		const CSR * Ak = &(A->sub[i]);
		checkMatUnit(Ak);
	}
	assert (nnz == A->nnz);
	assert (row == A->row);

	/* assert( A->sub_nb == 1); */
#endif

}

/* FREE */

void freeMatDense(DNS * A)
{
	free(A->ptr);
}

void freeMatUnit(CSR * A)
{
	free(A->start);
	free(A->colid);
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
	dimen_t i;
	for (i=0 ; i < A->sub_nb ; ++i) {
		freeMatUnit(&(A->sub[i])) ;
	}
	free(A->sub);
}

/* CONVERSION */

void convert_CSR_2_DNS(DNS * D, const GBMatrix_t * S )
{
	D->row = S->row ;
	D->col = S->col ;
	D->ld = ALIGN(D->col) ;
	D->mod = S->mod ;
	SAFE_CALLOC(D->ptr,(index_t)D->row * (index_t)D->ld,elemt_t);

	dimen_t k ;
	for (k = 0 ; k < S->sub_nb ; ++k) {
		CSR * S_k = &(S->sub[k]) ;
		dimen_t i ;
		index_t i_off = k * MAT_ROW_BLK ;
		for (i = 0 ; i < S_k->row ; ++i) {
			index_t jz ;
			for ( jz = S_k->start[i] ; jz < S_k->start[i+1] ; ++jz) {
				dimen_t j = S_k->colid[jz] ;
				D->ptr[(i_off+i)*D->ld+j] = S_k->data[jz] ;
			}
		}
	}


}


#ifdef BLOCK_CSR
void convert_CSR_2_CSR_block(GBMatrix_t * B, const GBMatrix_t * S )
{
	B->row = S->row;
	B->col = S->col;
	B->nnz = S->nnz;
	B->mod = S->mod;
	dimen_t k ;
	for (k = 0 ; k < S->sub_nb ; ++k) {
		const CSR * B_k = &(S->sub[k]) ;
		appendMatrix(B);
		CSR * Bd = getLastMatrix(B);
		Bd->row = B_k->row ;
		Bd->col = B_k->col ;
		Bd->nnz = B_k->nnz ;
		SAFE_REALLOC(Bd->start,Bd->row+1,index_t);
		dimen_t i ;
		for (i = 0 ; i <= Bd->row ; ++i) {
			Bd->start[i] = 0 ;
		}
		SAFE_CALLOC(Bd->data,UNRL*ALIGN(Bd->nnz),elemt_t);
		SAFE_MALLOC(Bd->colid,ALIGN(Bd->nnz),dimen_t);
		/* index_t there = (index_t)-1 ; */

		SAFE_CALLOC_DECL(b_there,(Bd->row+1),index_t);
		b_there[0] = (index_t)-1 ;
		for (i = 0 ; i < B_k->row ; ++i) {
			index_t jz ;
			dimen_t last_j = (dimen_t) -1 ;
			for (jz = B_k->start[i] ; jz < B_k->start[i+1] ; ++jz) {
				dimen_t j = B_k->colid[jz] ;
				if (j/UNRL  != last_j) {
					last_j = j/UNRL ;
					b_there[i+1] += 1 ;
				}
			}
		}
		for (i = 0 ; i < B_k->row ; ++i) {
			b_there[i+1] += b_there[i] ;
		}

#ifdef _OPENMP
#pragma omp parallel for
#endif
		for (i = 0 ; i < B_k->row ; ++i) {
			index_t there = b_there[i] ;
			index_t jz ;
			dimen_t last_j = (dimen_t) -1 ;
			for (jz = B_k->start[i] ; jz < B_k->start[i+1] ; ++jz) {
				dimen_t j = B_k->colid[jz] ;
				if (j/UNRL  == last_j) {
					Bd->data[UNRL*there+j%UNRL] = B_k->data[jz] ;
				}
				else {
					last_j = j/UNRL ;
					++there ;
					Bd->colid[there] = last_j*UNRL ;
					Bd->data[UNRL*there+j%UNRL] = B_k->data[jz] ;
					Bd->start[i+1] += 1 ;
					assert(there<Bd->nnz);
				}
			}
#ifndef NDEBUG
			if (B_k->start[i+1]-B_k->start[i]) { assert(Bd->start[i+1]) ; }
#endif

			/* Bd->start[i] += 1 ; */
			/* ++there ; */
		}

		index_t there = b_there[Bd->row] + 1 ;

		free(b_there);

		for ( i = 0 ; i < Bd->row ; ++i) {
			/* assert(Bd->start[i+1]); */
			Bd->start[i+1] += Bd->start[i];
		}
		assert(Bd->start[Bd->row] == there);
		SAFE_REALLOC(Bd->data,(UNRL*there),elemt_t);
		SAFE_REALLOC(Bd->colid,there,dimen_t);
	}
	return;

}
#endif


#endif /* __GB_matrix_H */
/* vim: set ft=c: */
