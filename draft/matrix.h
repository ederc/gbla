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


#ifndef __GBLA_matrix_H
#define __GBLA_matrix_H

#include "config.h"
#include "types.h"
#include "macros.h"
#include "printer.h"


typedef struct CSR_pol {
	dimen_t nb ;
	index_t  * start_pol ;
	elemt_t  * data_pol ;
} CSR_pol;

typedef struct CSR {
	dimen_t row;
	dimen_t col;
	index_t  nnz ; /* this is start[row] */
	index_t  * start ;
	dimen_t * colid ;
	dimen_t * map_zo_pol ; /* == NULL iff data != NULL */
	elemt_t   * data ;
} CSR;

typedef struct GBMatrix_t {

	dimen_t   row ;
	dimen_t   col ;
	index_t    nnz ;
	elemt_t     mod ;
	dimen_t   sub_row ;  /* there are sub_row x sub_col matrices sub */
	dimen_t   sub_col ;  /* there are sub_row x sub_col matrices sub */
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


/* INIT */

static void initSparseUnit(CSR * mat)
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

static void initSparse(GBMatrix_t * A)
{
	if (!A) return;
	A->row = 0 ;
	A->col = 0 ;
	A->nnz = 0 ;
	A->sub_row = 0 ;
	A->sub_col = 0 ;
	A->sub = NULL ;
}

static void initDenseUnit (DNS * A)
{
	A->row = 0 ;
	A->col = 0 ;
	A->ld  = 0 ;
	A->mod = 0 ;
	A->nnz = 0 ;
	A->ptr = NULL ;
}

/* GROW */


static void setRow(CSR * mat
		, dimen_t i
		, dimen_t * colid
		, index_t size
		, dimen_t pol
	   )
{
	dimen_t old = mat->start[i];
	dimen_t * mat_colid = mat->colid + old ;
	index_t j  ;
	for ( j = 0 ; j < size ; ++j) {
		mat_colid[j] = colid[j] ;
	}
	mat->start[i+1] = old + size ;
	mat->map_zo_pol[i] = pol;
}



/* PRINT */

static void printMatUnit(CSR * A)
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

static void printMat(GBMatrix_t * A)
{
	fprintf(stderr,"matrix %u x %u - %lu\n",A->row, A->col,A->nnz);
	fprintf(stderr,"mod ");
	Mjoin(print,elemt_t)(A->mod);
	fprintf(stderr,"\n");
	dimen_t i,j  ;
	fprintf(stderr," subdivision : %u x %u\n",A->sub_row , A->sub_col );
	for ( i=0 ; i < A->sub_row ; ++i ) {
		for ( j=0 ; j < A->sub_col ; ++j ) {
			printMatUnit(A->sub + i*A->sub_col+j);
		}
	}
}

static void printMatDense(DNS * A)
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

static void printPoly(CSR_pol * P)
{
	fprintf(stderr,"polys (%u)\n",P->nb);
	fprintf(stderr,"start\n");
	dimen_t i = 0 ;
	for ( ; i < P->nb+1 ; ++i) {
		fprintf(stderr,"%lu ", P->start_pol[i]);
	}
	fprintf(stderr,"\ndata:\n");
	i = 0 ;
	for ( ; i < P->nb ; ++i) {
		index_t j = P->start_pol[i] ;
		for ( ; j < P->start_pol[i+1] ; ++j) {
			Mjoin(print,elemt_t)(P->data_pol[j]);
			fprintf(stderr," ");
		}
		fprintf(stderr,"\n");
	}
}

/* STATS */

static index_t occupancySparse(GBMatrix_t * A)
{
	dimen_t k,j ;
	index_t  acc = 0;
	for ( k = 0 ; k < A->sub_row  ; ++k ) {
		for ( j = 0 ; j < A->sub_col  ; ++j ) {
			CSR * Ad = A->sub + k * A->sub_col + j;
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
	}
	acc *= UNRL ;

	return acc ;
}

static index_t occupancyDense ( DNS * D)
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

static void checkMatUnit(const CSR *Ak)
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

static void checkMat(const GBMatrix_t *A)
{
	if (A == NULL) exit(-1);
#ifndef NDEBUG
	dimen_t row = 0 ;
	index_t nnz = 0 ;
	dimen_t i,j  ;
	for ( i = 0 ; i < A->sub_row ; ++i) {
		for ( j = 0 ; j < A->sub_col ; ++j) {
			const CSR * Ak = A->sub + i * A->sub_col + j ;
			row += Ak->row;
			nnz += Ak->nnz;
			checkMatUnit(Ak);
		}
	}
	assert (nnz == A->nnz);
	assert (row == A->row);

#endif

}

/* FREE */

static void freeMatDense(DNS * A)
{
	free(A->ptr);
}

static void freeMatUnit(CSR * A)
{
	free(A->start);
	free(A->colid);
	free(A->data);
	free(A->map_zo_pol);
	/* free(A); */
}

static void freePol( CSR_pol * A)
{
	free (A->start_pol);
	free (A->data_pol);
}


static void freeMat(GBMatrix_t * A)
{
	dimen_t i;
	for (i=0 ; i < A->sub_row*A->sub_col ; ++i) {
		freeMatUnit(A->sub+i) ;
	}
	free(A->sub);
}

/* CONVERSION */

static void convert_CSR_2_DNS(DNS * D, const GBMatrix_t * S )
{
	D->row = S->row ;
	D->col = S->col ;
	D->ld = ALIGN(D->col) ;
	D->mod = S->mod ;
	SAFE_CALLOC(D->ptr,(index_t)D->row * (index_t)D->ld,elemt_t);

	dimen_t k,j ;
	for (k = 0 ; k < S->sub_row ; ++k) {
		for (j = 0 ; j < S->sub_col ; ++j) {
			CSR * S_k = S->sub+k*S->sub_col+j ;
			dimen_t i ;
			index_t i_off = k * MAT_ROW_BLK ;
			index_t j_off = j * MAT_COL_BLK ;
			elemt_t * D_off = D->ptr + i_off*D->ld+j_off ;
			for (i = 0 ; i < S_k->row ; ++i) {
				index_t jz ;
				for ( jz = S_k->start[i] ; jz < S_k->start[i+1] ; ++jz) {
					dimen_t j = S_k->colid[jz] ;
					D_off[i*D->ld + j] = S_k->data[jz] ;
				}
			}
		}
	}


}

static void convert_CSR_2_CSR_block(GBMatrix_t * B, const GBMatrix_t * S )
{
	B->row = S->row;
	B->col = S->col;
	B->nnz = S->nnz;
	B->mod = S->mod;
	B->sub_row = S->sub_row ;
	B->sub_col = S->sub_col ;
	SAFE_MALLOC(B->sub,B->sub_row*B->sub_col,CSR);

	dimen_t k,l ;
	for (k = 0 ; k < S->sub_row ; ++k) {
		for (l = 0 ; l < S->sub_col ; ++l) {
			const CSR * B_k = S->sub+k*S->sub_col+l ;
			CSR * Bd = B->sub + k *B->sub_col+l;
			Bd->row = B_k->row ;
			Bd->col = B_k->col ;
			Bd->nnz = B_k->nnz ;
			Bd->map_zo_pol = NULL ;

			SAFE_CALLOC(Bd->start,(Bd->row+1),index_t);
			SAFE_CALLOC(Bd->data,UNRL*ALIGN(Bd->nnz),elemt_t);
			SAFE_MALLOC(Bd->colid,ALIGN(Bd->nnz),dimen_t);

			SAFE_CALLOC_DECL(b_there,(Bd->row+1),index_t);
			b_there[0] = (index_t)-1 ;

			dimen_t i ;
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

			}

			index_t there = b_there[Bd->row] + 1 ;

			free(b_there);

			for ( i = 0 ; i < Bd->row ; ++i) {
				Bd->start[i+1] += Bd->start[i];
			}
			assert(Bd->start[Bd->row] == there);
			SAFE_REALLOC(Bd->data,(UNRL*there),elemt_t);
			SAFE_REALLOC(Bd->colid,there,dimen_t);
		}
	}
	return;

}

/* SUBMATRIX SPLIT */

static void createSubmatricesUnit(CSR * Ad, dimen_t col, const CSR * Ah)
{

	/* get sizes of submatrices */
	SAFE_CALLOC_DECL(sizes,col,dimen_t);
	dimen_t i ;
	for (i = 0 ; i < Ah->row ; ++i) {
		index_t jz ;
		for ( jz = Ah->start[i] ; jz < Ah->start[i+1] ; ++jz) {
			dimen_t macol = Ah->colid[jz];
			assert(macol/MAT_COL_BLK < col);
			sizes[macol/MAT_COL_BLK] += 1 ;
		}
	}

	/* allocate submatrices */
	dimen_t j ;
	for ( j = 0 ; j < col ; ++j) {
		CSR * Aj = Ad + j ;
		if (sizes[j] == 0)
			Aj = NULL ;
		else {
			Aj->row = Ad->row ;
			Aj->col = 0 ;
			Aj->nnz = sizes[j] ;
			SAFE_CALLOC(Aj->start,Aj->row+1,index_t); /* could be dimen_t */
			SAFE_MALLOC(Aj->colid,Aj->nnz,dimen_t);
			SAFE_MALLOC(Aj->data, Aj->nnz,elemt_t);
			Aj->map_zo_pol = NULL ;
		}

	}

	/* populate submatrices */
	SAFE_CALLOC_DECL(here,col,index_t);
	for (i = 0 ; i < Ah->row ; ++i) {
		index_t jz ;
		for ( jz = Ah->start[i] ; jz < Ah->start[i+1] ; ++jz) {
			dimen_t macol = Ah->colid[jz];
			dimen_t mamat = macol / MAT_COL_BLK ;
			dimen_t colid = macol % MAT_COL_BLK ;
			CSR * Aloc = Ad + mamat ;
			Aloc->start[i+1] += 1 ;
			Aloc->colid[here[mamat]] = colid ;
			Aloc->data [here[mamat]] = Ah->data[jz] ;
			here[mamat] += 1 ;
		}
	}

	/* fixing start */
	for ( j = 0 ; j < col ; ++j) {
		CSR * Aj = Ad + j ;
		if (Aj != NULL) {
			for (i = 0 ; i < Aj->row ; ++i) {
				Aj->start[i+1] += Aj->start[i] ;
			}
			/* could reduce row if start stagnates at end */
		}
	}

	/* fixing col */
	for ( j = 0 ; j < col ; ++j) {
		CSR * Aj = Ad + j ;
		if (Aj != NULL) {
			for (i = 0 ; i < Aj->row ; ++i) {
				/* if (Aj->start[i+1] > Aj->start[i]) { |+ non empty row +| */
				if (Aj->start[i+1] > 0) {
					Aj->col = max(Aj->col,Aj->colid[Aj->start[i+1]-1]);
				}
				/* } */
			}
			assert(Aj->col < MAT_COL_BLK);
		}
	}
	return ;
}

static void createSubmatrices(GBMatrix_t *A, const GBMatrix_t * AH)
{
	A->row = AH->row;
	A->col = AH->col;
	A->nnz = AH->nnz;
	A->mod = AH->mod;

	A->sub_row = AH->sub_row;
	A->sub_col = DIVIDE_INTO(A->col,MAT_COL_BLK);

	SAFE_MALLOC(A->sub,A->sub_row*A->sub_col,CSR);

	/* if (A->sub_col == 1) return ; */

	dimen_t i ;
	for (i = 0 ; i < A->sub_row ; ++i) {
		const CSR * Ah = AH->sub + i ;
		CSR * Ak =  A->sub + i * A->sub_col ;
		createSubmatricesUnit( Ak, A->sub_col, Ah) ;
	}
	return ;
}

#endif /* __GBLA_matrix_H */
/* vim: set ft=c: */
