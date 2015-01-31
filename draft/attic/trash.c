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


void do_permute_columns_up(
		uint32_t k,
		GBMatrix_t * A_init,
		GBMatrix_t * A,
		DenseMatrix_t * B
		/* GBMatrix_t * Bt */
		/* , uint32_t * perm_i */
		/* , uint32_t * perm_j */
		/* , uint32_t perm_size */
		, CSR_pol * polys
		)
{z
	/* uint32_t k = A_init->row; */

	A->row = B->row = k ;
	A->mod = B->mod = A_init->mod;
	A->col = k ;
	B->col = A_init->col - k ;

	SAFE_CALLOC(B->data,B->row*B->col,elem_t);

	uint32_t blk = 0;
	for ( ; blk < A_init->matrix_nb ; ++blk) {
		/* uint32_t i = 0 ; */

		CSR_zo * A_k = &(A_init->matrix_zo[blk]); /* A_k band of A_init */

		SAFE_MALLOC_DECL(start_b,A_k->row,uint32_t);

		/*  permute columns and create A  */

		permuteCSR(A_k,A,start_b,k
				/* ,perm_i,perm_j,perm_size */
				,polys);


		/* get B */

		uint32_t row_offset = blk*MAT_ROW_BLK ;
		elem_t * data = B->data + (row_offset * B->col) ;

		permuteDNS(A_k,A,start_b,data);

		A->nnz += A_k->nnz ;
		/* Bt->nnz += Bd->nnz ; */
	}

	/* fprintf(stderr,"A:"); */
	/* printMat(A); */

	/* fprintf(stderr,"Bt:"); */
	/* printMatDense(B); */
	/* printMat(Bt); */
}


void do_permute_columns_lo(
		uint32_t k,
		GBMatrix_t * C_init,
		GBMatrix_t * C,
		DenseMatrix_t * D
		/* , uint32_t * perm_i */
		/* , uint32_t * perm_j */
		/* , uint32_t perm_size */
		, CSR_pol * polys
		)
{
	C->row = D->row = C_init->row ;
	C->mod = D->mod = C_init->mod;
	C->col = k ;
	D->col = C_init->col - k ;

	SAFE_CALLOC(D->data,D->row*D->col,elem_t);

	uint32_t blk = 0;
	for ( ; blk < C_init->matrix_nb ; ++blk) {
		/* uint32_t i = 0 ; */

		CSR_zo * C_k = &(C_init->matrix_zo[blk]); /* A_k band of A_init */
		SAFE_MALLOC_DECL(start_b,C_k->row,uint32_t);

		/* C */
		permuteCSR(C_k,C,start_b,k
				/* ,perm_i,perm_j,perm_size */
				,polys);

		/*  D */
		uint32_t row_offset = blk*MAT_ROW_BLK ;
		elem_t * data = D->data + (row_offset * D->col) ;

		permuteDNS(C_k,C,start_b,data);

		C->nnz += C_k->nnz ;
	}

	/* fprintf(stderr,"C:"); */
	/* printMat(C); */

	/* fprintf(stderr,"D:"); */
	/* printMatDense(D); */


}


#if 0
/* transforms perm where perm[i] = j into two lists such that ther permutations
 * are (perm_i[k],perm_j[k]) for perm and perm^(-1).  perm_i is increasing.
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
#endif

void split_columns(
		GBMatrix_t * A_init
		, GBMatrix_t * C_init
		, CSR_pol * polys
		, GBMatrix_t * A
		, DenseMatrix_t * B
		/* , GBMatrix_t * Bt */
		, GBMatrix_t * C
		, DenseMatrix_t * D
		)
{
	/* search for columns to permute to make A sparser */
	/* SAFE_MALLOC_DECL(perm,A->row,uint32_t); |+ what columns the first row ones should be exchanged with +| */
	/* uint32_t i = 0; */
	/* for ( ; i < A_init->row ; ++i) */
		/* perm[i] = i ; */
	/* uint32_t trans = get_permute_columns(A_init,perm); */

	/* fprintf(stderr,"--------------\n"); */
	/* i =  0 ; for( ; i < A_init->row ; ++i) fprintf(stderr, "%u ", perm[i]); fprintf(stderr,"\n"); i = 0; */
	/* fprintf(stderr,"--------------\n"); */


	/* SAFE_MALLOC_DECL(perm_i,2*trans,uint32_t); */
	/* SAFE_MALLOC_DECL(perm_j,2*trans,uint32_t); */
	/* split_permutations(perm,A_init->row,perm_i,perm_j); */
	/* i =  0 ; for( ; i < 2*trans ; ++i) fprintf(stderr, "%u ", perm_i[i]); fprintf(stderr,"\n"); i = 0; */
	/* i =  0 ; for( ; i < 2*trans ; ++i) fprintf(stderr, "%u ", perm_j[i]); fprintf(stderr,"\n"); i = 0; */

	/* copy last columns of A,C to B,D in proper format */

	do_permute_columns_up(A_init->row,A_init,A,B
			/* ,perm_i,perm_j,2*trans */
			,polys);
	do_permute_columns_lo(A_init->row,C_init,C,D
			/* ,perm_i,perm_j ,2*trans*/
			,polys);

	checkMat(A);
	checkMat(C);
	fprintf(stderr," A is %u x %u (sparsity : %f)\n",A->row,A->col,(double)(A->nnz)/(double)(A->row*A->col));
	fprintf(stderr," B is %u x %u (sparsity : %f)\n",B->row,B->col,(double)(A_init->nnz-A->nnz)/(double)(B->row*B->col));
	fprintf(stderr," C is %u x %u (sparsity : %f)\n",C->row,C->col,(double)(C->nnz)/(double)(C->row*C->col));
	fprintf(stderr," D is %u x %u (sparsity : %f)\n",D->row,D->col,(double)(C_init->nnz-C->nnz)/(double)(D->row*D->col));

	return ;
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
					last_elt[A_k->colid_zo[j]-A->row] = (A->matrix_nb-1)*MAT_ROW_BLK+i ; /* zero based, but 0 is an index if column not empty :-) */
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
		/* , uint32_t * perm_i */
		/* , uint32_t * perm_j */
		/* , uint32_t perm_size */
		, CSR_pol * polys
	       )
{
	/* apply permutations (the rows keep the same length : in place) */
	uint32_t i = 0 ;
#if 0
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
#endif
	assert(A_k->start_zo[A_k->row] == A_k->nnz);
	SAFE_MALLOC(A_k->data,A_k->nnz,elem_t);
	/* sort rows for we have messed them up */
	i = 0;
	for ( ; i < A_k->row ; ++i) {
		uint32_t j0 = A_k->start_zo[i] ;
		uint32_t j1 = A_k->start_zo[i+1] ;
		uint32_t start_idx = polys->start_pol[ A_k->map_zo_pol[i] ] ;
		/* fprintf(stderr,"row %u, poly %u, start %u\n",i,A_k->map_zo_pol[i],start_idx); */
		MEMCPY(A_k->data+j0,polys->data_pol+start_idx,j1-j0);
		/* uint32_t jj = 0 ; for ( ; jj < j1-j0 ; ++jj) fprintf(stderr,"%u ",A_k->data[j0+jj]) ; fprintf(stderr,"\n"); */
		/* insert_sort(&A_k->colid_zo[j0], j1-j0) ; */
		/* insert_sort_duo_data(A_k->colid_zo+j0, j1-j0, A_k->data+j0) ; */
		/* jj = 0 ; for ( ; jj < j1-j0 ; ++jj) fprintf(stderr,"%u ",A_k->data[j0+jj]) ; fprintf(stderr,"\n"); */
	}

	/* fprintf(stderr,"permuted : "); */
	/* i =  0 ; for( ; i < A_k->nnz ; ++i) fprintf(stderr, "%u ", A_k->colid_zo[i]); fprintf(stderr,"\n"); i = 0; */

	/* where B starts in A_init */
	i = 0;
	for ( ; i < A_k->row ; ++i) {
		uint32_t j0 = A_k->start_zo[i] ;
		uint32_t j1 = A_k->start_zo[i+1] ;
		while(j0 < j1  && A_k->colid_zo[j0] < k)
			++j0;
		start_b[i] = j0 ;
	}
	/* fprintf(stderr,"col start: "); */
	/* i =  0 ; for( ; i < A_k->row ; ++i) fprintf(stderr, "%u ", start_b[i]); fprintf(stderr,"\n"); i = 0; */

	/* copy A */

	appendMatrix(A);
	CSR_zo * Ad = getLastMatrix(A);
	/* XXX BUG : MAT_ROW_BLK */
	SAFE_REALLOC(Ad->start_zo,A_k->row+1,uint64_t);
	Ad->start_zo[0] = 0 ;
	Ad->row = A_k->row ;
	Ad->col = k ;
	SAFE_MALLOC(Ad->map_zo_pol ,A_k->row,uint32_t);
	MEMCPY(Ad->map_zo_pol,A_k->map_zo_pol,A_k->row);

	i = 0;
	for ( ; i < A_k->row ; ++i) {
		Ad->nnz += start_b[i]-A_k->start_zo[i];
	}
	SAFE_MALLOC(Ad->colid_zo,Ad->nnz,uint32_t);
	SAFE_MALLOC(Ad->data,Ad->nnz,elem_t);

	i = 0;
	for ( ; i < A_k->row ; ++i) {
		uint32_t b = start_b[i]-A_k->start_zo[i];
		/* fprintf(stderr,"b : %u\n",b); */
		Ad->start_zo[i+1]= Ad->start_zo[i]+b;
		/* fprintf(stderr,"start %u+1=%lu\n",i,Ad->start_zo[i+1]); */
		MEMCPY(Ad->colid_zo+Ad->start_zo[i], A_k->colid_zo+A_k->start_zo[i], b);
		MEMCPY(Ad->data+Ad->start_zo[i], A_k->data+A_k->start_zo[i], b);
	}
	assert(Ad->start_zo[Ad->row] == Ad->nnz);


	/* fprintf(stderr,"Ad:"); */
	/* printMatUnit(Ad); */
}

#if 0

void permuteCSRT(void) {
	appendMatrix(Bt);
	CSR_zo * Bd = getLastMatrix(Bt);

	Bd->col = A_k->row ;
	Bd->row = A_init->col - A_init->row ;
	SAFE_CALLOC(Bd->start_zo,Bd->row+1,uint64_t);
	uint64_t other_nnz = getLastMatrix(A)->nnz ;
	Bd->nnz = A_k->nnz - other_nnz ;
	SAFE_MALLOC(Bd->colid_zo,Bd->nnz,uint32_t);
	SAFE_MALLOC(Bd->map_zo_pol ,A_k->row,uint32_t);
	MEMCPY(Bd->map_zo_pol,A_k->map_zo_pol,A_k->row); /* useless */

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
	SAFE_MALLOC(Bd->data, Bd->nnz, elem_t);

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
}
#endif

void permuteDNS(CSR_zo * C_k, GBMatrix_t * C, uint32_t * start_b, elem_t * data)
{

	/* uint64_t other_nnz = getLastMatrix(C)->nnz ; */

	/* C->nnz += other_nnz  ; */

	uint32_t k =  C->col;

	uint32_t i = 0 ;
	uint32_t ldd = C_k->col - k ;
	for ( ; i < C_k->row ; ++i) {
		uint32_t j = start_b[i] ;
		for (; j < C_k->start_zo[i+1] ; ++j) {
			uint32_t jj =  C_k->colid_zo[j]-k;
			data[ldd*i+jj] = C_k->data[j];
		}
	}
}


/* convert from dense */
#ifdef USE_B_SPARSE
	/* XXX convert B to sparse */
	SAFE_MALLOC_DECL(Bsp,1,GBMatrix_t);
	initSparse(Bsp);

	Bsp->mod = B->mod;
	Bsp->col = B->col;

#ifndef BLOCK_CSR
	SAFE_MALLOC_DECL(col_buf,B->col,taille_t);
#else
	SAFE_CALLOC_DECL(col_buf,DIVIDE_INTO(B->col,UNRL),taille_t);
#endif
	SAFE_CALLOC_DECL(dat_buf,B->ld,elem_t);

	for ( i = 0 ; i < B ->row ; ++i) {
		taille_t j ;
		taille_t length = 0 ;
#ifndef BLOCK_CSR
		for (j = 0 ; j < B->col ; ++j) {
			if (Bd[i*ldb+j] != 0) {
				col_buf[length]   = j ;
				dat_buf[length++] = Bd[(index_t)i*(index_t)ldb+(index_t)j];
			}
		}
		appendRowData(Bsp,col_buf,length,dat_buf);
#else
		taille_t last_j = (taille_t)-1 ;
		taille_t nnz = 0 ;
		taille_t start = 0 ;
		for (j = 0 ; j < B->col ; ++j) {
			if (Bd[i*ldb+j] != 0) {
				nnz ++ ;
				if (last_j == (taille_t)-1) { /* first chunk not started */
					col_buf[length] = j ;
					last_j = j ;
					dat_buf[start] = Bd[i*ldb+j];
				} else {
					if (last_j + UNRL > j) {
						dat_buf[start+j-last_j] =  Bd[i*ldb+j];
					}
					else {
						start += UNRL ;
						++ length ;
						col_buf[length] = j ;
						last_j = j ;
						dat_buf[start] = Bd[i*ldb+j];
					}
				}
			}
			assert(4*length == start);
		}
		appendRowData(Bsp,col_buf,length,dat_buf,nnz);
#endif
	}

	assert(Bsp->nnz == B->nnz);
	assert(Bsp->row == B->row);
	assert(Bsp->col == B->col);

	free(col_buf);
	free(dat_buf);

#ifndef NDEBUG
	checkMat(Bsp);
#endif
#else



#endif
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

void insert_sort_duo_data(uint32_t * liste, uint32_t  size, elem_t * copain)
{
	uint32_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && liste[d] < liste[d-1]) {
			SWAP(liste[d],liste[d-1]);
			SWAP(copain[d],copain[d-1]);
			d--;
		}
	}
}

void insert_sort_duo_data_rev(uint32_t * liste, uint32_t  size, elem_t * copain)
{
	uint32_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && liste[d] > liste[d-1]) {
			SWAP(liste[d],liste[d-1]);
			SWAP(copain[d],copain[d-1]);
			d--;
		}
	}
}


#ifdef BLOCK_CSR
void appendRowDataUnit_block(CSR * mat
		, dimen_t * colid
		, index_t size
		, elemt_t * data
		, index_t nnz
		)
{
	index_t new_colsize ;
	if (size > 0) {

		index_t old = mat->nnz ;
		index_t old_colsize = mat->start[mat->row]  ;
		new_colsize = old_colsize + size ;
		mat->nnz = old + nnz;
		/* XXX this may be slow */
		SAFE_REALLOC(mat->colid, new_colsize ,dimen_t);
		SAFE_REALLOC(mat->data , (UNRL * new_colsize) ,elemt_t);
		index_t i ;
		for ( i = 0 ; i < size ; ++i) {
			mat->colid[old_colsize+i] = colid[i] ;
		}
		for ( i = 0 ; i < UNRL*size ; ++i) {
			mat->data[UNRL*old_colsize+i] = data[i] ;
		}
	}

	mat->row ++ ;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->start,mat->row+1,index_t);
	mat->start[mat->row] = new_colsize ;
}
#endif

void appendRowData(GBMatrix_t * A
		, dimen_t * colid
		, index_t size
		, elemt_t * data
		)
{
	A->row += 1;
	A->nnz += size ;

	if ( ( A->sub_nb == 0) || ((A->sub[A->sub_nb-1]).row == MAT_ROW_BLK) ){
		appendMatrix(A);
	}

	appendRowDataUnit(&(A->sub[A->sub_nb-1]),colid,size,data);
}

void appendRowDataUnit(CSR * mat
		, dimen_t * colid
		, index_t size
		, elemt_t * data
		)
{
	if (size > 0) {

		index_t old = mat->start[mat->row] ;
		mat->nnz = old + size;
		/* XXX this may be slow */
		SAFE_REALLOC(mat->colid,mat->nnz,dimen_t);
		SAFE_REALLOC(mat->data ,mat->nnz,elemt_t);
		index_t i ;
		for ( i = 0 ; i < size ; ++i) {
			mat->colid[old+i] = colid[i] ;
		}
		for ( i = 0 ; i < size ; ++i) {
			mat->data[old+i] = data[i] ;
		}
	}

	mat->row ++ ;
	/* XXX this may be slow */
	SAFE_REALLOC(mat->start,mat->row+1,index_t);
	mat->start[mat->row] = mat->nnz  ;
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
		CSR * C_k = C->sub + k ;
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
				CSR * A_k =  A->sub + kk ;
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
					CSR * A_k = A->sub + kk ;
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
					CSR * A_k = A->sub + kk ;
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
				CSR * A_k = A->sub + kk ;
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
				CSR * A_k = A->sub + kk ;
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
				CSR * A_k = A->sub + kk ;
				dimen_t sz = (dimen_t)(A_k->start[jj+1]-A_k->start[jj]);
				assert(kk*MAT_ROW_BLK+jj == j);
				spaxpyn(diag,A_k->data+A_k->start[jj],
						sz,
						A_k->colid+A_k->start[jj]
						,temp_C,jump,nbnz,ldc);

				/* temp_D -= temp_C[j] * B[j] */

				CSR * B_k =  B->sub + kk ;
				sz = (dimen_t)(B_k->start[jj+1]-B_k->start[jj]) ;
				spaxpyn(diag,B_k->data+B_k->start[jj],
						sz,
						B_k->colid+B_k->start[jj],
						temp_D,jump,nbnz,ldd);
			}
}

