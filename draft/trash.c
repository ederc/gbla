
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

		uint32_t row_offset = blk*MAT_ROW_BLOCK ;
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
		uint32_t row_offset = blk*MAT_ROW_BLOCK ;
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
		MEMCPY(A_k->data+j0,polys->vals_pol+start_idx,j1-j0);
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
	/* XXX BUG : MAT_ROW_BLOCK */
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



