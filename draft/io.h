#ifndef __GB_io_H
#define __GB_io_H

#define ERROR_READING \
	if (verbose > 0) { \
		printf("Error while reading file '%s'\n",fn); \
		fclose(fh); \
		return -1 ; \
	}

#define SAFE_MALLOC(ptr,size,elt) \
	elt * ptr = (elt *) malloc(size*sizeof(elt)); \
	assert(ptr)

#define SAFE_CALLOC(ptr,size,elt) \
	elt * ptr = (elt *) calloc(size*sizeof(elt)); \
	assert(ptr)

#define SAFE_REALLOC(ptr,size,elt) \
	ptr = realloc(ptr,size*sizeof(elt)); \
	assert(ptr)

#define SAFE_READ_v(val,elt,file) \
	assert(fread(&(val),sizeof(elt),1,file)==1)

#define SAFE_READ_V(val,elt,file) \
	elt val ; \
	SAFE_READ_v(val,elt,file)

#define SAFE_READ_p(val,size,elt,file) \
	assert(fread(val,sizeof(elt),size,file)==size)

#define SAFE_READ_P(val,size,elt,file) \
	SAFE_MALLOC(val,size,elt); \
	SAFE_READ_p(val,size,elt,file)



int prepare_split(uint32_t * buffer
		, uint32_t buffer_rows
		, uint32_t * pivots
		, uint32_t * sparse
		, uint32_t * discard
		, uint32_t * discarded
		, uint32_t * start
		, uint32_t i
		)
{
	for (int j = 0 ; j < buffer_row_size ; ++j) {
		pivots[j] = buffer[start[i+j]] ; // first column in each row
		sparse[j] = start[i+j+1]-start[i]; // length of each row
		repet[pivots[j]] += 1 ;
	}
	for (int k = 0 ; k < buffer_row_size ; ++k) { // check that we are increasing
		if (pivots[k+1] < pivots[k]) exit(-2);
	}
	for (int k = 0 ; k < buffer_row_size ; ++k) {
		if (repet[k] == 0) {
			break;
		}
		if (repet[k] > 1) {
			uint32_t idx = k; // find the sparsest
			for (int l = 1 ; l < repet[k] ; ++l)
				if (sparse[l] < sparse[idx]) {
					idx = l ;
				}
			for (int l = 0 ; l < repet[k] ; ++l) {
				if (l != idx) {
					discard[discarded++] = l ;
				}
			}
		}


	}
}

void insert_row(CSR_zo * A, uint32_t * colid, uint32_t nnz, uint32_t pol)
{
	A->row += 1 ;
	A->map_zo_pol = (uint32_t*) realloc(A->map_zo_pol,A->row*sizeof(uint32_t));
	A->map_zo_pol[A->row-1] = pol ;
	uint32_t last = A->nnz();
	A->start_zo = (uint32_t *) realloc(A->start_zo,(A->row+1)*sizeof(uint32_t));
	A->start_zo[A->row] = last+nnz ;
	A->colid_zo = (uint32_t *) realloc(A->colid_zo, A->nnz() *sizeof(uint32_t));
	for (uint32_t i = 0 ; i < nnz ; ++i) {
		A->colid_zo[last+i] = colid[i] ;
	}
}

void fix_last_row(
		GBMatrix_t * A
		GBMatrix_t * C
		)
{
	CSR_zo * A_mat = A->matrix_zo[A->matrix_nb-1] ;
	CSR_zo * B_mat = C->matrix_zo[C->matrix_nb-1] ;
	if (B_mat->row == MAT_ROW_BLOCK) {
		B_mat->matrix_nb += 1 ;
		B_mat->matrix_zo = (CSR_zo  *) realloc(B_mat->matrix_zo, B_mat->matrix_nb * sizeof(CSR_zo));
		B_mat.init();
	}
	uint32_t last_row_size = A_mat->nnz() ;

	insert_row(B_mat,A_mat->getRow(A_mat->row),last_row_size, A_mat->map_zo_pol[A_mat->row]);

	A_mat->row -= 1 ;
}

int append_rows_in_buffer(
		GBMatrix_t * A
		, GBMatrix_t * C
		, int32_t  * buffer
		, int buffer_row_size
		, uint32_t * start_zo
		, uint32_t row0
		, uint32_t * discard
		, uint32_t * discarded
		, uint32_t * map_zo_pol
		)
{
	int i = 0 ;
	for (i = 0 ; i < buffer_size ; ++i)  {
		// until next non discarded
		populate C ;
		// is discarded
		populate A
	}
}

// matrix reader and row splitter
int read_file(GBMatrix * A_init, GBMatrix * B_init
		, GBpolynomials * polys
		, FILE * fh
		, int verbose = 0)
{
	uint32_t i = 0 ;

	uint32_t m, n , mod;
	int64_t nnz ;

	// sanity
	assert(fh);

	// READ in row col nnz mod
	SAFE_READ_V(m,uint32_t,fh);
	SAFE_READ_V(n,uint32_t,fh);
	SAFE_READ_V(mod,uint32_t,fh);
	assert(mod > 1);
	SAFE_READ_V(nnz,uint32_t,fh);


	SAFE_MALLOC(start_zo, m+1, uint32_t);


	// READ in ZO start
	SAFE_READ_V(start_zo,m,uint32_t,fh);

	// largest row
	uint32_t big_row = 0 ;
	for (i = 0 ; i < m ; ++i)
		big_row = max(big_row,start_zo[i]);

	// pol/zo correspondance
	SAFE_READ_V(map_zo_pol,m,uint32_t,fh);

	SAFE_MALLOC(map_pol_zo_A,m,uint32_t);
	uint32_t map_pol_zo_A_size = 0 ;
	SAFE_MALLOC(map_pol_zo_B,m,uint32_t);
	uint32_t map_pol_zo_B_size = 0 ;


	// colid in ZO
	SAFE_MALLOC(buffer,MAT_ROW_BLOCK*big_row,uint32_t);
	SAFE_READ_V(colid_size,uint32_t,fh);
	uint32_t last_start = 0 ;
	i = 0 ;
	// FIXME
	if (MAT_ROW_BLOCK <= m) {
		uint32_t * pivots = (uint32_t) malloc(MAT_ROW_BLOCK*sizeof(uint32_t));
		uint32_t * sparse = (uint32_t) malloc(MAT_ROW_BLOCK*sizeof(uint32_t));
		uint32_t * repet  = (uint32_t) calloc(MAT_ROW_BLOCK*sizeof(uint32_t));
		uint32_t * discard = (uint32_t *) malloc(MAT_ROW_BLOCK*sizeof(int32_t));
		uint32_t discarded = 0 ;

		// read in MAT_ROW_BLOCK rows
		for (; i < m ; i += MAT_ROW_BLOCK ) {
			uint32_t nnz_buffer = start[i+MAT_ROW_BLOCK+1]-start[i];
			SAFE_READ_v(buffer,nnz_buffer,uint32_t,fh);
			// split buffer
			prepare_split(buffer,MAT_ROW_BLOCK,pivots,sparse,discard,&discarded,start_zo,i);


			// for the first one, check that it is not in the previous batch
			if (i > 0) {
				if (A->last_pivot() == pivots[0]) {
					if (A->sparsity(A->row) > sparse[pivots[0]]) {
						fix_last_row(A,C);
					}
				}
			}
			// GROW A,C with seleted pivots
			append_rows_in_buffer(A,C,buffer,READ_MAT_ROW_BLOCK,start_zo,i,discard,discarded,map_pol_zo);
			// ADD POLYS

		}
	}


	{
		uint32_t nnz_buffer = start[i+MAT_ROW_BLOCK+1]-start[i];
		// same
		// split buffer
	}

	free(buffer);

	// size of nnz for pols:
	SAFE_READ_P(size_pol,uint64_t,fh);

	// create GBpolynomials shared by A_init and B_init
	SAFE_READ_V(start_pol,m,uint32_t,fh);
	// populate A and C

	SAFE_READ_V(vals_pol,size_pol,uint32_t,fh);
	// populate A and C

	return 0 ;
}

void permute_columns(
		GBMatrix_t * A
		// , GBMatrix_t * C
		, uint32_t *  perm
		)
{
	SAFE_MALLOC(col_size,A->col,uint32_t);
	SAFE_MALLOC(last_elt,A->row,uint32_t);
	// copy last columns of A,C to B,D

	for (uint32_t k = 0 ; k < A->matrix_nb ; ++k ) {
		CSR_zo * A_k = A->matrix_zo[k] ;

		for (uint32_t i = 0 ; i < A_k->row ; ++i) {
			for (uint32_t j = A_k->start_zo[i] ; j < A_k->start_zo[i+1] ; ++j) {
				col_size[A_k->colid_zo[j]] += 1 ;
				last_elt[A_k->colid_zo[j]] = i ;
			}
		}
	}
	for (uint32_t j = A->row ; j < A->col ; ++j) {
		uint32_t k = last_elt[j] ;
		if (perm[k] != 0 && (col_size[perm[k]] > col_size[j]));
		perm[k] = j ;
	}
}



int split_columns(
		GBMatrix_t * A
		, GBMatrix_t * B
		, GBMatrix_t * C
		, DenseMatrix_t * D
		)
{
	// search for columns to permute to make A sparser
	SAFE_MALLOC(perm,A->row,uint32_t);
	for (uint32_t i = 0 ; i < A->row ; ++i)
		perm[i] = i ;
	permute_columns(A,perm);
	// get col size of C
	uint32_t kC =0;
	for (uint32_t i = 0 ; i < A->row ; ++i) {
		if (perm[i] != i) {
			k1 ++ ;
		}
	}
	// copy last columns of A,C to B,D in proper format
}

#endif // __GB_io_H
