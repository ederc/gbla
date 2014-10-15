#ifndef __GB_io_H
#define __GB_io_H

#define ERROR_READING \
{ \ 	if (verbose > 0) \
		printf("Error while reading file '%s'\n",fn); \
	fclose(fh); \
	return -1 ; \
}

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

void fix_last_row(
		GBMatrix_t * A
		GBMatrix_t * B
		)
{
	// in A reduce  row, reduce nnz
	// move last to last
	// update start
	// update zo_map
	// in B augment row, add nnz
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
	if (fh == NULL) {
		ERROR_READING ;
	}

	// READ in row col nnz mod
	if (fread(&m, sizeof(uint32_t), 1, fh) != 1) {
		ERROR_READING ;
	}

	if (fread(&n, sizeof(uint32_t), 1, fh) != 1) {
		ERROR_READING ;
	}

	if (fread(&nnz, sizeof(uint64_t), 1, fh) != 1) {
		ERROR_READING ;
	}

	if ((fread(&mod, sizeof(uint32_t), 1, fh) != 1) || (mod ==1)) {
		ERROR_READING ;
	}

	uint32_t * start_zo = (int32_t *) malloc((m+1)*sizeof(uint32_t));
	assert(start_zo);


	// READ in ZO start
	if (fread(&start_zo, sizeof(uint32_t),m+1,fh) != 1) {
		free(start_zo);
		ERROR_READING ;
	}

	// largest row
	int32_t big_row = 0 ;
	for (i = 0 ; i < m ; ++i)
		big_row = max(start_zo[i+1]-start_zo[i]);


	// pol/zo correspondance
	uint32_t * map_pol_zo = (uint32_t *) malloc((m)*sizeof(uint32_t));
	if (fread(&map_pol_zo, sizeof(uint32_t),m,fh) != 1) {
		free(start_zo);
		free(map_pol_zo);
		ERROR_READING ;
	}

	uint32_t * map_pol_zo_A = (uint32_t *) malloc((m)*sizeof(uint32_t));
	uint32_t map_pol_zo_A_size = 0 ;
	uint32_t * map_pol_zo_B = (uint32_t *) malloc((m)*sizeof(uint32_t));
	uint32_t map_pol_zo_B_size = 0 ;


	// colid in ZO
	// uint64_t header_seek_offset = 3*sizeof(uint32_t) + sizeof(uint64_t);
	// uint64_t start_zo_seek_offset = (m+1)*sizeof(uint32_t) ;
	// uint64_t map_pol_zo_seek_offset = (m)*sizeof(uint32_t) ;
	// uint64_t colid_zo_seek_offset = (nnz)*sizeof(int32_t) ;

	int32_t * buffer = (int32_t *) malloc(MAT_ROW_BLOCK*big_row*sizeof(int32_t));
	uint32_t last_start = 0 ;
	i = 0 ;
	if (MAT_ROW_BLOCK <= m) {
		uint32_t * pivots = (uint32_t) malloc(MAT_ROW_BLOCK*sizeof(uint32_t));
		uint32_t * sparse = (uint32_t) malloc(MAT_ROW_BLOCK*sizeof(uint32_t));
		uint32_t * repet  = (uint32_t) calloc(MAT_ROW_BLOCK*sizeof(uint32_t));
		uint32_t * discard = (uint32_t *) malloc(MAT_ROW_BLOCK*sizeof(int32_t));
		uint32_t discarded = 0 ;

		// read in MAT_ROW_BLOCK rows
		for (; i < m ; i += MAT_ROW_BLOCK ) {
			if (fread(&buffer, sizeof(int32_t),start[i+MAT_ROW_BLOCK+1]-start[i]) != 1) {
				free(start_zo);
				free(map_pol_zo);
				free(buffer);
				ERROR_READING ;
			}
			// split buffer
			prepare_split(buffer,MAT_ROW_BLOCK,pivots,sparse,discard,&discarded,start_zo,i);


			// for the first one, check that it is not in the previous batch
			if (i > 0) {
				if (A->last_pivot() == pivots[0]) {
					if (A->sparsity(A->row) > sparse[pivots[0]])
						fix_last_row(A,B);
				}
			}
			// GROW A with seleted pivots
			append_rows_A();
			// GROW B with discarded
			append_rows_B();
			// ADD POLYS
			map_pol_zo_A()  ;
			map_pol_zo_B() ;

		}
	}
	int append_rows_A(
			GBMatrix_t * A
			, GBMatrix_t * B
			, int32_t  * buffer
			, int buffer_row_size
			, uint32_t * start_zo
			, uint32_t row0
			, uint32_t * discard
			, uint32_t * discarded
			)
	{
		int i = 0 ;
		for (i = 0 ; i < buffer_size ; ++i)  {
			// until next non discarded
			populate B ;
			// is discarded
			populate A
		}
	}

	{
		if (fread(&buffer, sizeof(int32_t),start[m+1]-start[i]) != 1) {
			free(start_zo);
			free(map_pol_zo);
			free(buffer);
			ERROR_READING ;
		}
		// split buffer
	}

	free(buffer);

	// create GBpolynomials shared by A_init and B_init
	uint32_t * start_pol =  (uint32_t*) malloc((m+1)*sizeof(uint32_t));
	if (fread(&start_pol, sizeof(int32_t),(m+1)) != 1) {
		free(start_pol);
		// free(start_zo);
		// free(map_pol_zo);
		ERROR_READING;
	}
	// populate A and B

	uint32_t * vals_pol =  (uint32_t*) malloc(start_pol[m+1]*sizeof(uint32_t));
	if (fread(&vals_pol, sizeof(int32_t),m) != 1) {
		free(vals_pol);
		// free(start_zo);
		// free(map_pol_zo);
		ERROR_READING;
	}
	// populate A and B

	return 0 ;
}

#endif // __GB_io_H
