#ifndef __GB_io_H
#define __GB_io_H

#define ERROR_READING \
{ \ 	if (verbose > 0) \
		printf("Error while reading file '%s'\n",fn); \
	fclose(fh); \
	return -1 ; \
}


#define FREE_MATRIX(m) \
	free(m)


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
	uint32_t * map_pol_zo = (uint32_t *) malloc((m+1)*sizeof(uint32_t));
	if (fread(&map_pol_zo, sizeof(uint32_t),m,fh) != 1) {
		free(start_zo);
		free(map_pol_zo);
		ERROR_READING ;
	}


	// colid in ZO
	// uint64_t header_seek_offset = 3*sizeof(uint32_t) + sizeof(uint64_t);
	// uint64_t start_zo_seek_offset = (m+1)*sizeof(uint32_t) ;
	// uint64_t map_pol_zo_seek_offset = (m)*sizeof(uint32_t) ;
	// uint64_t colid_zo_seek_offset = (nnz)*sizeof(int32_t) ;

	int32_t * buffer = (int32_t *) malloc(MAT_ROW_BLOCK*big_row*sizeof(int32_t));
	uint32_t last_start = 0 ;
	i = 0 ;
	if (MAT_ROW_BLOCK <= m) {
		for (; i < m ; i += MAT_ROW_BLOCK ) {
			if (fread(&buffer, sizeof(int32_t),start[i+MAT_ROW_BLOCK+1]-start[i]) != 1) {
				free(start_zo);
				free(map_pol_zo);
				free(buffer);
				ERROR_READING ;
			}
			// split buffer
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


	while( doing_zo) {
		// read in MAT_ROW_BLOCK rows
		tmp_rows[MAT_ROW_BLOCK] ;
		while (reading) {
			in >> tmp_rows[i] ;
			// get the first nz in this row (colid) and density (start[i+1]-start[i])
			get_first_nz(tmp_rows[i],map_first,map_sparse);
		}
		// select the rows in map_first/map_sparse
		while (update_a_init) {
			if (current_line >= MAT_ROW_BLOCK) {
				matrix_nb ++ ;
				realloc(positions,matrix_nb);
				realloc(map_zo_pol,matrix_nb);
				// create new positions
			}
			// update current position[matrix_nb]
			copy(A_init.positions[matrix_nb].colid_zo,tmp_rows[j]);
			// update the start_zo in the current matrix
			update(A_init.positions[matrix_nb].start_zo,start_zo,j);
			// update map_zo_pol[matrix_nb] for this one
			update(A_init.map_zo_pol[matrix_nb],map_zo_pol,j);
		}
		// update_b_init the same
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
	// populate

	uint32_t * vals_pol =  (uint32_t*) malloc(start_pol[m+1]*sizeof(uint32_t));
	if (fread(&vals_pol, sizeof(int32_t),m) != 1) {
		free(vals_pol);
		// free(start_zo);
		// free(map_pol_zo);
		ERROR_READING;
	}
	// populate

	return 0 ;
}

#endif // __GB_io_H
