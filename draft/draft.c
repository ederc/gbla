#define MAT_ROW_BLOCK 128      // read matrix MAT_ROW_BLOCK by MAT_ROW_BLOCK
#define Element       int32_t  // Element representation mod p
#define index_t       int32_t  // indexing elements

struct matrixFile {

	index_t row, col, nnz ; // 1
	Element mod ;           // 1
	index_t * start_zo ;    // 2
	index_t * map_zo_pol ;  // 2
	index_t * colid_zo ;    // 3
	index_t * start_pol;    // 4
	Element * vals_pol;   // 4
} ;

struct CSR_POL {
	int row, int col ;
	index_t * start_pol ;
	Element * vals_pol ;
} ;

struct CSR_ZO {
	index_t row; index_t col;
	index_t * start_zo ;
	index_t * colid_zo ;
} ;


struct GBMatrix {
	
	index_t block_size = MAT_ROW_BLOCK  ;
	index_t row, col, nnz ; 
	index_t mod ;           
	index_t * matrix_nb ; // nb of 0/1 matrices
	CSR_ZO  * positions ; // 0/1 matrices reprensenting positions
	index_t ** map_zo_pol ; // maps rows in positions to polynomials
} ;

struct GBpolynomials {
	CSR_POL * polynomials ; // polynomials value list
} ;


/*
 * Splitting the rows
 *
 *
 * A'
 * --
 * C'
 * 
 * - read 256 x 256 and construct on the fly
 * A' and C' have mat_file format
 * - select those rows with sub diagnal pivots and put them in C'
 * - proceed to next 256 ; allocate if necessary
 * - create map_zo_pol for C' as soon as a row is put in C', create the one in A' too
 * - vals_pol, start_pol are shared by A' and B'
 *
 */



// matrix reader and row splitter
int read_file(GBMatrix * A_init, GBMatrix * B_init
		, GBpolynomials * polys
		, FILE * in) 
{
	// first lines
	in >> row >> col >> nnz ;
	in >> mod ;
	// start in ZO
	in >> start_zo ;
	// pol/zo correspondance
	in >> map_zo_pol ;
	// colid in ZO
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
	// create GBpolynomials shared by A_init and B_init
	in >> start_pol ;
	in >> vals_pol ;

	return status ;
}


 /* Splitting the columns
 * 
 * A | B
 * -----
 * C | D
 *
 *
 * A is now upper triangular
 * B is stored column wise
 * D is dense
 * C is stored row wise
 * A, C are stored sparse.
 *
 *  Permution of columns :
 * - for each column in A' with same last non zero entry select the sparsest.
 * - scanning for sparsity of column with last non zero element is easy in CSC
 * - looking for the sparsity of each column in CSR can be done by recording the first non visited column in a start vector
 * - need to insert/remove in a vector : possibly easier with ** pointers
 * - removing columns in a CSR can be done starting from the end towards the beginning and memove
 * - convert from CSR to CSC just one column is done by insertions/remembering the places.
 * - do the column permutations in 2 steps : scanning then swapping. may be parallelisable.

 */

// now splitting columns
int split_columns(GBMatrix * A_init, GBMatrix * B_init
		, GBMatrix * C, DenseMatrix * D
		, GBpolynomials * polys
		) 
{
	index_t k =  A_init.row ;
	// copy n-k last columns from A_init to  C, in col major order
	// A_init is shrunk to save space
	splitA(A_init,C);
	// copy n-k last columns from B_init to  D, in row major order, and D is dense.
	// C_init is shrunk to save space
	splitC(C_init,D);
	// prepare for swapping
	index_t sparsest_cols[k]; // sparsest_cols[i] is the col id of sparsest column with last 1 at i

	// swap columns in  A_init and C
	swap_cols(A_init,C,sparsest_cols);
	// swap columns in  B_init and D
	swap_cols(B_init,D,sparsest_cols);
	// update polys
	swap_cols(polys, sparsest_cols);
}

 /* Operation C(A^-1 B) - D
 *
 * - TRSM using ATL_dreftrsmLUTN.c
 * - block by block (of full rows of B and full columns of C), convert to dense, do dense fgemm
 * - need the coeffients of Ainv B (as a dense matrix) and of C (as a dense matrix)
 * - convert C (CSR) to dense: use an extra start that remembers where we left at the previous conversion
 *
 *
 */

int gauss_elim(GBMatrix * A, GBMatrix * Bt, GBMatrix *C, DenseMatrix *D)
{
	index_t k1 = B.cols  ;
	index_t d  = A.rows ;
	index_t k2 = C.rows ;
	while (blocks_of_B) {
		DenseMatrix B_tmp(B,d,k1);
		// perform A^-1 B_tmp  (gauss elim)
		do_elim(A,B_tmp);
		DenseMatrix C_tmp(C,k2,d);
		// perform dense D -= C_tmp B_tmp  ;
		fgemm(C_tmp,B_tmp,1,D,1);
	}
}
