

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
