


struct mat_file {

	int row, col, nnz ; // 1
	int mod ;           // 1
	int * start_zo ;    // 2
	int * map_zo_pol ;  // 2
	int ** colid_zo ;    // 3
	int * start_pol;    // 4
	short * vals_pol;   // 4
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
 * 
 * Splitting the columns
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
 *
 * Operation C(A^-1 B) - D
 *
 * - TRSM using ATL_dreftrsmLUTN.c
 * - block by block (of full rows of B and full columns of C), convert to dense, do dense fgemm
 * - need the coeffients of Ainv B (as a dense matrix) and of C (as a dense matrix)
 * - convert C (CSR) to dense: use an extra start that remembers where we left at the previous conversion
 *
 *
 * Permution of columns :
 * - for each column in A' with same last non zero entry select the sparsest.
 * - scanning for sparsity of column with last non zero element is easy in CSC
 * - looking for the sparsity of each column in CSR can be done by recording the first non visited column in a start vector
 * - need to insert/remove in a vector : possibly easier with ** pointers
 * - removing columns in a CSR can be done starting from the end towards the beginning and memove
 * - convert from CSR to CSC just one column is done by insertions/remembering the places.
 * - do the column permutations in 2 steps : scanning then swapping. may be parallelisable.
 */
