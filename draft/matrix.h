#ifndef __GB_matrix_H
#define __GB_matrix_H


struct matrixFile_t {

	index_t row ;
	index_t col ;
	index_t nnz ;           // 1
	Element mod ;           // 1
	index_t * start_zo ;    // 2
	index_t * map_zo_pol ;  // 2
	index_t * colid_zo ;    // 3
	index_t * start_pol;    // 4
	Element * vals_pol;   // 4
} ;

struct CSR_pol {
	index_t row ;
	index_t col ;
	index_t * start_pol ;
	Element * vals_pol ;
} ;

struct CSR_zo {
	index_t row;
	// index_t col; // useless
	index_t nnz;
	index_t * start_zo ; // size row+1
	index_t * colid_zo ; // size nnz
	index_t * map_zo_pol ; // size row
} ;


struct GBMatrix_t {

	index_t block   ; // size of blocks
	index_t row ;
	index_t col ;
	index_t nnz ;
	index_t mod ;
	index_t * matrix_nb ; // nb of 0/1 matrices
	CSR_zo  * matrix_zo ; // 0/1 matrices reprensenting positions
} ;

struct GBpolynomials {
	index_t   poly_nb ;  // size of polynomial values
	CSR_POL * poly_val ; // polynomials value list
} ;

#endif // __GB_matrix_H
