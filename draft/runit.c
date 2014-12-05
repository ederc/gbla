#include <stdint.h>


#include "io.h"


int main(int ac, char **av) {

	if (ac != 2) {
		fprintf(stderr,"usage %s nom_fichier\n",av[0]);
		return -1;
	}


	FILE * fh = fopen(av[1],"r") ;
	SAFE_MALLOC_DECL(AA,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CC,1,GBMatrix_t);
	init(AA);
	init(CC);
	SAFE_MALLOC_DECL(POL,1,CSR_pol);

	/* fprintf(stderr,"READ FILE \n"); */
	read_file(AA,CC,POL,fh);

	/* fprintf(stderr,">>>>**************\n"); */
	/* fprintf(stderr,"FIRST SPLIT\n"); */

	/* printMat(AA); */

	/* fprintf(stderr,"--------------\n"); */

	/* printMat(CC); */

	/* fprintf(stderr,"--------------\n"); */

	/* printPoly(POL); */
	/* fprintf(stderr,"<<<<<**************\n"); */


	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	/* SAFE_MALLOC_DECL(Bt,1,GBMatrix_t); */
	SAFE_MALLOC_DECL(B,1,DenseMatrix_t);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DenseMatrix_t);

	init(A);
	/* init(Bt); */
	init(C);
	/* init(D); */

	/* fprintf(stderr,"SPLIT COLUMNS\n"); */
	split_columns(AA,CC,POL,A,B,C,D);

	/* fprintf(stderr,">>>>>**************\n"); */
	/* fprintf(stderr,"SECOND SPLIT\n"); */
	/* printMat(A); */

	/* fprintf(stderr,"--------------\n"); */

	/* printMat(Bt); */
	/* printMatDense(B); */

	/* fprintf(stderr,"--------------\n"); */

	/* printMat(C); */

	/* fprintf(stderr,"--------------\n"); */

	/* printMatDense(D); */

	/* fprintf(stderr,"--------------\n"); */

	/* printPoly(POL); */
	/* fprintf(stderr,"<<<<**************\n"); */

	/* fprintf(stderr,"SOLVING\n"); */
	reduce(A,B,C,D);
	/* fprintf(stderr,">>>>>**************\n"); */
	/* printMatDense(B); */
	/* fprintf(stderr,"--------------\n"); */
	/* printMatDense(D); */
	/* fprintf(stderr,"<<<<**************\n"); */
	return 0;
}
