#include "io.h"

int main() {


	FILE * fh = fopen("test.gb","r") ;
	SAFE_MALLOC_DECL(AA,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CC,1,GBMatrix_t);
	init(AA);
	init(CC);
	SAFE_MALLOC_DECL(POL,1,CSR_pol);

	read_file(AA,CC,POL,fh);

	fprintf(stderr,"FIRST SPLIT\n");

	printMat(AA);

	fprintf(stderr,"--------------\n");

	printMat(CC);

	fprintf(stderr,"--------------\n");

	printPoly(POL);


	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	SAFE_MALLOC_DECL(Bt,1,GBMatrix_t);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DenseMatrix_t);

	init(A);
	init(Bt);
	init(C);
	/* init(D); */

	split_columns(AA,CC,A,Bt,C,D);

	fprintf(stderr,"SECOND SPLIT\n");
	printMat(A);

	fprintf(stderr,"--------------\n");

	printMat(Bt);

	fprintf(stderr,"--------------\n");

	printPoly(POL);

	return 0;
}
