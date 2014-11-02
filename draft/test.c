#include "io.h"

int main() {


	FILE * fh = fopen("test_new.gb","r") ;
	SAFE_MALLOC_DECL(AA,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CC,1,GBMatrix_t);
	init(AA);
	init(CC);
	SAFE_MALLOC_DECL(POL,1,CSR_pol);

	fprintf(stderr,"READ FILE \n");
	read_file(AA,CC,POL,fh);

	fprintf(stderr,">>>>**************\n");
	fprintf(stderr,"FIRST SPLIT\n");

	printMat(AA);

	fprintf(stderr,"--------------\n");

	printMat(CC);

	fprintf(stderr,"--------------\n");

	printPoly(POL);
	fprintf(stderr,"<<<<<**************\n");


	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	SAFE_MALLOC_DECL(Bt,1,GBMatrix_t);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DenseMatrix_t);

	init(A);
	init(Bt);
	init(C);
	/* init(D); */

	fprintf(stderr,"SPLIT COLUMNS\n");
	split_columns(AA,CC,POL,A,Bt,C,D);

	fprintf(stderr,">>>>>**************\n");
	fprintf(stderr,"SECOND SPLIT\n");
	printMat(A);

	fprintf(stderr,"--------------\n");

	printMat(Bt);

	fprintf(stderr,"--------------\n");

	printPoly(POL);
	fprintf(stderr,"<<<<**************\n");

	return 0;
}
