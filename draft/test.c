#include "io.h"

int main() {


	FILE * fh = fopen("test.gb","r") ;
	SAFE_MALLOC_DECL(AA,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CC,1,GBMatrix_t);
	init(AA);
	init(CC);
	SAFE_MALLOC_DECL(POL,1,CSR_pol);

	read_file(AA,CC,POL,fh);

	printMat(AA);

	fprintf(stderr,"--------------\n");

	printMat(CC);

	fprintf(stderr,"--------------\n");

	printPoly(POL);

	return 0;
}
