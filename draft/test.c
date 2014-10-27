#include "io.h"

int main() {


	FILE * fh = fopen("test.gb","r") ;
	SAFE_MALLOC_DECL(AA,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CC,1,GBMatrix_t);

	read_file(AA,CC,fh);

	printMat(AA);

	fprintf(stderr,"--------------");
	printMat(CC);

	return 0;
}
