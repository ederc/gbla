#include <stdint.h>


#include "io.h"
#include <sys/time.h>


int main(int ac, char **av) {

	if (ac != 2) {
		fprintf(stderr,"usage %s nom_fichier\n",av[0]);
		return -1;
	}

	struct timeval start,end ;
	gettimeofday(&start,NULL);

	FILE * fh = fopen(av[1],"r") ;
	SAFE_MALLOC_DECL(AA,1,GBMatrix_t);
	SAFE_MALLOC_DECL(CC,1,GBMatrix_t);
	init(AA);
	init(CC);
	SAFE_MALLOC_DECL(POL,1,CSR_pol);

	/* fprintf(stderr,"READ FILE \n"); */
	read_file(AA,CC,POL,fh);

	gettimeofday(&end,NULL);

	fprintf(stderr," LOAD time : %.3f\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec) - start.tv_usec)/1e6);

	/* fprintf(stderr,">>>>**************\n"); */
	/* fprintf(stderr,"FIRST SPLIT\n"); */

	/* printMat(AA); */

	/* fprintf(stderr,"--------------\n"); */

	/* printMat(CC); */

	/* fprintf(stderr,"--------------\n"); */

	/* printPoly(POL); */
	/* fprintf(stderr,"<<<<<**************\n"); */

	gettimeofday(&start,NULL);

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

	gettimeofday(&end,NULL);

	fprintf(stderr," SPLIT 4 time : %.3f\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec) - start.tv_usec)/1e6);

	gettimeofday(&start,NULL);
	/* fprintf(stderr,"SOLVING\n"); */
	reduce(A,B,C,D);
	/* fprintf(stderr,">>>>>**************\n"); */
	/* printMatDense(B); */
	/* fprintf(stderr,"--------------\n"); */
	/* printMatDense(D); */
	/* fprintf(stderr,"<<<<**************\n"); */

	gettimeofday(&end,NULL);

	fprintf(stderr," REDUCE  time : %.3f\n", (double)((end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec) - start.tv_usec)/1e6);

	return 0;
}
