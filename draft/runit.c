#include <stdint.h>


#include "io.h"
#include <sys/time.h>


int main(int ac, char **av) {

	if (ac != 2) {
		fprintf(stderr,"usage %s nom_fichier\n",av[0]);
		return -1;
	}

	struct timeval start,end ;

	/* READ and SPLIT  */
	gettimeofday(&start,NULL);

	FILE * fh = fopen(av[1],"r") ;
	uint32_t * col_perm;

	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	SAFE_MALLOC_DECL(B,1,DenseMatrix_t);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DenseMatrix_t);

	init(A);
	/* init(Bt); */
	init(C);
	/* init(D); */


	col_perm = readFileSplit(A,B,C,D,fh);

	gettimeofday(&end,NULL);

	fprintf(stderr," LOAD time : %.3f\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec - start.tv_usec)/1e6));


	/* REDUCE */

	gettimeofday(&start,NULL);

	reduce(A,B,C,D);

	gettimeofday(&end,NULL);

	fprintf(stderr," REDUCE  time : %.3f\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec - start.tv_usec)/1e6));

	gettimeofday(&start,NULL);

	/* ECHELON */

	echelonD(A,D);

	gettimeofday(&end,NULL);

	fprintf(stderr," ECHELON time : %.3f\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec - start.tv_usec)/1e6));


	return 0;
}
