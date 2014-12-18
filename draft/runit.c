#include <stdint.h>


#include "io.h"
#include <sys/time.h>


int main(int ac, char **av) {

	if (ac < 2 || ac > 3) {
		fprintf(stderr,"usage %s nom_fichier in new rev sorted format [1|0] (1 if reduce)]\n",av[0]);
		return -1;
	}

	struct timeval start,end ;

	/* READ and SPLIT  */
	gettimeofday(&start,NULL);

	FILE * fh = fopen(av[1],"r") ;
	int red = 0 ;
	if (ac == 3) {
		red = atoi(av[2]);
		assert(red == 1 || red == 0);
	}
#if 0
	uint32_t * col_perm;
#endif

	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	SAFE_MALLOC_DECL(B,1,DenseMatrix_t);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DenseMatrix_t);

	init(A);
	/* init(Bt); */
	init(C);
	/* init(D); */

#if 0
	col_perm =
#endif
		readFileSplit(A,B,C,D,fh);

	gettimeofday(&end,NULL);

	fprintf(stderr," LOAD time : %.3f\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec - start.tv_usec)/1e6));


	/* REDUCE */

	gettimeofday(&start,NULL);

	if (red == 1) {
		reduce(A,B,C,D);
	}
	else {
		reduce_fast(A,B,C,D);
	}

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
