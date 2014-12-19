#include <stdint.h>


#include <sys/time.h>
#include <string.h>
#include "io.h"
/* #include <omp.h> */


void usage(char * nom) {
		fprintf(stderr,"usage %s [-r] (if reduce) nom_fichier in new rev sorted format  ]\n",nom);
}

FILE* ouvrir(char* a,char* b)
{
	if ( (strcmp(a,"-") == 0) && (strcmp(b,"r") == 0) )
		return stdin;
	else
	{
		FILE *f=fopen(a,b);
		if (!f)
		{
			fprintf(stderr,"Can't open %s\n",a);
			exit(1);
		}
		return f;
	}
}



int main(int ac, char **av) {

	/* omp_set_numthreads(1); */
	if (ac < 2 ) {
		usage(av[0]);
		return -1;
	}

	int red = 0 ;

	if (ac > 1 && ( (strcmp(av[1],"-h") == 0) ||(strcmp(av[1],"--help") == 0) || (strcmp(av[1],"-?") == 0) ) ) {
		usage(av[0]);
		return -1;
	}

	if (ac > 1 &&  (strcmp(av[1],"-r") == 0)) {
		red = 1 ;
		ac-- ;
		av++ ;
	}

	if (ac < 1) {
		usage(av[0]);
		return -1 ;
	}

	char * fic = av[1] ;

	FILE * fh = ouvrir(fic,"r");

	struct timeval start,end ;
	struct timeval aa ;

	/* READ and SPLIT  */
	gettimeofday(&start,NULL);
	gettimeofday(&aa,NULL);

#if 0
	uint32_t * col_perm;
#endif

	fprintf(stderr," reducing ? %u\n",(red==1));

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

	fprintf(stderr," LOAD    time         : %.3f s\n", ((double)(end.tv_sec - start.tv_sec)
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

	fprintf(stderr," REDUCE  time         : %.3f s\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec - start.tv_usec)/1e6));

	gettimeofday(&start,NULL);

	/* ECHELON */

	uint32_t r = echelonD(A,D);

	gettimeofday(&end,NULL);

	fprintf(stderr," ECHELON time         : %.3f s\n", ((double)(end.tv_sec - start.tv_sec)
				           +(double)(end.tv_usec - start.tv_usec)/1e6));

	fprintf(stderr,"  -- result           : %u\n",r);

	fprintf(stderr," TOTAL   time         : %.3f s\n", ((double)(end.tv_sec - aa.tv_sec)
				           +(double)(end.tv_usec - aa.tv_usec)/1e6));


	return 0;
}
