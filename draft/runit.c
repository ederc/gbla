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

	struct timeval tic,tac ;
	struct timeval aa ;

	/* READ and SPLIT  */
	gettimeofday(&tic,NULL);
	gettimeofday(&aa,NULL);

	taille_t * col_perm;

	fprintf(stderr," reducing ? %u\n",(red==1));

	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	SAFE_MALLOC_DECL(B,1,DNS);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DNS);

	initSparse(A);
	initDenseUnit(B);
	initSparse(C);
	initDenseUnit(D);

	col_perm = readFileSplit(A,B,C,D,fh);

	fclose(fh);
	gettimeofday(&tac,NULL);

	fprintf(stderr," LOAD    time         : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				           +(double)(tac.tv_usec - tic.tv_usec)/1e6));
	fprintf(stderr,"   -- sparsity of A   : %.3f%% (%u x %u - %lu)\n",(double)A->nnz/(double)A->row/(double)A->col*100.,A->row,A->col,A->nnz);
	fprintf(stderr,"   -- sparsity of B   : %.3f%% (%u x %u - %lu)\n",(double)B->nnz/(double)B->row/(double)B->col*100.,B->row,B->col,B->nnz);
	fprintf(stderr,"   -- sparsity of C   : %.3f%% (%u x %u - %lu)\n",(double)C->nnz/(double)C->row/(double)C->col*100.,C->row,C->col,C->nnz);
	fprintf(stderr,"   -- sparsity of D   : %.3f%% (%u x %u - %lu)\n",(double)D->nnz/(double)D->row/(double)D->col*100.,D->row,D->col,D->nnz);


	/* REDUCE */

	gettimeofday(&tic,NULL);

	if (red == 1) {
		reduce(A,B,C,D);
	}
	else {
		reduce_fast(A,B,C,D);
	}

	gettimeofday(&tac,NULL);

	fprintf(stderr," REDUCE  time         : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				           +(double)(tac.tv_usec - tic.tv_usec)/1e6));

	gettimeofday(&tic,NULL);

	/* ECHELON */

	uint32_t r = echelonD(A,D);

	gettimeofday(&tac,NULL);

	fprintf(stderr," ECHELON time         : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				           +(double)(tac.tv_usec - tic.tv_usec)/1e6));

	fprintf(stderr,"  -- result           : %u\n",r);

	fprintf(stderr," TOTAL   time         : %.3f s\n", ((double)(tac.tv_sec - aa.tv_sec)
				           +(double)(tac.tv_usec - aa.tv_usec)/1e6));

	free(col_perm);
	freeMat(A);
	free(A);
	freeMat(C);
	free(C);
	freeMatDense(B);
	free(B);
	freeMatDense(D);
	free(D);

	return 0;
}
