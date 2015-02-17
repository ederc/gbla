/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Brice Boyer <brice.boyer@lip6.fr>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdint.h>


#include <sys/time.h>
#include <string.h>
#include "io.h"
#ifdef _OPENMP
#include "reduce_B_omp.h"
#include "reduce_C_omp.h"
#else
#include "reduce_B.h"
#include "reduce_C.h"
#endif
#include "echelon.h"

#include "ouvrir.h"
#include <omp.h>


void usage(char * nom) {
	fprintf(stderr,"usage %s [options]  nom_fichier \n",nom);
	fprintf(stderr,"         nom_fichier in new rev sorted format\n");
	fprintf(stderr,"  -r     full row echelon form\n");
	fprintf(stderr,"  -t k   use k threads\n");

}


int main(int ac, char **av) {

	if (ac < 2 ) {
		usage(av[0]);
		return -1;
	}

	int red = 0 ;
	int nth = 1;

#ifdef _OPENMP
	nth = omp_get_num_threads() ;
#endif

	while (ac > 2) {

		if (ac > 1 && ( (strcmp(av[1],"-h") == 0) ||(strcmp(av[1],"--help") == 0) || (strcmp(av[1],"-?") == 0) ) ) {
			usage(av[0]);
			return -1;
		}

		if (ac > 1 && ( (strcmp(av[1],"-t") == 0)  ) ) {
			nth = atoi(av[2]) ;
			ac -=2;
			av +=2;
#ifdef _OPENMP
			if (nth > 1)
				omp_set_num_threads(nth);
#endif

		}
#ifndef _OPENMP
		assert(nth == 1);
#endif

		if (ac > 1 &&  (strcmp(av[1],"-r") == 0)) {
			red = 1 ;
			ac-- ;
			av++ ;
		}
	}

#ifdef _OPENMP
	if (nth == 1) {
		omp_set_num_threads(nth);
	}
#endif

	fprintf(stderr,"using  %d thread(s)\n",nth);

	if (ac != 2) {
		usage(av[0]);
		return -1 ;
	}


	char * fic = av[1] ;

	FILE * fh = ouvrir(fic,"r");

	struct timeval tic,tac,toc,tuc ;

	/* READ and SPLIT  */
	gettimeofday(&tic,NULL);
	gettimeofday(&tuc,NULL);

	dimen_t * col_perm;

	fprintf(stderr," reducing ? %u\n",(red==1));

	SAFE_MALLOC_DECL(A,1,GBMatrix_t);
	SAFE_MALLOC_DECL(B,1,GBMatrix_t);
	SAFE_MALLOC_DECL(C,1,GBMatrix_t);
	SAFE_MALLOC_DECL(D,1,DNS);

	initSparse(A);
	initSparse(B);
	initSparse(C);
	initDenseUnit(D);


	col_perm = readFileSplit(A,B,C,D,fh,&toc);

	fclose(fh);
	gettimeofday(&tac,NULL);

	fprintf(stderr," LOAD    time         : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

#ifndef ONLY_FFLAS
	/* REDUCE */

	gettimeofday(&tic,NULL);

	if (red == 1) {
		SAFE_MALLOC_DECL(Bd,1,DNS);
		initDenseUnit(Bd);
		convert_CSR_2_DNS(Bd,B);
		reduce_B(A,Bd,C,D,nth);
		freeMatDense(Bd);
	}
	else {
		int conv_a = 0 , conv_b = 0 , conv_c = 0 , algo_red = 1 ;
#ifdef CONV_A
		conv_a = 1 ;
#endif
#ifdef CONV_B
		conv_b = 1 ;
#endif
#ifdef CONV_C
		conv_c = 1 ;
#endif
#ifdef USE_SAXPY2
		algo_red = 2;
#endif
		reduce_C(A,conv_a,B,conv_b,C,conv_c,D,algo_red,nth);
	}

	gettimeofday(&tac,NULL);

	fprintf(stderr," REDUCE  time         : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

#endif /* ONLY_FFLAS */

	gettimeofday(&tic,NULL);

	/* ECHELON */

	uint32_t r = echelonD(A,D,nth);

	gettimeofday(&tac,NULL);

	fprintf(stderr," ECHELON time         : %.3f s\n", ((double)(tac.tv_sec - tic.tv_sec)
				+(double)(tac.tv_usec - tic.tv_usec)/1e6));

	fprintf(stderr,"  -- result           : %u\n",r);

	tuc.tv_sec  += toc.tv_sec ;
	tuc.tv_usec += toc.tv_usec ;
	fprintf(stderr," TOTAL   time         : %.3f s\n", ((double)(tac.tv_sec - tuc.tv_sec)
				+(double)(tac.tv_usec - tuc.tv_usec)/1e6));

	free(col_perm);
	freeMat(A);
	free(A);
	freeMat(C);
	free(C);
	freeMat(B);
	free(B);
	freeMatDense(D);
	free(D);

	return 0;
}

