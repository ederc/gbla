/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Brice Boyer <brice.boyer@lip6.fr>
 *                    Jean-Charles Faugère <jean.charles-faugere@inria.fr>
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


#define true 1
#define false 0
#define nil 0
#ifdef __GNUC__
#if __GNUC__<3 || (! defined(__GNUG__))
#define and &&
#define or ||
#define not !
#endif /*  __GNUC__<3 || (! defined(__GNUG__)) */
#endif
#define NEQ !=
#define EQ ==
#define loop while(true)

/* Style et type */
typedef short I16;
typedef long unsigned UI32_or_64;
#ifndef Visual_cpp
typedef long long unsigned UI64;
#else
typedef long long unsigned UI64;
#define not !
#define and &&
#define or ||
#endif /* ndef Visual_cpp */
typedef long I32_or_64;

#ifndef vis_String
typedef char* String;
typedef char const * const Cste_String;
#endif /* ndef vis_String */

typedef char unsigned UI8;
typedef float F32;
typedef double F64;

typedef unsigned UI32;
typedef short unsigned UI16;
typedef int I32;
typedef char I8;

/* Style et Basic */
#ifndef vis_Boolean
typedef I32 Boolean;
#endif /* ndef vis_Boolean */

typedef void* Any;

#include <stdlib.h>
#include <string.h> /* strcmp */
#include <stdio.h>
#include <assert.h>

#include <stdint.h> /* uint32_t et cie. */

#include "types.h"
#include "macros.h"
#include "ouvrir.h"


static void dump_matrix_old(char *fic,int all,int strict,int magma)
{
	FILE* f=ouvrir(fic,"r");
	/* short int *nz; */
	/* unsigned int       *pos; */
	/* unsigned int       *row; */
	/* unsigned int        *sz; */
	unsigned int n;
	unsigned int m;
	unsigned int  mod;
	unsigned long long  nb;
	fprintf(stderr,"current format printing\n");

	SAFE_READ_V(n,unsigned int,f);
	SAFE_READ_V(m,unsigned int,f);
	SAFE_READ_V(mod,unsigned int,f);
	SAFE_READ_V(nb,unsigned long long,f);


	fprintf(stderr,"%u x %u matrix\n",n,m);
	fprintf(stderr,"mod %u\n",mod);
	{
		double Nz=(double)(n)*(double)(m);
		Nz=(double)(nb)/Nz;
		Nz*=100.0;
		fprintf(stderr,"Nb of Nz elements %Lu (density %.2f)\n",nb,Nz);
	}

	if (all)
	{
		unsigned int i;
		SAFE_READ_DECL_P(nz,  nb, short int, f);
		SAFE_READ_DECL_P(pos, nb, unsigned int, f);
		SAFE_READ_DECL_P(sz,  n,  unsigned int, f);
		if (magma)
		{
#if 1
			printf("K:=GF(%d);\n",mod);
			printf("sz:=[0 : i in [1..%u*%u]];\n",n,m);
			printf("A:=Matrix(K,%u,%u,sz);\n",n,m);
#else
			printf("A:=matrix(%u,%u):\n",n,m);
#endif /*  1 */
		}
		else
			if (strict)
				printf("%u %u M\n",n,m);
			else
				printf("%u\n%u\n",n,m);

		for(i=0;i<n;i++)
		{
			const unsigned int szi=sz[i];
			unsigned int j;
			/* fprintf(stderr,"<%u>",szi); */
			if (magma)
				for(j=0;j<szi;j++)
					printf("A[%u,%u]:=%d;\n",i+1,pos[j]+1,(int)(nz[j]));
			else
				for(j=0;j<szi;j++)
					printf("%u %u %d\n",i+1,pos[j]+1,(int)(nz[j]));
			nz+=szi;
			pos+=szi;
		}
		fprintf(stderr,"\n");
	}

	if (magma)
		printf("// ----------------------------------------\n");
	else
		if (strict)
			printf("0 0 0\n");
	fclose(f);
}

#include "print_helper.h"

#include "tools.h"

#define TEMPL_TYPE int8_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE int16_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE int32_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE int64_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE uint8_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE uint16_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE uint32_t
#include "dump_matrix.h"
#undef TEMPL_TYPE
#define TEMPL_TYPE uint64_t
#include "dump_matrix.h"
#undef TEMPL_TYPE

/* matrix file is :
 * (everything in uint32_t unless specified)
 * bitseq :
 * * bit 0   | 0 iff unsigned
 * * bit 1,2 | 0 : 8 bit ; 1 : 16 bit ; 2 : 32bit ; 3 : 64 bit (this is 8*2^k)
 * m n mod nnz (u64) everything is 0 based.
 * start_zo (u64) of size m+1 s.t. start_zo[i] points at the beginning of row[i] and start_zo[m] = nnz
 * map_zo_pol of size m that tells what polynomial is used in row i
 * colid_zo_size (u64) size of compressed colid_zo
 * colid_zo first one of each sequence followed by repetition, of size colid_sz
 * size_pol  : number of polynomials
 * start_pol : of size size_pol+1
 * vals_pol (s32) of size start_pol[size_pol]
 * idx_pol : indexation for polyomials (size size_pol) ie. unique id for each pol.
 */

static void dump_matrix_new(char *fic,int all,int strict,int magma)
{
	FILE* fh=ouvrir(fic,"r");

	fprintf(stderr,"proposed format printing\n");

	SAFE_READ_DECL_V(b,  uint32_t,fh);

	if ((b & VERMASK) != VERMASK)  {
		fprintf(stderr,"error : file %s has bad version number \n",fic);
	}

	b = b ^ VERMASK ;

	switch(b&1)
	{
		case (0) :
			switch (b>>1)
		{
			case  0 :
				dump_uint8_t(fh,all,strict,magma);
				break;
			case  1 :
				dump_uint16_t(fh,all,strict,magma);
				break;
			case  2 :
				dump_uint32_t(fh,all,strict,magma);
				break;
			case  3 :
				dump_uint64_t(fh,all,strict,magma);
				break;
			default :
				exit(-2);
		}
			break ;
		case (1) :
			switch (b>>1)
		{
			case  0 :
				dump_int8_t(fh,all,strict,magma);
				break;
			case  1 :
				dump_int16_t(fh,all,strict,magma);
				break;
			case  2 :
				dump_int32_t(fh,all,strict,magma);
				break;
			case  3 :
				dump_int64_t(fh,all,strict,magma);
				break;
			default :
				exit(-2);
		}
	}

	if (magma)
		printf("// ----------------------------------------\n");
	else
		if (strict)
			printf("0 0 0\n");
	fclose(fh);
}

void help(char * prog) {
	printf("usage : %s [options] mat \n",prog);
	printf("usage : %s [options] -\n\n",prog);

	printf(" [default] prints the header of the matrix if no option is given\n");
	printf(" options can be one of :\n");
	printf("     -h : print this help\n");
	printf("     -l : dump the whole matrix in SMS format\n");
	printf("     -m : dump the whole matrix in magma format\n");
	printf("     -s : dump the whole matrix in strict SMS format\n");
	printf(" options can be additionally:\n");
	printf("     -n : uses new format [default : old format]\n");
	printf(" filename argument: \n");
	printf("    mat : represents the matrix is the correct format.\n");
	printf("    -   : matrix read from stdin as in zcat mat.gz | %s -\n",prog);
	printf("warning: filename has to be in the correct format\n");
}

int main(int nargs,char** argv)
{
	int all   =0;
	int strict=0;
	int magma =0;
	int new   =0;

	if (nargs < 2 || nargs > 4) {
		fprintf(stderr,"erreur : wrong number of arguments\n");
		help(argv[0]);
		exit(-1);
	}

	while (nargs > 2) {

		if ( (strcmp(argv[1],"-h") EQ 0) or (strcmp(argv[1],"-?") EQ 0) or (strcmp(argv[1],"--help") EQ 0) )
		{
			help(argv[0]);
			return 0;
		}

		else if (strcmp(argv[1],"-l") EQ 0)
		{
			all=1;
			strict=0;
			nargs--;
			argv++;
		}

		else if (strcmp(argv[1],"-n") EQ 0)
		{
			new=1;
			nargs--;
			argv++;
		}

		else if (strcmp(argv[1],"-m") EQ 0)
		{
			all=1;
			strict=0;
			magma=1;
			nargs--;
			argv++;
		}

		/* SMS format */
		else if (strcmp(argv[1],"-s") EQ 0)
		{
			all=1;
			strict=1;
			nargs--;
			argv++;
		}

		else {
			fprintf(stderr,"erreur : bad argument %s\n",argv[1]);
			help(argv[0]);
			exit(-1);
		}
	}

	if ( (strcmp(argv[1],"-h") EQ 0) or (strcmp(argv[1],"-?") EQ 0) or (strcmp(argv[1],"--help") EQ 0) )
	{
		help(argv[0]);
		return 0;
	}

	char* fic= argv[1] ;
	if (new == 1)
		dump_matrix_new(fic,all,strict,magma);
	else {
		printf("old, %u %u %u\n",all,strict,magma);
		dump_matrix_old(fic,all,strict,magma);
	}

	return 0;
}
