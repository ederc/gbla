/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Brice Boyer <brice.boyer@lip6.fr>
 *                    Jean-Charles Faugère <jean-charles.faugere@inria.fr>
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
#include <stdio.h>
#include <assert.h>

FILE* ouvrir(char* a,char* b)
{
  if ((strcmp(a,"-") EQ 0) and (strcmp(b,"r") EQ 0))
    return stdin;
  else
    {
      FILE *f=fopen(a,b);
      if (not f)
	{
	  fprintf(stderr,"Can't open %s\n",a);
	  exit(1);
	}
      return f;
    }
}

static void dump_matrix(char *fic,int all,int strict,int magma)
{
  FILE* f=ouvrir(fic,"r");
  unsigned short int *nz;
  unsigned int       *pos;
  unsigned int       *row;
  unsigned int        *sz;
  unsigned int n;
  unsigned int m;
  unsigned int  mod;
  unsigned long long  nb;
  size_t n_inp;

  assert(n_inp=fread(&n,sizeof(unsigned int),        1,f)==1);
  assert(n_inp=fread(&m,sizeof(unsigned int),       1,f)==1);
  assert(n_inp=fread(&mod,sizeof(unsigned int),     1,f)==1);
  assert(n_inp=fread(&nb,sizeof(unsigned long long),1,f)==1);


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
      nz  = malloc(nb*sizeof(short    int));
      pos = malloc(nb*sizeof(unsigned int));
      sz  = malloc(n *sizeof(unsigned int));
      assert(fread(nz, sizeof(short int),   nb,f)==nb);
      assert(fread(pos,sizeof(unsigned int),nb,f)==nb);
      assert(fread(sz, sizeof(unsigned int),n, f)==n );
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
	      printf("A[%u,%u]:=%u;\n",i+1,pos[j]+1,(unsigned int)(nz[j]));
	  else
	    for(j=0;j<szi;j++)
	      printf("%u %u %u\n",i+1,pos[j]+1,(unsigned int)(nz[j]));
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

int main(int nargs,char** argv)
{
  int all=0;
  int strict=0;
  int magma=0;
  if (nargs>1 and (strcmp(argv[1],"-l") EQ 0))
    {
      all=1;
      strict=0;
      nargs--;
      argv++;
    }

  if (nargs>1 and (strcmp(argv[1],"-m") EQ 0))
    {
      all=1;
      strict=0;
      magma=1;
      nargs--;
      argv++;
    }

  /* SMS format */
  if (nargs>1 and (strcmp(argv[1],"-s") EQ 0))
    {
      all=1;
      strict=1;
      nargs--;
      argv++;
    }

  {
    char* fic=(nargs>1 ? argv[1] : "mat1");
    dump_matrix(fic,all,strict,magma);
  }

  return 0;
}
