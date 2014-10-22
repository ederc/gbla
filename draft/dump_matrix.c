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

#ifndef _PROPOSED_FORMAT
static void dump_matrix(char *fic,int all,int strict,int magma)
{
  FILE* f=ouvrir(fic,"r");
  unsigned short int *nz;
  unsigned int       *pos;
  /* unsigned int       *row; */
  unsigned int        *sz;
  unsigned int n;
  unsigned int m;
  unsigned int  mod;
  unsigned long long  nb;
  size_t n_inp;
  fprintf(stderr,"current format printing");

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
#else



#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })

#define NEGMASK  (1U<<31)

#define SAFE_MALLOC(ptr,size,elt) \
	        elt * ptr = (elt *) malloc(size*sizeof(elt)); \
	        assert(ptr)

#define SAFE_CALLOC(ptr,size,elt) \
	        elt * ptr = (elt *) calloc(size*sizeof(elt)); \
	        assert(ptr)

#define SAFE_REALLOC(ptr,size,elt) \
	        ptr = realloc(ptr,size*sizeof(elt)); \
	        assert(ptr)

#define SAFE_READ_V(val,elt,file) \
	elt val ; \
	assert(fread(&(val),sizeof(elt),1,file)==1)


#define SAFE_READ_P(val,size,elt,file) \
	SAFE_MALLOC(val,size,elt); \
	assert(fread(val,sizeof(elt),size,file)==size)

#include "printer.h"


#define TYPE int8_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE int16_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE int32_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE int64_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE uint8_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE uint16_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE uint32_t
#include "dump_matrix.h"
#undef TYPE
#define TYPE uint64_t
#include "dump_matrix.h"
#undef TYPE



/* matrix file is :
 * (everything in uint32_t unless specified)
 * bitseq :
 * * bit 0   | 0 iff unsigned
 * * bit 1,2 | 0 : 8 bit ; 1 : 16 bit ; 2 : 32bit ; 3 : 64 bit (this is 8*2^k)
 * m n mod nnz (u64) everything is 0 based.
 * start_zo (u64) of size m+1 s.t. start_zo[i] points at the beginning of row[i] and start_zo[m] = nnz
 * map_pol_zo of size m that tells what polynomial is used in row i
 * colid_zo_size (u64) size of compressed colid_zo
 * colid_zo first one of each sequence followed by repetition, of size colid_sz
 * size_pol  : number of polynomials
 * start_pol : of size size_pol+1
 * vals_pol (s32) of size start_pol[size_pol]
 * idx_pol : indexation for polyomials (size size_pol) ie. unique id for each pol.
 */

static void dump_matrix(char *fic,int all,int strict,int magma)
{
	FILE* fh=ouvrir(fic,"r");

	fprintf(stderr,"proposed format printing");

	SAFE_READ_V(b,  uint32_t,fh);

	switch(b&0)
	{
		case (0) :
			switch (b<<1)
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
			switch (b<<1)
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
	}

	if (magma)
		printf("// ----------------------------------------\n");
	else
		if (strict)
			printf("0 0 0\n");
	fclose(fh);
}

#endif

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
