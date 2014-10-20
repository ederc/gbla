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

#ifndef PROPOSED_FORMAT
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

/* matrix file is :
 * (everything in uint32_t unless specified
 * m n mod nnz (u64)
 * start_zo (u64) of size m+1 s.t. start_zo[i] points at the beginning of row[i] and start_zo[m] = nnz
 * map_pol_zo of size m that tells what polynomial is used in row i
 * colid_zo_size (u64) size of compressed colid_zo
 * colid_zo first one of each sequence followed by repetition, of size colid_sz
 * size_pol (u64) : number of polynomials
 * start_pol (u64) of size size_pol
 * vals_pol (s32) of size start_pol[size_pol]
 */

static void dump_matrix(char *fic,int all,int strict,int magma)
{
	FILE* fh=ouvrir(fic,"r");

	fprintf(stderr,"proposed format printing");

	SAFE_READ_V(m,  uint32_t,fh);
	SAFE_READ_V(n,  uint32_t,fh);
	SAFE_READ_V(mod,uint32_t,fh);
	SAFE_READ_V(nnz,uint64_t,fh);


	fprintf(stderr,"%u x %u matrix\n",n,m);
	fprintf(stderr,"mod %u\n",mod);
	{
		double Nz=(double)(n)*(double)(m);
		Nz=(double)(nnz)/Nz;
		Nz*=100.0;
		fprintf(stderr,"Nb of Nz elements %lu (density %.2f)\n",nnz,Nz);
	}

	if (all)
	{
		uint32_t i;

		/* start_zo */
		SAFE_READ_P(start_zo, m+1, uint64_t, fh); /* last element is an uint64_t. I don't beleive it is necessary to compress (as in row_length on uint16_t for instance but saving 3m/4, really ?) */

		/* map_pol_zo correspondance */
		SAFE_READ_P(map_pol_zo,m,uint32_t,fh); /* rows are numbered on 32 bits */

		/* if compressed, we need to know the size of colid_zo */
		SAFE_READ_V(colid_zo_size,uint32_t,fh); /* the size may be large (but < nnz/2) */

		/* colid_zo */
		SAFE_READ_P(colid_zo,colid_zo_size,uint32_t,fh); /*   columns are numbered on 32 bits */

		/* nnz_pol */
		SAFE_READ_V(nnz_pol,uint64_t,fh); /* number of coefficients may be large */

		/* start_pol */
		SAFE_READ_P(start_pol,m,uint64_t,fh); /* last number is uint64_t. Can be shrunk to save some space. at most m/4 i guess. */

		/* vals_pol */
		SAFE_READ_P(vals_pol,nnz_pol,int32_t,fh); /* elements are int32_t (or anything else) */


		if (magma)
		{
#if 1
			printf("K:=GF(%d);\n",mod);
			printf("sz:=[0 : i in [1..%u*%u]];\n",m,n);
			printf("A:=Matrix(K,%u,%u,sz);\n",m,n);
#else
			printf("A:=matrix(%u,%u):\n",n,m);
#endif /*  1 */
		}
		else
			if (strict)
				printf("%u %u M\n",m,n);
			else
				printf("%u\n%u\n",m,n);

		for(i=0;i<m;i++) {
			/* pointer to the values of the polynomial corresponding to row i  */
			int32_t * vals_pol_begin = vals_pol + start_pol[map_pol_zo[i]];
			uint32_t v ; /* just the value */
			uint32_t j = start_zo[i] ; /* C99 inside for :-( */
			for (  ; j < start_zo [i+1] ; ) {
				/* NEGMASK flag (last bit set) says that column
				 * is next column for that line has 0 ;
				 * otherwise, next element is the number of
				 * consecutive columns with non zeros at this
				 * row. */
				uint32_t first = colid_zo[j++] ;
				uint32_t repet = ((first & NEGMASK) == NEGMASK)?1:(colid_zo[j++]);
				first ^= NEGMASK ; /* get the actual first column by unmasking */
				uint32_t k = 0 ; /* C99 for this inside for :-( */
				for ( ; k < repet ; ++k) {
					v = vals_pol_begin[k]; /* consecutive values */

					/* fprintf(stderr,"<%u>",szi); */
					if (magma)
						printf("A[%u,%u]:=%u;\n",i+1,first+k+1,v);
					else
						printf("%u %u %d\n",i+1,first+k+1,v);

				}
				/* assert something sur first */
				vals_pol_begin += repet ; /* jump to the next "first" */
			}
		}
		fprintf(stderr,"\n");
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
