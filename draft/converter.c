#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <search.h>

#include "types.h"

#define OLD_TYPE uint16_t

#include "tools.h"

#include "selecter.h"


/* #define REVERT */
#ifdef REVERT
#warning "reverting rows"
#endif /* REVERT */
#ifdef SORT
#warning "reverting rows"
#endif /* SORT */



/* ./convert toto.gb */
int main( int ac, char ** av)
{
	assert(ac == 2);
	FILE * titi =fopen(av[1],"r"); /* in */
	char out[1024]; /* not too large the path... */
	strcpy(out,av[1]);
	strcat(out,"_new");
#ifdef REVERT
	strcat(out,"_rev");
#endif
	FILE * toto =fopen(out,"wb"); /* out */
	uint32_t un = Mjoin(select,elem_s)();
	fwrite(&un,sizeof(uint32_t),1,toto);

	SAFE_READ_DECL_V(m,uint32_t,titi);
	fwrite(&m,sizeof(uint32_t),1,toto);
	SAFE_READ_DECL_V(n,uint32_t,titi);
	fwrite(&n,sizeof(uint32_t),1,toto);
	SAFE_READ_DECL_V(mod,uint32_t,titi);
	fwrite(&mod,sizeof(elem_s),1,toto);
	SAFE_READ_DECL_V(nnz,uint64_t,titi);
	fwrite(&nnz,sizeof(uint64_t),1,toto);

	SAFE_READ_DECL_P(data,nnz,OLD_TYPE,titi);
	SAFE_READ_DECL_P(cols,nnz,uint32_t,titi);
	SAFE_READ_DECL_P(rows,m  ,uint32_t,titi);

	fclose(titi);

	/* start */
	SAFE_MALLOC_DECL(start_zo,m+1,uint64_t);
	start_zo[0] =  0 ;
	uint32_t i = 0 ;
	for ( ; i<m ; ++i)
	{
		start_zo[i+1] = start_zo[i] + rows[i] ;
	}

	assert(start_zo[m] == nnz);
	/* i = 0 ; for ( ; i < m+1 ; ++i) fprintf(stderr,"%lu ",start_zo[i]) ; fprintf(stderr,"\n"); */

	SAFE_MALLOC_DECL(map_zo_pol,m,uint32_t);

	SAFE_MALLOC_DECL(hash_row_pol,m,uint64_t);

	/* uint64_t comp = 0 ; */

	SAFE_MALLOC_DECL(colid_zo,nnz,uint32_t);


	/* compress colid */
	/* uint64_t cons = 0 ; */
	uint64_t here = 0 ;

#ifdef SORT
	SAFE_MALLOC_DECL(pivots,m,uint32_t);
	SAFE_MALLOC_DECL(permut,m,uint32_t);
	i = 0 ;
	for ( ; i < m ;++i) {
		pivots[i] = cols[start_zo[i]] ;
		permut[i] = i ;
	}
	insert_sort_duo_rev(pivots,m,permut);

#endif /* SORT */

#ifndef REVERT
	i = 0 ;
	for ( ; i < m ; ++i)
#else
	i = m ;
	for ( ; i-- ; )
#endif /* REVERT */
	{
#ifdef SORT
		uint64_t k = permut[i] ;
#else
		uint64_t k = i ;
#endif
		uint64_t j = start_zo[k] ;
		if ( j == start_zo[k+1] ) {  /* zero element */
			exit(-15);
		}
		/* saving element. will either be masked or next one is number
		 * of consecutive column indexes (>=2 then) */
		colid_zo[here] = cols[j] ;
		/* fprintf(stderr,"first in row %u\n", cols[j]); */
		if ( j + 1 == start_zo[k+1] ) { /* just one element */
			colid_zo[here++] |= NEGMASK ;
			/* fprintf(stderr,"was only in row, masking : %u\n",colid_zo[here]); */
			continue ;
		}
		++ j ;
		uint64_t cons = 0 ;
		for ( ; j < start_zo[k+1] ; ++j) { /* at least 2 elements */
			if (cols[j] == cols[j-1]+1) {
				/* fprintf(stderr,"next %u (%u)\n", cols[j], cols[j-1]); */
				cons += 1 ;
			}
			else { /* not consecutive */
				/* fprintf(stderr,"not consecutive %u (%u)\n", cols[j], cols[j-1]); */
				if (cons == 0) { /* last element was unit */
					/* fprintf(stderr,"was single : %u\n",colid_zo[here]); */
					colid_zo[here] |= NEGMASK;
					/* fprintf(stderr,"was single, masking : %u\n",colid_zo[here]); */
				}
				else { /* last element was last of sequence */
					colid_zo[++here] = cons+1 ;
					/* fprintf(stderr,"was not single : %u x %u\n",colid_zo[here-1],colid_zo[here]); */
				}
				/* next element to consider */
				cons = 0 ;
				colid_zo[++here] = cols[j];
			}
			if (j + 1 == start_zo[k+1]) { /* last element in row */
				/* fprintf(stderr,"end of row"); */
				if (cons == 0) {
					colid_zo[here] |= NEGMASK;
				}
				else {
					colid_zo[++here] = cons+1 ;
				}
				++ here ; /* prepare to write the next element */
				break;
			}
		}
	}

	fprintf(stderr,"saved %lu / %lu\n",here,nnz);
	SAFE_REALLOC(colid_zo,here,uint32_t);

	/* i = 0 ; for ( ; i < here ; ++i) fprintf(stderr,"%u ",colid_zo[i]) ; fprintf(stderr,"\n"); */

	ENTRY item;
	ENTRY *result;
	hcreate(m);

	typedef struct row { uint32_t i ; } row ;

	row * d ;

	/* map_zo_pol */
	uint32_t pol_nb = 0 ;
	i = 0 ;
	for ( ; i < m ; ++i) {
#ifdef SORT
		uint64_t k = permut[i];
#else
		uint64_t k = i ;
#endif

		uint64_t j0 = start_zo[k] ;
		uint64_t j1 = start_zo[k+1] ;
		char * seq = (char*) (data+j0);
		uint64_t length = (j1-j0)*sizeof(OLD_TYPE)/sizeof(char) ;

		/* uint64_t ii = 0 ; for ( ; ii < length ; ++ii) fprintf(stderr,"%u ",seq[ii]) ; fprintf(stderr,"\n"); */
		/* fprintf(stderr,"length %lu\n",length); */
		uint64_t key = JOAAT_hash(seq,length);
		/* fprintf(stderr,"key %lu\n",key); */
		char key_char[64] ;
		sprintf(key_char, "%lu", key);
		key_char[63]='\0';
		item.key =  key_char ;
		SAFE_MALLOC(d,1,row);
		d->i = pol_nb ;
		item.data = (char*) d;
		result = hsearch(item, FIND);
		if (result == NULL) {
			/* fprintf(stderr,"new\n"); */
			hsearch(item,ENTER);
#ifndef REVERT
			map_zo_pol[k] = pol_nb;
			hash_row_pol[pol_nb] = k ;
#else
			map_zo_pol[m-k-1] = pol_nb;
			hash_row_pol[pol_nb] = k ;
#endif
			pol_nb++;
		}
		else {
			/* fprintf(stderr,"old %u\n", ((row*)result->data)->i); */
#ifndef REVERT
			map_zo_pol[k] = ((row*)result->data)->i;
#else
			map_zo_pol[m-k-1] = ((row*)result->data)->i;
#endif
		}
	}
	hdestroy();
	/* uint64_t ii = 0 ; for ( ; ii <  m; ++ii) fprintf(stderr,"%u ",map_zo_pol[ii]) ; fprintf(stderr,"\n"); */
	/* ii = 0 ; for ( ; ii <  pol_nb; ++ii) fprintf(stderr,"%u ",hash_row_pol[ii]) ; fprintf(stderr,"\n"); */

	SAFE_MALLOC_DECL(pol_start,pol_nb+1,uint32_t);
	/* pol_start  */
	pol_start[0] = 0 ;
	i = 0 ;
	for ( ; i < pol_nb ; ++i) {
		uint64_t j0 = start_zo[hash_row_pol[i]] ;
		uint64_t j1 = start_zo[hash_row_pol[i]+1] ;
		pol_start[i+1] = pol_start[i]+(j1-j0);
	}

	uint32_t pol_nnz = pol_start[pol_nb];
	SAFE_MALLOC_DECL(pol_data,pol_nnz,elem_s);

	/* pol_data */
	i =  0 ;
	for ( ; i < pol_nb ; ++i) {
		uint64_t j0 = start_zo[hash_row_pol[i]] ;
		uint64_t j1 = start_zo[hash_row_pol[i]+1] ;
		/* fprintf(stderr,"%lu to %lu\n",j0,j1); */
		uint32_t k = 0 ;
		for ( ; k< j1-j0 ; ++k)
			pol_data[pol_start[i]+k]=data[j0+k];
	}
	/* i = 0 ; for ( ; i < pol_nnz ; ++i) fprintf(stderr,"%u ",pol_data[i]) ; fprintf(stderr,"\n"); */

#ifdef SORT
	i = 0 ;
	for ( ; i < m ; ++i) {
		uint64_t k = permut[i];
		start_zo[i+1] = start_zo[i] + rows[k] ;
	}

#endif /* SORT */

#ifdef REVERT
	i = 0 ;
	for ( ; i < m ; ++i) {
		start_zo[i+1] = start_zo[i] + rows[m-i-1] ;
	}
#endif /* REVERT */

	fwrite(start_zo,sizeof(uint64_t),m+1,toto);
       	fwrite(map_zo_pol,sizeof(uint32_t),m,toto);
	/* i = 0 ; for ( ; i < m ; ++i) fprintf(stderr,"%u ",map_zo_pol[i]) ; fprintf(stderr,"\n"); */

	fwrite(&here,sizeof(uint64_t),1,toto);
	fwrite(colid_zo,sizeof(uint32_t),here,toto);

	fwrite(&pol_nb,sizeof(uint32_t),1,toto);


	/* i = 0 ; for ( ; i < pol_nb+1 ; ++i) fprintf(stderr,"%u ",pol_start[i]) ; fprintf(stderr,"\n"); */
	fwrite(pol_start,sizeof(uint32_t),pol_nb+1,toto);


	fwrite(pol_data,sizeof(elem_s),pol_nnz,toto);



	fclose(toto);

	return 0;
}

