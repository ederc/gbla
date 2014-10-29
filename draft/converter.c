#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <search.h>

#include "types.h"

#define MASK (1U<<31)

uint64_t JOAAT_hash(char *key, size_t len)
{
	uint64_t hash, i;
	for(hash = i = 0; i < len; ++i)
	{
		hash += key[i];
		hash += (hash << 10);
		hash ^= (hash >> 6);
	}
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	return hash;
}

/* ./convert toto.gb */
int main( int ac, char ** av)
{
	assert(ac == 2);
	FILE * titi =fopen(av[1],"wb"); /* in */
	char out[80];
	strcpy(out,av[1]);
	strcpy(out,"_new");
	FILE * toto =fopen(out,"r"); /* out */
	uint32_t un = 1 ;
	fwrite(&un,sizeof(uint32_t),1,toto);


	SAFE_READ_DECL_V(m,uint32_t,titi);
	fwrite(&m,sizeof(uint32_t),1,toto);
	SAFE_READ_DECL_V(n,uint32_t,titi);
	fwrite(&n,sizeof(uint32_t),1,toto);
	SAFE_READ_DECL_V(mod,uint32_t,titi);
	fwrite(&mod,sizeof(uint32_t),1,toto);
	SAFE_READ_DECL_V(nnz,uint32_t,titi);
	fwrite(&nnz,sizeof(uint32_t),1,toto);

	SAFE_READ_DECL_P(data,nnz,uint16_t,titi);
	SAFE_READ_DECL_P(cols,nnz,uint32_t,titi);
	SAFE_READ_DECL_P(rows,m  ,uint32_t,titi);

	fclose(titi);

	/* start */
	SAFE_MALLOC_DECL(start_zo,m+1,uint64_t);
	start_zo[0] =  0 ;
	uint32_t i = 0 ;
	for ( ; i<m ; ++i) {
		start_zo[i+1] = start_zo[i] + rows[i] ;
	}
	fwrite(start_zo,sizeof(uint64_t),m+1,toto);

	SAFE_MALLOC_DECL(map_zo_pol,m,uint32_t);

	SAFE_MALLOC_DECL(hash_row_pol,m,uint64_t);

	uint64_t comp = 0 ;

	SAFE_MALLOC_DECL(colid_zo,nnz,uint32_t);


	/* compress colid */
	uint64_t cons = 0 ;
	uint64_t here = 0 ;

	i = 0 ;
	for ( ; i < m ; ++i) {
		uint64_t j = start_zo[i] ;
		if ( j == start_zo[i+1] ) {  /* zero element */
			exit(-15);
			continue ;
		}
		colid_zo[here] = data[cols[j]] ;
		if ( j + 1 == start_zo[i+1] ) { /* just one element */
			colid_zo[here] |= MASK ;
			continue ;
		}
		++ j ;
		uint64_t cons = 0 ;
		for ( ; j < start_zo[i+1] ; ++j) { /* at least 2 elements */
			if (cols[j] == cols[j-1]) {
				cons += 1 ;
			}
			else { /* not consecutive */
				if (cons == 0) { /* last element was unit */
					colid_zo[here] |= MASK;
				}
				else { /* last element was last of sequence */
					colid_zo[++here] = cons ;
				}
				cons = 0 ;
				++ here ;
			}
			if (j + 1 == start_zo[i+1]) { /* last element in row */
				if (cons == 0) {
					colid_zo[here] |= MASK;
				}
				else {
					colid_zo[++here] = cons ;
				}
				++ here ;
				break;
			}
		}
	}

	SAFE_REALLOC(colid_zo,here,uint32_t);

	ENTRY item;
	ENTRY *result;
	hcreate(m);

	/* map_zo_pol */
	uint32_t pol_nb = 0 ;
	i = 0 ;
	for ( ; i < m ; ++i) {
		uint64_t j0 = start_zo[i] ;
		uint64_t j1 =  start_zo[i+1] ;
		char * data = (char*) (data+j0);
		char   length = (j1-j0)*sizeof(uint64_t)/sizeof(char) ;
		uint64_t key = JOAAT_hash(data,length);
		char key_char[64] ;
		sprintf(key_char, "%lu", key);
		key_char[63]='\0';
		item.key =  key_char ;
		item.data = (uint32_t*) &i ;
		result = hsearch(item, FIND);
		if (result == NULL) {
			hsearch(item,ENTER);
			hash_row_pol[pol_nb] = i ;
			pol_nb++;
		}
		else {
			map_zo_pol[i] = (uint32_t) result->data ;
		}
	}
	hdestroy();

       	fwrite(map_zo_pol,sizeof(uint32_t),m,toto);

	fwrite(&here,sizeof(uint64_t),1,toto);
	fwrite(colid_zo,sizeof(uint32_t),here,toto);

	fwrite(&pol_nb,sizeof(uint32_t),1,toto);

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
	SAFE_MALLOC_DECL(pol_data,pol_nnz,uint32_t);

	/* pol_data */
	for ( ; i < pol_nb ; ++i) {
		uint64_t j0 = start_zo[hash_row_pol[i]] ;
		uint64_t j1 = start_zo[hash_row_pol[i]+1] ;
		uint32_t k = 0 ;
		for ( ; k< j1-j0 ; ++k)
			pol_data[pol_start[i]+k]=data[j0+k];
	}


	fwrite(pol_start,sizeof(uint32_t),pol_nb+1,toto);
	fwrite(pol_data,sizeof(uint32_t),pol_nnz,toto);



	fclose(toto);
}

