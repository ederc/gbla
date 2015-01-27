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


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <search.h>

#include "types.h"

#include "tools.h"

#include "selecter.h"

#include "ouvrir.h"


void usage(char * av)
{
	printf(" usage : \n");
	printf(" (1)   %s [option] mat\n",av);
	printf(" (2)   %s [option] mat_name - \n\n",av);

	printf(" options can be:\n");
	printf("   -n : convert from new format to old format. [default: old to new]\n");
	printf("   -s : sorts the row by weight (smallest first) [default : enabled, new format only]\n");
	printf("   -S : don't sort the rows by weight (smallest first) \n");
	printf("   -r : reverts order [default : enabled, new format only] \n");
	printf("   -R : don't revert order \n");
	printf(" (1) will produce a file named mat.gbm in the new format from matrix mat\n");
	printf(" (2) will produce a file named mat_name.gbm in the new format from  stdin \n");
	printf(" this would happen like : zcat mat.gz | %s mat - \n",av);
	printf(" Warning : the input matrix has to be a file in the old format\n");
	printf(" Warning : new matrix format is suffixed by .gbm\n");
}

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


void insert_sort_duo_rev(dimen_t * liste, dimen_t size, dimen_t * copain)
{
	dimen_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && (liste)[d] > (liste)[d-1]) {
			SWAP((liste)[d],(liste)[d-1]);
			SWAP((copain)[d],(copain)[d-1]);
			d--;
		}
	}
}


void convert_old2new(char * out, FILE * titi, int rev, int sor)
{
	strcat(out,".gbm");

	FILE * toto = fopen(out,"wb"); /* out */
	dimen_t un = Mjoin(select,elemt_s)();
	un = un | VERMASK ;
	fwrite(&un,sizeof(dimen_t),1,toto);

	SAFE_READ_DECL_V(m,  stor_t,titi);
	SAFE_READ_DECL_V(n,  stor_t,titi);
	SAFE_READ_DECL_V(mod,stor_t,titi);
	SAFE_READ_DECL_V(nnz,larg_t,titi);

	fwrite(&m,  sizeof(dimen_t),1,toto);
	fwrite(&n,  sizeof(dimen_t),1,toto);
	fwrite(&mod,sizeof(elemt_m),1,toto);
	fwrite(&nnz,sizeof(index_t),1,toto);

	SAFE_READ_DECL_P(data,nnz,elem_o,titi);
	SAFE_READ_DECL_P(cols,nnz,stor_t,titi);
	SAFE_READ_DECL_P(rows,m  ,stor_t,titi);

	fclose(titi);

	/* start */
	SAFE_MALLOC_DECL(start,(m+1),index_t);
	start[0] =  0 ;
	dimen_t i  ;
	for ( i=0; i<m ; ++i)
	{
		start[i+1] = start[i] + rows[i] ;
	}

	assert(start[m] == nnz);

	SAFE_MALLOC_DECL(map_zo_pol,m,dimen_t);

	SAFE_MALLOC_DECL(hash_row_pol,m,index_t);


	dimen_t * pivots = NULL ;
	dimen_t * permut = NULL ;

	if (sor) {
		SAFE_MALLOC(pivots,m,dimen_t);
		SAFE_MALLOC(permut,m,dimen_t);
		for ( i=0 ; i < m ;++i) {
			pivots[i] = cols[start[i]] ;
			permut[i] = i ;
		}
		insert_sort_duo_rev(pivots,m,permut);
		free(pivots);
	} /* SORT */

	dimen_t * cols_reord;
	index_t here = 0 ;
	if (rev || sor) {
		SAFE_MALLOC(cols_reord,nnz,dimen_t);
		if (rev) {
			for ( i = m ; i-- ;)
			{
				dimen_t k ;
				if (sor)
					k = permut[i] ;
				else
					k = i ;
				index_t j  ;
				for ( j = start[k] ; j < start[k+1] ; ++j) {
					cols_reord[here++] = cols[j];
				}
			}
		}
		else {
			for ( i=0 ; i < m ; ++i)
			{
				dimen_t k ;
				if (sor)
					k = permut[i] ;
				else
					k = i ;
				index_t j  ;
				for ( j = start[k] ; j < start[k+1] ; ++j) {
					cols_reord[here++] = cols[j];
				}
			}
		}
		free(cols);
	}
	else
		cols_reord = cols ;

	SAFE_MALLOC_DECL(colid,nnz,dimen_t);


	/* compress colid */
	here = 0 ;
	index_t j = 0 ;
	colid[here] = cols_reord[j++] ;
	int cons = 0;
	assert(nnz > 1);
	for ( ; j < nnz-1 ; ++ j) {
		if (cols_reord[j] == cols_reord[j-1]+1) {
			++cons;
		}
		else {
			if (cons == 0) {
				colid[here] |= NEGMASK ;
			}
			else {
				colid[++here] = cons+1 ;
			}
			cons = 0 ;
			colid[++here] = cols_reord[j];
		}
	}
	if (cols_reord[j] != cols_reord[j-1]+1) { /* last one */
		if (cons == 0) {
			colid[here] |= NEGMASK;
		}
		else {
			colid[++here] = cons+1 ;
		}
		colid[++here] = cols_reord[j] | NEGMASK;

	}
	else {
		colid[++here] = cons+2 ;
	}
	++here;

	if (! sor && ! rev)
		free(cols);

	fprintf(stderr,"saved %lu / %lu\n",here,(uint64_t)nnz);
	SAFE_REALLOC(colid,here,dimen_t);


	ENTRY item;
	ENTRY *result;
	hcreate(m);

	typedef struct row { dimen_t i ; } row ;

	row * d ;
	/* map_zo_pol */
	dimen_t pol_nb = 0 ;
	for ( i = 0 ; i < m ; ++i) {
		dimen_t k;
		if (sor)
			k = permut[i];
		else
			k = i ;

		index_t j0 = start[k] ;
		index_t j1 = start[k+1] ;
		char * seq = (char*) (data+j0);
		index_t length = (j1-j0)*sizeof(elem_o)/sizeof(char) ;

		uint64_t key = JOAAT_hash(seq,length);
		char key_char[64] ;
		sprintf(key_char, "%lu", key);
		key_char[63]='\0';
		item.key =  key_char ;
		SAFE_MALLOC(d,1,row); /* XXX this is a memory leak */
		d->i = pol_nb ;
		item.data = (char*) d;
		/* row d ; */
		/* d.i = pol_nb ; */
		/* item.data = (char*) &d; */
		result = hsearch(item, FIND);
		if (result == NULL) {
			hsearch(item,ENTER);
			if (!rev)
				map_zo_pol[k] = pol_nb;
			else
				map_zo_pol[m-k-1] = pol_nb;
			hash_row_pol[pol_nb] = k ;
			pol_nb++;
		}
		else {
			if (!rev)
				map_zo_pol[k] = ((row*)result->data)->i;
			else
				map_zo_pol[m-k-1] = ((row*)result->data)->i;
		}
	}
	hdestroy();

	SAFE_MALLOC_DECL(pol_start,(pol_nb+1),dimen_t);
	/* pol_start  */
	pol_start[0] = 0 ;
	for ( i=0 ; i < pol_nb ; ++i) {
		index_t j0 = start[hash_row_pol[i]] ;
		index_t j1 = start[hash_row_pol[i]+1] ;
		pol_start[i+1] = pol_start[i]+(j1-j0);
	}

	index_t pol_nnz = pol_start[pol_nb];
	SAFE_MALLOC_DECL(pol_data,pol_nnz,elemt_s);

	/* pol_data */
	for ( i=0 ; i < pol_nb ; ++i) {
		index_t j0 = start[hash_row_pol[i]] ;
		index_t j1 = start[hash_row_pol[i]+1] ;
		dimen_t k  ;
		for ( k=0; k< j1-j0 ; ++k)
			pol_data[pol_start[i]+k]=data[j0+k];
	}

	free(data);
	free(hash_row_pol);

	if (sor) {
		for ( i=0 ; i < m ; ++i) {
			dimen_t k = permut[i];
			if (rev)
				start[i+1] = start[i] + rows[m-k-1] ;
			else
				start[i+1] = start[i] + rows[k] ;
		}

		free(permut);

	}
	else {
		if (rev) {
			for (i=0 ; i < m ; ++i) {
				start[i+1] = start[i] + rows[m-i-1] ;
			}
		}
		else {
			/* nothing */
		}
	}

	/* XXX we don't need to write this, we could write rows */
	for ( i = 0 ; i < m ; ++i)
		rows[i] = start[i+1]-start[i] ; /* odd :) */

	fwrite(rows,sizeof(dimen_t),m,toto);
	free(rows);
	free(start);

	fwrite(map_zo_pol,sizeof(dimen_t),m,toto);
	free(map_zo_pol);

	fwrite(&here,sizeof(index_t),1,toto);
	fwrite(colid,sizeof(dimen_t),here,toto);
	free(colid);

	fwrite(&pol_nb,sizeof(dimen_t),1,toto);
	fwrite(&pol_nnz,sizeof(index_t),1,toto);

	SAFE_MALLOC_DECL(pol_rows,pol_nb,dimen_t);
	for ( i = 0 ; i < pol_nb ; ++i)
		pol_rows[i] = pol_start[i+1]-pol_start[i] ;

	fwrite(pol_rows,sizeof(dimen_t),pol_nb,toto);

	free(pol_rows);
	free(pol_start);

	fwrite(pol_data,sizeof(elemt_s),pol_nnz,toto);

	free(pol_data);


	fclose(toto);

	printf("created file %s\n",out);
}

void convert_new2old(char * out, FILE * fh)
{

	uint32_t lout = strlen(out);
	char nouv [1024];
	strncpy(nouv,out,lout-4);

	if (strcmp(strrchr(out,'.'),".gbm") != 0) {
		fprintf(stderr,"erreur : file %s has incorrect extension. Expects .gbm\n",out);
		fclose(fh);
		exit(-1);
	}
	nouv[lout-4]='\0';

	/* format */
	SAFE_READ_DECL_V(b,dimen_t,fh);
	/* XXX set elemt_s here and C++-ise*/
	if((b ^ VERMASK) != Mjoin(select,elemt_s)()) {
		fprintf(stderr,"bad format for %s. Expected %u, got %u\n",out,Mjoin(select,elemt_s)(),(b ^ VERMASK));
		exit(-1);
	}



	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,dimen_t,fh);
	/* assert(m < MAT_ROW_BLK); */
	SAFE_READ_DECL_V(n,dimen_t,fh);
	SAFE_READ_DECL_V(mod,elemt_m,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,index_t,fh);

	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %.3f%%) mod %lu\n",m,n,nnz,(double)(nnz)/((double)m*(double)n)*100.,(int64_t)mod);

	/* READ in ZO start */
	SAFE_READ_DECL_P(rows,m,dimen_t,fh);

	/* pol/zo correspondance */
	SAFE_READ_DECL_P(map_zo_pol,m,dimen_t,fh);


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,index_t,fh);
	SAFE_READ_DECL_P(buffer,colid_size,dimen_t,fh); /* buffer has the matrix */

	/* size of nnz for pols: */
	SAFE_READ_DECL_V(pol_nb,dimen_t,fh);

	SAFE_READ_DECL_V(pol_nnz,index_t,fh);

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_DECL_P(pol_rows,pol_nb,dimen_t,fh);

	/* XXX what if elemt_s == elemt_t ??? */
	SAFE_READ_DECL_P(data_pol,pol_nnz,elemt_s,fh);


	fclose(fh);

	FILE * toto =fopen(nouv,"wb"); /* out */

	stor_t m_ = m ;
	stor_t n_ = n ;
	stor_t p_ = mod;
	larg_t z_ = nnz ;
	fwrite(&m_,sizeof(stor_t),1,toto);
	fwrite(&n_,sizeof(stor_t),1,toto);
	fwrite(&p_,sizeof(stor_t),1,toto);
	fwrite(&z_,sizeof(larg_t),1,toto);

	SAFE_MALLOC_DECL(start,m+1,index_t);
	start[0] = 0 ;
	dimen_t i ;
	for ( i = 0 ; i < m ; ++i)
		start[i+1] = start[i] + rows[i];

	SAFE_MALLOC_DECL(start_pol,pol_nb+1,dimen_t);
	start_pol[0] = 0 ;
	for ( i = 0 ; i < pol_nb ; ++i)
		start_pol[i+1] = start_pol[i]+pol_rows[i];
	free(pol_rows);

	assert(start_pol[pol_nb] == pol_nnz);

	SAFE_MALLOC_DECL(colid,nnz,dimen_t); /* colid expands buffer */


	expandColid(buffer,colid_size,colid
#ifndef NDEBUG
			,nnz,n
#endif
		   );
	free(buffer);

	SAFE_MALLOC_DECL(data,nnz,elem_o);

	for (i = 0 ; i < m ; ++i) {
		dimen_t start_p = start_pol[ map_zo_pol[i] ] ;
		elemt_s * d = data_pol+start_p ;
		index_t jz ;
		for (jz = start[i] ; jz < start[i+1] ; ++jz) {
			data[jz] = (elem_o) *d ;
			++d ;
		}
	}

	free(map_zo_pol);

	SAFE_MALLOC_DECL(sz, m_,stor_t);
	SAFE_MALLOC_DECL(pos,z_,stor_t);
	SAFE_MALLOC_DECL(nz, z_,elem_o);

	for (i = 0 ; i < m_ ; ++i) sz[i] = rows[i] ;
	for (i = 0 ; i < z_ ; ++i) pos[i] = colid[i] ;
	for (i = 0 ; i < z_ ; ++i) nz[i] = data[i] ;

	fwrite(nz, sizeof(elem_o),z_,toto);
	fwrite(pos,sizeof(stor_t),z_,toto);
	fwrite(sz, sizeof(stor_t),m_,toto);

	fclose(toto);

	free(sz);
	free(pos);
	free(nz);

	free(rows);
	free(data);
	free(colid);
	free(start);
	free(data_pol);
	free(start_pol);
}

/* ./convert toto.gb */
int main( int ac, char ** av)
{
	int new = 0 ;
	int sor = 1 ;
	int rev = 1 ;

	if (ac < 2 || ac > 7) {
		usage(av[0]);
		exit(-1);
	}

	int options = 1 ;
	while ( options ) {
		if ( (strcmp(av[1],"-h") == 0) || (strcmp(av[1],"-?") == 0) || (strcmp(av[1],"--help") == 0) ) {
			usage(av[0]);
			return(0);
		}

		if (strcmp(av[1],"-n") == 0) {
			new = 1 ;
			av ++ ; ac -- ;
		}
		else if (strcmp(av[1],"-s") == 0) {
			sor = 1 ;
			av ++ ; ac -- ;
		} else if (strcmp(av[1],"-S") == 0) {
			sor = 0 ;
			av ++ ; ac -- ;
		} else if (strcmp(av[1],"-r") == 0) {
			rev = 1 ;
			av ++ ; ac -- ;
		} else if (strcmp(av[1],"-R") == 0) {
			rev = 0 ;
			av ++ ; ac -- ;
		} else
			options = 0 ;
	}

	if (! rev || ! sor) {
		fprintf(stderr, " warning : not sorting (%d) or not reverting (%d) rows !\n",sor,rev);
	}
	FILE * titi ;
	char out[1024]; /* not too large the path... */
	/* out[0] = '\0'; */

	if (ac == 3) {
		if ( strcmp(av[2],"-") != 0)  {
			usage(av[0]);
			exit(-1);
		}
	}

	titi = ouvrir(av[ac-1],"r"); /* in */

	strcpy(out,av[1]);

	if (new ==0) {
		printf("converting from old to new\n");
		convert_old2new(out,titi,rev,sor);
	}
	else {
		printf("converting from new to old\n");
		convert_new2old(out,titi);
	}

	return 0;
}

