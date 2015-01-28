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
	printf(" (2)   %s [option] - \n\n",av);

	printf(" options can be:\n");
	printf("   -n : convert from new format to old format. [default: old to new]\n");
	printf("   -s : sorts the row by weight (smallest first) [default : enabled, new format only]\n");
	printf("   -S : don't sort the rows by weight (smallest first) \n");
	printf("   -r : reverts order [default : enabled, new format only] \n");
	printf("   -R : don't revert order \n");
	printf(" (1) will produce a file on stdout from matrix mat\n");
	printf(" (2) will produce a file on stdout from  stdin \n");
	printf(" this would happen like : zcat mat.gz | %s mat - > mat.gbm\n",av);
	printf(" Warning : the input matrix has to be in the correct format\n");
	printf(" Warning : new matrix format is suffixed by .gbm (but not necessary)\n");
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

void convert_old2new( FILE * titi, int rev, int sor)
{
	FILE * toto = stdout; /* out */
	dimen_t un = Mjoin(select,elemt_s)();
	un = un | VERMASK ;
	fwrite(&un,sizeof(dimen_t),1,toto);

	SAFE_READ_DECL_V(m,  stor_t,titi);
	SAFE_READ_DECL_V(n,  stor_t,titi);
	SAFE_READ_DECL_V(mod,stor_t,titi);
	SAFE_READ_DECL_V(nnz,larg_t,titi);


	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %.3f%%) mod %lu\n Reading...\n",m,n,nnz,(double)(nnz)/((double)m*(double)n)*100.,(int64_t)mod);

	fwrite(&m,  sizeof(dimen_t),1,toto);
	fwrite(&n,  sizeof(dimen_t),1,toto);
	fwrite(&mod,sizeof(elemt_m),1,toto);
	fwrite(&nnz,sizeof(index_t),1,toto);

	SAFE_READ_DECL_P(data,nnz,elem_o,titi);
	SAFE_READ_DECL_P(cols,nnz,stor_t,titi);
	SAFE_READ_DECL_P(rows,m  ,stor_t,titi);

	fclose(titi);
	fprintf(stderr," ...Finished reading\n");

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


	dimen_t * permut = NULL ;

	if (sor) {
		SAFE_MALLOC_DECL(pivots,m,dimen_t);
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
	here = compressColid(cols_reord,nnz,colid);

	if (! sor && ! rev)
		free(cols);

	fprintf(stderr,"saved %lu / %lu\n",here,(uint64_t)nnz);
	SAFE_REALLOC(colid,here,dimen_t);


	ENTRY item;
	ENTRY *result;
	hcreate(m);

	typedef struct row { dimen_t id ; } row ;

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
		d->id = pol_nb ;
		item.data = (char*) d;
		/* row d ; */
		/* d.id = pol_nb ; */
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
		else { /* found */
			{ /* check hash a little */

				dimen_t o = hash_row_pol[((row*)result->data)->id];
				index_t k0 = start[o];
				index_t k1 = start[o+1];
				int vrai = 1 ;
				if ( (k1-k0) != (j1-j0) ) {
					vrai = 0 ;
				}

				if (vrai) {
					index_t kk  ; /* first is suppose to be 1 */
					for (kk = 0 ; vrai && kk < (k1-k0)/10 ; ++kk)
						if (data[k0+1+kk*10] != data[k0+1+kk*10]) {
							vrai = 0;
						}
#if 0
					index_t zz,yy ;
					for (yy = j0, zz=k0 ; zz < k1 ; ++yy,++zz) {
						if (data[yy] != data[zz]) {
							fprintf(stderr,"bad hash %u %u\n",o,k);
							exit(-3);
						}
					}
#endif
				}

				if (!vrai) {
					if (!rev)
						map_zo_pol[k] = pol_nb;
					else
						map_zo_pol[m-k-1] = pol_nb;
					hash_row_pol[pol_nb] = k ;
					pol_nb++;
					continue ;
				}

			}
			if (!rev)
				map_zo_pol[k] = ((row*)result->data)->id;
			else
				map_zo_pol[m-k-1] = ((row*)result->data)->id;
		}
	}
	hdestroy();

	SAFE_MALLOC_DECL(pol_start,(pol_nb+1),index_t);
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
		for ( k=0; k< (dimen_t)(j1-j0) ; ++k)
			pol_data[pol_start[i]+(index_t)k]=data[j0+(index_t)k];
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
		pol_rows[i] = (dimen_t)(pol_start[i+1]-pol_start[i]) ;

	fwrite(pol_rows,sizeof(dimen_t),pol_nb,toto);

	free(pol_rows);
	free(pol_start);

	fwrite(pol_data,sizeof(elemt_s),pol_nnz,toto);

	free(pol_data);


	fclose(toto);

}

void convert_new2old( FILE * fh)
{


	/* format */
	SAFE_READ_DECL_V(b,dimen_t,fh);
	/* XXX set elemt_s here and C++-ise*/
	if((b ^ VERMASK) != Mjoin(select,elemt_s)()) {
		fprintf(stderr,"bad format. Expected %u, got %u\n",Mjoin(select,elemt_s)(),(b ^ VERMASK));
		exit(-1);
	}



	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,dimen_t,fh);
	/* assert(m < MAT_ROW_BLK); */
	SAFE_READ_DECL_V(n,dimen_t,fh);
	SAFE_READ_DECL_V(mod,elemt_m,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,index_t,fh);

	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %.3f%%) mod %lu\n Reading...\n",m,n,nnz,(double)(nnz)/((double)m*(double)n)*100.,(int64_t)mod);

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
	SAFE_READ_DECL_P(pol_data,pol_nnz,elemt_s,fh);


	fclose(fh);

	fprintf(stderr," ...Finished reading\n");

	FILE * toto = stdout ; /* out */

	stor_t m_ = m ;
	stor_t n_ = n ;
	stor_t p_ = mod;
	larg_t z_ = nnz ;
	fwrite(&m_,sizeof(stor_t),1,toto);
	fwrite(&n_,sizeof(stor_t),1,toto);
	fwrite(&p_,sizeof(stor_t),1,toto);
	fwrite(&z_,sizeof(larg_t),1,toto);

	SAFE_MALLOC_DECL(start,(m+1),index_t);
	start[0] = 0 ;
	dimen_t i ;
	for ( i = 0 ; i < m ; ++i)
		start[i+1] = start[i] + (index_t)rows[i];

	SAFE_MALLOC_DECL(pol_start,(pol_nb+1),index_t);
	pol_start[0] = 0 ;
	for ( i = 0 ; i < pol_nb ; ++i)
		pol_start[i+1] = pol_start[i]+(index_t)pol_rows[i];
	free(pol_rows);

	assert(pol_start[pol_nb] == pol_nnz);

	SAFE_MALLOC_DECL(colid,nnz,dimen_t); /* colid expands buffer */


	expandColid(buffer,colid_size,colid
#ifndef NDEBUG
			,nnz,n
#endif
		   );
	free(buffer);

	SAFE_MALLOC_DECL(data,nnz,elem_o);

	for (i = 0 ; i < m ; ++i) {
		index_t start_p = pol_start[ map_zo_pol[i] ] ;
		elemt_s * d = pol_data+start_p ;
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
	free(pol_data);
	free(pol_start);
}

/* ./convert toto.gb */
int main( int ac, char ** av)
{
	int new = 0 ;
	int sor = 1 ;
	int rev = 1 ;

	if (ac < 2 || ac > 5) {
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

	if (ac != 2) {
		usage(av[0]);
		exit(-1);
	}

	char * in = av[ac-1];
	titi = ouvrir(in,"r"); /* in */

	fprintf(stderr,"converting matrix %s",in);
	if (new ==0) {
		fprintf(stderr," (from old to new)\n");
		convert_old2new(titi,rev,sor);
	}
	else {
		fprintf(stderr," (from new to old)\n");
		if ( strcmp(in,"-") != 0) {
			char * found = strrchr(in,'.');
			if (!found || strcmp(found,".gbm") != 0) {
				fprintf(stderr,"warning : file %s has incorrect extension. Expects .gbm\n",in);
			}
		}
		convert_new2old(titi);
	}

	fprintf(stderr,"created file from %s\n",in);
	return 0;
}

