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

#define OLD_TYPE short int

#include "tools.h"

#include "selecter.h"

#include "ouvrir.h"

#ifndef REVERT
#warning "not reverting rows"
#endif /* REVERT */
#ifndef SORT
#warning "not sorting rows"
#endif /* SORT */


void usage(char * av) {
	printf(" usage : \n");
	printf(" (1)   %s [option] mat\n",av);
	printf(" (2)   %s [option] mat_name - \n\n",av);

	printf(" option can be:\n");
	printf("   -n : convert from new format to old format. Default: old to new\n");
	printf(" (1) will produce a file named mat.gbm in the new format from matrix mat\n");
	printf(" (2) will produce a file named mat_name.gbm in the new format from  stdin \n");
	printf(" this would happen like : zcat mat.gz | %s mat - \n",av);
	printf(" Warning : the input matrix has to be a file in the old format\n");
	printf(" Warning : %s will sort the matrix by denser to sparser rows\n",av);
	printf("           compile without -DREVERT or -DSORT to change that (see Makefile).\n");
	printf(" Warning : new matrix format is suffixed by .gbm\n");
}


void convert_old2new(char * out, FILE * titi) {
#if defined(REVERT) && defined(SORT)
	strcat(out,".gbm");
#else
	strcat(out,"_new");
#ifdef REVERT
	strcat(out,"_rev");
#endif

#ifdef SORT
	strcat(out,"_srt");
#endif
#endif

	FILE * toto = fopen(out,"wb"); /* out */
	uint32_t un = Mjoin(select,elemt_s)();
	printf("un : %u\n",un);
	un = un | VERMASK ;
	fwrite(&un,sizeof(uint32_t),1,toto);

	SAFE_READ_DECL_V(m,unsigned int,titi);
	SAFE_READ_DECL_V(n,unsigned int,titi);
	SAFE_READ_DECL_V(mod,unsigned int,titi);
	SAFE_READ_DECL_V(nnz,unsigned long long int,titi);

	fwrite(&m,sizeof(uint32_t),1,toto);
	fwrite(&n,sizeof(uint32_t),1,toto);
	fwrite(&mod,sizeof(elemt_s),1,toto);
	/* XXX we don't need to write this */
	fwrite(&nnz,sizeof(uint64_t),1,toto);

	SAFE_READ_DECL_P(data,nnz,OLD_TYPE,titi);
	SAFE_READ_DECL_P(cols,nnz,unsigned int,titi);
	SAFE_READ_DECL_P(rows,m  ,unsigned int,titi);

	fclose(titi);

	/* start */
	SAFE_MALLOC_DECL(start,(m+1),uint64_t);
	start[0] =  0 ;
	uint32_t i  ;
	for ( i=0; i<m ; ++i)
	{
		start[i+1] = start[i] + rows[i] ;
	}

	assert(start[m] == nnz);

	SAFE_MALLOC_DECL(map_zo_pol,m,uint32_t);

	SAFE_MALLOC_DECL(hash_row_pol,m,uint64_t);

	/* uint64_t comp = 0 ; */


#ifdef SORT
	SAFE_MALLOC_DECL(pivots,m,uint32_t);
	SAFE_MALLOC_DECL(permut,m,uint32_t);
	for ( i=0 ; i < m ;++i) {
		pivots[i] = cols[start[i]] ;
		permut[i] = i ;
	}
	insert_sort_duo_rev(pivots,m,permut);
	free(pivots);
#endif /* SORT */

#if 0
	SAFE_MALLOC_DECL(colid,nnz,uint32_t);
	uint64_t here = 0 ;

#ifndef REVERT
        for ( i=0 ; i < m ; ++i)
#else
        for ( i = m ; i-- ; )
#endif /* REVERT */
        {
#ifdef SORT
                uint64_t k = permut[i] ;
#else
                uint64_t k = i ;
#endif
                uint64_t j = start[k] ;
                if ( j == start[k+1] ) {  /* zero element */
                        exit(-15);
                }
                /* saving element. will either be masked or next one is number
                 * of consecutive column indexes (>=2 then) */
                colid[here] = cols[j] ;
                assert(colid[here] < n);
                if ( j + 1 == start[k+1] ) { /* just one element */
                        colid[here++] |= NEGMASK ;
                        continue ;
                }
                ++ j ;
                uint64_t cons = 0 ;
                for ( ; j < start[k+1] ; ++j) { /* at least 2 elements */
                        if (cols[j] == cols[j-1]+1) {
                                cons += 1 ;
                        }
                        else { /* not consecutive */
                                if (cons == 0) { /* last element was unit */
                                        colid[here] |= NEGMASK;
                                }
                                else { /* last element was last of sequence */
                                        colid[++here] = cons+1 ;
                                }
                                /* next element to consider */
                                cons = 0 ;
                                colid[++here] = cols[j];
                                assert(colid[here] < n);
                        }
                        if (j + 1 == start[k+1]) { /* last element in row */
                                if (cons == 0) {
                                        colid[here] |= NEGMASK;
                                }
                                else {
                                        colid[++here] = cons+1 ;
                                }
                                ++ here ; /* prepare to write the next element */
                                break;
                        }
                }

        }
        free(cols);

#else
	uint32_t * cols_reord;
	uint64_t here = 0 ;
#if defined(REVERT) ||  defined(SORT)
SAFE_MALLOC(cols_reord,nnz,uint32_t);
	for (
#ifndef REVERT
			i=0 ; i < m ; ++i
#else
			i = m ; i-- ;
#endif /* REVERT */
	    )
	{
#ifdef SORT
		uint64_t k = permut[i] ;
#else
		uint64_t k = i ;
#endif
		uint64_t j  ;
		for ( j = start[k] ; j < start[k+1] ; ++j) {
			cols_reord[here++] = cols[j];
		}
	}
	free(cols);
#else
	cols_reord = cols ;
#endif

	SAFE_MALLOC_DECL(colid,nnz,uint32_t);


	/* compress colid */
	/* uint64_t cons = 0 ; */
	here = 0 ;
	for ( i=0 ; i < m ; ++i) {
		uint64_t j = start[i] ;
		if ( j == start[i+1] ) {  /* zero element */
			exit(-15);
		}
		/* saving element. will either be masked or next one is number
		 * of consecutive column indexes (>=2 then) */
		colid[here] = cols_reord[j] ;
		assert(colid[here] < n);
		if ( j + 1 == start[i+1] ) { /* just one element */
			colid[here++] |= NEGMASK ;
			continue ;
		}
		++ j ;
		uint64_t cons = 0 ;
		for ( ; j < start[i+1] ; ++j) { /* at least 2 elements */
			if (cols_reord[j] == cols_reord[j-1]+1) {
				cons += 1 ;
			}
			else { /* not consecutive */
				if (cons == 0) { /* last element was unit */
					colid[here] |= NEGMASK;
				}
				else { /* last element was last of sequence */
					colid[++here] = cons+1 ;
				}
				/* next element to consider */
				cons = 0 ;
				colid[++here] = cols_reord[j];
				assert(colid[here] < n);
			}
			if (j + 1 == start[i+1]) { /* last element in row */
				if (cons == 0) {
					colid[here] |= NEGMASK;
				}
				else {
					colid[++here] = cons+1 ;
				}
				++ here ; /* prepare to write the next element */
				break;
			}
		}

	}

#if !defined(SORT) && !defined(REVERT)
	free(cols);
#endif
#endif

	fprintf(stderr,"saved %lu / %llu\n",here,nnz);
	SAFE_REALLOC(colid,here,uint32_t);


	ENTRY item;
	ENTRY *result;
	hcreate(m);

	typedef struct row { uint32_t i ; } row ;

	row * d ;
	/* map_zo_pol */
	uint32_t pol_nb = 0 ;
	for ( i = 0 ; i < m ; ++i) {
#ifdef SORT
		uint64_t k = permut[i];
#else
		uint64_t k = i ;
#endif

		uint64_t j0 = start[k] ;
		uint64_t j1 = start[k+1] ;
		char * seq = (char*) (data+j0);
		uint64_t length = (j1-j0)*sizeof(OLD_TYPE)/sizeof(char) ;

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
#ifndef REVERT
			map_zo_pol[k] = pol_nb;
#else
			map_zo_pol[m-k-1] = pol_nb;
#endif
			hash_row_pol[pol_nb] = k ;
			pol_nb++;
		}
		else {
#ifndef REVERT
			map_zo_pol[k] = ((row*)result->data)->i;
#else
			map_zo_pol[m-k-1] = ((row*)result->data)->i;
#endif
		}
	}
	hdestroy();

	SAFE_MALLOC_DECL(pol_start,(pol_nb+1),uint32_t);
	/* pol_start  */
	pol_start[0] = 0 ;
	for ( i=0 ; i < pol_nb ; ++i) {
		uint64_t j0 = start[hash_row_pol[i]] ;
		uint64_t j1 = start[hash_row_pol[i]+1] ;
		pol_start[i+1] = pol_start[i]+(j1-j0);
	}

	uint64_t pol_nnz = pol_start[pol_nb];
	SAFE_MALLOC_DECL(pol_data,pol_nnz,elemt_s);

	/* pol_data */
	for ( i=0 ; i < pol_nb ; ++i) {
		uint64_t j0 = start[hash_row_pol[i]] ;
		uint64_t j1 = start[hash_row_pol[i]+1] ;
		uint32_t k  ;
		for ( k=0; k< j1-j0 ; ++k)
			pol_data[pol_start[i]+k]=data[j0+k];
	}

	free(data);
	free(hash_row_pol);

#ifdef SORT
	for ( i=0 ; i < m ; ++i) {
		uint64_t k = permut[i];
		start[i+1] = start[i] + rows[k] ;
	}

	free(permut);

#endif /* SORT */

#ifdef REVERT
	for (i=0 ; i < m ; ++i) {
		start[i+1] = start[i] + rows[m-i-1] ;
	}
#endif /* REVERT */

	/* XXX we don't need to write this, we could write rows */
	for ( i = 0 ; i < m ; ++i)
		rows[i] = start[i+1]-start[i] ;

	fwrite(rows,sizeof(uint32_t),m,toto);
	free(rows);
	free(start);

	fwrite(map_zo_pol,sizeof(uint32_t),m,toto);
	free(map_zo_pol);

	fwrite(&here,sizeof(uint64_t),1,toto);
	fwrite(colid,sizeof(uint32_t),here,toto);
	free(colid);

	fwrite(&pol_nb,sizeof(uint32_t),1,toto);
	fwrite(&pol_nnz,sizeof(uint64_t),1,toto);

	SAFE_MALLOC_DECL(pol_rows,pol_nb,uint32_t);
	for ( i = 0 ; i < pol_nb ; ++i)
		pol_rows[i] = pol_start[i+1]-pol_start[i] ;

	fwrite(pol_rows,sizeof(uint32_t),pol_nb,toto);

	free(pol_rows);
	free(pol_start);

	fwrite(pol_data,sizeof(elemt_s),pol_nnz,toto);

	free(pol_data);


	fclose(toto);

	printf("created file %s\n",out);
}


void expandColid(
	       	const dimen_t * compress
		, index_t          size_compressed
		, dimen_t       * expand
#ifndef NDEBUG
		, index_t          size_expand
		, dimen_t         n
#endif
		)
{
	uint32_t mask = (1<<31);
	index_t i = 0 ;
	index_t j = 0 ;
	dimen_t col ;
	for ( ; i < size_compressed ;) {
		col = compress[i++] ;
		if (col & mask) {
			expand[j++] = col ^ mask ;
		}
		else {
			dimen_t k = 0 ;
			for (; k < compress[i] ;++k) {
				expand[j++] = col + k;
			}
			++i;
		}
	}
	assert(j == size_expand);
	assert(i == size_compressed);
#ifndef NDEBUG
	for ( i = 0 ; i < size_expand ; ++i) {
		assert(expand[i] < n);
	}
#endif
}


void convert_new2old(char * out, FILE * fh) {

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
	SAFE_READ_DECL_V(b,uint32_t,fh);
	/* XXX set elemt_s here and C++-ise*/
	if((b ^ VERMASK) != Mjoin(select,elemt_s)()) {
		fprintf(stderr,"bad format for %s\n",out);
		exit(-1);
	}

	assert( sizeof( short int ) == sizeof( int16_t ) );
	assert( sizeof( unsigned int ) == sizeof( uint32_t ) );
	assert( sizeof( unsigned long long int ) == sizeof( uint64_t ) );



	/* READ in row col nnz mod */
	SAFE_READ_DECL_V(m,uint32_t,fh);
	/* assert(m < MAT_ROW_BLK); */
	SAFE_READ_DECL_V(n,uint32_t,fh);
	SAFE_READ_DECL_V(mod,elemt_s,fh);
	assert(mod > 1);
	SAFE_READ_DECL_V(nnz,uint64_t,fh);

	fprintf(stderr," Mat is %u x %u - %lu (sparsity : %.3f%%) mod %lu\n",m,n,nnz,(double)(nnz)/((double)m*(double)n)*100.,(int64_t)mod);

	/* READ in ZO start */
	SAFE_READ_DECL_P(rows,m,uint32_t,fh);

	/* pol/zo correspondance */
	SAFE_READ_DECL_P(map_zo_pol,m,uint32_t,fh);


	/* colid in ZO */
	SAFE_READ_DECL_V(colid_size,uint64_t,fh);
	SAFE_READ_DECL_P(buffer,colid_size,uint32_t,fh); /* buffer has the matrix */

	/* size of nnz for pols: */
	SAFE_READ_DECL_V(pol_nb,uint32_t,fh);

	SAFE_READ_DECL_V(pol_nnz,uint64_t,fh);

	/* create GBpolynomials shared by A_init and B_init */
	SAFE_READ_DECL_P(pol_rows,pol_nb,uint32_t,fh);

	/* XXX what if elemt_s == elemt_t ??? */
	SAFE_READ_DECL_P(data_pol,pol_nnz,elemt_s,fh);


	fclose(fh);

	FILE * toto =fopen(nouv,"wb"); /* out */

	unsigned int m_ = m ;
	unsigned int n_ = n ;
	unsigned int p_ = mod;
	unsigned long long z_ = nnz ;
	fwrite(&m_,sizeof(unsigned int),1,toto);
	fwrite(&n_,sizeof(unsigned int),1,toto);
	fwrite(&p_,sizeof(unsigned int),1,toto);
	fwrite(&z_,sizeof(unsigned long long int),1,toto);

	SAFE_MALLOC_DECL(start,m+1,uint64_t);
	start[0] = 0 ;
	dimen_t i ;
	for ( i = 0 ; i < m ; ++i)
		start[i+1] = start[i] + rows[i];

	SAFE_MALLOC_DECL(start_pol,pol_nb+1,uint32_t);
	start_pol[0] = 0 ;
	for ( i = 0 ; i < pol_nb ; ++i)
		start_pol[i+1] = start_pol[i]+pol_rows[i];
	free(pol_rows);

	assert(start_pol[pol_nb] == pol_nnz);

	SAFE_MALLOC_DECL(colid,nnz,unsigned int); /* colid expands buffer */


	expandColid(buffer,colid_size,colid
#ifndef NDEBUG
			,nnz,n
#endif
		   );
	free(buffer);

	SAFE_MALLOC_DECL(data,nnz,short int);

	for (i = 0 ; i < m ; ++i) {
		dimen_t start_p = start_pol[ map_zo_pol[i] ] ;
		elemt_s * d = data_pol+start_p ;
		index_t jz ;
		for (jz = start[i] ; jz < start[i+1] ; ++jz) {
			data[jz] = (int16_t) *d ;
			++d ;
		}
	}

	free(map_zo_pol);

	SAFE_MALLOC_DECL(sz,m_,unsigned int);
	SAFE_MALLOC_DECL(pos,z_,unsigned int);
	SAFE_MALLOC_DECL(nz,z_,short int);

	for (i = 0 ; i < m_ ; ++i) sz[i] = rows[i] ;
	for (i = 0 ; i < z_ ; ++i) pos[i] = colid[i] ;
	for (i = 0 ; i < z_ ; ++i) nz[i] = data[i] ;

	fwrite(nz,sizeof(short int),z_,toto);
	fwrite(pos,sizeof(unsigned int),z_,toto);
	fwrite(sz,sizeof(unsigned int),m_,toto);

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

	if (ac < 2 || ac > 5) {
		usage(av[0]);
		exit(-1);
	}

	if ( (strcmp(av[1],"-h") == 0) || (strcmp(av[1],"-?") == 0) || (strcmp(av[1],"--help") == 0) ) {
		usage(av[0]);
		return(0);
	}

	if (strcmp(av[1],"-n") == 0) {
		new = 1 ;
	}

	FILE * titi ;
	char out[1024]; /* not too large the path... */


	if (ac-new == 3) {
		if ( strcmp(av[new+2],"-") != 0)  {
			usage(av[0]);
			exit(-1);
		}
	}

	titi = ouvrir(av[ac-1],"r"); /* in */

	strcpy(out,av[new+1]);

	if (new ==0) {
		printf("converting from old to new\n");
		convert_old2new(out,titi);
	}
	else {
		printf("converting from new to old\n");
		convert_new2old(out,titi);
	}

	return 0;
}

