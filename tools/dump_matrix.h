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

void Mjoin(dump,TEMPL_TYPE)(FILE * fh,int all,int strict,int magma)
{
	SAFE_READ_DECL_V(m,  uint32_t,fh);
	SAFE_READ_DECL_V(n,  uint32_t,fh);
	SAFE_READ_DECL_V(mod,TEMPL_TYPE,fh);

	SAFE_READ_DECL_V(nnz,uint64_t,fh);


	fprintf(stderr,"%u x %u matrix\n",m,n);
	/* fprintf(stderr,"mod %u\n",mod); */
	{
		double Nz=(double)(n)*(double)(m);
		Nz=(double)(nnz)/Nz;
		Nz*=100.0;
		fprintf(stderr,"Nb of Nz elements %lu (density %.2f)\n",nnz,Nz);
	}

	if (all)
	{
		uint32_t i;

		/* rows */
		SAFE_READ_DECL_P(rows, m, uint32_t, fh); /* length of each row on 32 bits*/

		SAFE_MALLOC_DECL(start,m+1,uint64_t);
		start[0] = 0 ;
		for ( i = 0 ; i < m ; ++i) {
			start[i+1] = start[i] + rows[i] ;
		}
		free(rows);

		/* map_zo_pol correspondance */
		SAFE_READ_DECL_P(map_zo_pol,m,uint32_t,fh); /* rows are numbered on 32 bits */

		/* if compressed, we need to know the size of colid */
		SAFE_READ_DECL_V(colid_size,uint64_t,fh); /* the size may be large (but < nnz/2) */

		fprintf(stderr,"%lu\n",colid_size);

		/* colid */
		SAFE_READ_DECL_P(colid,colid_size,uint32_t,fh); /*   columns are numbered on 32 bits */

		/* pol_nb */
		SAFE_READ_DECL_V(pol_nb,uint32_t,fh); /* number of polynomials */

		/* pol_nnz */
		SAFE_READ_DECL_V(pol_nnz,uint64_t,fh); /* size of polynomial data */

		/* pol_start */
		SAFE_READ_DECL_P(pol_rows,pol_nb,uint32_t,fh); /* length of each polynomial. less than number of rows */

		SAFE_MALLOC_DECL(pol_start,pol_nb+1,uint64_t);
		pol_start[0]=0;

		for ( i = 0 ; i < pol_nb ; ++i) {
			pol_start[i+1] = pol_start[i] + pol_rows[i] ;
		}
		free(pol_rows);

		assert(pol_nnz = pol_start[pol_nb]);

		/* pol_vals */
		SAFE_READ_DECL_P(pol_vals,pol_nnz,TEMPL_TYPE,fh); /* elements are int32_t (or anything else) */


		if (magma)
		{
#if 1
			Mjoin(print_mod,TEMPL_TYPE)(mod);
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

		uint32_t here = 0;
		for(i=0;i<m;i++) {
			/* pointer to the values of the polynomial corresponding to row i  */
			TEMPL_TYPE * pol_vals_begin = pol_vals + pol_start[map_zo_pol[i]];
			uint32_t v ; /* just the value */
			uint32_t j ;
			for ( j = start[i]  ; j < start[i+1] ; ) {
				/* NEGMASK flag (last bit set) says that column
				 * is next column for that line has 0 ;
				 * otherwise, next element is the number of
				 * consecutive columns with non zeros at this
				 * row. */
				uint32_t first = colid[here++] ;
				uint32_t repet = 1 ;
				if ((first & NEGMASK) == NEGMASK) {
					first ^= NEGMASK ; /* get the actual first column by unmasking */
				}
				else  {
					repet = colid[here++];
				}
				assert(first < n);
				assert(repet < n);
				uint32_t k = 0 ; /* C99 for this inside for :-( */
				for ( ; k < repet ; ++k) {
					v = pol_vals_begin[k]; /* consecutive values */

					Mjoin(print_line,TEMPL_TYPE)(i+1,first+k+1,v,magma);

				}
				j += repet ;
				/* assert something sur first */
				pol_vals_begin += repet ; /* jump to the next "first" */
			}
		}
		fprintf(stderr,"\n");

		free(pol_start);
		free(pol_vals);
		free(map_zo_pol);
		free(colid);
		free(start);
	}
	else {
		printf("%u  %u %lu %ld \n",m,n,nnz,(int64_t)mod);
	}

}

/* vim: set ft=c: */
