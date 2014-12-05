
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

		/* start_zo */
		SAFE_READ_DECL_P(start_zo, m+1, uint64_t, fh); /* last element is an uint64_t. I don't beleive it is necessary to compress (as in row_length on uint16_t for instance but saving 3m/4, really ?) */

		/* map_zo_pol correspondance */
		SAFE_READ_DECL_P(map_zo_pol,m,uint32_t,fh); /* rows are numbered on 32 bits */

		/* if compressed, we need to know the size of colid_zo */
		SAFE_READ_DECL_V(colid_zo_size,uint64_t,fh); /* the size may be large (but < nnz/2) */

		/* colid_zo */
		SAFE_READ_DECL_P(colid_zo,colid_zo_size,uint32_t,fh); /*   columns are numbered on 32 bits */

		/* size_pol */
		SAFE_READ_DECL_V(size_pol,uint32_t,fh); /* number of coefficients may be large */

		/* start_pol */
		SAFE_READ_DECL_P(start_pol,size_pol+1,uint32_t,fh); /* last number is uint32_t.  */

		/* vals_pol */
		SAFE_READ_DECL_P(vals_pol,start_pol[size_pol],TEMPL_TYPE,fh); /* elements are int32_t (or anything else) */


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
			TEMPL_TYPE * vals_pol_begin = vals_pol + start_pol[map_zo_pol[i]];
			uint32_t v ; /* just the value */
			uint32_t j = start_zo[i] ; /* C99 inside for :-( */
			for (  ; j < start_zo[i+1] ; ) {
				/* NEGMASK flag (last bit set) says that column
				 * is next column for that line has 0 ;
				 * otherwise, next element is the number of
				 * consecutive columns with non zeros at this
				 * row. */
				uint32_t first = colid_zo[here++] ;
				uint32_t repet = 1 ;
				if ((first & NEGMASK) == NEGMASK) {
					first ^= NEGMASK ; /* get the actual first column by unmasking */
				}
				else  {
					repet = colid_zo[here++];
				}
				assert(first < n);
				assert(repet < n);
				uint32_t k = 0 ; /* C99 for this inside for :-( */
				for ( ; k < repet ; ++k) {
					v = vals_pol_begin[k]; /* consecutive values */

					Mjoin(print_line,TEMPL_TYPE)(i+1,first+k+1,v,magma);

				}
				j += repet ;
				/* assert something sur first */
				vals_pol_begin += repet ; /* jump to the next "first" */
			}
		}
		fprintf(stderr,"\n");
	}

}

/* vim: set ft=c: */
