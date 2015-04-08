/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 *                    Brice Boyer <brice.boyer@lip6.fr>
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




#include "io.h"
#include "tools/selecter.h"
#include "tools/types.h"
#include "tools/macros.h"
#include "tools/ouvrir.h"
#include "tools/tools.h"


/*  ========== TIMINGS ========== */
double walltime(struct timeval t_start) {
	struct timeval t_end;
	gettimeofday(&t_end, NULL);
	return (double)((t_end.tv_sec - t_start.tv_sec) * 1000000 + t_end.tv_usec - t_start.tv_usec);
}

/* #define USE_SEEK */

/*  ========== READING ========== */
sm_t *load_jcf_matrix(const char *fn, int verbose, int new_format) {

	/*  meta information of matrix */
	double density;
	double fs;
	char *fsu;
	/*  timing structs */
	struct timeval t_load_start;

	if (verbose > 1) {
		gettimeofday(&t_load_start, NULL);
	}

	/*  start loading the matrix */
	ri_t m;
	ci_t n;
	mod_t     mod;
	nnz_t     nnz;
	uint64_t  fl = 0;

	/*  open in binary mode first to get file size with fseek */

#ifdef USE_SEEK
	{
	  FILE*	fh        = fopen(fn,"rb");
	  if (fh == NULL) {
	    if (verbose > 0)
	      printf("File not found!\n");
	    return NULL;
	  }
	  else {
	    fseek(fh, 0L, SEEK_END);
	    fl  = ftell(fh);
	    fclose(fh);
	  }
	}
#endif /* USE_SEEK */

	/*  now read data from file */
	FILE* fh = ouvrir(fn,"r");

	if (new_format == 1) {
		uint32_t b;
		if (fread(&b, sizeof(uint32_t), 1, fh) != 1) {
			if (verbose > 0)
				printf("Error while reading file '%s'\n",fn);
			fclose(fh);
			return NULL;
		}

		if (( b & VERMASK ) != VERMASK) {
			if (verbose > 0)
				printf("Error while reading file '%s' (bad version)\n",fn);
			fclose(fh);
			return NULL;
		}

		b = b ^ VERMASK ;

		if ( b != Mjoin(select,elem_o)() || (sizeof(re_t) != sizeof(elem_o) || (re_t)(-1) < 0 ) ){
			if (verbose >0)
				printf("Error the data is not in expected representation\n");
			fclose(fh);
			return NULL;
		}

		fl += sizeof(uint32_t);
	}

	/*  get columns */
	if (fread(&m, sizeof(uint32_t), 1, fh) != 1) {
		if (verbose > 0)
			printf("Error while reading file '%s'\n",fn);
		fclose(fh);
		return NULL;
	}
	fl += sizeof(uint32_t);

	/*  get rows */
	if (fread(&n, sizeof(uint32_t), 1, fh) != 1) {
		if (verbose > 0)
			printf("Error while reading file '%s'\n",fn);
		fclose(fh);
		return NULL;
	}
	fl += sizeof(uint32_t);

	/*  get modulus */
	if (new_format == 1) {
		elemt_m mode;
		if ((fread(&mode, sizeof(elemt_m), 1, fh) != 1) || (mode == 1)) {
			if (verbose > 0)
				printf("Error while reading file '%s' (modulo)\n",fn);
			fclose(fh);
			return NULL;
		}
		mod = (mod_t) mode ; /* this is wrong in general */
		fl += sizeof(elemt_m);
	}
	else {
		mod_s mode ;
		if ((fread(&mode, sizeof(mod_s), 1, fh) != 1) || (mode == 1)) {
			if (verbose > 0)
				printf("Error while reading file '%s' (modulo)\n",fn);
			fclose(fh);
			return NULL;
		}
		mod = (mod_t) mode ; /* this is wrong in general */
		fl += sizeof(mod_s);
	}
	/*  get number of nonzero elements */
	if (fread(&nnz, sizeof(uint64_t), 1, fh) != 1) {
		if (verbose > 0)
			printf("Error while reading file '%s' (nnz)\n",fn);
		fclose(fh);
		return NULL;
	}
	fl += sizeof(uint64_t);

	/*  density of matrix */
	density =   (double) n * (double) m;
	density =   (double) (nnz) / density;
	density *=  100.0;
	/*  read entries from file */
	sm_t *M   = (sm_t *) malloc(sizeof(sm_t));
	M->rows   = (re_t **)malloc(m*sizeof(re_t *));
	M->pos    = (ci_t **)malloc(m*sizeof(ci_t *));
	M->rwidth = (ci_t *) malloc(m*sizeof(ci_t));

	/*  get meta data */
	M->nrows    = m;
	M->ncols    = n;
	M->nnz      = nnz;
	M->mod      = mod;
	M->density  = (float)density;

	if (new_format == 0) {

#ifndef USE_SEEK
		/*  maximal possible nonzero entries per row is n*sizeof(entry_t) */
		re_s *nze = (re_s *)malloc(nnz * sizeof(re_s));
		ci_t *pos = (ci_t *)malloc(nnz * sizeof(ci_t));
		ci_t *sz  = (ci_t *)malloc(m   * sizeof(ci_t));

		if (fread(nze, sizeof(re_s), nnz, fh) != nnz) {
			if (verbose > 0)
				printf("Error while reading file '%s' (data)\n",fn);
			free(M); /* XXX and many others */
			fclose(fh);
			return NULL ;
		}

		fl +=  nnz * sizeof(re_t) ;

		if (fread(pos,sizeof(ci_t),nnz,fh) != nnz) {
			if (verbose > 0)
				printf("Error while reading file '%s' (cols)\n",fn);
			free(M); /* XXX and many others */
			fclose(fh);
			return NULL ;
		}

		fl +=   nnz *sizeof(ci_t);

		if (fread(sz,sizeof(ci_t),m,fh) != m) {
			if (verbose > 0)
				printf("Error while reading file '%s' (rows)\n",fn);
			free(M); /* XXX and many others */
			fclose(fh);
			return NULL ;
		}

		fl +=  m *sizeof(ci_t);

		ri_t i;
		ci_t j;
		nnz_t here = 0 ;
		for (i = 0 ; i < m ; ++i) {
			M->rows[i] = (re_t *)malloc(sz[i] * sizeof(re_t));
			M->pos[i]  = (ci_t *)malloc(sz[i] * sizeof(ci_t));
			for (j = 0; j < sz[i]; ++j) {
				M->rows[i][j] = nze[here+(nnz_t)j];
				M->pos[i][j]  = pos[here+(nnz_t)j];
			}
			M->rwidth[i]  = sz[i];
			here += (nnz_t)sz[i] ;
		}
		free(sz);

#else /* USE SEEK */

		re_s *nze = (re_s *)malloc(n * sizeof(re_s));
		ci_t *pos = (ci_t *)malloc(n * sizeof(ci_t));


		/*  store header size of file */
		/*  size of m, n, mod and nb */
		uint32_t hs  = 3 * sizeof(uint32_t) + sizeof(uint64_t);

		/*  offsets for file handling */
		uint64_t row_size_offset  = nnz * sizeof(re_s) + nnz * sizeof(uint32_t) + hs;
		uint64_t row_val_offset   = hs;
		uint64_t row_pos_offset   = nnz * sizeof(re_s) + hs;

		ri_t i;
		ci_t j;
		ci_t sz;

		for (i = 0; i < m; ++i) {
			fseek(fh, row_size_offset, SEEK_SET);
			if (fread(&sz, sizeof(ci_t), 1, fh) != 1) {
				if (verbose > 0)
					printf("Error while reading file '%s' (nnz)\n",fn);
				free(M);
				fclose(fh);
				return NULL;
			}

			row_size_offset +=  sizeof(ci_t);

			fseek(fh, row_val_offset, SEEK_SET);
			if (fread(nze, sizeof(re_s), sz, fh) != sz) {
				if (verbose > 0)
					printf("Error while reading file '%s' (data)\n",fn);
				free(M);
				fclose(fh);
				return NULL;
			}

			row_val_offset +=  sz * sizeof(re_s);

			fseek(fh, row_pos_offset, SEEK_SET);
			if (fread(pos, sizeof(ci_t), sz, fh) != sz) {
				if (verbose > 0)
					printf("Error while reading file '%s' (cols)\n",fn);
				free(M);
				fclose(fh);
				return NULL;
			}

			row_pos_offset +=  sz * sizeof(ci_t);

			/*  reserve memory in matrix M for rows[i] */
			M->rows[i] = (re_t *)malloc(sz * sizeof(re_t));
			M->pos[i]  = (ci_t *)malloc(sz * sizeof(ci_t));
			for (j = 0; j < sz; ++j) {
				M->rows[i][j] = nze[j];
				M->pos[i][j]  = pos[j];
			}
			M->rwidth[i]  = sz;
		}

#endif /* USE_SEEK */

		free(nze);
		free(pos);
	}
	else { /* new_format == 1 */

		/* if ((sizeof(ci_t) != sizeof(uint32_t)) ||  ((ci_t)-1 < 0)) exit(-1); */

		dimen_t *row = (dimen_t *)malloc((m) * sizeof(dimen_t));
		if (fread(row, sizeof(dimen_t), m , fh) != m) {
			if (verbose > 0)
				printf("Error while reading file '%s' (rows)\n",fn);
			fclose(fh);
			return NULL;
		}
		fl +=  m*sizeof(dimen_t) ;

		dimen_t *mzp = (dimen_t*)malloc(m * sizeof(dimen_t));
		if (fread(mzp, sizeof(dimen_t), m , fh) != m) {
			if (verbose > 0)
				printf("Error while reading file '%s' (mat_zo_pol)\n",fn);
			fclose(fh);
			return NULL;
		}

		fl +=  m*sizeof(dimen_t) ;

		index_t czs;
		if (fread(&czs, sizeof(index_t), 1, fh) != 1) {
			if (verbose > 0)
				printf("Error while reading file '%s' (col size))\n",fn);
			fclose(fh);
			return NULL;
		}
		fl += sizeof(index_t) ;

		dimen_t * cz = (dimen_t*)malloc(czs * sizeof(dimen_t));
		if (fread(cz, sizeof(dimen_t), czs, fh) != czs) {
			if (verbose > 0)
				printf("Error while reading file '%s' (cols)\n",fn);
			fclose(fh);
			return NULL;
		}

		fl +=  sizeof(dimen_t)*(czs);

		dimen_t np;
		if (fread(&np, sizeof(dimen_t), 1, fh) != 1) {
			if (verbose > 0)
				printf("Error while reading file '%s' (nb pols)\n",fn);
			fclose(fh);
			return NULL;
		}
		fl += sizeof(dimen_t);

		index_t zp;
		if (fread(&zp, sizeof(index_t), 1, fh) != 1) {
			if (verbose > 0)
				printf("Error while reading file '%s' (nb pol data)\n",fn);
			fclose(fh);
			return NULL;
		}
		fl += sizeof(index_t);

		dimen_t * rp = (dimen_t*)malloc((np) * sizeof(dimen_t)); /* row length */
		if (fread(rp, sizeof(dimen_t), np, fh) != (np)) {
			if (verbose > 0)
				printf("Error while reading file '%s' (pol rows)\n",fn);
			fclose(fh);
			return NULL;
		}
		fl += sizeof(dimen_t)*(np);

		ri_t i;
		index_t * sp = (index_t*)malloc((np+1) * sizeof(index_t)); /* row pointers */
		sp[0] = 0 ;
		for ( i = 0 ; i < np ; ++i) {
			sp[i+1] = sp[i] + rp[i] ;
		}
		free(rp);


		/* BUG re_s depends on flag in first word  */

		re_s * vp = (re_s*)malloc((zp) * sizeof(re_s)) ;
		if (fread(vp, sizeof(re_s), zp, fh) != zp) {
			if (verbose > 0)
				printf("Error while reading file '%s' (pol data)\n",fn);
			fclose(fh);
			return NULL;
		}

		fl += sizeof(re_s)*(zp);

		dimen_t * pos = (dimen_t*)malloc(nnz * sizeof(dimen_t));

		expandColid( cz , czs , pos
#ifndef NDEBUG
		 , nnz, n
#endif
		);


		ci_t j;
		nnz_t here = 0 ;
		re_s *nze;
		for (i = 0 ; i < m ; ++i) {
			ci_t sz     = row[i];
			M->rows[i]  = (re_t *)malloc(sz * sizeof(re_t));
			M->pos[i]   = (ci_t *)malloc(sz * sizeof(ci_t));
			nze         = vp + sp[mzp[i]] ;
			for (j = 0; j < sz; ++j) {
				M->rows[i][j] = nze[j];
				M->pos[i][j]  = pos[here+(nnz_t)j];
			}
			M->rwidth[i]  = sz;
			here += (nnz_t)sz ;
		}

		/*  free data */
		free(pos);
		free(vp);
		free(sp);
		free(cz);
		free(mzp);
		free(row);
	}

	/*  file size of matrix */
	fs  = (double) fl / 1024 / 1024;
	fsu = "MB";
	if (fs > 1000) {
		fs  = fs / 1024;
		fsu = "GB";
	}

	M->fs       = (float)fs;
	M->fsu      = fsu;

	if (strcmp(fn,"-") !=0)
	  fclose(fh);
	return M;
}


/*  ========== WRITING ========== */

void write_jcf_matrix_to_file(sm_t *M, const char *fn, int verbose) {
}

void write_jcf_matrix_to_pbm(sm_t *M, const char *fn, int verbose) {
	char buffer[512];
	unsigned char out_byte  = 0;

	ri_t m = M->nrows;
	ci_t n = M->ncols;

	FILE *fh  = fopen(fn, "wb");

	/*  magic PBM header */
#ifdef __LP64__ /*  64bit machine */
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#else /*  32bit machine */
	sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), fh);

	ri_t i;
	ci_t j, k;
	/*  row width: number of nonzero elements in current row */
	ci_t sz;

	for (i = 0; i < m; ++i) {
		k   = 0;
		sz  = M->rwidth[i];
		for (j = 0; j < n; ++j) {
			if (k < sz && M->pos[i][k] == j) {
				out_byte  |=  (1 << (7 - (j % 8)));
				k++;
			} else {
				out_byte  &=  ~(1 << (7 - (j % 8)));
			}
			if (j % 8 == 7) {
				fwrite(&out_byte, sizeof(unsigned char), 1, fh);
				out_byte  = 0;
			}
		}
		if (j % 8 != 0)
			fwrite(&out_byte, sizeof(unsigned char), 1, fh);

		fflush(fh);
	}
	fclose(fh);
}

void print_mem_usage() {
	char    *unit = "KB";
	double  vms   = 0.0; /*  virtual memory size */
	double  rss   = 0.0; /*  resident set size */
	/*  possibly x86-64 is configured to use 2MB pages */
	long    page_size_kb  = sysconf(_SC_PAGE_SIZE) / 1024;

	unsigned long _vms = 0;
	long          _rss = 0;
	/*  get memory usage from 'file' stat which is not perfect, but gives most */
	/*  reliable information. */
	/*  Note: This corresponds to Martani's memory usage printing in his LELA */
	/*  implementation, thus it is used for comparison reasons. It might be changed */
	/*  in later versions of the gb. */
	const char *fn  = "/proc/self/stat";
	FILE *fh        = fopen(fn,"r");

	/*  dummy vars for leading entries in /proc/self/stat we are not interested in */
	char pid[1024] = "\0", comm[1024] ="\0", state[1024] ="\0", ppid[1024] ="\0", pgrp[1024] ="\0", session[1024] ="\0", tty_nr[1024] ="\0";
	char tpgid[1024] ="\0", flags[1024] ="\0", minflt[1024] ="\0", cminflt[1024] ="\0", majflt[1024] ="\0", cmajflt[1024] ="\0";
	char utime[1024] ="\0", stime[1024] ="\0", cutime[1024] ="\0", cstime[1024] ="\0", priority[1024] ="\0", nice[1024] ="\0";
	char nthrds[1024] ="\0", itrealvalue[1024] ="\0", starttime[1024] ="\0";

	/*  dummy reading of useless information */
	fscanf(fh, "%1023s", &pid[0]);
	fscanf(fh, "%1023s", &comm[0]);
	fscanf(fh, "%1023s", &state[0]);
	fscanf(fh, "%1023s", &ppid[0]);
	fscanf(fh, "%1023s", &pgrp[0]);
	fscanf(fh, "%1023s", &session[0]);
	fscanf(fh, "%1023s", &tty_nr[0]);
	fscanf(fh, "%1023s", &tpgid[0]);
	fscanf(fh, "%1023s", &flags[0]);
	fscanf(fh, "%1023s", &minflt[0]);
	fscanf(fh, "%1023s", &cminflt[0]);
	fscanf(fh, "%1023s", &majflt[0]);
	fscanf(fh, "%1023s", &cmajflt[0]);
	fscanf(fh, "%1023s", &utime[0]);
	fscanf(fh, "%1023s", &stime[0]);
	fscanf(fh, "%1023s", &cutime[0]);
	fscanf(fh, "%1023s", &cstime[0]);
	fscanf(fh, "%1023s", &priority[0]);
	fscanf(fh, "%1023s", &nice[0]);
	fscanf(fh, "%1023s", &nthrds[0]);
	fscanf(fh, "%1023s", &itrealvalue[0]);
	fscanf(fh, "%1023s", &starttime[0]);

	/*  get real memory information */
	fscanf(fh, "%lu", &_vms);
	fscanf(fh, "%ld", &_rss);

	/*  close file */
	fclose(fh);

	/*  TODO: How to read /proc/self/stat ??? */

	vms = _vms / 1024.0;
	rss = _rss * page_size_kb;

	/*  MB ? */
	if (vms > 1024) {
		vms   = vms/1024.0;
		rss   = rss/1024.0;
		unit  = "MB";
	}
	/*  GB ? */
	if (vms > 1024) {
		vms   = vms/1024.0;
		rss   = rss/1024.0;
		unit  = "GB";
	}
	/*  TB ? Just joking! */
	if (vms > 1024) {
		vms   = vms/1024.0;
		rss   = rss/1024.0;
		unit  = "TB";
	}
	printf("MMRY\tRSS - %.3f %s | VMS - %.3f %s\n", rss, unit, vms, unit);
}

/* STATIC STUFF */

 void compute_density_ml_submatrix(sm_fl_ml_t *A) {
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / __GB_NROWS_MULTILINE);
  uint32_t i, j;
  for (i=0; i<rlA; ++i) {
    if (A->ml[i].sz>0) {
      for (j=0; j<A->ml[i].sz; ++j) {
        if (A->ml[i].val[2*j] != 0) {
          A->nnz++;
        }
        if (A->ml[i].val[2*j+1] != 0) {
          A->nnz++;
        }
      }
    }
  }
 A->density  = compute_density(A->nnz, A->nrows, A->ncols);
}

 void compute_density_block_submatrix(sbm_fl_t *A) {
  const uint32_t rlA  = (uint32_t) ceil((float)A->nrows / A->bheight);
  const uint32_t clA  = (uint32_t) ceil((float)A->ncols / A->bwidth);
  uint32_t i, j, k, l;
  for (i=0; i<rlA; ++i) {
    for (j=0; j<clA; ++j) {
      if (A->blocks[i][j] != NULL) {
        for (k=0; k<A->bwidth/2; ++k) {
          for (l=0; l<A->blocks[i][j][k].sz; ++l) {
            if (A->blocks[i][j][k].val[2*l]) {
              A->nnz++;
            }
            if (A->blocks[i][j][k].val[2*l+1]) {
              A->nnz++;
            }
          }
        }
      }
    }
  }
 A->density  = compute_density(A->nnz, A->nrows, A->ncols);
}


/* vim:sts=2:sw=2:ts=2:
 */
