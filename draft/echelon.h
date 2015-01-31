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

#ifndef __GB_echelon_H
#define __GB_echelon_H

dimen_t RowReduce_int32_t ( int32_t p, int32_t * A, dimen_t m, dimen_t n, dimen_t lda) ;
dimen_t RowReduce_double  ( double p, double * A, dimen_t m, dimen_t n, dimen_t lda, uint32_t nt) ;

dimen_t echelonD(
		GBMatrix_t      * A
		, DNS * D
		, uint32_t nt)
{

	dimen_t r = Mjoin(RowReduce,elemt_t)(D->mod,D->ptr,D->row,D->col,D->ld, nt);

	fprintf(stderr,"  -- residual rank    : %u\n",r);

	return r + A->row;
}


#endif /*  __GB_echelon_H */

/* vim: set ft=c: */
