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



#ifndef __GB_tools_H
#define __GB_tools_H

#include "macros.h"


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
	index_t i = 0 ;
	index_t j = 0 ;
	dimen_t col ;
	for ( ; i < size_compressed ;) {
		col = compress[i++] ;
		if ((col & NEGMASK) == NEGMASK) {
			expand[j++] = col ^ NEGMASK ;
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

dimen_t compressColid(
		const dimen_t * cols_reord
		, const index_t nnz
		, dimen_t * colid
		)
{
	dimen_t here = 0 ;
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

	return here ;
}

#endif /* __GB_tools_H */
