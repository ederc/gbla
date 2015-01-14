#ifndef __GB_tools_H
#define __GB_tools_H

#include "macros.h"

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


void insert_sort_duo_rev(uint32_t * liste, uint32_t  size, uint32_t * copain)
{
	uint32_t d , c = 1 , t ;
	for ( ; c < size ; ++c) {
		d = c;
		while ( d > 0 && (liste)[d] > (liste)[d-1]) {
			SWAP((liste)[d],(liste)[d-1]);
			SWAP((copain)[d],(copain)[d-1]);
			d--;
		}
	}
}


#endif /* __GB_tools_H */
