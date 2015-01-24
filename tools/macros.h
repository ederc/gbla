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

#ifndef __GB_macros_H
#define __GB_macros_H

#define Mjoin(pre,nam) my_join(pre , nam)
#define my_join(pre, nam) pre ## _ ## nam

#define max(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a > _b ? _a : _b; })

#define min(a,b) \
	({ __typeof__ (a) _a = (a); \
	 __typeof__ (b) _b = (b); \
	 _a < _b ? _a : _b; })

#define SWAP(a,b)  \
	t = a ; \
	a = b ; \
	b = t

#define DIVIDE_INTO(x,y) (((x) + (y) - 1)/(y))
#define ALIGN(x) DIVIDE_INTO(x,UNRL)*UNRL


#ifndef NDEBUG
#define checkassert(a,b) assert(a==b);
#else
#define checkassert(a,b) a
#endif

#include <stdlib.h>
#include <stdalign.h>
#include <malloc.h>
#include <immintrin.h>


static void* new_malloc(size_t size)  {
    void *alignedPointer;
#ifndef NDEBUG
    int alignError = 0;

    alignError =
#endif
	    posix_memalign(&alignedPointer, 32, size);

#ifndef NDEBUG
    assert(!alignError);
#endif

    /* alignedPointer = aligned_alloc(32,size); */

    return alignedPointer;
}


#define erase(a,size,elt) \
{ \
/*	size_t i = 0 ; \
	for ( i = 0 ; i < size ; ++i) \
	    a[i] = (elt) 0 ; \
	    */ \
	memset(a,0,size*sizeof(elt)); \
}

#define SAFE_MALLOC(ptr,size,elt) \
	ptr = (elt *) new_malloc((size)*sizeof(elt)); \
	assert(ptr)


#define SAFE_CALLOC(ptr,size,elt) \
	ptr = (elt *) new_malloc((size)*sizeof(elt)); \
        erase((ptr),(size),elt); \
	assert(ptr)

#define SAFE_REALLOC(ptr,size,elt) \
	ptr = (elt *) realloc((ptr),(size)*sizeof(elt)); \
	assert(ptr)

	/* ptr = (elt *) new_realloc((ptr),(size)*sizeof(elt)); \ */
	 /* ptr = (elt *) realloc((ptr),(size)*sizeof(elt)); \ */
	/* assert(ptr && !( (uintptr_t)ptr % 32)) */

#define SAFE_MALLOC_DECL(ptr,size,elt) \
	elt * SAFE_MALLOC((ptr),(size),elt)

#define SAFE_CALLOC_DECL(ptr,size,elt) \
	elt * SAFE_CALLOC((ptr),(size),elt)

#define SAFE_READ_V(val,elt,file) \
	checkassert(fread(&(val),sizeof(elt),1,file),1)

#define SAFE_READ_DECL_V(val,elt,file) \
	elt val ; \
	SAFE_READ_V(val,elt,file)

#define SAFE_READ_P(val,size,elt,file) \
	SAFE_MALLOC(val,size,elt); \
	checkassert(fread(val,sizeof(elt),(size),(file)),(size))

#define SAFE_READ_DECL_P(val,size,elt,file) \
	SAFE_MALLOC_DECL(val,size,elt); \
	checkassert(fread(val,sizeof(elt),(size),(file)),(size))

#define MEMCPY(to,from,size) \
{ \
	uint32_t iii = 0 ; \
	for ( ; iii < (size) ; ++iii) {  \
		(to)[iii] = (from)[iii] ; \
	} \
}


#define MEMCPY_CVT(ptr_a,elm_a,ptr_b,nb) \
{  \
	uint32_t iiii = 0 ; \
	for ( ; iiii < (nb) ; ++iiii) { \
		(ptr_a)[iiii] = (elm_a) (ptr_b)[iiii]; \
	} \
}

#define SAFE_MEMCPY_CVT(ptr_a,elm_a,ptr_b,nb) \
	SAFE_MALLOC(ptr_a,(nb),elm_a); \
        MEMCPY_CVT(ptr_a,elm_a,ptr_b,(nb))

#ifdef SIMD

#ifdef AVX
#define SET1  _mm256_set1_pd
#define STORE _mm256_storeu_pd
#define LOAD  _mm256_loadu_pd
#define ADD   _mm256_add_pd
#define MUL   _mm256_mul_pd
#define ELEM  __m256d
#endif /* AVX */

#ifdef SSE
#define SET1  _mm_set1_pd
#define STORE _mm_store_pd
#define LOAD  _mm_load_pd
#define ADD   _mm_add_pd
#define MUL   _mm_mul_pd
#define ELEM  __m128d
#endif /* SSE */

#define SET1_SIMD(a,b) \
	ELEM a =  SET1(b)

#define AXPY_SIMD(y,a,x) \
	STORE((y),  ADD(LOAD((y)), MUL((a), LOAD((x)))))
/* can do better if FMA */

#define COPY_SIMD(y,x) \
	STORE(y, LOAD(x))

#endif /* SIMD */

#endif /* __GB_macros_H */

/* vim: set ft=c: */
