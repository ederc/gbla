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


#ifndef NDEBUG
#define checkassert(a,b) assert(a==b);
#else
#define checkassert(a,b) a
#endif

#define SAFE_MALLOC(ptr,size,elt) \
	ptr = (elt *) malloc((size)*sizeof(elt)); \
	assert(ptr)

#define SAFE_CALLOC(ptr,size,elt) \
	ptr = (elt *) calloc((size),sizeof(elt)); \
	assert(ptr)

#define SAFE_MALLOC_DECL(ptr,size,elt) \
	elt * SAFE_MALLOC(ptr,size,elt)

#define SAFE_CALLOC_DECL(ptr,size,elt) \
	elt * SAFE_CALLOC(ptr,size,elt)

#define SAFE_REALLOC(ptr,size,elt) \
	ptr = (elt *) realloc((ptr),(size)*sizeof(elt)); \
	assert(ptr)

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



#endif /* __GB_macros_H */

/* vim: set ft=c: */
