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

/* #define READ_MAT_ROW_BLK 128      |+ read matrix MAT_ROW_BLK by MAT_ROW_BLK +| */
/* #define MAT_ROW_BLK (1UL<<31)     |+ write matrix MAT_ROW_BLK by MAT_ROW_BLK +| */
#ifndef _OPENMP
#define MAT_SUB_BLK 16
#define MAT_ROW_BLK (MAT_SUB_BLK*8) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#else
#define MAT_SUB_BLK 16
#define MAT_ROW_BLK (MAT_SUB_BLK*32) /* write matrix MAT_ROW_BLK by MAT_ROW_BLK */
#endif
#define USE_B_SPARSE
#define DEROULE
#define UNRL 4
#undef USE_SAXPY
#define USE_SAXPY2
#undef USE_SAXPYn
/* #define GROW_REALLOC 16 |+ grow by GROW_REALLOC when using realloc +| */

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


