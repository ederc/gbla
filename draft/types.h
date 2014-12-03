#ifndef __GB_types_H
#define __GB_types_H


#ifndef TYPE
#define TYPE uint16_t
#endif

/* #define READ_MAT_ROW_BLOCK 128      |+ read matrix MAT_ROW_BLOCK by MAT_ROW_BLOCK +| */
#define MAT_ROW_BLOCK 128000      /* write matrix MAT_ROW_BLOCK by MAT_ROW_BLOCK */
#define GROW_REALLOC 16 /* grow by GROW_REALLOC when using realloc */
#define storage_t       int32_t  /* Element representation mod p on file */
#define index_t       int32_t  /* indexing elements */
#define element_t     int64_t /* Element representation mod p in memory */
#define integer_t     int32_t /* modulo representation/storage */


#define NEGMASK (1U<<31)

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
	assert(fread(&(val),sizeof(elt),1,file)==1)

#define SAFE_READ_DECL_V(val,elt,file) \
	elt val ; \
	SAFE_READ_V(val,elt,file)

#define SAFE_READ_P(val,size,elt,file) \
	SAFE_MALLOC(val,size,elt); \
	assert(fread(val,sizeof(elt),(size),(file))==(size))

#define SAFE_READ_DECL_P(val,size,elt,file) \
	SAFE_MALLOC_DECL(val,size,elt); \
	assert(fread(val,sizeof(elt),(size),(file))==(size))



#endif /* __GB_types_H */
/* vim: set ft=c: */
