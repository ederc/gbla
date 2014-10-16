#ifndef __GB_types_H
#define __GB_types_H

#define READ_MAT_ROW_BLOCK 128      // read matrix MAT_ROW_BLOCK by MAT_ROW_BLOCK
#define BLOCK_MAT_ROW 128      // write matrix MAT_ROW_BLOCK by MAT_ROW_BLOCK
#define GROW_REALLOC 16 // grow by GROW_REALLOC when using realloc
#define storage_t       int32_t  // Element representation mod p on file
#define index_t       int32_t  // indexing elements
#define element_t     int64_t // Element representation mod p in memory
#define integer_t     int32_t // modulo representation/storage

#define SAFE_MALLOC(ptr,size,elt) \
{ \
	elt * ptr = (elt *) malloc(size*sizeof(elt)); \
	assert(ptr); \
}

#define SAFE_CALLOC(ptr,size,elt) \
{ \
	elt * ptr = (elt *) calloc(size*sizeof(elt)); \
	assert(ptr); \
}



#define SAFE_REALLOC(ptr,size,elt) \
{ \
	ptr = realloc(ptr,size*sizeof(elt)); \
	assert(ptr); \
}


#endif // __GB_types_H
