#ifndef __GB_types_H
#define __GB_types_H


#ifndef TYPE
#define elem_s uint16_t /* field element storage */
#define elem_t double   /* field element type    */
#endif

#define storage_t       int32_t  /* Element representation mod p on file */
#define index_t       uint64_t  /* indexing elements */
#define taille_t      uint32_t /* size_t */
/* #define element_t     int64_t |+ Element representation mod p in memory +| */
/* #define integer_t     int32_t |+ modulo representation/storage +| */


#define NEGMASK (1U<<31)
#define VERMASK (1U<<31)





#endif /* __GB_types_H */
/* vim: set ft=c: */
