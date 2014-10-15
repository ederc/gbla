#ifndef __GB_field_H
#define __GB_field_H


// a + b
Element_t * add( Element_t * a, Element_t * b, Integer_t * p) ;
// a - b
Element_t * sub( Element_t * a, Element_t * b, Integer_t * p) ;
// a * b
Element_t * mul( Element_t * a, Element_t * b, Integer_t * p) ;
// c + a * b
Element_t * fma( Element_t * c, Element_t * a, Element_t * b, Integer_t * p) ;

#endif // __GB_field_H
