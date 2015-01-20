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
/* vim: set ft=c: */
