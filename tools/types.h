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


#ifndef __GB_types_H
#define __GB_types_H


#ifndef TYPE
#define elemt_s  int16_t /* field element storage */
#define elemt_t double   /* field element type    */
#endif

#define storage_t       int32_t  /* Element representation mod p on file */
#define index_t       uint64_t  /* indexing elements */
#define dimen_t       uint32_t /* size_t */
/* #define element_t     int64_t |+ Element representation mod p in memory +| */
/* #define integer_t     int32_t |+ modulo representation/storage +| */


#define NEGMASK (1U<<31)
#define VERMASK (1U<<31)





#endif /* __GB_types_H */
/* vim: set ft=c: */
