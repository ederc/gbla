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



#ifndef __GB_selecter_h
#define __GB_selecter_h

uint32_t select_int8_t()   { return 1 | (0<<1) ; }
uint32_t select_uint8_t()  { return 0 | (0<<1) ; }
uint32_t select_int16_t()  { return 1 | (1<<1) ; }
uint32_t select_uint16_t() { return 0 | (1<<1) ; }
uint32_t select_int32_t()  { return 1 | (2<<1) ; }
uint32_t select_uint32_t() { return 0 | (2<<1) ; }
uint32_t select_int64_t()  { return 1 | (3<<1) ; }
uint32_t select_uint64_t() { return 0 | (3<<1) ; }

#endif
/* vim: set ft=c: */
