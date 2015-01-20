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


#ifndef __GB_printer_h
#define __GB_printer_h

void print_double(double m)
{
	fprintf(stderr,"%f",m);
}
void print_float(float m)
{
	fprintf(stderr,"%f",m);
}
void print_int16_t(int16_t m)
{
	fprintf(stderr,"%d",m);
}
void print_uint16_t(uint16_t m)
{
	fprintf(stderr,"%u",m);
}
void print_int32_t(int32_t m)
{
	fprintf(stderr,"%d",m);
}
void print_uint32_t(uint32_t m)
{
	fprintf(stderr,"%u",m);
}
void print_int64_t(int64_t m)
{
	fprintf(stderr,"%ld",m);
}
void print_uint64_t(uint64_t m)
{
	fprintf(stderr,"%lu",m);
}


#endif /* __GB_printer_h */

/* vim: set ft=c: */

