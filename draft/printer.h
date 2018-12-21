/* gbla: Gröbner Basis Linear Algebra
 * This file is part of gbla.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA, or see  <http://www.gnu.org/licenses/>. */


#ifndef __GBLA_printer_h
#define __GBLA_printer_h

static void print_double(double m)
{
	fprintf(stderr,"%f",m);
}
static void print_float(float m)
{
	fprintf(stderr,"%f",m);
}
static void print_int16_t(int16_t m)
{
	fprintf(stderr,"%d",m);
}
static void print_uint16_t(uint16_t m)
{
	fprintf(stderr,"%u",m);
}
static void print_int32_t(int32_t m)
{
	fprintf(stderr,"%d",m);
}
static void print_uint32_t(uint32_t m)
{
	fprintf(stderr,"%u",m);
}
static void print_int64_t(int64_t m)
{
	fprintf(stderr,"%ld",m);
}
static void print_uint64_t(uint64_t m)
{
	fprintf(stderr,"%lu",m);
}


#endif /* __GBLA_printer_h */

/* vim: set ft=c: */

