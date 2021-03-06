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


#ifndef __GBLA_print_help_h
#define __GBLA_print_help_h

void print_line_float(uint32_t i, uint32_t j, float v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%f;\n",i,j,v);
	else
		printf("%u %u %f\n",i,j,v);
}
void print_line_double(uint32_t i, uint32_t j, double v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%f;\n",i,j,v);
	else
		printf("%u %u %f\n",i,j,v);
}

void print_line_int8_t(uint32_t i, uint32_t j, int8_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%d;\n",i,j,v);
	else
		printf("%u %u %d\n",i,j,v);
}
void print_line_uint8_t(uint32_t i, uint32_t j, uint8_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%d;\n",i,j,v);
	else
		printf("%u %u %d\n",i,j,v);
}

void print_line_int16_t(uint32_t i, uint32_t j, int16_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%d;\n",i,j,v);
	else
		printf("%u %u %d\n",i,j,v);
}
void print_line_uint16_t(uint32_t i, uint32_t j, uint16_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%d;\n",i,j,v);
	else
		printf("%u %u %d\n",i,j,v);
}

void print_line_int32_t(uint32_t i, uint32_t j, int32_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%d;\n",i,j,v);
	else
		printf("%u %u %d\n",i,j,v);
}
void print_line_uint32_t(uint32_t i, uint32_t j, uint32_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%d;\n",i,j,v);
	else
		printf("%u %u %d\n",i,j,v);
}

void print_line_int64_t(uint32_t i, uint32_t j, int64_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%ld;\n",i,j,v);
	else
		printf("%u %u %ld\n",i,j,v);
}
void print_line_uint64_t(uint32_t i, uint32_t j, uint64_t v, int magma)
{
	/* fprintf(stderr,"<%u>",szi); */
	if (magma)
		printf("A[%u,%u]:=%ld;\n",i,j,v);
	else
		printf("%u %u %ld\n",i,j,v);
}

void print_mod_int8_t(int8_t m)
{
	printf("K:=GF(%d);\n",m);
}
void print_mod_uint8_t(uint8_t m)
{
	printf("K:=GF(%u);\n",m);
}
void print_mod_int16_t(int16_t m)
{
	printf("K:=GF(%d);\n",m);
}
void print_mod_uint16_t(uint16_t m)
{
	printf("K:=GF(%u);\n",m);
}
void print_mod_int32_t(int32_t m)
{
	printf("K:=GF(%d);\n",m);
}
void print_mod_uint32_t(uint32_t m)
{
	printf("K:=GF(%u);\n",m);
}
void print_mod_int64_t(int64_t m)
{
	printf("K:=GF(%ld);\n",m);
}
void print_mod_uint64_t(uint64_t m)
{
	printf("K:=GF(%lu);\n",m);
}


#endif /* __GBLA_print_help_h */

/* vim: set ft=c: */
