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


#ifndef __GB_field_ops_H
#define __GB_field_ops_H

#include <math.h>

/* double */

double fmod_double(double a, double p)
{
	double t = fmod(a,p);
	if (t<0.) t += p ;
	return t;
}
void fmodin_double(double * a, double p)
{
	*a = fmod(*a,p);
	if (*a<0.) *a += p ;
	return;
}

double invert_double(double a, double p)
{
	int x_int, y_int, tx, ty;
	x_int = (int32_t) (p);
	y_int = (int32_t) (a);
	tx = 0;
	ty = 1;

	while (y_int != 0) {
		int q, temp;
		q = x_int / y_int; /* integer quotient */
		temp = y_int; y_int = x_int - q * y_int;
		x_int = temp;
		temp = ty; ty = tx - q * ty;
		tx = temp;
	}

	a = (double) tx ;
	if (a<0.) a += p ;
	return a ;
}

/* int32_t */

int32_t fmod_int32_t(int32_t a, int32_t p)
{
	return a%p;
}
void fmodin_int32_t(int32_t * a, int32_t p)
{
	 *a = *a%p;
	 assert(*a>=0);
	 return;
}

int32_t invert_int32_t(int32_t a, int32_t p)
{
	assert(a>0);
	if (a == 1 )
		return a;
	int32_t nouv = 1, old = 0, q = p,r,h ;
	int pos = 0 ;
	while(a > 0) {
		r = q%a ;
		q = q/a ;
		h = q*nouv+old ;
		old = nouv ;
		nouv = h;
		q = a ;
		a = r ;
		pos = !pos ;
	}
	return pos ? old : (p-old);
}


#endif /* __GB_field_ops_H */
