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
	int32_t new = 1, old = 0, q = p,r,h ;
	int pos = 0 ;
	while(a > 0) {
		r = q%a ;
		q = q/a ;
		h = q*new+old ;
		old = new ;
		new = h;
		q = a ;
		a = r ;
		pos = !pos ;
	}
	return pos ? old : (p-old);
}


#endif /* __GB_field_ops_H */
