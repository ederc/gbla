/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
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

/**
 * \file types.h
 * \brief General typedefs
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_TYPES_H
#define GB_TYPES_H

/* #define GBLA_USE_DOUBLE */
/* #define GBLA_USE_AVX */

/*  #include <stdint.h> */

/** index type */
typedef uint32_t  mli_t;
/** block index type */
typedef uint16_t  bi_t;

/** storage type for entries */
typedef uint16_t  re_s;
/** storage type for mod */
typedef uint32_t mod_s;

#ifdef GBLA_USE_DOUBLE
/** matrix row entry type */
typedef double re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef double re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef double re_m_t;
/** type of field characteristic */
typedef double mod_t;
#else
/** matrix row entry type */
typedef uint16_t  re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef uint64_t  re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef uint32_t  re_m_t;
/** type of field characteristic */
typedef uint32_t  mod_t;
#endif
/** row and column index types */
typedef uint32_t  ci_t;
typedef uint32_t  ri_t;
/** number of nonzero elements type */
typedef uint64_t  nnz_t;


#define ALIGNT 32

/* field_ops.h */
#ifdef GBLA_USE_DOUBLE
#define MODP(a,b) \
	fmod((a),(b))
	/* static double MODP(double a, double b) { assert(a>=0) ; assert(b>0) ; double c = fmod(a,b) ; assert(c >=0) ; return c; } */
#define CAST(a) \
	(double) (a)
#else
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x000000000000ffff
#endif


#endif /* GB_TYPES_H */

/* vim:sts=2:sw=2:ts=2:ft=c:
 */
