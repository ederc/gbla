/**
 * \file types.h
 * \brief General typedefs
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_TYPES_H
#define GB_TYPES_H

#include <stdint.h>

/// block index type
typedef uint16_t  bi_t;
/// matrix row entry type
typedef uint16_t  re_t;
/// row and column index types
typedef uint32_t  ci_t;
typedef uint32_t  ri_t;
/// number of nonzero elements type
typedef uint64_t  nnz_t;
/// tyoe of field characteristic
typedef uint32_t  mod_t;

#endif
