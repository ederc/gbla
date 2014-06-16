/**
 * \file mapping.h
 * \brief Interfaces for sparse matrix index maps used for subdividing the matrix
 * in Faugère-Lachartre style
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_MAPPING_H
#define GB_MAPPING_H

#include <gb_config.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <matrix.h>

/**
 * \brief Indexer for subdividing sparse matrix into 4 parts as described by
 * Faugère and Lachartre in http://dx.doi.org/10.1145/1837210.1837225.
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 */

typedef struct map_fl_t {
  ci_t *pc;               /*!<  map of pivot columns: from input matrix M to
                                submatrix A; has length M->ncols, maps non-pivot
                                columns to __GB_MINUS_ONE_32 */
  ci_t *npc;              /*!<  map of non-pivot columns: from input matrix M to
                                submatrix B; has length M->ncols, maps pivot
                                columns to __GB_MINUS_ONE_32 */
  ci_t *pc_rev;           /*!<  reverse map of pivot columns: from submatrices
                                A and B to input matrix M */
  ci_t *npc_rev;          /*!<  reverse map of non-pivot columns: from submatrices
                                A and B to input matrix M */
  ri_t *pr_idxs_by_entry; /*!<  has length M->nrows, maps pivot columns to
                                their corresponding row index, maps non-pivot
                                columns to __GB_MINUS_ONE_32 */
  ri_t *npr_idxs;         /*!<  indexes of non-pivot rows */

  int nthrds;             /*!<  number of threads to be used for the indexer */
} map_fl_t;


/**
 * \brief Initializes a map for the input matrix M with all entries set
 * to __GB_MINUS_ONE_8.
 *
 * \param input matrix M
 * \return initialized map
 */
static inline map_fl_t *init_fl_map(sm_t *M) {
  map_fl_t *map   = NULL;
  map             = (map_fl_t *)malloc(sizeof(map_fl_t));
  
  // initialize map arrays and 
  // set initial values to __GB_MINUS_ONE_8
  map->pc = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->pc, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));

  map->npc  = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->npc, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));
  
  map->pc_rev = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->pc_rev, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));
  
  map->npc_rev  = (ci_t *)malloc(M->ncols * sizeof(ci_t));
  memset(map->npc_rev, __GB_MINUS_ONE_8, M->ncols * sizeof(ci_t));
  
  map->pr_idxs_by_entry  = (ri_t *)malloc(M->nrows * sizeof(ri_t));
  memset(map->pr_idxs_by_entry, __GB_MINUS_ONE_8, M->nrows * sizeof(ri_t));
  
  map->npr_idxs = (ri_t *)malloc(M->nrows * sizeof(ri_t));
  memset(map->npr_idxs, __GB_MINUS_ONE_8, M->nrows * sizeof(ri_t));
}

/**
 * \brief Constructs an indexer map for a Faugère-Lachartre decomposition of the
 * input matrix M
 * 
 * \param input matrix M
 *
 * \return indexer map for M
 */
map_fl_t *construct_fl_map(sm_t *M);

/**
 * \brief Constructs index maps for the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *  \param input matrix M
 *
 *  \return corresponding index maps for ABCD decomposition
 */
void splice_fl_matrix(sm_t *M, sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D, map_fl_t *map);

#endif
