/**
 * \file maps.h
 * \brief Interfaces for sparse matrix index maps used for subdividing the matrix
 * in Faugère-Lachartre style
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_MAPS_H
#define GB_MAPS_H

#include <sm.h>

/**
 * \brief Indexer for subdividing sparse matrix into 4 parts as described by
 * Faugère and Lachartre in http://dx.doi.org/10.1145/1837210.1837225.
 *
 *                 A | B
 * M     ---->     --+--
 *                 C | D
 */

typedef struct maps_fl_t {
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

  int idx_maps_done;      /*!<  0 if indexer for M is not constructed
                                1 if indexer for M is constructed */
  int nthrds;             /*!<  number of threads to be used for the indexer */
} maps_fl_t;


/**
 * \brief Constructs index maps for the subdivision of M into ABCD in the
 * Faugère-Lachartre style
 *
 *  \param Input matrix M
 *
 *  \return Corresponding index maps for ABCD decomposition
 */
maps_fl_t *construct_fl_maps(sm_t *M);

#endif
