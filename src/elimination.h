/**
 * \file elimination.h
 * \brief Different Gaussian Elimination methods
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_ELIMINATION_H
#define GB_ELIMINATION_H

#include <mapping.h>

/**
 * \brief Elimination procedure which reduces the block submatrices A, B, C and
 * D correspondingly.
 *
 * \param block submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param block submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_block(sbm_fl_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D);

/**
 * \brief Elimination procedure which reduces the multiline submatrix A
 * and the block submatrices B, C and D.
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param block submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_ml_block(sm_fl_ml_t *A, sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D);

/**
 * \brief Elimination procedure which reduces the multiline submatrices A and C
 * and the block submatrices B and D.
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param multiline submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_ml(sm_fl_ml_t *A, sbm_fl_t *B, sm_fl_ml_t *C, sbm_fl_t *D);

#endif
