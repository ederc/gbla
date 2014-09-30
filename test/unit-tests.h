/**
 * \file unit-tests.h
 * \brief Unit test suite for gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_UNIT_TESTS_H
#define GB_UNIT_TESTS_H

#include <stdio.h>
#include <assert.h>

int tests_run = 0;
#define gb_assert(message, test) do { if (!(test)) return message; } while (0)
#define gb_test(test) do { char *message = test(); tests_run++; if (message) return message; } while (0)
#define KNRM "\x1B[0m"
#define KRED "\x1B[1;37;41m"
#define KGRN "\x1B[1;42;37m"

#endif
