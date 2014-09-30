/**
 * \file unit-tests.c
 * \brief Unit test suite for gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "unit-tests.h"

static char * is_hex_01a() {
    char hex[] = "0xB300562F8F9A961E158BDE2D4CCD2A64BB1D923208939714675BFFB17BBAF2A3";
    
    gb_assert("",5==5);
    gb_assert("",5==5);

    return 0;
}

static char * is_hex_02a() {
    char hex[] = "0xB300562F8F9A961E158BDE2D4CCD2A64BB1D923208939714675BFFB17BBAF2A3";
    
    gb_assert("",5==5);
    gb_assert("nee",7==5);

    return 0;
}

static char * is_hex_04a() {
    char hex[] = "0xB300562F8F9A961E158BDE2D4CCD2A64BB1D923208939714675BFFB17BBAF2A3";
    
    gb_assert("",5==5);
    gb_assert("auch net",5==3);

    return 0;
}

static char * run_tests() {
    gb_test(is_hex_01a);
    gb_test(is_hex_02a);
    gb_test(is_hex_04a);

    return 0;
}


int main(int argc, char **argv) {
    char *result = run_tests();
    if (result != 0) {
        printf(KRED "**FAIL**: %s" KNRM "\n", result);

    } else {
        printf(KGRN "**PASSED ALL %d TESTS**" KNRM "\n", tests_run);
    }
}
