/**
 * \file unit-tests.c
 * \brief Unit test suite for gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include <stdio.h>
#include <assert.h>

static char * is_hex_03a() {
    char hex[] = "0xB300562F8F9A961E158BDE2D4CCD2A64BB1D923208939714675BFFB17BBAF2A3";
    
    assert(5==5);
    assert(3==5);

    return 0;
}

int main(int argc, char **argv) {
  is_hex_03a();
  return 0;
}
