AUTOMAKE_OPTIONS = gnu subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libgb.la

libgb_la_SOURCES = \
	src/gb.h	src/gb.c

BUILT_SOURCES = src/gb_config.h

pkgincludesubdir = $(includedir)/m4ri
pkgincludesub_HEADERS = \
	src/gb.h

nodist_pkgincludesub_HEADERS = src/gb_config.h

EXTRA_DIST=src/Doxyfile autogen.sh

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gb.pc

libgb_la_LDFLAGS = -release $(RELEASE) -no-undefined
libgb_la_LIBADD = $(LIBPNG_LIBADD)

bin_PROGRAMS = gb

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
gb_SOURCES = \
	src/cli/gb.h	src/cli/gb.c

gb_LDADD = $(top_builddir)/libgb.la $(DEPS_LIBS)

check_PROGRAMS=test_multiplication test_elimination test_trsm test_ple test_solve test_kernel test_random test_smallops test_transpose test_colswap test_invert test_misc
test_multiplication_SOURCES=testsuite/test_multiplication.c
test_multiplication_LDFLAGS=-lm4ri -lm
test_multiplication_CFLAGS=$(AM_CFLAGS)

test_elimination_SOURCES=testsuite/test_elimination.c
test_elimination_LDFLAGS=-lm4ri -lm
test_elimination_CFLAGS=$(AM_CFLAGS)

test_trsm_SOURCES=testsuite/test_trsm.c
test_trsm_LDFLAGS=-lm4ri -lm
test_trsm_CFLAGS=$(AM_CFLAGS)

test_ple_SOURCES=testsuite/test_ple.c
test_ple_LDFLAGS=-lm4ri -lm
test_ple_CFLAGS=$(AM_CFLAGS)

test_solve_SOURCES=testsuite/test_solve.c
test_solve_LDFLAGS=-lm4ri -lm
test_solve_CFLAGS=$(AM_CFLAGS)

test_kernel_SOURCES=testsuite/test_kernel.c
test_kernel_LDFLAGS=-lm4ri -lm
test_kernel_CFLAGS=$(AM_CFLAGS)

test_random_SOURCES=testsuite/test_random.c
test_random_LDFLAGS=-lm4ri -lm
test_random_CFLAGS=$(AM_CFLAGS)

test_smallops_SOURCES=testsuite/test_smallops.c testsuite/testing.c testsuite/testing.h
test_smallops_LDFLAGS=-lm4ri -lm
test_smallops_CFLAGS=$(AM_CFLAGS)

test_transpose_SOURCES=testsuite/test_transpose.c
test_transpose_LDFLAGS=-lm4ri -lm
test_transpose_CFLAGS=$(AM_CFLAGS)

test_colswap_SOURCES=testsuite/test_colswap.c
test_colswap_LDFLAGS=-lm4ri -lm
test_colswap_CFLAGS=$(AM_CFLAGS)

test_invert_SOURCES=testsuite/test_invert.c
test_invert_LDFLAGS=-lm4ri -lm
test_invert_CFLAGS=$(AM_CFLAGS)

test_misc_SOURCES=testsuite/test_misc.c
test_misc_LDFLAGS=-lm4ri -lm
test_misc_CFLAGS=$(AM_CFLAGS)

TESTS = test_multiplication test_elimination test_trsm test_ple test_solve test_kernel test_random test_smallops test_transpose test_colswap test_invert test_misc

