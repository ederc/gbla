AUTOMAKE_OPTIONS = gnu subdir-objects parallel-tests
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS} -lm
	
TESTS = unit_tests unit_tests2
XFAIL_TESTS = unit_tests2
check_PROGRAMS = unit_tests unit_tests2

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
unit_tests_SOURCES = \
	unit-tests.h \
	unit-tests.c

unit_tests_LDADD = $(top_builddir)/libgb.la $(DEPS_LIBS)

unit_tests2_SOURCES = \
	unit-tests2.h \
	unit-tests2.c

unit_tests2_LDADD = $(top_builddir)/libgb.la $(DEPS_LIBS)
