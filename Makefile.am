AUTOMAKE_OPTIONS = gnu subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libgb.la

libgb_la_SOURCES = \
	src/sparse_matrix.h	src/sparse_matrix.c

BUILT_SOURCES = src/gb_config.h

pkgincludesubdir = $(includedir)/src
pkgincludesub_HEADERS = \
	src/sparse_matrix.h

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
	cli/gb.h	cli/gb.c \
	cli/io.h	cli/io.c 

gb_LDADD = $(top_builddir)/libgb.la $(DEPS_LIBS)
