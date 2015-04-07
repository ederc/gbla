AUTOMAKE_OPTIONS = gnu subdir-objects parallel-tests
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=${SIMD_FLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}

lib_LTLIBRARIES = libgbla.la

libgbla_la_SOURCES = \
	src/matrix.h	src/matrix.c \
	src/mapping.h	src/mapping.c \
	src/elimination.h	src/elimination.c

BUILT_SOURCES = src/gbla_config.h

pkgincludesubdir = $(includedir)/src
pkgincludesub_HEADERS = \
	src/types.h \
	src/matrix.h \
	src/mapping.h \
	src/elimination.h \
	tools/dump_matrix.h \
	tools/macros.h \
	tools/ouvrir.h  \
	tools/print_helper.h \
	tools/selecter.h \
	tools/tools.h \
	tools/types.h


nodist_pkgincludesub_HEADERS = src/gbla_config.h

EXTRA_DIST=autogen.sh tools/Makefile tools/converter.c tools/dump_matrix.c

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gbla.pc

libgbla_la_LDFLAGS = $(LDFLAGS) -release $(RELEASE) -no-undefined
libgbla_la_LIBADD = $(LIBPNG_LIBADD)

bin_PROGRAMS = gbla

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
gbla_SOURCES = \
	cli/gbla.h	cli/gbla.c	\
	cli/io.h	cli/io.c

gbla_LDADD = $(top_builddir)/libgbla.la $(DEPS_LIBS) -lm

SUBDIRS = test

check:
