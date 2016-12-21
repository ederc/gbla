AUTOMAKE_OPTIONS = gnu subdir-objects
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS	=	-std=gnu99 -march=native -Wextra -Wall -Wfloat-equal -Wundef \
						-Wpointer-arith -Wcast-align -Wstrict-prototypes \
						-Wmissing-prototypes -Wstrict-overflow=5 -Wwrite-strings \
						-Waggregate-return -Wswitch-default \
						-Wswitch-enum -Wconversion -Wunreachable-code \
						-Wno-unused-result -save-temps -fverbose-asm -masm=intel \
						-Werror-implicit-function-declaration -pedantic \
						-Wformat=2 -Wformat-y2k -Wno-format-extra-args \
						-Wno-format-zero-length -Wformat-nonliteral \
						-Wformat-security -Wno-unknown-pragmas -mpopcnt \
						-funroll-loops -finline-functions -fpic \
						${OPENMP_CFLAGS} ${CL_CFLAGS}

#AM_LIBS	=	${CL_LIBS}

lib_LTLIBRARIES = libgbla.la

libgbla_la_SOURCES = \
	src/matrix.h	src/matrix.c \
	src/mapping.h	src/mapping.c \
	src/elimination.h	src/elimination.c \
	src/cl-helper.h	src/cl-helper.c

BUILT_SOURCES = src/gbla_config.h

dev_stuff= \
	tools/dump_matrix.h \
	tools/macros.h \
	tools/ouvrir.h  \
	tools/print_helper.h \
	tools/selecter.h \
	tools/tools.h \
	tools/types.h\
	draft/matrix.h \
	draft/echelon.h


pkgincludesubdir = $(includedir)/gbla/src/
pkgincludesub_HEADERS = \
	src/types.h \
	src/matrix.h \
	src/mapping.h \
	src/elimination.h \
	src/cl-helper.h \
	src/cl-timing.h

pkginclude_HEADERS = config.h


EXTRA_DIST=autogen.sh tools/Makefile tools/converter.c tools/dump_matrix.c

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = gbla.pc

libgbla_la_LDFLAGS = $(LDFLAGS) -release $(RELEASE) -no-undefined

bin_PROGRAMS = gbla

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
gbla_SOURCES = \
	cli/gbla.h	cli/gbla.c	\
	cli/io.h	cli/io.c

gbla_LDADD = $(top_builddir)/libgbla.la ${CL_LIBS} $(DEPS_LIBS) -lm

SUBDIRS = 

TESTS = cli-tests poly-tests
