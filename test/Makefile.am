AUTOMAKE_OPTIONS = gnu subdir-objects parallel-tests
ACLOCAL_AMFLAGS = -I m4

AM_CFLAGS=-I. -I.. ${SIMD_CFLAGS} ${OPENMP_CFLAGS} ${DEBUG_FLAGS}
AM_LDFLAGS=-lm

TESTS =	\
	read-jcf-matrix-mod			\
	read-jcf-matrix-ncols		\
	read-jcf-matrix-nrows		\
	read-jcf-matrix-nnz			\
	read-jcf-matrix-density

check_PROGRAMS =	\
	read-jcf-matrix-mod			\
	read-jcf-matrix-ncols		\
	read-jcf-matrix-nrows		\
	read-jcf-matrix-nnz			\
	read-jcf-matrix-density

# set up the console program. Listing the headers in sources ensure that
# those files are included in distributions.
read_jcf_matrix_mod_SOURCES = \
	../cli/io.h	../cli/io.c		\
	read-jcf-matrix-mod.c
read_jcf_matrix_mod_LDADD = $(top_builddir)/libgbla.la $(DEPS_LIBS)

read_jcf_matrix_ncols_SOURCES = \
	../cli/io.h	../cli/io.c		\
	read-jcf-matrix-ncols.c
read_jcf_matrix_ncols_LDADD = $(top_builddir)/libgbla.la $(DEPS_LIBS)

read_jcf_matrix_nrows_SOURCES = \
	../cli/io.h	../cli/io.c		\
	read-jcf-matrix-nrows.c
read_jcf_matrix_nrows_LDADD = $(top_builddir)/libgbla.la $(DEPS_LIBS)

read_jcf_matrix_nnz_SOURCES = \
	../cli/io.h	../cli/io.c		\
	read-jcf-matrix-nnz.c
read_jcf_matrix_nnz_LDADD = $(top_builddir)/libgbla.la $(DEPS_LIBS)

read_jcf_matrix_density_SOURCES = \
	../cli/io.h	../cli/io.c		\
	read-jcf-matrix-density.c
read_jcf_matrix_density_LDADD = $(top_builddir)/libgbla.la $(DEPS_LIBS)
